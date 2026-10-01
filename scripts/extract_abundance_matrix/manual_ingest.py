#!/usr/bin/env python3
"""Feed a MANUALLY downloaded file (the ones NCBI PMC's proof-of-work
wall or any other block prevented stage7 from fetching itself) back into
the pipeline through the exact same sniff/normalise logic stage7 uses,
so the result lands in the same abundance-matrix / studies.tsv /
sample.tsv outputs -- no hand-copying values, no separate ad-hoc format.

Typical flow for a paper in blocked_manual_download.tsv:
  1. Open its `landing_url` (or the PMC article page, for europepmc_supp)
     in a real browser and download the file(s) named in `attempted_files`.
  2. Run this script once per downloaded file -- it re-uses stage6's
     sniff_frame/iter_tables/to_long_rows on that local file, and if it
     finds an abundance table, does everything stage7's automated path
     would have done: assigns/reuses the paper's study_id, writes (or
     APPENDS to, if the paper already has partial data from an earlier
     manually- or auto-ingested file) the wide matrix, and replaces the
     paper's "blocked" record in progress.jsonl with an "ok" one -- so
     the next `stage7 ... build` run picks it up into studies.tsv/
     sample.tsv automatically, exactly like an automated success would.

Usage:
    python3 manual_ingest.py \\
        --abundance-ready Database/results/stage5_final/abundance_ready.csv \\
        --datasets        Database/results/stage5_final/dataset_candidates_final.csv \\
        --progress        Database/results/stage7_extract/progress.jsonl \\
        --data-dir        Database/data \\
        --paper-id        "10.1002/ece3.11617" \\
        --dataset-id      "europepmc_supp:PMC12927239" \\
        --file            /path/to/your/manually_downloaded_Dataset_1.xlsx

If the file has no abundance-shaped table, the script says so and does
NOT touch progress.jsonl -- it never claims success it didn't verify.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import threading
from pathlib import Path

csv.field_size_limit(sys.maxsize)

sys.path.insert(0, str((Path(__file__).resolve().parent.parent / "find_papers")))
from stage6_verify_abundance import (          # noqa: E402
    iter_tables, sniff_frame, to_long_rows, retain_raw_file, append_raw_manifest,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage7_extract_abundance_matrix import (         # noqa: E402
    build_table_matrix, fetch_first_author_year, assign_study_id, MatrixTooLargeError,
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--abundance-ready", required=True)
    ap.add_argument("--datasets", required=True)
    ap.add_argument("--progress", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--raw-store", required=True,
                    help="the file you downloaded by hand is COPIED here (byte-for-byte, SHA256 "
                         "recorded) -- same permanent audit trail as an automated download.")
    ap.add_argument("--paper-id", required=True)
    ap.add_argument("--dataset-id", required=True, help="e.g. europepmc_supp:PMC12927239 "
                    "(must already exist as a row in --datasets, for provenance)")
    ap.add_argument("--file", required=True, help="path to the file you downloaded by hand")
    args = ap.parse_args()

    with open(args.abundance_ready, newline="", encoding="utf-8") as f:
        papers_by_id = {r["paper_id"]: r for r in csv.DictReader(f)}
    paper = papers_by_id.get(args.paper_id)
    if paper is None:
        print(f"FATAL: paper_id {args.paper_id!r} not found in {args.abundance_ready}", file=sys.stderr)
        sys.exit(1)

    local = Path(args.file)
    if not local.exists():
        print(f"FATAL: {local} does not exist", file=sys.stderr)
        sys.exit(1)

    repo = args.dataset_id.split(":", 1)[0] if ":" in args.dataset_id else ""
    acc = args.dataset_id.split(":", 1)[1] if ":" in args.dataset_id else ""
    manifest_row = retain_raw_file(local, Path(args.raw_store), args.paper_id, args.dataset_id,
                                   repo, acc, source_url="manually_downloaded_by_user",
                                   original_filename=local.name)
    append_raw_manifest(Path(args.raw_store) / "raw_files_manifest.csv", manifest_row)
    local = Path(manifest_row["retained_path"])
    print(f"retained a permanent copy at {local} (sha256={manifest_row['sha256'][:12]}...)")

    print(f"opening {local} (full read, not the row-capped preview) ...")
    import pandas as pd
    existing_records = []
    if Path(args.progress).exists():
        with open(args.progress, encoding="utf-8") as f:
            existing_records = [json.loads(l) for l in f if l.strip()]
    prior = next((r for r in existing_records if r["paper_id"] == args.paper_id), None)
    taken = {r.get("study_id") for r in existing_records if r.get("study_id")}
    if prior and prior.get("study_id"):
        study_id = prior["study_id"]
    else:
        surname, year, err = fetch_first_author_year(paper.get("pmid", ""), paper.get("doi", ""))
        if surname and year:
            study_id = assign_study_id(surname, year, taken)
        else:
            study_id = f"paper{args.paper_id.replace('/', '_').replace(':', '_')[:40]}"
        if err:
            print(f"note: author lookup failed ({err}) -- study_id falls back to {study_id!r}")

    data_dir = Path(args.data_dir) / study_id
    data_dir.mkdir(parents=True, exist_ok=True)
    n_checked = 0
    new_tables = []
    for sub_id, _hr, d in iter_tables(local, preview=False):
        n_checked += 1
        v = sniff_frame(d)
        if not v["is_abundance"]:
            continue
        source_file = f"{local.name}!{sub_id}"
        if v["needs_review"]:
            print(f"  {source_file}: weak taxonomic signal ({v.get('reason')}) -- "
                  "NOT auto-included; review it yourself before deciding.")
            continue
        rows = to_long_rows(d, v, args.paper_id, args.dataset_id, source_file)
        if rows.empty:
            continue
        try:
            wide, dup = build_table_matrix(rows)
        except MatrixTooLargeError as e:
            print(f"  {source_file}: refusing to build matrix -- {e}")
            continue
        if wide.empty:
            continue
        safe_name = __import__("re").sub(r"[^A-Za-z0-9._-]", "_", f"{args.dataset_id}__{source_file}")[:150]
        matrix_path = data_dir / f"{safe_name}.matrix.tsv"
        long_path = data_dir / f"{safe_name}.long.tsv.gz"
        wide.to_csv(matrix_path, sep="\t", index=False)
        rows.to_csv(long_path, sep="\t", index=False, compression="gzip")
        if not dup.empty:
            dup.to_csv(data_dir / f"{safe_name}.duplicates.tsv", sep="\t", index=False)
        blank_ids = sorted(rows.loc[rows["sample_flag"] == "blank_or_control", "sample_id"].unique().tolist())
        seq_hints = sorted({h for h in rows["sequencing_type_hint"] if h})
        new_tables.append({
            "dataset_id": args.dataset_id, "source_file": source_file,
            "matrix_path": str(matrix_path), "long_path": str(long_path),
            "n_sample_ids_raw": int(rows["sample_id"].nunique()),
            "n_verified_biological_samples": "",
            "n_taxa": int(wide.shape[1] - 2),
            "domains": ";".join(sorted(rows["domain"].unique())),
            "sequencing_type_hints": ";".join(seq_hints),
            "blank_or_control_sample_ids": ";".join(blank_ids),
            "n_duplicate_cells": int(dup.shape[0]) if not dup.empty else 0,
            "raw_file_sha256": manifest_row["sha256"],
        })
        print(f"  {source_file}: {wide.shape[0]} sample_ids x {wide.shape[1]-2} taxa -> {matrix_path}")

    if not new_tables:
        print(f"\nno abundance-shaped table accepted from {local.name} ({n_checked} sheet(s)/table(s) checked). "
              "progress.jsonl was NOT modified -- this file doesn't get counted as a success.")
        return

    # Merge into the paper's EXISTING record (never replace its other
    # already-extracted tables) -- this file ADDS tables, it never removes
    # or overwrites ones from a prior automated or manual run.
    if prior is not None:
        merged_tables = prior.get("tables", []) + new_tables
        new_rec = {**prior, "status": "ok", "study_id": study_id, "tables": merged_tables}
    else:
        new_rec = {"paper_id": args.paper_id, "status": "ok", "study_id": study_id,
                   "tables": new_tables, "needs_review": [], "blocked_files": [],
                   "author_lookup_error": ""}

    kept = [r for r in existing_records if r["paper_id"] != args.paper_id]
    kept.append(new_rec)
    with open(args.progress, "w", encoding="utf-8") as f:
        for r in kept:
            f.write(json.dumps(r) + "\n")

    print(f"\ndone. {len(new_tables)} new table(s) added for {args.paper_id} (study_id={study_id!r}).")
    print("run the build phase next to get this into studies.tsv/sample.tsv/table_manifest.tsv.")


if __name__ == "__main__":
    main()
