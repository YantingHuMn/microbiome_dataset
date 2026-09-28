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
from stage6_verify_abundance import iter_tables, sniff_frame, to_long_rows  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage7_extract_abundance_matrix import (         # noqa: E402
    build_wide_matrix, fetch_first_author_year, assign_study_id,
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--abundance-ready", required=True)
    ap.add_argument("--datasets", required=True)
    ap.add_argument("--progress", required=True)
    ap.add_argument("--data-dir", required=True)
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

    print(f"opening {local} ...")
    all_long = []
    n_checked = 0
    for _sub_id, _hr, d in iter_tables(local):
        n_checked += 1
        v = sniff_frame(d)
        if not v["is_abundance"]:
            continue
        rows = to_long_rows(d, v, args.paper_id, args.dataset_id, local.name)
        if not rows.empty:
            all_long.append(rows)

    if not all_long:
        print(f"no abundance-shaped table found in {local.name} ({n_checked} sheet(s)/table(s) checked). "
              "progress.jsonl was NOT modified -- this file doesn't get counted as a success.")
        return

    import pandas as pd
    long_df = pd.concat(all_long, ignore_index=True)

    # if this paper already has a partial "ok" record (some datasets already
    # succeeded automatically, this one didn't), merge with its prior long
    # rows so the manual file ADDS to, rather than replaces, existing data.
    existing_records = []
    if Path(args.progress).exists():
        with open(args.progress, encoding="utf-8") as f:
            existing_records = [json.loads(l) for l in f if l.strip()]
    prior = next((r for r in existing_records if r["paper_id"] == args.paper_id and r["status"] == "ok"), None)
    if prior and prior.get("matrix_path") and Path(prior["matrix_path"]).exists():
        old_wide = pd.read_csv(prior["matrix_path"], sep="\t")
        prior_path = prior["matrix_path"]
        print(f"paper already has a matrix from an earlier run ({prior_path}) "
              "-- merging this file's data into it rather than starting over.")
        study_id = prior["study_id"]
    else:
        old_wide = None
        taken = {r.get("study_id") for r in existing_records if r.get("study_id")}
        surname, year, err = fetch_first_author_year(paper.get("pmid", ""), paper.get("doi", ""))
        if surname and year:
            study_id = assign_study_id(surname, year, taken)
        else:
            study_id = f"paper{args.paper_id.replace('/', '_').replace(':', '_')[:40]}"
        if err:
            print(f"note: author lookup failed ({err}) -- study_id falls back to {study_id!r}")

    wide = build_wide_matrix(long_df)
    if old_wide is not None:
        wide = pd.concat([old_wide, wide], ignore_index=True).fillna(0)
        # a sample appearing in both (re-ingesting the same file twice) should
        # not be duplicated -- keep the newest values for any repeated sample_id.
        wide = wide.drop_duplicates(subset="sample_id", keep="last")

    data_dir = Path(args.data_dir); data_dir.mkdir(parents=True, exist_ok=True)
    out_path = data_dir / f"{study_id}_abundance_matrix.tsv"
    wide.to_csv(out_path, sep="\t", index=False)

    new_rec = {"paper_id": args.paper_id, "status": "ok", "study_id": study_id,
               "n_samples": int(wide.shape[0]), "n_taxa": int(wide.shape[1] - 2),
               "matrix_path": str(out_path), "blocked_files": [],
               "author_lookup_error": "", "domains": sorted(long_df["domain"].unique().tolist()),
               "notes": f"manually ingested from {local.name} via manual_ingest.py"}

    kept = [r for r in existing_records if r["paper_id"] != args.paper_id]
    kept.append(new_rec)
    with open(args.progress, "w", encoding="utf-8") as f:
        for r in kept:
            f.write(json.dumps(r) + "\n")

    print(f"\ndone. {wide.shape[0]} samples x {wide.shape[1]-2} taxa written to {out_path}")
    print(f"progress.jsonl updated: {args.paper_id} is now status=ok, study_id={study_id!r}")
    print("run the build phase next to get this into studies.tsv/sample.tsv.")


if __name__ == "__main__":
    main()
