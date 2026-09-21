#!/usr/bin/env python3
"""classify paper/dataset candidates before any data download."""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

ABUND = re.compile(r"abundance|taxonomic.profile|feature.table|otu.table|asv.table|\.biom(?:\.|$)|\.qza(?:\.|$)|bracken|kraken|metaphlan|centrifuge|kaiju|species.profile|genus.profile|viral.profile", re.I)
READS = re.compile(r"\.f(?:ast)?q(?:\.gz)?$|fastq_ftp|\.sra$", re.I)
SEQUENCE_REPO = {"bioproject", "sra_study", "sra_run", "sra_experiment", "mgnify_study", "mgnify_analysis", "gsa"}
COMMUNITY_ASSAY = re.compile(r"16S_amplicon|shotgun_metagenomics|viral_metagenomics|metatranscriptomics", re.I)
WRONG = re.compile(r"isolate_or_single_genome|host_only_assay", re.I)


def write(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)


def main() -> None:
    base = Path(__file__).resolve().parent / "results"
    ap = argparse.ArgumentParser()
    ap.add_argument("--papers", default=str(base / "stage3_links/papers_with_data_text.csv"))
    ap.add_argument("--datasets", default=str(base / "stage4_datasets/datasets_master.csv"))
    ap.add_argument("--outdir", required=True, help="Directory for stage 5 outputs")
    args = ap.parse_args()
    datasets = list(csv.DictReader(open(args.datasets, newline="", encoding="utf-8")))
    bypaper = defaultdict(list)
    drows = []
    for d in datasets:
        filetext = d.get("listed_files", "")
        abundance = bool(ABUND.search(filetext))
        reads = bool(READS.search(filetext)) or d["repository"] in SEQUENCE_REPO
        if abundance:
            status, reason = "abundance_ready", "processed_abundance_filename"
        elif reads:
            status, reason = "raw_reads_needs_assay_check", "sequence_repository_or_fastq"
        elif d.get("resolve_status") == "error":
            status, reason = "manual_review", "repository_metadata_query_failed"
        else:
            status, reason = "manual_review", "no_decisive_filename_evidence"
        d.update({"has_abundance_hint": abundance, "has_raw_reads_hint": reads,
                  "candidate_status": status, "classification_reason": reason})
        drows.append(d)
        for p in d.get("paper_ids", "").split(";"):
            if p: bypaper[p].append(d)
    prows = []
    for p in csv.DictReader(open(args.papers, newline="", encoding="utf-8")):
        ds = bypaper.get(p["paper_id"], [])
        has_abund = any(x["candidate_status"] == "abundance_ready" for x in ds)
        has_reads = any(x["has_raw_reads_hint"] == "True" or x["has_raw_reads_hint"] is True for x in ds)
        assay_ok = bool(COMMUNITY_ASSAY.search(p.get("assay_prediction", "")))
        wrong = bool(WRONG.search(p.get("screen_reason", "")))
        if has_abund:
            status, reason = "abundance_ready", "repository_file_listing_has_abundance_product"
        elif has_reads and assay_ok and not wrong:
            status, reason = "raw_reads_ready", "community_assay_and_public_sequence_accession"
        elif has_reads:
            status, reason = "raw_reads_needs_assay_check", "reads_found_but_assay_not_confirmed"
        elif ds:
            status, reason = "manual_review", "data_link_found_but_type_unresolved"
        else:
            status, reason = "no_data_link_found", "no_accession_or_data_url_extracted"
        p.update({"n_linked_datasets": len(ds), "final_status": status, "final_reason": reason,
                  "dataset_ids": ";".join(x["dataset_id"] for x in ds)})
        prows.append(p)
    out = Path(args.outdir)
    df = list(drows[0]) if drows else []
    pf = list(prows[0]) if prows else []
    write(out / "dataset_candidates_final.csv", drows, df)
    write(out / "paper_candidates_final.csv", prows, pf)
    write(out / "abundance_ready.csv", [r for r in prows if r["final_status"] == "abundance_ready"], pf)
    write(out / "raw_reads_ready.csv", [r for r in prows if r["final_status"] == "raw_reads_ready"], pf)
    write(out / "manual_review.csv", [r for r in prows if r["final_status"] in {"manual_review", "raw_reads_needs_assay_check"}], pf)
    print(f"classified {len(prows)} papers and {len(drows)} datasets -> {out}")


if __name__ == "__main__":
    main()