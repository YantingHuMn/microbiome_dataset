#!/usr/bin/env python3
"""classify paper/dataset candidates before any data download.

`abundance_ready` here means "filename evidence suggests a processed
abundance product" -- it is NOT content-verified. A filename can be
generic (Supplementary_Table_1.xlsx) and still hold a real OTU table, or
can match ABUND and still be something else (an LEfSe/alpha-diversity
export). Rows that carry a table-shaped file but no decisive filename
match are routed to `needs_content_check` rather than silently merged
into `manual_review`, so stage6_verify_abundance.py has an explicit,
narrower queue to open and content-check instead of re-scanning everything.
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
import sys
csv.field_size_limit(sys.maxsize)

ABUND = re.compile(r"abundance|relative.abund|taxonomic.profile|taxonomy.*table|feature.table|otu.table|asv.table|zotu|\.biom(?:\.|$)|\.qza(?:\.|$)|bracken|kraken|metaphlan|centrifuge|kaiju|species.profile|genus.profile|phylum.profile|viral.profile|community.composition", re.I)
READS = re.compile(r"\.f(?:ast)?q(?:\.gz)?$|fastq_ftp|\.sra$", re.I)
TABLE_LIKE = re.compile(r"\.(?:xlsx|xls|csv|tsv|txt|tab)(?:!|$)", re.I)
SEQUENCE_REPO = {"bioproject", "sra_study", "sra_run", "sra_experiment", "mgnify_study", "mgnify_analysis", "gsa"}
PROCESSED_REPO = {"figshare", "zenodo", "dryad", "osf", "europepmc_supp", "mendeley"}
COMMUNITY_ASSAY = re.compile(r"16S_amplicon|shotgun_metagenomics|viral_metagenomics|metatranscriptomics", re.I)
WRONG = re.compile(r"isolate_or_single_genome|host_only_assay", re.I)

# Every downstream record ends up in exactly one of two buckets: a matrix
# someone already produced (processed), or reads that still need OUR
# pipeline (raw_reads). For the processed bucket we also want to know WHICH
# pipeline the original authors used, so a downstream harmonisation step can
# tell a QIIME2 relative-abundance table apart from a Kraken2 count table
# without re-deriving that from scratch.
DATA_CATEGORY = {
    "abundance_ready": "processed", "needs_content_check": "processed",
    "raw_reads_ready": "raw_reads", "raw_reads_needs_assay_check": "raw_reads",
    "manual_review": "unknown", "no_data_link_found": "unknown",
}
PIPELINE_PATTERNS = [
    ("QIIME2", re.compile(r"\bqiime\s?2\b|\.qza(?:\.|$)|\.qzv(?:\.|$)", re.I)),
    ("QIIME1", re.compile(r"\bqiime\b(?!\s?2)", re.I)),
    ("mothur", re.compile(r"\bmothur\b", re.I)),
    ("DADA2", re.compile(r"\bdada2\b", re.I)),
    ("USEARCH/UPARSE", re.compile(r"\busearch\b|\buparse\b", re.I)),
    ("VSEARCH", re.compile(r"\bvsearch\b", re.I)),
    ("MetaPhlAn", re.compile(r"\bmetaphlan\d?\b", re.I)),
    ("Kraken/Bracken", re.compile(r"\bkraken\d?\b|\bbracken\b", re.I)),
    ("Kaiju", re.compile(r"\bkaiju\b", re.I)),
    ("Centrifuge", re.compile(r"\bcentrifuge\b", re.I)),
    ("HUMAnN", re.compile(r"\bhumann\d?\b", re.I)),
    ("PICRUSt2", re.compile(r"\bpicrust\d?\b", re.I)),
    ("LotuS2", re.compile(r"\blotus\d?\b", re.I)),
    ("MEGAN", re.compile(r"\bmegan\b", re.I)),
    ("MG-RAST", re.compile(r"\bmg-?rast\b", re.I)),
    ("phyloseq", re.compile(r"\bphyloseq\b", re.I)),
]


def detect_pipelines(text: str) -> str:
    if not text:
        return ""
    hits = [name for name, rx in PIPELINE_PATTERNS if rx.search(text)]
    return ";".join(hits)


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
        table_like = bool(TABLE_LIKE.search(filetext))
        if abundance:
            status, reason = "abundance_ready", "processed_abundance_filename"
        elif reads:
            status, reason = "raw_reads_needs_assay_check", "sequence_repository_or_fastq"
        elif d.get("resolve_status") == "error":
            status, reason = "manual_review", "repository_metadata_query_failed"
        elif table_like and d["repository"] in PROCESSED_REPO:
            # A table-shaped file sits in a repository known to host processed
            # data, but its name alone doesn't confirm an abundance matrix.
            # Queue for content verification instead of dropping into the
            # same bucket as datasets with no lead at all.
            status, reason = "needs_content_check", "table_file_present_filename_inconclusive"
        else:
            status, reason = "manual_review", "no_decisive_filename_evidence"
        d.update({"has_abundance_hint": abundance, "has_raw_reads_hint": reads,
                  "has_table_like_file": table_like,
                  "candidate_status": status, "classification_reason": reason})
        drows.append(d)
        for p in d.get("paper_ids", "").split(";"):
            if p: bypaper[p].append(d)
    prows = []
    for p in csv.DictReader(open(args.papers, newline="", encoding="utf-8")):
        ds = bypaper.get(p["paper_id"], [])
        has_abund = any(x["candidate_status"] == "abundance_ready" for x in ds)
        has_reads = any(x["has_raw_reads_hint"] == "True" or x["has_raw_reads_hint"] is True for x in ds)
        needs_check = any(x["candidate_status"] == "needs_content_check" for x in ds)
        assay_ok = bool(COMMUNITY_ASSAY.search(p.get("assay_prediction", "")))
        wrong = bool(WRONG.search(p.get("screen_reason", "")))
        if has_abund:
            status, reason = "abundance_ready", "repository_file_listing_has_abundance_product"
        elif needs_check:
            status, reason = "needs_content_check", "table_file_present_filename_inconclusive"
        elif has_reads and assay_ok and not wrong:
            status, reason = "raw_reads_ready", "community_assay_and_public_sequence_accession"
        elif has_reads:
            status, reason = "raw_reads_needs_assay_check", "reads_found_but_assay_not_confirmed"
        elif ds:
            status, reason = "manual_review", "data_link_found_but_type_unresolved"
        else:
            status, reason = "no_data_link_found", "no_accession_or_data_url_extracted"
        pipelines = detect_pipelines(p.get("data_availability_text", ""))
        p.update({"n_linked_datasets": len(ds), "final_status": status, "final_reason": reason,
                  "dataset_ids": ";".join(x["dataset_id"] for x in ds),
                  "data_category": DATA_CATEGORY.get(status, "unknown"),
                  "source_pipeline_hint": pipelines})
        prows.append(p)
    out = Path(args.outdir)
    df = list(drows[0]) if drows else []
    pf = list(prows[0]) if prows else []
    write(out / "dataset_candidates_final.csv", drows, df)
    write(out / "paper_candidates_final.csv", prows, pf)
    write(out / "abundance_ready.csv", [r for r in prows if r["final_status"] == "abundance_ready"], pf)
    write(out / "raw_reads_ready.csv", [r for r in prows if r["final_status"] == "raw_reads_ready"], pf)
    write(out / "needs_content_check.csv", [r for r in prows if r["final_status"] == "needs_content_check"], pf)
    write(out / "manual_review.csv", [r for r in prows if r["final_status"] in {"manual_review", "raw_reads_needs_assay_check"}], pf)
    print(f"classified {len(prows)} papers and {len(drows)} datasets -> {out}")
    print(f"  abundance_ready={sum(r['final_status'] == 'abundance_ready' for r in prows)}  "
          f"needs_content_check={sum(r['final_status'] == 'needs_content_check' for r in prows)}  "
          f"raw_reads_ready={sum(r['final_status'] == 'raw_reads_ready' for r in prows)}")


if __name__ == "__main__":
    main()