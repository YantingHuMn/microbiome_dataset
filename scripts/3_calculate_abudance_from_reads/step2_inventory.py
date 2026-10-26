#!/usr/bin/env python3
"""Inventory the downloaded raw reads, one row per BIOLOGICAL sample.

Reads only what is already on disk plus the ENA run report, and writes
results/3_download_raw_reads/raw_reads_inventory.csv.

Two traps this is built to avoid:
  * a run is not a sample - Castro-Nallar2015 submits each library as two
    single-end lane runs, so 28 runs are 14 samples;
  * R1/R2 are one sample, not two.
Both are resolved through the ENA sample_accession, not file names.

Assay type, host and environment come from the ENA run report. The amplicon
target is measured from the reads themselves (primer match in the first few
thousand reads) rather than assumed, and the evidence is recorded per row.

Usage: step2_inventory.py [--out-dir DIR] [--paper NAME]
"""

import argparse
import collections
import csv
import gzip
import io
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

SAVE_DIR = Path(
    "/hickory/proj/didonglab/dataset/virus/Database/data/abundance_matrix_precise_v2")
STATUS_CSV = Path(
    "/hickory/proj/didonglab/dataset/virus/Database/data/abudance_integration/"
    "1_2_status_raw_reads_only.csv")
OUT_DIR = Path(
    "/hickory/proj/didonglab/dataset/virus/Database/results/3_download_raw_reads")

ENA_FIELDS = ("run_accession,sample_accession,experiment_accession,study_accession,"
              "library_strategy,library_source,library_selection,library_layout,"
              "instrument_model,read_count,scientific_name,host,isolation_source,"
              "environment_material,sample_title")

# 5' primer prefixes; a hit in >50% of sampled reads is taken as that target
# still carrying its primer.
PRIMERS = [
    ("16S V1-V3 (27F)", "AGAGTTTGATC"),
    ("16S V3-V4 (341F)", "CCTACGGG"),
    ("16S V4 (515F)", "GTGYCAGC"),
    ("ITS1 (ITS1F)", "CTTGGTCATTTAG"),
    ("ITS2 (ITS3)", "GCATCGATGAAGAACGCAGC"),
]


def ena_runs(study: str) -> list[dict]:
    q = urllib.parse.urlencode({"accession": study, "result": "read_run",
                                "fields": ENA_FIELDS, "format": "tsv", "limit": "0"})
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + q
    with urllib.request.urlopen(url, timeout=90) as r:
        return list(csv.DictReader(io.StringIO(r.read().decode()), delimiter="\t"))


def primer_scan(path: Path, n: int = 3000) -> tuple[str, str]:
    """-> (target, evidence). Measured, not assumed."""
    seqs = []
    try:
        with gzip.open(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    seqs.append(line.strip())
                if len(seqs) >= n:
                    break
    except OSError as exc:
        return "unknown", f"could not read {path.name}: {exc}"
    if not seqs:
        return "unknown", "no reads"
    for name, prefix in PRIMERS:
        hits = sum(1 for s in seqs if s.startswith(prefix[:8]))
        if hits > len(seqs) * 0.5:
            return name, f"{100*hits//len(seqs)}% of {len(seqs)} reads start with {prefix[:8]}"
    common = collections.Counter(s[:10] for s in seqs).most_common(1)[0]
    return ("no primer detected",
            f"most common 5' 10-mer {common[0]} in {100*common[1]//len(seqs)}% of reads "
            f"(primers already trimmed, or not an amplicon library)")


def analysis_targets(strategy: str, source: str, target: str) -> tuple[str, str]:
    """What abundance this library can and cannot support."""
    s, src = strategy.upper(), source.upper()
    if s == "AMPLICON":
        if target.startswith("ITS"):
            return "fungal ITS ASV abundance", "no bacteria, no viruses (ITS is fungal-only)"
        return ("bacterial/archaeal 16S ASV abundance",
                "NO viral abundance - 16S does not amplify viruses; no fungi")
    if "METATRANSCRIPTOMIC" in src:
        return ("expressed taxa / RNA viruses (needs rRNA-depletion-aware handling)",
                "not comparable with DNA abundance; DNA phages only if transcribed")
    if s in ("WGS", "WGA") and "METAGENOMIC" in src:
        return ("bacterial + viral (vOTU) abundance from DNA shotgun",
                "viral fraction is bulk-derived, prophage-biased")
    return "unknown", "assay not resolved"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--save-dir", type=Path, default=SAVE_DIR)
    ap.add_argument("--status-csv", type=Path, default=STATUS_CSV)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--paper", action="append")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    status = {r["paper_folder"].strip(): r for r in
              csv.DictReader(args.status_csv.open(encoding="utf-8-sig"))}
    completed = [r["paper_folder"] for r in csv.DictReader(
        (args.save_dir / "download_completed.csv").open(encoding="utf-8-sig"))]
    if args.paper:
        completed = [p for p in completed if p in set(args.paper)]

    rows = []
    for paper in completed:
        manifest = args.save_dir / paper / "raw_reads" / "download_manifest.json"
        if not manifest.exists():
            print(f"  ! {paper}: no manifest, skipped", file=sys.stderr)
            continue
        man = json.loads(manifest.read_text(encoding="utf-8"))
        files = man["files"]
        study = sorted({f["study_accession"] for f in files})[0]
        try:
            meta = {r["run_accession"]: r for r in ena_runs(study)}
        except Exception as exc:
            print(f"  ! {paper}: ENA lookup failed ({exc})", file=sys.stderr)
            meta = {}

        first = args.save_dir / paper / "raw_reads" / files[0]["run_accession"] / \
            files[0]["url"].rsplit("/", 1)[-1]
        target, evidence = primer_scan(first)

        by_sample = collections.defaultdict(list)
        for f in files:
            by_sample[f["sample_accession"]].append(f)

        for sample, fs in sorted(by_sample.items()):
            runs = sorted({f["run_accession"] for f in fs})
            m = meta.get(runs[0], {})
            names = [f["url"].rsplit("/", 1)[-1] for f in fs]
            layout = (m.get("library_layout") or fs[0]["library_layout"] or "").upper()
            per_run = collections.Counter(f["run_accession"] for f in fs)
            if layout == "PAIRED":
                paired = ("complete" if all(v == 2 for v in per_run.values())
                          else f"INCOMPLETE: files per run {dict(per_run)}")
            else:
                paired = f"single-end, {len(runs)} run(s) per sample"
            on_disk = sum(
                1 for f in fs
                if (args.save_dir / paper / "raw_reads" / f["run_accession"] /
                    f["url"].rsplit("/", 1)[-1]).exists())
            can, cannot = analysis_targets(m.get("library_strategy", ""),
                                           m.get("library_source", ""), target)
            todo = []
            if len(runs) > 1:
                todo.append(f"{len(runs)} runs for one BioSample - concatenate, "
                            f"do not treat as separate samples")
            if not m:
                todo.append("ENA metadata unavailable; assay/host unverified")
            if (m.get("isolation_source") or "").lower() in ("missing", ""):
                todo.append("isolation_source not given in ENA")
            rows.append({
                "paper_folder": paper,
                "doi": status.get(paper, {}).get("doi", ""),
                "title": status.get(paper, {}).get("title", "")[:150],
                "study_accession": study,
                "biosample_accession": sample,
                "sample_title": m.get("sample_title", ""),
                "run_accessions": ";".join(runs),
                "n_runs_for_this_sample": len(runs),
                "n_fastq_files": len(fs),
                "n_fastq_on_disk": on_disk,
                "library_strategy": m.get("library_strategy", ""),
                "library_source": m.get("library_source", ""),
                "library_selection": m.get("library_selection", ""),
                "library_layout": layout,
                "instrument_model": m.get("instrument_model", ""),
                "read_count_run1": m.get("read_count", ""),
                "amplicon_target_measured": target,
                "target_evidence": evidence,
                "host": m.get("host", "") or m.get("scientific_name", ""),
                "isolation_source": m.get("isolation_source", ""),
                "environment_material": m.get("environment_material", ""),
                "fastq_dir": str(args.save_dir / paper / "raw_reads"),
                "fastq_files": ";".join(sorted(names)),
                "pairing_status": paired,
                "download_status": "complete" if on_disk == len(fs) else "INCOMPLETE",
                "existing_preprocessing": (
                    "none found on disk; primers still present"
                    if target not in ("no primer detected", "unknown")
                    else "none found on disk; primers absent (pre-trimmed or not amplicon)"),
                "usable_for": can,
                "not_usable_for": cannot,
                "info_source": "ENA filereport + local manifest + primer scan of reads",
                "open_questions": "; ".join(todo),
            })
        print(f"  {paper:<20} {len(by_sample):>4} samples, {len(files):>4} files, {target}",
              flush=True)

    out = args.out_dir / "raw_reads_inventory.csv"
    with out.open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n{len(rows)} biological samples from {len(completed)} papers -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
