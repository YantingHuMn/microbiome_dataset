#!/usr/bin/env python3
"""resolve accessions to landing/download locations; never download data."""
from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

UA = {"User-Agent": "microbiome-dataset-resolver/1.0 (academic research)"}
TEMPLATES = {
    "bioproject": "https://www.ncbi.nlm.nih.gov/bioproject/{a}",
    "sra_study": "https://www.ncbi.nlm.nih.gov/sra/?term={a}", "sra_run": "https://www.ncbi.nlm.nih.gov/sra/{a}",
    "sra_experiment": "https://www.ncbi.nlm.nih.gov/sra/?term={a}", "biosample": "https://www.ncbi.nlm.nih.gov/biosample/{a}",
    "mgnify_study": "https://www.ebi.ac.uk/metagenomics/studies/{a}", "mgnify_analysis": "https://www.ebi.ac.uk/metagenomics/analyses/{a}",
    "arrayexpress": "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{a}", "biostudies": "https://www.ebi.ac.uk/biostudies/studies/{a}",
    "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={a}", "gsa": "https://ngdc.cncb.ac.cn/gsa/browse/{a}",
    "omix": "https://ngdc.cncb.ac.cn/omix/release/{a}", "zenodo": "https://zenodo.org/records/{a}",
    "figshare": "https://figshare.com/articles/dataset/_/{a}", "dryad": "https://datadryad.org/dataset/doi:10.5061/dryad.{a}",
    "mendeley": "https://data.mendeley.com/datasets/{a}",
}
ENA_TYPES = {"bioproject", "sra_study", "sra_run", "sra_experiment"}


def get_json(url: str) -> tuple[dict | list | None, str]:
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
            return json.loads(r.read()), ""
    except Exception as e: return None, str(e)


def list_record(repo: str, acc: str) -> tuple[list[str], str, str]:
    """Return metadata file names, API endpoint, error. No deposited file content."""
    if repo == "zenodo":
        api = f"https://zenodo.org/api/records/{acc}"; d, err = get_json(api)
        return [x.get("key", "") for x in (d or {}).get("files", [])], api, err
    if repo == "figshare":
        api = f"https://api.figshare.com/v2/articles/{acc}"; d, err = get_json(api)
        return [x.get("name", "") for x in (d or {}).get("files", [])], api, err
    if repo == "mgnify_study":
        api = f"https://www.ebi.ac.uk/metagenomics/api/v1/studies/{acc}"; d, err = get_json(api)
        return [], api, err
    if repo in ENA_TYPES:
        result = {"bioproject": "read_run", "sra_study": "read_run", "sra_run": "read_run", "sra_experiment": "read_run"}[repo]
        fields = "study_accession,sample_accession,experiment_accession,run_accession,library_strategy,library_source,library_selection,instrument_platform,fastq_ftp,fastq_md5,fastq_bytes"
        api = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode(
            {"accession": acc, "result": result, "fields": fields, "format": "json", "limit": 0})
        d, err = get_json(api)
        names = []
        for x in d or []:
            names.extend([p.rsplit("/", 1)[-1] for p in (x.get("fastq_ftp") or "").split(";") if p])
        return names, api, err
    return [], "", ""


def main() -> None:
    base = Path(__file__).resolve().parent / "results"
    ap = argparse.ArgumentParser()
    ap.add_argument("--links", default=str(base / "stage3_links/paper_data_links.csv"))
    ap.add_argument("--outdir", required=True, help="Directory for stage 4 outputs")
    ap.add_argument("--sleep", type=float, default=0.2)
    args = ap.parse_args()
    grouped = {}
    for r in csv.DictReader(open(args.links, newline="", encoding="utf-8")):
        repo, acc = r["repository"], r["accession"]
        key = (repo, acc or r["raw_url"])
        x = grouped.setdefault(key, {"repository": repo, "accession": acc, "paper_ids": set(),
                                    "pmids": set(), "dois": set(), "raw_url": r["raw_url"]})
        for src, dst in (("paper_id", "paper_ids"), ("pmid", "pmids"), ("doi", "dois")):
            if r.get(src): x[dst].add(r[src])
    rows = []
    for i, ((repo, key), x) in enumerate(grouped.items(), 1):
        acc = x["accession"]
        landing = TEMPLATES.get(repo, "").format(a=acc) if acc else x["raw_url"]
        files, api, err = list_record(repo, acc) if acc else ([], "", "")
        rows.append({"dataset_id": f"{repo}:{key}", "repository": repo, "accession": acc,
                     "paper_ids": ";".join(sorted(x["paper_ids"])), "pmids": ";".join(sorted(x["pmids"])),
                     "dois": ";".join(sorted(x["dois"])), "landing_url": landing,
                     "metadata_api_url": api, "n_listed_files": len(files),
                     "listed_files": ";".join(files), "resolve_status": "error" if err else "ok",
                     "resolve_error": err[:500], "data_downloaded": False})
        if api: time.sleep(args.sleep)
        if i % 100 == 0: print(f"resolved {i}/{len(grouped)}", flush=True)
    out = Path(args.outdir) / "datasets_master.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["dataset_id", "repository", "accession", "paper_ids", "pmids", "dois", "landing_url", "metadata_api_url", "n_listed_files", "listed_files", "resolve_status", "resolve_error", "data_downloaded"]
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)
    print(f"wrote {len(rows)} dataset records -> {out}")


if __name__ == "__main__":
    main()
