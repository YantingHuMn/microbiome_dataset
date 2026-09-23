#!/usr/bin/env python3
"""resolve accessions to landing/download locations; never persist deposited data.

Note on europepmc_supp: this stage fetches the PMC supplementary-files ZIP
into memory to read its file *names* (needed to score abundance likelihood
in stage5), then discards the bytes. It writes nothing to disk. Actually
downloading and opening file *contents* is a separate, explicit stage
(stage6_verify_abundance.py) so the "never persist deposited data" contract
of stages 0-5 still holds.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import time
import urllib.parse
import urllib.request
import zipfile
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
    "osf": "https://osf.io/{a}/", "europepmc_supp": "https://pmc.ncbi.nlm.nih.gov/articles/{a}/",
}
ENA_TYPES = {"bioproject", "sra_study", "sra_run", "sra_experiment"}
MAX_SUPP_ZIP_MB = 300  # abandon a listing if the supplementary bundle is absurdly large


def get_json(url: str) -> tuple[dict | list | None, str]:
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
            return json.loads(r.read()), ""
    except Exception as e: return None, str(e)


def standardize_bioproject(acc: str) -> tuple[str, str, str]:
    if not acc.isdigit(): return acc, "", ""
    api = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(
        {"db": "bioproject", "id": acc, "retmode": "json"})
    d, err = get_json(api)
    uids = (d or {}).get("result", {}).get("uids", [])
    standard = (d or {}).get("result", {}).get(uids[0], {}).get("project_acc", "") if uids else ""
    return (standard.upper() if standard else acc), api, err or ("BioProject accession not found" if not standard else "")


def list_zip_names(url: str, max_mb: float = MAX_SUPP_ZIP_MB) -> tuple[list[str], str]:
    """Fetch a zip fully into memory, list member names (recursing one level
    into nested zips -- MDPI/OUP ship zip-inside-zip), then drop the bytes.
    Nothing is written to disk."""
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=180) as r:
            cl = r.headers.get("Content-Length")
            if cl and int(cl) > max_mb * 1e6:
                return [], f"supplementary bundle too large ({int(cl)} bytes)"
            data = r.read(int(max_mb * 1e6) + 1)
        if len(data) > max_mb * 1e6:
            return [], "supplementary bundle exceeded size cap while streaming"
        names = []
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            for info in z.infolist():
                if info.is_dir():
                    continue
                if info.filename.lower().endswith(".zip"):
                    try:
                        with z.open(info) as nested_fh:
                            nested_bytes = nested_fh.read()
                        with zipfile.ZipFile(io.BytesIO(nested_bytes)) as z2:
                            names.extend(f"{info.filename}!{n.filename}"
                                        for n in z2.infolist() if not n.is_dir())
                    except zipfile.BadZipFile:
                        names.append(info.filename)
                else:
                    names.append(info.filename)
        return names, ""
    except Exception as e:
        return [], str(e)


def list_record(repo: str, acc: str) -> tuple[list[str], str, str, dict[str, str]]:
    """Return metadata file names, API endpoint, error. No deposited file content."""
    if repo == "zenodo":
        api = f"https://zenodo.org/api/records/{acc}"; d, err = get_json(api)
        return [x.get("key", "") for x in (d or {}).get("files", [])], api, err, {}
    if repo == "figshare":
        api = f"https://api.figshare.com/v2/articles/{acc}"; d, err = get_json(api)
        return [x.get("name", "") for x in (d or {}).get("files", [])], api, err, {}
    if repo == "europepmc_supp":
        api = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{acc}/supplementaryFiles"
        names, err = list_zip_names(api)
        return names, api, err, {}
    if repo == "dryad":
        doi_path = f"doi:10.5061/dryad.{acc}".replace("/", "%2F").replace(":", "%3A")
        api = f"https://datadryad.org/api/v2/datasets/{doi_path}"
        d, err = get_json(api)
        ver = ((d or {}).get("_links", {}).get("stash:version", {}) or {}).get("href", "")
        if not ver:
            return [], api, err or "no dataset version link", {}
        files_api = f"https://datadryad.org{ver}/files"
        jf, err2 = get_json(files_api)
        names = [it.get("path", "") for it in (jf or {}).get("_embedded", {}).get("stash:files", [])]
        return names, f"{api};{files_api}", err or err2, {}
    if repo == "osf":
        api = f"https://api.osf.io/v2/nodes/{acc}/files/osfstorage/"; d, err = get_json(api)
        names = [it.get("attributes", {}).get("name", "") for it in (d or {}).get("data", [])
                 if it.get("attributes", {}).get("kind") == "file"]
        return names, api, err, {}
    if repo == "mgnify_study":
        api = f"https://www.ebi.ac.uk/metagenomics/api/v1/studies/{acc}"; d, err = get_json(api)
        return [], api, err, {}
    if repo in ENA_TYPES:
        result = {"bioproject": "read_run", "sra_study": "read_run", "sra_run": "read_run", "sra_experiment": "read_run"}[repo]
        fields = "study_accession,sample_accession,experiment_accession,run_accession,library_strategy,library_source,library_selection,library_layout,instrument_platform,fastq_ftp,fastq_md5,fastq_bytes"
        api = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode(
            {"accession": acc, "result": result, "fields": fields, "format": "json", "limit": 0})
        d, err = get_json(api)
        names = []
        for x in d or []:
            names.extend([p.rsplit("/", 1)[-1] for p in (x.get("fastq_ftp") or "").split(";") if p])
        metadata = {field: ";".join(sorted({x.get(field, "") for x in d or [] if x.get(field, "")}))
                    for field in ("library_strategy", "library_source", "library_layout", "instrument_platform")}
        return names, api, err, metadata
    return [], "", "", {}


def main() -> None:
    base = Path(__file__).resolve().parent / "results"
    ap = argparse.ArgumentParser()
    ap.add_argument("--links", default=str(base / "stage3_links/paper_data_links.csv"))
    ap.add_argument("--outdir", required=True, help="Directory for stage 4 outputs")
    ap.add_argument("--sleep", type=float, default=0.2)
    args = ap.parse_args()
    grouped = {}
    bioproject_cache = {}
    for r in csv.DictReader(open(args.links, newline="", encoding="utf-8")):
        repo, acc = r["repository"], r["accession"]
        normalization_api = normalization_error = ""
        if repo == "bioproject" and acc.isdigit():
            if acc not in bioproject_cache: bioproject_cache[acc] = standardize_bioproject(acc)
            acc, normalization_api, normalization_error = bioproject_cache[acc]
        key = (repo, acc or r["raw_url"])
        x = grouped.setdefault(key, {"repository": repo, "accession": acc, "paper_ids": set(),
                                    "pmids": set(), "dois": set(), "raw_url": r["raw_url"],
                                    "normalization_api": normalization_api, "normalization_error": normalization_error})
        if normalization_api: x["normalization_api"] = normalization_api
        if normalization_error: x["normalization_error"] = normalization_error
        for src, dst in (("paper_id", "paper_ids"), ("pmid", "pmids"), ("doi", "dois")):
            if r.get(src): x[dst].add(r[src])
    rows = []
    for i, ((repo, key), x) in enumerate(grouped.items(), 1):
        acc = x["accession"]
        landing = TEMPLATES.get(repo, "").format(a=acc) if acc else x["raw_url"]
        files, api, err, metadata = list_record(repo, acc) if acc else ([], "", "", {})
        api = ";".join(v for v in (x["normalization_api"], api) if v)
        err = "; ".join(v for v in (x["normalization_error"], err) if v)
        rows.append({"dataset_id": f"{repo}:{acc or key}", "repository": repo, "accession": acc,
                     "paper_ids": ";".join(sorted(x["paper_ids"])), "pmids": ";".join(sorted(x["pmids"])),
                     "dois": ";".join(sorted(x["dois"])), "landing_url": landing,
                     "metadata_api_url": api, "n_listed_files": len(files),
                     "listed_files": ";".join(files), "library_strategy": metadata.get("library_strategy", ""),
                     "library_source": metadata.get("library_source", ""), "library_layout": metadata.get("library_layout", ""),
                     "instrument_platform": metadata.get("instrument_platform", ""), "resolve_status": "error" if err else "ok",
                     "resolve_error": err[:500], "data_downloaded": False})
        if api: time.sleep(args.sleep)
        if i % 100 == 0: print(f"resolved {i}/{len(grouped)}", flush=True)

    # MGnify shortcut: a BioProject that only has raw reads may already have
    # a processed taxonomy-abundance table on MGnify, avoiding a pipeline
    # re-run. This never overwrites the original bioproject row -- it adds a
    # sibling dataset record so stage5 can credit either source.
    seen_bioprojects = {r["accession"] for r in rows if r["repository"] == "bioproject" and r["accession"]}
    for i, bp in enumerate(sorted(seen_bioprojects), 1):
        api = f"https://www.ebi.ac.uk/metagenomics/api/v1/studies?bioproject={bp}"
        d, err = get_json(api)
        for it in (d or {}).get("data", []):
            mgys = it.get("id", "")
            if not mgys:
                continue
            src_row = next(r for r in rows if r["repository"] == "bioproject" and r["accession"] == bp)
            files, mapi, merr, _ = list_record("mgnify_study", mgys)
            rows.append({"dataset_id": f"mgnify_from_bioproject:{mgys}", "repository": "mgnify_study",
                        "accession": mgys, "paper_ids": src_row["paper_ids"], "pmids": src_row["pmids"],
                        "dois": src_row["dois"], "landing_url": TEMPLATES["mgnify_study"].format(a=mgys),
                        "metadata_api_url": f"{api};{mapi}", "n_listed_files": len(files),
                        "listed_files": ";".join(files), "library_strategy": "", "library_source": "",
                        "library_layout": "", "instrument_platform": "",
                        "resolve_status": "error" if (err or merr) else "ok",
                        "resolve_error": "; ".join(v for v in (err, merr) if v)[:500],
                        "data_downloaded": False})
        time.sleep(args.sleep)
        if i % 50 == 0: print(f"mgnify-checked {i}/{len(seen_bioprojects)} bioprojects", flush=True)

    out = Path(args.outdir) / "datasets_master.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["dataset_id", "repository", "accession", "paper_ids", "pmids", "dois", "landing_url", "metadata_api_url", "n_listed_files", "listed_files", "library_strategy", "library_source", "library_layout", "instrument_platform", "resolve_status", "resolve_error", "data_downloaded"]
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)
    print(f"wrote {len(rows)} dataset records -> {out}")


if __name__ == "__main__":
    main()
