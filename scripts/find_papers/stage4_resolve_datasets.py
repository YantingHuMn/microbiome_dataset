#!/usr/bin/env python3
"""resolve accessions to landing/download locations; never persist deposited data.

Note on europepmc_supp: this stage fetches the PMC supplementary-files ZIP
into memory to read its file *names* (needed to score abundance likelihood
in stage5), then discards the bytes. It writes nothing to disk. Actually
downloading and opening file *contents* is a separate, explicit stage
(stage6_verify_abundance.py) so the "never persist deposited data" contract
of stages 0-5 still holds.

--workers > 1 resolves multiple accessions concurrently; see _netutil.py for
the shared per-host Throttle that keeps this safe -- more workers overlaps
latency across DIFFERENT hosts, it does not send more requests/sec to any
ONE host than a single worker would.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import threading
import time
import urllib.parse
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _netutil import GLOBAL_THROTTLE  # noqa: E402

csv.field_size_limit(sys.maxsize)  # listed_files can exceed the 131072-byte default

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
        GLOBAL_THROTTLE.wait(url)
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=20) as r:
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
        GLOBAL_THROTTLE.wait(url)
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


PMC_BIN_RE = re.compile(r'/articles/instance/(\d+)/bin/([^"\']+)')


def list_pmc_associated_data(pmcid: str) -> list[str]:
    """Filenames scraped from the live PMC article page's Associated Data
    section. Europe PMC's supplementaryFiles zip mirrors ONE NIHMS
    manuscript submission and can miss files from a LATER revision that
    the live PMC page already shows -- observed directly: for PMC6342642,
    the zip contained NIHMS80310's 14 figure/reporting-summary files,
    while the live page's Associated Data listed NIHMS1510763's 8 xlsx
    datasets + 2 PDFs, a completely different manuscript submission. Used
    here to widen listed_files so stage5 doesn't miss these; stage6 has
    the matching fetch-side fix to actually download them.
    """
    try:
        GLOBAL_THROTTLE.wait(f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/")
        req = urllib.request.Request(f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/", headers=UA)
        with urllib.request.urlopen(req, timeout=30) as r:
            html_text = r.read().decode("utf-8", "replace")
    except Exception:
        return []
    seen, names = set(), []
    for _num, fname in PMC_BIN_RE.findall(html_text):
        if fname not in seen:
            seen.add(fname); names.append(fname)
    return names


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
        extra = list_pmc_associated_data(acc)
        merged = list(dict.fromkeys(names + extra))  # union, de-duped, order-preserving
        page_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{acc}/"
        return merged, (f"{api};{page_url}" if extra else api), err, {}
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


def resolve_one(repo: str, key, x: dict) -> dict:
    acc = x["accession"]
    landing = TEMPLATES.get(repo, "").format(a=acc) if acc else x["raw_url"]
    files, api, err, metadata = list_record(repo, acc) if acc else ([], "", "", {})
    api = ";".join(v for v in (x["normalization_api"], api) if v)
    err = "; ".join(v for v in (x["normalization_error"], err) if v)
    return {"dataset_id": f"{repo}:{acc or key}", "repository": repo, "accession": acc,
            "paper_ids": ";".join(sorted(x["paper_ids"])), "pmids": ";".join(sorted(x["pmids"])),
            "dois": ";".join(sorted(x["dois"])), "landing_url": landing,
            "metadata_api_url": api, "n_listed_files": len(files),
            "listed_files": ";".join(files), "library_strategy": metadata.get("library_strategy", ""),
            "library_source": metadata.get("library_source", ""), "library_layout": metadata.get("library_layout", ""),
            "instrument_platform": metadata.get("instrument_platform", ""), "resolve_status": "error" if err else "ok",
            "resolve_error": err[:500], "data_downloaded": False}


MGNIFY_MAX_MATCHES = 10  # a popular BioProject can link 25+ MGnify studies; listing
                         # files for each one is a separate serial network call inside
                         # this single worker slot -- one such BioProject observed taking
                         # 20s (25 matches) against ~2.5s for a BioProject with none. This
                         # lookup is a best-effort shortcut, not required for correctness,
                         # so cap it rather than let a handful of popular BioProjects blow
                         # up the whole pass's tail latency.


def mgnify_lookup_one(bp: str, src_row: dict) -> list[dict]:
    api = f"https://www.ebi.ac.uk/metagenomics/api/v1/studies?bioproject={bp}"
    d, err = get_json(api)
    out = []
    matches = (d or {}).get("data", [])
    for it in matches[:MGNIFY_MAX_MATCHES]:
        mgys = it.get("id", "")
        if not mgys:
            continue
        files, mapi, merr, _ = list_record("mgnify_study", mgys)
        out.append({"dataset_id": f"mgnify_from_bioproject:{mgys}", "repository": "mgnify_study",
                    "accession": mgys, "paper_ids": src_row["paper_ids"], "pmids": src_row["pmids"],
                    "dois": src_row["dois"], "landing_url": TEMPLATES["mgnify_study"].format(a=mgys),
                    "metadata_api_url": f"{api};{mapi}", "n_listed_files": len(files),
                    "listed_files": ";".join(files), "library_strategy": "", "library_source": "",
                    "library_layout": "", "instrument_platform": "",
                    "resolve_status": "error" if (err or merr) else "ok",
                    "resolve_error": "; ".join(v for v in (err, merr) if v)[:500],
                    "data_downloaded": False})
    return out


FIELDS = ["dataset_id", "repository", "accession", "paper_ids", "pmids", "dois", "landing_url",
         "metadata_api_url", "n_listed_files", "listed_files", "library_strategy", "library_source",
         "library_layout", "instrument_platform", "resolve_status", "resolve_error", "data_downloaded"]


def main() -> None:
    base = Path(__file__).resolve().parent / "results"
    ap = argparse.ArgumentParser()
    ap.add_argument("--links", default=str(base / "stage3_links/paper_data_links.csv"))
    ap.add_argument("--outdir", required=True, help="Directory for stage 4 outputs")
    ap.add_argument("--sleep", type=float, default=0.2,
                    help="ignored when --workers > 1 -- the shared Throttle paces requests instead")
    ap.add_argument("--workers", type=int, default=1,
                    help="parallel accessions resolved at once. Safe to raise (8-16 is reasonable) "
                         "-- see _netutil.Throttle for why this doesn't hit any one API harder.")
    ap.add_argument("--skip-mgnify", action="store_true",
                    help="skip the MGnify-from-BioProject reverse lookup entirely. This pass is a "
                         "best-effort shortcut (it can convert some raw_reads-only BioProjects into "
                         "abundance_ready for free), NOT required for correctness -- skipped "
                         "BioProjects simply stay classified as raw_reads/needs-own-pipeline, which "
                         "is already a valid outcome. Even at ~5.5s/item post-fix, ~24k BioProjects "
                         "is still ~35-40h; use this to get the main pipeline result first and run "
                         "the MGnify pass separately later, on its own time budget.")
    args = ap.parse_args()

    out = Path(args.outdir) / "datasets_master.csv"
    out.parent.mkdir(parents=True, exist_ok=True)

    # Resumable like stage6: every resolved dataset_id is appended to disk
    # immediately, not accumulated in memory until the end. A job that hits
    # its walltime limit partway through can just be resubmitted with the
    # same --outdir -- already-resolved dataset_ids are skipped, not redone.
    done_ids = set()
    if out.exists():
        with out.open(newline="", encoding="utf-8") as fh:
            done_ids = {r["dataset_id"] for r in csv.DictReader(fh)}
        if done_ids:
            print(f"resuming: {len(done_ids)} dataset records already resolved in {out}", flush=True)
    fh = out.open("a", newline="", encoding="utf-8")
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if out.stat().st_size == 0:
        w.writeheader()

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
        for s, dst in (("paper_id", "paper_ids"), ("pmid", "pmids"), ("doi", "dois")):
            if r.get(s): x[dst].add(r[s])
    pending = {key: x for key, x in grouped.items()
              if f"{key[0]}:{x['accession'] or key[1]}" not in done_ids}
    print(f"{len(grouped)} grouped datasets, {len(pending)} pending after resume", flush=True)

    bp_rows: dict[str, dict] = {}  # accession -> row, for the mgnify pass below
    def emit(row: dict) -> None:
        w.writerow({k: row.get(k, "") for k in FIELDS}); fh.flush()
        if row["repository"] == "bioproject" and row["accession"]:
            bp_rows[row["accession"]] = row

    if args.workers <= 1:
        for i, ((repo, key), x) in enumerate(pending.items(), 1):
            emit(resolve_one(repo, key, x))
            if x.get("normalization_api"): time.sleep(args.sleep)
            if i % 100 == 0: print(f"resolved {i}/{len(pending)}", flush=True)
    else:
        with ThreadPoolExecutor(args.workers) as ex:
            futs = {ex.submit(resolve_one, repo, key, x): 1 for (repo, key), x in pending.items()}
            for i, fu in enumerate(as_completed(futs), 1):
                emit(fu.result())
                if i % 100 == 0: print(f"resolved {i}/{len(pending)}", flush=True)

    # For the MGnify pass, bp_rows only has bioprojects resolved THIS run --
    # on a resumed run we also need bioprojects that were already resolved
    # in a prior (interrupted) run, so re-read them from disk.
    if done_ids:
        with out.open(newline="", encoding="utf-8") as f2:
            for r in csv.DictReader(f2):
                if r["repository"] == "bioproject" and r["accession"]:
                    bp_rows.setdefault(r["accession"], r)

    # MGnify shortcut: a BioProject that only has raw reads may already have
    # a processed taxonomy-abundance table on MGnify, avoiding a pipeline
    # re-run. This never overwrites the original bioproject row -- it adds a
    # sibling dataset record so stage5 can credit either source. Optional
    # (--skip-mgnify): skipped BioProjects simply stay classified as
    # raw_reads/needs-own-pipeline downstream, which is already valid.
    if args.skip_mgnify:
        print("skipping MGnify-from-BioProject lookup (--skip-mgnify)", flush=True)
        fh.close()
        with out.open(newline="", encoding="utf-8") as f2:
            n_total = sum(1 for _ in csv.DictReader(f2))
        print(f"done. {n_total} dataset records total -> {out}")
        return

    # Resuming this pass by dataset_id alone does NOT work: a checked
    # bioproject with ZERO MGnify matches writes no output row at all, so
    # there is nothing in datasets_master.csv to prove it was ever checked
    # -- a naive resume would re-query every bioproject from scratch even
    # if 24,000 of them were already checked and simply came back empty.
    # Track "checked" explicitly in its own ledger, independent of whether
    # a match was found.
    mgnify_checked_path = Path(args.outdir) / "mgnify_checked_bioprojects.txt"
    already_checked_bp = set()
    if mgnify_checked_path.exists():
        already_checked_bp = {l.strip() for l in mgnify_checked_path.open(encoding="utf-8") if l.strip()}
        if already_checked_bp:
            print(f"resuming mgnify pass: {len(already_checked_bp)} bioprojects already checked", flush=True)
    mgnify_fh = mgnify_checked_path.open("a", encoding="utf-8")
    mgnify_lock = threading.Lock()

    def mark_checked(bp: str) -> None:
        with mgnify_lock:
            mgnify_fh.write(bp + "\n"); mgnify_fh.flush()

    seen_bioprojects = sorted(bp_rows)
    pending_bp = [bp for bp in seen_bioprojects if bp not in already_checked_bp]
    print(f"{len(seen_bioprojects)} bioprojects total, {len(pending_bp)} pending after resume", flush=True)
    if args.workers <= 1:
        for i, bp in enumerate(pending_bp, 1):
            for row in mgnify_lookup_one(bp, bp_rows[bp]):
                if row["dataset_id"] not in done_ids:
                    emit(row)
            mark_checked(bp)
            time.sleep(args.sleep)
            if i % 50 == 0: print(f"mgnify-checked {i}/{len(pending_bp)} bioprojects", flush=True)
    else:
        with ThreadPoolExecutor(args.workers) as ex:
            futs = {ex.submit(mgnify_lookup_one, bp, bp_rows[bp]): bp for bp in pending_bp}
            for i, fu in enumerate(as_completed(futs), 1):
                bp = futs[fu]
                for row in fu.result():
                    if row["dataset_id"] not in done_ids:
                        emit(row)
                mark_checked(bp)
                if i % 50 == 0: print(f"mgnify-checked {i}/{len(pending_bp)} bioprojects", flush=True)
    mgnify_fh.close()

    fh.close()
    with out.open(newline="", encoding="utf-8") as f2:
        n_total = sum(1 for _ in csv.DictReader(f2))
    print(f"done. {n_total} dataset records total -> {out}")


if __name__ == "__main__":
    main()
