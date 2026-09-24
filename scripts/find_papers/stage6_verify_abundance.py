#!/usr/bin/env python3
"""Open the actual files stage4/5 only listed, and judge abundance-matrix
content by STRUCTURE, not filename. This is the one stage in the pipeline
that downloads deposited data -- stages 0-5 stay listing-only by design.

Input: dataset_candidates_final.csv (stage5) filtered to candidate_status in
{abundance_ready, needs_content_check} -- i.e. exactly the datasets stage5
could not, or should not, resolve from filenames alone.

Memory/disk discipline (this is the point of splitting stage6 out):
  - one dataset is downloaded, sniffed, normalised, and its raw bytes
    deleted before that dataset's task ends -- peak disk usage per unit of
    work is bounded by ONE study's supplementary/deposit size, not the sum
    of all of them.
  - long-format rows are appended to disk immediately after each dataset,
    never accumulated in memory across the whole run.
  - a --max-file-mb and --max-study-mb budget skips anything oversized
    rather than filling the disk.
  - gc.collect() after every dataset to drop pandas/zipfile intermediates.

Speed (--workers): this stage is network-latency bound, not CPU bound --
most of the wall-clock time is spent waiting on a download, not parsing
it. --workers N runs N datasets concurrently in a thread pool; every
worker still follows the exact discipline above for ITS OWN dataset, so
raising --workers scales peak disk/memory roughly as N x --max-study-mb,
not unboundedly -- pick N so that stays under your node's actual free
scratch space and RAM (e.g. --max-study-mb 500 --workers 8 wants a scratch
volume that can hold ~4 GB comfortably, not the whole corpus). Concurrent
requests to any ONE host (Zenodo, figshare, EBI, ...) are still paced to
the same per-host rate as a single worker -- see _netutil.Throttle --
so --workers shortens wall-clock time by overlapping DIFFERENT hosts'
latency, it does not hit any one API harder. Writes to the shared
abundance_long.tsv.gz are serialised with a lock; each worker's own
downloaded files never collide since they live under
scratch/<dataset_id>/.

Usage
-----
    python stage6_verify_abundance.py \\
        --datasets  Database/results/stage5_final/dataset_candidates_final.csv \\
        --papers    Database/results/stage5_final/paper_candidates_final.csv \\
        --outdir    Database/results/stage6_verified \\
        --scratch   /scratch/$USER/mb_stage6 \\
        --workers   8   # fast/local scratch disk with room for ~8 x --max-study-mb/mb_stage6   # fast local disk, NOT $HOME

Outputs (outdir)
-----------------
    dataset_verification.csv   one row per dataset actually opened: verdict,
                                axis, value_type, n_taxa, n_samples, sha256,
                                note -- the audit trail for every abundance_ready
                                claim.
    paper_verification.csv     full audit trail: paper_candidates_final.csv +
                                content_verified, data_category, source_pipeline,
                                next_action -- every column, every paper.
    abundance_final.csv        DELIVERABLE 1: papers with >=1 content-verified
                                abundance matrix. Nothing left to do.
    run_own_pipeline.csv        DELIVERABLE 2: papers for which only raw reads
                                are available (next_action=run_own_pipeline).
    manual_content_review.csv   DELIVERABLE 3: ambiguous cases that still need
                                a human to open and inspect the deposited file
                                (next_action=manual_content_review).
    abundance_long.tsv.gz      full long-format table (paper_id, dataset_id,
                                sample_id, taxon, kingdom..species, value,
                                value_type, domain) -- append-only during the run.
    matrices/<paper_id>__<dataset_id>.tsv   one taxa x sample matrix per
                                verified dataset (written at the end from the
                                long table, not held in memory meanwhile).

Re-running skips datasets already present in dataset_verification.csv, so a
killed cluster job resumes instead of restarting.
"""
from __future__ import annotations

import argparse
import csv
import gc
import gzip
import io
import json
import re
import shutil
import sys
import threading
import time
import urllib.request
import zipfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _netutil import GLOBAL_THROTTLE  # noqa: E402

csv.field_size_limit(sys.maxsize)
UA = {"User-Agent": "microbiome-dataset-verifier/1.0 (academic research)"}

TABLE_EXT = {".xlsx", ".xls", ".csv", ".tsv", ".txt", ".tab", ".biom", ".qza"}

LIST_APIS = {
    "zenodo": lambda acc: f"https://zenodo.org/api/records/{acc}",
    "figshare": lambda acc: f"https://api.figshare.com/v2/articles/{acc}",
    "europepmc_supp": lambda acc: f"https://www.ebi.ac.uk/europepmc/webservices/rest/{acc}/supplementaryFiles",
}

TAXA_RE = re.compile(
    r"(^|[;|,\s])[kpcofgsd]__"
    r"|;\s*[A-Z][a-z]{2,}"
    r"|(aceae|ales|bacteria|bacteriota|archaea|mycota|virales|viridae|"
    r"coccus|bacillus|monas|bacter|vibrio|spirillum|clostridium|prevotella|"
    r"bacteroides|lactobacillus|firmicutes|proteobacteria|actinobacteri|"
    r"verrucomicrobi|bifidobacterium|akkermansia|faecalibacterium)"
    r"|^(otu|asv|zotu|sv)[_ -]?\d+", re.I)
RANK_WORDS = {"kingdom", "domain", "superkingdom", "phylum", "class", "order",
              "family", "genus", "species", "taxon", "taxonomy", "otu", "asv",
              "otu_id", "asv_id", "feature id", "featureid", "#otu id",
              "lineage", "clade", "taxa"}
ABUND_RE = re.compile(r"(abundance|relative[_ ]?abund|count|reads|freq|"
                      r"rpkm|tpm|cpm|value|proportion)", re.I)
SAMPLE_RE = re.compile(r"(sample|subject|specimen|run|library|host|individual)", re.I)
RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
RANK_PREFIX = {"d": "kingdom", "k": "kingdom", "p": "phylum", "c": "class",
               "o": "order", "f": "family", "g": "genus", "s": "species"}
MIN_TAXA_FRACTION = 0.50
REVIEW_TAXA_FRACTION = 0.30

# Every accepted dataset is either a matrix someone already produced
# (processed) or reads still needing OUR pipeline (raw_reads is never
# routed here -- stage5 keeps those out of --statuses). For processed
# datasets we record WHICH pipeline made the matrix, read off the actual
# filename/column evidence we just opened -- more reliable than the
# abstract-text guess stage5 makes, since this is the file itself.
PIPELINE_PATTERNS = [
    ("QIIME2", re.compile(r"\bqiime\s?2\b|\.qza(?:!|\.|$)|\.qzv(?:!|\.|$)|feature-table", re.I)),
    ("QIIME1", re.compile(r"\bqiime\b(?!\s?2)", re.I)),
    ("mothur", re.compile(r"\bmothur\b", re.I)),
    ("DADA2", re.compile(r"\bdada2\b|seqtab", re.I)),
    ("USEARCH/UPARSE", re.compile(r"\busearch\b|\buparse\b|\bzotu", re.I)),
    ("VSEARCH", re.compile(r"\bvsearch\b", re.I)),
    ("MetaPhlAn", re.compile(r"\bmetaphlan\d?\b", re.I)),
    ("Kraken/Bracken", re.compile(r"\bkraken\d?\b|\bbracken\b", re.I)),
    ("Kaiju", re.compile(r"\bkaiju\b", re.I)),
    ("Centrifuge", re.compile(r"\bcentrifuge\b", re.I)),
    ("HUMAnN", re.compile(r"\bhumann\d?\b", re.I)),
    ("PICRUSt2", re.compile(r"\bpicrust\d?\b", re.I)),
    ("LotuS2", re.compile(r"\blotus\d?\b", re.I)),
    ("MEGAN", re.compile(r"\bmegan\b", re.I)),
    ("phyloseq", re.compile(r"\bphyloseq\b", re.I)),
    (".biom format", re.compile(r"\.biom(?:!|\.|$)", re.I)),
]


def detect_pipeline(*texts: str) -> str:
    joined = " ".join(t for t in texts if t)
    if not joined:
        return ""
    return ";".join(name for name, rx in PIPELINE_PATTERNS if rx.search(joined))


# --------------------------------------------------------------------------- #
# download -- one file at a time, streamed, size-capped, never re-read twice
# --------------------------------------------------------------------------- #

def stream_download(url: str, dest: Path, max_mb: float) -> tuple[bool, str]:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        GLOBAL_THROTTLE.wait(url)
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=300) as r:
            cl = r.headers.get("Content-Length")
            if cl and int(cl) > max_mb * 1e6:
                return False, f"too_large:{cl}"
            written = 0
            with open(tmp, "wb") as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    written += len(chunk)
                    if written > max_mb * 1e6:
                        fh.close()
                        tmp.unlink(missing_ok=True)
                        return False, "too_large:stream"
                    fh.write(chunk)
        tmp.replace(dest)
        return True, "ok"
    except Exception as e:                          # noqa: BLE001
        tmp.unlink(missing_ok=True)
        return False, f"ERR:{type(e).__name__}:{e}"[:200]


PMC_BIN_RE = re.compile(r'/articles/instance/(\d+)/bin/([^"\']+)')


def list_pmc_associated_data(pmcid: str) -> list[tuple[str, str]]:
    """(filename, download_url) pairs scraped from the live PMC article
    page's Associated Data section. Europe PMC's supplementaryFiles zip
    mirrors ONE NIHMS manuscript submission and can miss files from a
    LATER revision the live PMC page already shows -- observed directly:
    for PMC6342642, the zip contained NIHMS80310's 14 figure/reporting-
    summary files, while the live page's Associated Data listed
    NIHMS1510763's 8 xlsx datasets + 2 PDFs, a completely different
    manuscript submission. Without this, those 8 xlsx files -- exactly
    the kind of file this whole stage exists to open -- are never even
    downloaded, let alone content-checked.
    """
    import urllib.request as ur
    url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"
    try:
        GLOBAL_THROTTLE.wait(url)
        with ur.urlopen(ur.Request(url, headers=UA), timeout=30) as r:
            html_text = r.read().decode("utf-8", "replace")
    except Exception:                                # noqa: BLE001
        return []
    seen, out = set(), []
    for num, fname in PMC_BIN_RE.findall(html_text):
        dl_url = f"https://pmc.ncbi.nlm.nih.gov/articles/instance/{num}/bin/{fname}"
        if fname not in seen:
            seen.add(fname); out.append((fname, dl_url))
    return out


def list_deposit_files(repo: str, acc: str) -> list[tuple[str, str, int]]:
    """(name, download_url, size) for one deposit. europepmc_supp returns a
    single zip URL -- its members are discovered after download -- PLUS
    any files found only on the live PMC page (see list_pmc_associated_data)."""
    import urllib.request as ur
    if repo == "europepmc_supp":
        files = [(f"{acc}_supplementary.zip", LIST_APIS["europepmc_supp"](acc), 0)]
        seen_names = set()
        for fname, dl_url in list_pmc_associated_data(acc):
            if fname not in seen_names:
                seen_names.add(fname); files.append((fname, dl_url, 0))
        return files
    if repo not in LIST_APIS:
        return []
    api_url = LIST_APIS[repo](acc)
    try:
        GLOBAL_THROTTLE.wait(api_url)
        with ur.urlopen(ur.Request(api_url, headers=UA), timeout=90) as r:
            d = json.loads(r.read())
    except Exception:                                # noqa: BLE001
        return []
    if repo == "zenodo":
        return [(f.get("key", ""), (f.get("links") or {}).get("self", ""), int(f.get("size") or 0))
                for f in d.get("files", [])]
    if repo == "figshare":
        return [(f.get("name", ""), f.get("download_url", ""), int(f.get("size") or 0))
                for f in d.get("files", [])]
    return []


# --------------------------------------------------------------------------- #
# sniff -- structural verdict, content only, filename is never trusted alone
# --------------------------------------------------------------------------- #

def taxa_fraction(values) -> float:
    vals = [str(v).strip() for v in values
            if v is not None and str(v).strip() and str(v).lower() != "nan"]
    if not vals:
        return 0.0
    return sum(bool(TAXA_RE.search(v)) for v in vals) / len(vals)


def iter_tables(path: Path, depth: int = 0):
    """Yield (sub_id, header_row, DataFrame). Recurses into zip/gz, including
    nested zips (MDPI/OUP) and QIIME2 .qza artifacts."""
    import pandas as pd  # local import: keeps stage0-5 dependency-free
    if depth > 3:
        return
    suf = path.suffix.lower()
    try:
        if suf == ".zip":
            with zipfile.ZipFile(path) as z:
                for info in z.infolist():
                    if info.is_dir():
                        continue
                    inner = Path(info.filename)
                    if inner.suffix.lower() not in TABLE_EXT | {".zip", ".gz"}:
                        continue
                    tmp = path.parent / f"__x{depth}_{re.sub(r'[^A-Za-z0-9._-]', '_', inner.name)}"
                    with z.open(info) as src, open(tmp, "wb") as out:
                        shutil.copyfileobj(src, out)
                    try:
                        for sid, hr, d in iter_tables(tmp, depth + 1):
                            yield f"{inner.name}!{sid}", hr, d
                    finally:
                        tmp.unlink(missing_ok=True)
        elif suf == ".qza":
            with zipfile.ZipFile(path) as z:
                for info in z.infolist():
                    if info.filename.endswith(("feature-table.biom", ".tsv", ".csv")):
                        tmp = path.parent / f"__q_{Path(info.filename).name}"
                        with z.open(info) as src, open(tmp, "wb") as out:
                            shutil.copyfileobj(src, out)
                        try:
                            for sid, hr, d in iter_tables(tmp, depth + 1):
                                yield f"qza!{sid}", hr, d
                        finally:
                            tmp.unlink(missing_ok=True)
        elif suf == ".gz":
            inner = path.with_suffix("")
            with gzip.open(path, "rb") as fi, open(inner, "wb") as fo:
                shutil.copyfileobj(fi, fo)
            try:
                yield from iter_tables(inner, depth + 1)
            finally:
                inner.unlink(missing_ok=True)
        elif suf in (".xlsx", ".xls"):
            xl = pd.ExcelFile(path)
            for sh in xl.sheet_names:
                for hr in range(0, 3):
                    try:
                        d = xl.parse(sh, header=hr, nrows=2000)
                    except Exception:                # noqa: BLE001
                        break
                    if d.shape[1] >= 2 and d.notna().sum().sum() > 0:
                        yield sh, hr, d
                    if hr == 0 and d.shape[1] >= 3:
                        break
        else:
            for sep in (",", "\t", ";"):
                try:
                    d = pd.read_csv(path, sep=sep, nrows=2000, engine="python",
                                    on_bad_lines="skip")
                except Exception:                    # noqa: BLE001
                    continue
                if d.shape[1] > 1:
                    yield path.name, 0, d
                    break
    except Exception as e:                           # noqa: BLE001
        print(f"    iter_tables error on {path.name}: {e}", file=sys.stderr)


def sniff_frame(d) -> dict:
    import numpy as np
    import pandas as pd
    v = {"is_abundance": False, "kind": "none", "axis": "", "label_col": None,
        "value_type": "unknown", "score": 0.0, "n_rows": int(d.shape[0]),
        "n_cols": int(d.shape[1]), "needs_review": False, "reason": ""}
    if d.shape[0] < 5 or d.shape[1] < 2:
        v["reason"] = "too small"
        return v
    d = d.dropna(axis=1, how="all").dropna(axis=0, how="all")
    cols = [str(c).strip() for c in d.columns]
    low = [c.lower() for c in cols]
    numeric = d.apply(pd.to_numeric, errors="coerce")
    num_frac = numeric.notna().mean()
    num_cols = [c for c, f in zip(d.columns, num_frac) if f > 0.80]
    lab_cols = [c for c in d.columns if c not in num_cols]

    has_tax = any(lc in RANK_WORDS for lc in low)
    has_abund = any(ABUND_RE.search(c) for c in cols)
    has_samp = any(SAMPLE_RE.search(c) for c in cols)
    if has_tax and has_abund and has_samp and d.shape[0] >= 10:
        v.update(is_abundance=True, kind="long", score=0.9, reason="sample/taxon/abundance columns")
        return v

    if len(num_cols) < 3:
        v["reason"] = f"only {len(num_cols)} numeric columns"
        return v
    row_best, row_col = 0.0, None
    for c in lab_cols[:3]:
        f = taxa_fraction(d[c].dropna().astype(str).head(300))
        if f > row_best:
            row_best, row_col = f, c
    col_best = taxa_fraction(cols)
    if row_best >= col_best:
        axis, frac, label = "rows", row_best, row_col
    else:
        axis, frac, label = "cols", col_best, None
    v.update(axis=axis, score=float(frac), label_col=(str(label) if label is not None else None))

    arr = numeric[num_cols].to_numpy(dtype=float, na_value=np.nan)
    finite = arr[np.isfinite(arr)]
    if finite.size:
        col_sums = np.nansum(arr, axis=0)
        col_sums = col_sums[np.isfinite(col_sums) & (col_sums > 0)]
        if col_sums.size and np.median(np.abs(col_sums - 1.0)) < 0.02:
            v["value_type"] = "relative_fraction"
        elif col_sums.size and np.median(np.abs(col_sums - 100.0)) < 1.0:
            v["value_type"] = "relative_percent"
        elif np.allclose(finite, np.round(finite)) and finite.max() > 10:
            v["value_type"] = "count"

    if frac >= MIN_TAXA_FRACTION:
        v.update(is_abundance=True, kind="wide", reason=f"{frac:.0%} of {axis} look taxonomic")
    elif frac >= REVIEW_TAXA_FRACTION:
        v.update(is_abundance=True, kind="wide", needs_review=True,
                 reason=f"weak taxonomic signal ({frac:.0%})")
    else:
        v["reason"] = f"no taxonomic signal ({frac:.0%})"
    return v


def split_taxonomy(label: str) -> dict[str, str]:
    out = {r: "" for r in RANKS}
    parts = re.split(r"[;|]", str(label or ""))
    unprefixed = []
    for part in parts:
        part = part.strip()
        m = re.match(r"^([dkpcofgs])__\s*(.*)$", part, re.I)
        if m:
            rank = RANK_PREFIX.get(m.group(1).lower())
            if rank in out:
                out[rank] = m.group(2).strip()
        elif part:
            unprefixed.append(part)
    if unprefixed and not any(out.values()):
        for rank, val in zip(RANKS, unprefixed):
            out[rank] = val
    return out


def to_long_rows(d, verdict: dict, paper_id: str, dataset_id: str, source_file: str):
    import pandas as pd
    import numpy as np
    if verdict["kind"] == "long":
        colmap = {str(c).lower(): c for c in d.columns}
        samp = next((colmap[c] for c in colmap if SAMPLE_RE.search(c)), None)
        tax = next((colmap[c] for c in colmap if c in RANK_WORDS), None)
        val = next((colmap[c] for c in colmap if ABUND_RE.search(c)), None)
        if not all([samp, tax, val]):
            return pd.DataFrame()
        out = d[[samp, tax, val]].copy()
        out.columns = ["sample_id", "taxon", "value"]
    else:
        if verdict["axis"] == "rows":
            label = verdict.get("label_col") or str(d.columns[0])
            if label not in d.columns:
                label = str(d.columns[0])
            num = d.drop(columns=[label]).apply(pd.to_numeric, errors="coerce")
            num = num.loc[:, num.notna().mean() > 0.5]
            out = num.copy()
            out.insert(0, "taxon", d[label].astype(str).values)
            out = out.melt(id_vars="taxon", var_name="sample_id", value_name="value")
        else:
            idx = d.columns[0]
            num = d.drop(columns=[idx]).apply(pd.to_numeric, errors="coerce")
            num = num.loc[:, num.notna().mean() > 0.5]
            out = num.copy()
            out.insert(0, "sample_id", d[idx].astype(str).values)
            out = out.melt(id_vars="sample_id", var_name="taxon", value_name="value")
    out = out.dropna(subset=["value"])
    out["value"] = pd.to_numeric(out["value"], errors="coerce")
    out = out.dropna(subset=["value"])
    if out.empty:
        return out
    out["sample_id"] = out["sample_id"].astype(str).str.strip()
    out["taxon"] = out["taxon"].astype(str).str.strip()
    ranks = out["taxon"].map(split_taxonomy).apply(pd.Series)
    out = pd.concat([out.reset_index(drop=True), ranks.reset_index(drop=True)], axis=1)
    out["paper_id"] = paper_id
    out["dataset_id"] = dataset_id
    out["value_type"] = verdict.get("value_type", "unknown")
    out["source_file"] = source_file
    out["domain"] = np.where(
        out["taxon"].str.contains(r"(?:virus|viridae|virales|phage)", case=False,
                                  regex=True, na=False), "virus", "prokaryote")
    return out


# --------------------------------------------------------------------------- #
# per-dataset pipeline: download -> sniff -> normalise -> flush -> delete
# --------------------------------------------------------------------------- #

def verify_dataset(row: dict, scratch: Path, max_file_mb: float, max_study_mb: float,
                   long_path: Path, matrices_dir: Path,
                   long_lock: threading.Lock | None = None) -> dict:
    repo, acc, dsid = row["repository"], row["accession"], row["dataset_id"]
    paper_ids = row.get("paper_ids", "")
    study_dir = scratch / re.sub(r"[^A-Za-z0-9._-]", "_", dsid)[:100]
    verdicts = []
    budget_used = 0.0
    n_rows_written = 0
    try:
        files = list_deposit_files(repo, acc)
        if not files:
            return {"dataset_id": dsid, "content_verified": False,
                    "n_tables_checked": 0, "n_tables_accepted": 0,
                    "note": "no files listed at verification time"}
        for name, url, size in files:
            if budget_used > max_study_mb * 1e6:
                verdicts.append({"file": name, "verdict": "skipped:study_budget"})
                continue
            if not url:
                continue
            local = study_dir / re.sub(r"[^A-Za-z0-9._-]", "_", name)[:100]
            ok, status = stream_download(url, local, max_file_mb)
            if not ok:
                verdicts.append({"file": name, "verdict": f"download_failed:{status}"})
                continue
            budget_used += local.stat().st_size
            n_checked = n_accepted = 0
            pipelines_found = set()
            for sub_id, hr, d in iter_tables(local):
                n_checked += 1
                v = sniff_frame(d)
                if not v["is_abundance"]:
                    continue
                n_accepted += 1
                pipe = detect_pipeline(name, sub_id, " ".join(str(c) for c in d.columns))
                if pipe:
                    pipelines_found.update(pipe.split(";"))
                long_rows = to_long_rows(d, v, paper_ids.split(";")[0] if paper_ids else "",
                                         dsid, f"{name}!{sub_id}")
                if not long_rows.empty:
                    long_rows["pipeline_source"] = pipe
                    # check-exists + append must be atomic across worker
                    # threads, or two datasets writing "first ever" rows at
                    # the same moment can both see "no file yet" and both
                    # write a header row into the same gzip stream.
                    lock_ctx = long_lock if long_lock is not None else threading.Lock()
                    with lock_ctx:
                        write_header = not long_path.exists()
                        long_rows.to_csv(long_path, sep="\t", index=False, mode="a",
                                         header=write_header, compression="gzip" if long_path.suffix == ".gz" else None)
                    n_rows_written += len(long_rows)
                verdicts.append({"file": f"{name}!{sub_id}", "verdict": "accepted" if not v["needs_review"] else "needs_review",
                                 "axis": v.get("axis", ""), "value_type": v.get("value_type", ""),
                                 "pipeline": pipe, "score": round(v.get("score", 0), 3)})
            local.unlink(missing_ok=True)  # always drop the raw file, verified or not
            gc.collect()
        accepted = sum(1 for v in verdicts if v["verdict"] in ("accepted", "needs_review"))
        all_pipelines = sorted({p for v in verdicts for p in v.get("pipeline", "").split(";") if p})
        return {"dataset_id": dsid, "content_verified": accepted > 0,
                "n_tables_checked": sum(1 for v in verdicts if "download_failed" not in v["verdict"] and "skipped" not in v["verdict"]),
                "n_tables_accepted": accepted, "n_rows_written": n_rows_written,
                "data_category": "processed" if accepted > 0 else "unknown",
                "pipeline_source": ";".join(all_pipelines),
                "note": json.dumps(verdicts)[:2000]}
    finally:
        if study_dir.exists():
            shutil.rmtree(study_dir, ignore_errors=True)
        gc.collect()


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #

def check_dependencies() -> None:
    """Fail fast and loud, before touching a single dataset, if a table
    format this stage needs to read can't actually be opened. Without this,
    a missing package degrades silently: iter_tables()'s broad except just
    logs "iter_tables error: ... openpyxl failed" per FILE to stderr and
    yields zero tables for it, so every .xlsx in the run gets silently
    scored as "no abundance matrix found" instead of "never actually
    opened" -- indistinguishable from a real negative in the output CSVs,
    and easy to miss in a run over tens of thousands of files. This was
    observed for real: an entire cluster run logged nothing but per-file
    openpyxl ImportErrors, meaning none of that run's .xlsx supplementary
    tables (the majority format encountered in practice) were ever
    content-checked at all.
    """
    missing = []
    try:
        import openpyxl  # noqa: F401  -- reads .xlsx
    except ImportError:
        missing.append(("openpyxl", ".xlsx"))
    try:
        import xlrd  # noqa: F401  -- reads legacy .xls
    except ImportError:
        missing.append(("xlrd", ".xls"))
    if missing:
        pkgs = " ".join(p for p, _ in missing)
        fmts = ", ".join(f"{p} ({fmt})" for p, fmt in missing)
        print(f"FATAL: missing required package(s) for reading spreadsheet tables: {fmts}\n"
             f"Install into the active environment before rerunning, e.g.:\n"
             f"    pip install {pkgs}\n"
             f"or:\n"
             f"    conda install -n <env> {pkgs}\n"
             f"Running without these does not fail loudly -- every affected file is silently "
             f"scored as \"no abundance matrix found\" instead of \"never opened\".",
             file=sys.stderr)
        sys.exit(1)


def main() -> None:
    check_dependencies()
    base = Path(__file__).resolve().parent / "results"
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", default=str(base / "stage5_final/dataset_candidates_final.csv"))
    ap.add_argument("--papers", default=str(base / "stage5_final/paper_candidates_final.csv"))
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--scratch", default=None,
                    help="fast local scratch dir for downloads; defaults to <outdir>/_scratch "
                         "but a cluster /scratch path is strongly preferred over $HOME")
    ap.add_argument("--statuses", default="abundance_ready,needs_content_check")
    ap.add_argument("--max-file-mb", type=float, default=200)
    ap.add_argument("--max-study-mb", type=float, default=500)
    ap.add_argument("--sleep", type=float, default=0.2,
                    help="ignored when --workers > 1 -- the shared Throttle paces requests instead")
    ap.add_argument("--workers", type=int, default=1,
                    help="parallel datasets downloaded/verified at once. Each worker still follows "
                         "the same download-one-file -> sniff -> delete -> gc.collect() discipline "
                         "per dataset, so peak disk/memory scales as roughly workers x max-study-mb "
                         "rather than accumulating across the whole run. 4-8 is reasonable; do not "
                         "set this anywhere near --max-study-mb x workers > your node's actual free "
                         "disk under --scratch.")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    scratch = Path(args.scratch) if args.scratch else outdir / "_scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    matrices_dir = outdir / "matrices"
    long_path = outdir / "abundance_long.tsv.gz"
    verif_path = outdir / "dataset_verification.csv"
    long_lock = threading.Lock()  # guards the shared abundance_long.tsv.gz across worker threads

    wanted = set(args.statuses.split(","))
    datasets = [r for r in csv.DictReader(open(args.datasets, newline="", encoding="utf-8"))
               if r.get("candidate_status") in wanted]
    print(f"{len(datasets)} datasets queued for content verification")

    done_ids = set()
    if verif_path.exists():
        done_ids = {r["dataset_id"] for r in csv.DictReader(open(verif_path, newline="", encoding="utf-8"))}
        print(f"resuming: {len(done_ids)} datasets already verified")
    pending = [row for row in datasets if row["dataset_id"] not in done_ids]

    fieldnames = ["dataset_id", "content_verified", "n_tables_checked",
                 "n_tables_accepted", "n_rows_written", "data_category",
                 "pipeline_source", "note"]
    fh = open(verif_path, "a", newline="", encoding="utf-8")
    w = csv.DictWriter(fh, fieldnames=fieldnames)
    if verif_path.stat().st_size == 0:
        w.writeheader()

    if args.workers <= 1:
        for i, row in enumerate(pending, 1):
            result = verify_dataset(row, scratch, args.max_file_mb, args.max_study_mb,
                                    long_path, matrices_dir, long_lock)
            w.writerow({k: result.get(k, "") for k in fieldnames})
            fh.flush()
            if i % 20 == 0:
                print(f"verified {i}/{len(pending)} (disk under {scratch} should be near-empty between rows)",
                     flush=True)
            time.sleep(args.sleep)
    else:
        # Every worker downloads its OWN dataset into scratch/<dataset_id>/
        # (verify_dataset derives that subdirectory from dataset_id, so
        # concurrent workers never share a path) and deletes it before
        # returning -- peak scratch usage is bounded by
        # (in-flight workers) x max-study-mb, not the whole queue.
        with ThreadPoolExecutor(args.workers) as ex:
            futs = {ex.submit(verify_dataset, row, scratch, args.max_file_mb, args.max_study_mb,
                              long_path, matrices_dir, long_lock): row["dataset_id"]
                   for row in pending}
            for i, fu in enumerate(as_completed(futs), 1):
                result = fu.result()
                w.writerow({k: result.get(k, "") for k in fieldnames})
                fh.flush()
                if i % 20 == 0:
                    print(f"verified {i}/{len(pending)} ({args.workers} workers, "
                         f"disk under {scratch} bounded by ~{args.workers}x --max-study-mb)",
                         flush=True)
    fh.close()

    # roll up to paper level: a paper is "processed" if ANY of its datasets
    # verified as an abundance matrix; otherwise it stays whatever stage5
    # called it (raw_reads / unknown) -- stage6 never downgrades a status,
    # it only confirms or leaves it unconfirmed.
    verif_rows = list(csv.DictReader(open(verif_path, newline="", encoding="utf-8")))
    verified = {r["dataset_id"]: r["content_verified"] in ("True", True) for r in verif_rows}
    pipeline_of = {r["dataset_id"]: r.get("pipeline_source", "") for r in verif_rows}
    dsets = list(csv.DictReader(open(args.datasets, newline="", encoding="utf-8")))
    verified_by_paper = defaultdict(bool)
    pipelines_by_paper = defaultdict(set)
    for d in dsets:
        ok = verified.get(d["dataset_id"], False)
        for p in d.get("paper_ids", "").split(";"):
            if p:
                verified_by_paper[p] = verified_by_paper[p] or ok
                if ok and pipeline_of.get(d["dataset_id"]):
                    pipelines_by_paper[p].update(pipeline_of[d["dataset_id"]].split(";"))

    papers = list(csv.DictReader(open(args.papers, newline="", encoding="utf-8")))
    for p in papers:
        cv = verified_by_paper.get(p["paper_id"], False)
        p["content_verified"] = cv
        p["verified_reason"] = "content_matches_abundance_matrix" if cv else "not_yet_content_verified"
        p["data_category"] = "processed" if cv else p.get("data_category", "unknown")
        p["source_pipeline"] = ";".join(sorted(pipelines_by_paper.get(p["paper_id"], set()))) or p.get("source_pipeline_hint", "")
        # next_action tells you what to DO with a paper that isn't a
        # confirmed abundance matrix yet: run your own pipeline on public
        # reads, or open the file yourself because content-checking could
        # not settle it automatically.
        if cv:
            p["next_action"] = "none_have_matrix"
        elif p.get("final_status") in ("raw_reads_ready", "raw_reads_needs_assay_check"):
            p["next_action"] = "run_own_pipeline"
        elif p.get("final_status") == "needs_content_check":
            p["next_action"] = "manual_content_review"
        elif p.get("final_status") == "manual_review":
            p["next_action"] = "manual_content_review"
        else:
            p["next_action"] = "no_data_found"
    pfields = list(papers[0]) if papers else []

    # Full audit trail (every paper, every column) -- keep this even though
    # the three files below are the actual deliverables, since it's the only
    # place n_linked_datasets / dataset_ids / classification_reason survive.
    with open(outdir / "paper_verification.csv", "w", newline="", encoding="utf-8") as f:
        wcsv = csv.DictWriter(f, fieldnames=pfields); wcsv.writeheader(); wcsv.writerows(papers)

    # THE THREE DELIVERABLE FILES. The two unresolved actions are deliberately
    # separate so downstream work does not need to filter a mixed queue.
    # 1. abundance_final.csv: content-verified matrix in hand, nothing left to do.
    with open(outdir / "abundance_final.csv", "w", newline="", encoding="utf-8") as f:
        wcsv = csv.DictWriter(f, fieldnames=pfields); wcsv.writeheader()
        wcsv.writerows([p for p in papers if p["content_verified"]])
    # 2. run_own_pipeline.csv: public raw reads that need our pipeline.
    with open(outdir / "run_own_pipeline.csv", "w", newline="", encoding="utf-8") as f:
        wcsv = csv.DictWriter(f, fieldnames=pfields); wcsv.writeheader()
        wcsv.writerows([p for p in papers if p["next_action"] == "run_own_pipeline"])
    # 3. manual_content_review.csv: unresolved table/file content for a human.
    with open(outdir / "manual_content_review.csv", "w", newline="", encoding="utf-8") as f:
        wcsv = csv.DictWriter(f, fieldnames=pfields); wcsv.writeheader()
        wcsv.writerows([p for p in papers if p["next_action"] == "manual_content_review"])

    # per-dataset matrices, built from the long table without holding it all
    # in memory: read the tsv.gz in chunks, spool one accumulator per dataset.
    if long_path.exists():
        import pandas as pd
        matrices_dir.mkdir(exist_ok=True)
        acc: dict[str, list] = defaultdict(list)
        for chunk in pd.read_csv(long_path, sep="\t", compression="gzip", chunksize=200_000):
            for dsid, g in chunk.groupby("dataset_id"):
                acc[dsid].append(g[["taxon", "sample_id", "value"]])
        for dsid, parts in acc.items():
            g = pd.concat(parts, ignore_index=True)
            mat = g.pivot_table(index="taxon", columns="sample_id", values="value", aggfunc="sum")
            safe = re.sub(r"[^A-Za-z0-9._-]", "_", dsid)[:120]
            mat.to_csv(matrices_dir / f"{safe}__matrix.tsv", sep="\t")
            del g, mat
        gc.collect()

    shutil.rmtree(scratch, ignore_errors=True)
    n_final = sum(p["content_verified"] for p in papers)
    n_pipeline = sum(p["next_action"] == "run_own_pipeline" for p in papers)
    n_review = sum(p["next_action"] == "manual_content_review" for p in papers)
    print(f"done. {n_final}/{len(papers)} papers -> abundance_final.csv; "
          f"{n_pipeline} -> run_own_pipeline.csv; "
          f"{n_review} -> manual_content_review.csv")


if __name__ == "__main__":
    main()
