#!/usr/bin/env python3
"""stage7: for every abundance_ready paper, download the actual file and
extract the abundance matrix -- reusing stage6's already-tested download/
sniff/normalise machinery (list_deposit_files, stream_download,
iter_tables, sniff_frame, to_long_rows) rather than reimplementing it.

Scope: ONLY candidate_status == "abundance_ready" datasets belonging to
final_status == "abundance_ready" papers (per user decision -- not
needs_content_check; that subset is stage6's job if/when it's run).

Outputs, in the exact schema the user's own Google Sheet template uses:
  - one wide matrix per study: <data_dir>/<study_id>_abundance_matrix.tsv
        sample_id, source, then one column per taxon, prefixed
        virus_/prok_ per stage6's domain classification. `source` is the
        dataset_id (repo:accession) the row's values came from.
  - studies.tsv   (appended, matching the user's studies.tsv template columns)
  - sample.tsv    (appended, matching the user's sample.tsv template columns --
                   only fields resolvable from what stage3/4/6 already carry
                   are filled; clinical/demographic fields are left blank
                   per the user's explicit choice, never guessed)
  - blocked_manual_download.tsv  -- abundance_ready papers where every
        attempted file failed to download/parse, WITH the reason (most
        commonly: NCBI PMC's proof-of-work anti-bot challenge -- see
        commit 736a9e6's discovery). These need a human to open a real
        browser and download manually.

study_id = "{FirstAuthorSurname}_{PubYear}" (e.g. "Ramirez_2026"), matching
the user's template exactly. First author comes from a lightweight Europe
PMC lookup per paper (this pipeline has never captured author info before
now) -- verified against a template row: PMID 41725012 ->
authorString "Ramirez AL, ..." + pubYear 2026 -> "Ramirez_2026". Collisions
(two papers, same surname+year) get a/b/c suffixes, assigned once and
recorded in the progress ledger so a resumed run keeps prior assignments
stable rather than recomputing them from scratch.

Two independently resumable phases, same pattern as
augment_europepmc_supp.py:
  extract  -- network phase. One JSON line per paper appended to
              --progress as it resolves (safe to interrupt/rerun; already-
              recorded paper_ids are skipped). Writes matrix files to
              --data-dir as papers succeed.
  build    -- local-only. Reads --progress + --abundance-ready and writes
              studies.tsv / sample.tsv / blocked_manual_download.tsv.
              Rerun freely; it is a pure, fast, idempotent merge.

Usage
-----
    python3 stage7_extract_abundance_matrix.py extract \\
        --abundance-ready Database/results/stage5_final/abundance_ready.csv \\
        --datasets        Database/results/stage5_final/dataset_candidates_final.csv \\
        --data-dir        Database/data \\
        --progress        Database/results/stage7_extract/progress.jsonl \\
        --scratch         /scratch/$USER/mb_stage7 \\
        --workers 8

    python3 stage7_extract_abundance_matrix.py build \\
        --abundance-ready Database/results/stage5_final/abundance_ready.csv \\
        --progress        Database/results/stage7_extract/progress.jsonl \\
        --studies-tsv     Database/metadata/studies.tsv \\
        --sample-tsv      Database/metadata/sample.tsv \\
        --blocked-tsv     Database/metadata/blocked_manual_download.tsv
"""
from __future__ import annotations

import argparse
import csv
import gc
import json
import re
import shutil
import sys
import threading
import time
import unicodedata
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

csv.field_size_limit(sys.maxsize)

# Reuse stage6's already-tested download/parse machinery -- no reimplementation.
sys.path.insert(0, str((Path(__file__).resolve().parent.parent / "find_papers")))
from stage6_verify_abundance import (          # noqa: E402
    list_deposit_files, stream_download, iter_tables, sniff_frame,
    to_long_rows, check_dependencies,
)
from _netutil import GLOBAL_THROTTLE           # noqa: E402

UA = {"User-Agent": "microbiome-dataset-resolver/1.0 (academic research)"}

STUDIES_FIELDS = ["study_id", "title", "doi", "pmid", "year", "disease", "body_site",
                   "sample_type", "n_samples", "data_available", "sequencing_type",
                   "sequencing_platform", "target_type", "available_abundance",
                   "profiling_possible", "metadata_available", "repository",
                   "submission_accession", "url", "include", "notes"]
SAMPLE_FIELDS = ["study_id", "submission_accession", "sample_id", "participant_id",
                  "group", "disease_status", "age", "sex", "country", "city_location",
                  "specimen_type", "run_accession", "collection_year",
                  "sequencing_platform", "antibiotic_use", "treatment_medication",
                  "download_url", "notes"]
BLOCKED_FIELDS = ["paper_id", "pmid", "pmcid", "doi", "title", "dataset_id", "repository",
                   "accession", "landing_url", "attempted_files", "reason", "notes"]

# --------------------------------------------------------------------------- #
# 1. study_id = FirstAuthorSurname_Year, via a lightweight Europe PMC lookup
# --------------------------------------------------------------------------- #

def _ascii_surname(raw: str) -> str:
    """Strip accents/diacritics to plain ASCII letters (Ramirez, not Ram\u00edrez)."""
    norm = unicodedata.normalize("NFKD", raw)
    ascii_only = "".join(c for c in norm if not unicodedata.combining(c))
    return re.sub(r"[^A-Za-z]", "", ascii_only)


def fetch_first_author_year(pmid: str, doi: str) -> tuple[str, str, str]:
    """(surname, year, error). Never fabricates -- returns ("","",err) on
    any failure rather than guessing."""
    query = f"EXT_ID:{pmid}" if pmid else (f"DOI:{doi}" if doi else "")
    if not query:
        return "", "", "no pmid or doi to look up"
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={query}&format=json&resultType=core"
    try:
        GLOBAL_THROTTLE.wait(url)
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=30) as r:
            d = json.loads(r.read())
    except Exception as e:                          # noqa: BLE001
        return "", "", f"{type(e).__name__}: {e}"[:300]
    res = d.get("resultList", {}).get("result", [])
    if not res:
        return "", "", "no Europe PMC record found"
    rec = res[0]
    author_string = rec.get("authorString", "")
    year = str(rec.get("pubYear", "")).strip()
    first = author_string.split(",")[0].strip() if author_string else ""
    # authorString looks like "Ram\u00edrez AL, P\u00e1ez L, ..." -- normalise accents
    # to ASCII BEFORE splitting, so "Ram\u00edrez AL" doesn't fail an ASCII-only regex
    # and get concatenated into "RamirezAL" by the accent-stripping fallback.
    first_ascii = unicodedata.normalize("NFKD", first)
    first_ascii = "".join(c for c in first_ascii if not unicodedata.combining(c))
    tokens = first_ascii.split()
    # last token is the initials (e.g. "AL", "JD") if it's short and all-caps;
    # everything before it is the surname (may itself be multi-word).
    if tokens and re.match(r"^[A-Z]{1,3}$", tokens[-1]):
        surname_raw = " ".join(tokens[:-1])
    else:
        surname_raw = " ".join(tokens)
    surname = re.sub(r"[^A-Za-z]", "", surname_raw)
    if not surname or not year:
        return "", "", f"could not parse author/year from authorString={author_string!r}"
    return surname, year, ""


def assign_study_id(surname: str, year: str, taken: set[str]) -> str:
    base = f"{surname}_{year}"
    if base not in taken:
        return base
    for suffix in "abcdefghijklmnopqrstuvwxyz":
        cand = f"{base}{suffix}"
        if cand not in taken:
            return cand
    # exceedingly unlikely (27+ papers, same first author, same year)
    return f"{base}_{len(taken)}"


# --------------------------------------------------------------------------- #
# 2. anti-bot / challenge-page detection (stage6 doesn't distinguish this --
#    a "download succeeded" that's actually an HTML challenge page just
#    fails later at iter_tables with a confusing "can't determine engine"
#    error; this makes the failure mode explicit for the blocked-list)
# --------------------------------------------------------------------------- #

CHALLENGE_SIGNATURES = (b"<html", b"<!DOCTYPE", b"cloudpmc", b"Preparing to download")
BINARY_TABLE_EXT = {".xlsx", ".xls", ".biom", ".qza"}


def looks_like_challenge_page(path: Path) -> bool:
    if path.suffix.lower() not in BINARY_TABLE_EXT:
        return False  # csv/tsv/txt are plausibly real even if small/textual
    try:
        head = path.read_bytes()[:512]
    except Exception:                                # noqa: BLE001
        return False
    return any(sig in head for sig in CHALLENGE_SIGNATURES)


# --------------------------------------------------------------------------- #
# 3. wide matrix: pivot to_long_rows' output into sample_id x taxon, with
#    virus_/prok_-prefixed column names and a `source` column
# --------------------------------------------------------------------------- #

def build_wide_matrix(long_df):
    """long_df: concatenated to_long_rows() output across all of a study's
    datasets. Returns a wide pandas DataFrame: sample_id, source, then one
    column per taxon (most-specific non-empty rank name), prefixed by
    domain (virus_/prok_)."""
    import pandas as pd
    if long_df.empty:
        return pd.DataFrame()
    df = long_df.copy()
    rank_cols = ["genus", "family", "order", "class", "phylum", "kingdom"]

    def taxon_label(row):
        genus, species = row.get("genus", ""), row.get("species", "")
        if genus and species:
            return f"{genus}_{species}"  # avoid collisions: two genera can share a species epithet
        for r in rank_cols:
            v = row.get(r, "")
            if v:
                return v
        return row.get("taxon", "unknown")

    df["taxon_label"] = df.apply(taxon_label, axis=1)
    prefix = df["domain"].map({"virus": "virus_", "prokaryote": "prok_"}).fillna("prok_")
    df["col_name"] = prefix + df["taxon_label"].astype(str).str.replace(r"[^\w]+", "_", regex=True)
    # a (sample, taxon) pair can appear more than once across multiple
    # files/datasets for the same study -- sum rather than silently drop.
    piv = df.pivot_table(index="sample_id", columns="col_name", values="value", aggfunc="sum", fill_value=0)
    piv = piv.reset_index()
    # `source` records which dataset(s) contributed to each sample's row.
    src = df.groupby("sample_id")["dataset_id"].agg(lambda s: ";".join(sorted(set(s)))).rename("source")
    piv = piv.merge(src, on="sample_id", how="left")
    cols = ["sample_id", "source"] + [c for c in piv.columns if c not in ("sample_id", "source")]
    return piv[cols]


# --------------------------------------------------------------------------- #
# 4. extract phase: one paper at a time -- download every abundance_ready
#    dataset it has, accumulate long rows, write the wide matrix, record
#    the outcome (including which files were anti-bot-blocked, if any)
# --------------------------------------------------------------------------- #

def process_paper(paper: dict, datasets: list[dict], scratch: Path,
                   max_file_mb: float, max_study_mb: float, data_dir: Path,
                   taken_study_ids: set[str], lock: threading.Lock) -> dict:
    paper_id = paper["paper_id"]
    study_dir = scratch / re.sub(r"[^A-Za-z0-9._-]", "_", paper_id)[:100]
    all_long = []
    blocked_files = []   # [(dataset_id, filename, reason)]
    budget_used = 0.0

    for ds in datasets:
        repo, acc, dsid = ds["repository"], ds["accession"], ds["dataset_id"]
        try:
            files = list_deposit_files(repo, acc)
        except Exception as e:                      # noqa: BLE001
            blocked_files.append((dsid, "<listing>", f"list_deposit_files error: {e}"))
            continue
        for name, url, _size in files:
            if budget_used > max_study_mb * 1e6 or not url:
                continue
            local = study_dir / re.sub(r"[^A-Za-z0-9._-]", "_", name)[:100]
            ok, status = stream_download(url, local, max_file_mb)
            if not ok:
                blocked_files.append((dsid, name, f"download_failed:{status}"))
                continue
            budget_used += local.stat().st_size
            if looks_like_challenge_page(local):
                blocked_files.append((dsid, name, "anti_bot_challenge_page"))
                local.unlink(missing_ok=True)
                continue
            try:
                for _sub_id, _hr, d in iter_tables(local):
                    v = sniff_frame(d)
                    if not v["is_abundance"]:
                        continue
                    rows = to_long_rows(d, v, paper_id, dsid, name)
                    if not rows.empty:
                        all_long.append(rows)
            except Exception as e:                   # noqa: BLE001
                blocked_files.append((dsid, name, f"parse_error:{type(e).__name__}"))
            local.unlink(missing_ok=True)
        gc.collect()

    shutil.rmtree(study_dir, ignore_errors=True)

    if not all_long:
        return {"paper_id": paper_id, "status": "blocked", "study_id": "",
                "n_samples": 0, "n_taxa": 0, "matrix_path": "",
                "blocked_files": blocked_files, "error": ""}

    import pandas as pd
    long_df = pd.concat(all_long, ignore_index=True)
    wide = build_wide_matrix(long_df)
    if wide.empty:
        return {"paper_id": paper_id, "status": "blocked", "study_id": "",
                "n_samples": 0, "n_taxa": 0, "matrix_path": "",
                "blocked_files": blocked_files, "error": "pivot produced empty matrix"}

    surname, year, err = fetch_first_author_year(paper.get("pmid", ""), paper.get("doi", ""))
    with lock:
        if surname and year:
            study_id = assign_study_id(surname, year, taken_study_ids)
        else:
            study_id = f"paper{paper_id.replace('/', '_').replace(':', '_')[:40]}"
        taken_study_ids.add(study_id)

    data_dir.mkdir(parents=True, exist_ok=True)
    out_path = data_dir / f"{study_id}_abundance_matrix.tsv"
    wide.to_csv(out_path, sep="\t", index=False)

    return {"paper_id": paper_id, "status": "ok", "study_id": study_id,
            "n_samples": int(wide.shape[0]), "n_taxa": int(wide.shape[1] - 2),
            "matrix_path": str(out_path), "blocked_files": blocked_files,
            "author_lookup_error": err, "domains": sorted(long_df["domain"].unique().tolist())}


def cmd_extract(args: argparse.Namespace) -> None:
    check_dependencies()
    with open(args.abundance_ready, newline="", encoding="utf-8") as f:
        papers = [r for r in csv.DictReader(f) if r.get("final_status") == "abundance_ready"]
    with open(args.datasets, newline="", encoding="utf-8") as f:
        all_ds = list(csv.DictReader(f))
    by_ds_id = {d["dataset_id"]: d for d in all_ds}

    def datasets_for(p: dict) -> list[dict]:
        ids = [x for x in p.get("dataset_ids", "").split(";") if x]
        return [by_ds_id[i] for i in ids if i in by_ds_id and by_ds_id[i].get("candidate_status") == "abundance_ready"]

    papers = sorted(papers, key=lambda p: p["paper_id"])  # deterministic order across resumed runs
    print(f"{len(papers)} abundance_ready papers")

    done_ids, taken_study_ids = set(), set()
    if Path(args.progress).exists():
        with open(args.progress, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                done_ids.add(rec["paper_id"])
                if rec.get("study_id"):
                    taken_study_ids.add(rec["study_id"])
    pending = [p for p in papers if p["paper_id"] not in done_ids]
    print(f"{len(done_ids)} already processed, {len(pending)} pending")
    if not pending:
        print("nothing to do -- run the build phase.")
        return

    scratch = Path(args.scratch); scratch.mkdir(parents=True, exist_ok=True)
    data_dir = Path(args.data_dir)
    out_fh = open(args.progress, "a", encoding="utf-8")
    write_lock = threading.Lock()
    study_id_lock = threading.Lock()
    n_done = [0]
    t0 = time.time()

    def one(p: dict) -> None:
        ds = datasets_for(p)
        if not ds:
            rec = {"paper_id": p["paper_id"], "status": "blocked", "study_id": "",
                   "n_samples": 0, "n_taxa": 0, "matrix_path": "",
                   "blocked_files": [], "error": "no abundance_ready dataset resolved for this paper"}
        else:
            rec = process_paper(p, ds, scratch, args.max_file_mb, args.max_study_mb,
                                data_dir, taken_study_ids, study_id_lock)
        with write_lock:
            out_fh.write(json.dumps(rec) + "\n")
            out_fh.flush()
            n_done[0] += 1
            if n_done[0] % 100 == 0:
                elapsed = time.time() - t0
                rate = n_done[0] / elapsed if elapsed > 0 else 0
                eta_min = (len(pending) - n_done[0]) / rate / 60 if rate > 0 else float("inf")
                print(f"  processed {n_done[0]}/{len(pending)}  ({rate:.2f}/s, ETA {eta_min:.0f} min)", flush=True)

    with ThreadPoolExecutor(args.workers) as ex:
        futs = [ex.submit(one, p) for p in pending]
        for fu in as_completed(futs):
            fu.result()
    out_fh.close()
    print(f"done. {n_done[0]} papers processed this run.")


# --------------------------------------------------------------------------- #
# 5. build phase: local-only merge into studies.tsv / sample.tsv / blocked
# --------------------------------------------------------------------------- #

def append_tsv(path: Path, fieldnames: list[str], rows: list[dict], key_fields: tuple[str, ...]) -> int:
    """Append rows not already present (by key_fields) -- idempotent on rerun.
    Returns the number of rows actually written."""
    path.parent.mkdir(parents=True, exist_ok=True)
    existing_keys = set()
    write_header = not path.exists()
    if path.exists():
        with open(path, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                existing_keys.add(tuple(r.get(k, "") for k in key_fields))
    new_rows = [row for row in rows if tuple(str(row.get(k, "")) for k in key_fields) not in existing_keys]
    with open(path, "a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        if write_header:
            w.writeheader()
        for row in new_rows:
            w.writerow({k: row.get(k, "") for k in fieldnames})
    return len(new_rows)


def cmd_build(args: argparse.Namespace) -> None:
    with open(args.abundance_ready, newline="", encoding="utf-8") as f:
        papers_by_id = {r["paper_id"]: r for r in csv.DictReader(f)}

    if not Path(args.progress).exists():
        print(f"FATAL: {args.progress} does not exist -- run 'extract' first", file=sys.stderr)
        sys.exit(1)

    studies_rows, sample_rows, blocked_rows = [], [], []
    with open(args.progress, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            paper = papers_by_id.get(rec["paper_id"])
            if paper is None:
                continue

            if rec["status"] == "ok":
                repo_accs = sorted({(d.split(":", 1)[0], d.split(":", 1)[1])
                                    for d in (paper.get("dataset_ids") or "").split(";") if ":" in d})
                studies_rows.append({
                    "study_id": rec["study_id"], "title": paper.get("title", ""),
                    "doi": paper.get("doi", ""), "pmid": paper.get("pmid", ""),
                    "year": paper.get("publication_year", ""), "disease": "",
                    "body_site": "", "sample_type": "", "n_samples": rec["n_samples"],
                    "data_available": "processed",
                    "sequencing_type": "", "sequencing_platform": "",
                    "target_type": "+".join(rec.get("domains", [])),
                    "available_abundance": "yes", "profiling_possible": "from_matrix",
                    "metadata_available": "no",
                    "repository": ";".join(sorted({r for r, _ in repo_accs})),
                    "submission_accession": ";".join(a for _, a in repo_accs),
                    "url": f"https://pmc.ncbi.nlm.nih.gov/articles/{paper.get('pmcid','')}/" if paper.get("pmcid") else "",
                    "include": "yes",
                    "notes": (f"extracted {rec['n_samples']} samples x {rec['n_taxa']} taxa "
                              f"from {rec['matrix_path']}."
                              + (f" author lookup failed ({rec.get('author_lookup_error')}) "
                                 "-- study_id falls back to paper_id." if rec.get("author_lookup_error") else "")),
                })
                import csv as _csv  # local alias, avoids shadowing concerns
                try:
                    with open(rec["matrix_path"], newline="", encoding="utf-8") as mf:
                        for row in _csv.DictReader(mf, delimiter="\t"):
                            sample_rows.append({
                                "study_id": rec["study_id"],
                                "submission_accession": ";".join(a for _, a in repo_accs),
                                "sample_id": row.get("sample_id", ""),
                                "participant_id": "", "group": "", "disease_status": "",
                                "age": "", "sex": "", "country": "", "city_location": "",
                                "specimen_type": "", "run_accession": "", "collection_year": "",
                                "sequencing_platform": "", "antibiotic_use": "",
                                "treatment_medication": "", "download_url": "",
                                "notes": "clinical/demographic fields not extracted at this stage "
                                        "-- only machine-derivable identifiers populated",
                            })
                except FileNotFoundError:
                    pass

            if rec["status"] == "blocked" or rec.get("blocked_files"):
                reasons = rec.get("blocked_files", []) or [("", "", rec.get("error", "no dataset resolved"))]
                for dsid, fname, reason in reasons:
                    ds_repo = dsid.split(":", 1)[0] if ":" in dsid else ""
                    ds_acc = dsid.split(":", 1)[1] if ":" in dsid else ""
                    blocked_rows.append({
                        "paper_id": rec["paper_id"], "pmid": paper.get("pmid", ""),
                        "pmcid": paper.get("pmcid", ""), "doi": paper.get("doi", ""),
                        "title": paper.get("title", ""), "dataset_id": dsid,
                        "repository": ds_repo, "accession": ds_acc, "landing_url": "",
                        "attempted_files": fname, "reason": reason,
                        "notes": "anti_bot_challenge_page = needs a real browser to download manually"
                                if reason == "anti_bot_challenge_page" else "",
                    })

    print(f"studies: {len(studies_rows)}   samples: {len(sample_rows)}   blocked entries: {len(blocked_rows)}")
    n1 = append_tsv(Path(args.studies_tsv), STUDIES_FIELDS, studies_rows, key_fields=("study_id",))
    n2 = append_tsv(Path(args.sample_tsv), SAMPLE_FIELDS, sample_rows, key_fields=("study_id", "sample_id"))
    n3 = append_tsv(Path(args.blocked_tsv), BLOCKED_FIELDS, blocked_rows,
                    key_fields=("paper_id", "dataset_id", "attempted_files"))
    print(f"wrote {n1} new study rows, {n2} new sample rows, {n3} new blocked entries "
          f"(rows already present from an earlier build run were skipped)")
    n_blocked_papers = len({r["paper_id"] for r in blocked_rows})
    print(f"\n{len(studies_rows)} studies got a matrix written to disk. "
          f"{n_blocked_papers} abundance_ready papers had every file blocked -- "
          f"see {args.blocked_tsv} for which ones and why (most commonly "
          f"anti_bot_challenge_page = manual browser download needed).")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    ep = sub.add_parser("extract", help="network phase: download + parse, write matrices + progress")
    ep.add_argument("--abundance-ready", required=True)
    ep.add_argument("--datasets", required=True)
    ep.add_argument("--data-dir", required=True)
    ep.add_argument("--progress", required=True)
    ep.add_argument("--scratch", required=True)
    ep.add_argument("--workers", type=int, default=8)
    ep.add_argument("--max-file-mb", type=float, default=200)
    ep.add_argument("--max-study-mb", type=float, default=500)

    bp = sub.add_parser("build", help="local phase: assemble studies.tsv/sample.tsv/blocked list")
    bp.add_argument("--abundance-ready", required=True)
    bp.add_argument("--progress", required=True)
    bp.add_argument("--studies-tsv", required=True)
    bp.add_argument("--sample-tsv", required=True)
    bp.add_argument("--blocked-tsv", required=True)

    args = ap.parse_args()
    if args.cmd == "extract":
        cmd_extract(args)
    else:
        cmd_build(args)


if __name__ == "__main__":
    main()
