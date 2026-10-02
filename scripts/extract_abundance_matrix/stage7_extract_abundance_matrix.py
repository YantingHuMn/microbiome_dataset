#!/usr/bin/env python3
"""stage7: for every abundance_ready paper, download the actual file and
extract the abundance matrix -- reusing stage6's already-tested download/
sniff/normalise machinery (list_deposit_files, stream_download,
iter_tables, sniff_frame, to_long_rows) rather than reimplementing it.

Scope: ONLY candidate_status == "abundance_ready" datasets belonging to
final_status == "abundance_ready" papers (per user decision -- not
needs_content_check; that subset is stage6's job if/when it's run).

Outputs, in the exact schema the user's own Google Sheet template uses:
  - one wide matrix per (dataset, source file, sub-table): <data_dir>/
        <study_id>/<safe_name>.matrix.tsv -- sample_id, source, then one
        column per taxon, using the BARE original name (no domain
        prefix -- see build_table_matrix/build_taxonomy_table). `source`
        is the dataset_id (repo:accession) the row's values came from.
  - a sibling <safe_name>.taxonomy.tsv for every matrix: one row per
        column, recording what that column name actually means -- the
        full raw taxon string as deposited, the per-rank breakdown, the
        domain call, and its evidence. This is the reference for "what is
        this column", kept separate from the matrix header on purpose: a
        wrong classification no longer means renaming a column
        invalidates anything already keyed off its name.
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

"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gc
import io
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
    to_long_rows, check_dependencies, retain_raw_file, append_raw_manifest, safe_pivot,
)
from _netutil import GLOBAL_THROTTLE           # noqa: E402

UA = {"User-Agent": "microbiome-dataset-resolver/1.0 (academic research)"}

STUDIES_FIELDS = ["study_id", "title", "doi", "pmid", "year", "disease", "body_site",
                   "sample_type", "n_tables", "n_sample_ids_raw_total", "n_verified_biological_samples",
                   "data_available", "sequencing_type", "sequencing_platform", "target_type",
                   "available_abundance", "profiling_possible", "metadata_available", "repository",
                   "submission_accession", "url", "include", "notes"]
SAMPLE_FIELDS = ["study_id", "dataset_id", "source_file", "submission_accession", "sample_id",
                  "sample_flag", "sequencing_type_hint", "participant_id",
                  "group", "disease_status", "age", "sex", "country", "city_location",
                  "specimen_type", "run_accession", "collection_year",
                  "sequencing_platform", "antibiotic_use", "treatment_medication",
                  "download_url", "notes"]
BLOCKED_FIELDS = ["paper_id", "pmid", "pmcid", "doi", "title", "dataset_id", "repository",
                   "accession", "landing_url", "attempted_files", "reason", "notes"]
# One row per (study, dataset, source_file) actually extracted -- the authoritative,
# never-merged per-table manifest. studies.tsv/sample.tsv are convenience rollups;
# THIS file is where n_sample_ids_raw/domains/flags/sha256 for any one table live.
TABLE_MANIFEST_FIELDS = ["study_id", "paper_id", "dataset_id", "source_file", "matrix_path",
                          "long_path", "n_sample_ids_raw", "n_verified_biological_samples",
                          "n_taxa", "domains", "sequencing_type_hints",
                          "blank_or_control_sample_ids", "n_duplicate_cells", "raw_file_sha256"]
NEEDS_REVIEW_FIELDS = ["paper_id", "dataset_id", "source_file", "reason", "score"]

# 1. study_id = FirstAuthorSurname_Year, via a lightweight Europe PMC lookup
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


# 2. anti-bot / challenge-page detection (stage6 doesn't distinguish this --
#    a "download succeeded" that's actually an HTML challenge page just
#    fails later at iter_tables with a confusing "can't determine engine"
#    error; this makes the failure mode explicit for the blocked-list)
CHALLENGE_SIGNATURES = (b"<html", b"<!DOCTYPE", b"cloudpmc", b"Preparing to download")


def looks_like_challenge_page(path: Path) -> bool:
    """Checked on EVERY downloaded file regardless of extension -- a bug
    fix from an earlier version that only checked a fixed set of "binary"
    extensions (.xlsx/.xls/.biom/.qza) and missed .zip entirely, letting
    real cases (mmc1.zip, mmc2.zip, FSN3-11-3154-s001.zip -- Elsevier/
    Wiley supplementary bundles scraped from a PMC page's Associated Data
    links, same anti-bot wall as commit 736a9e6) fall through: the
    "download" succeeds (HTTP 200, a small HTML challenge page) but
    zipfile.ZipFile() then raises "File is not a zip file" -- and because
    that happens inside stage6's iter_tables(), which SWALLOWS its own
    exceptions (prints to stderr, yields nothing), the caller never even
    sees an error: the file just silently looks like "no abundance data
    found here" instead of "blocked". The signature check itself is
    specific enough (exact HTML/challenge-page markers) that running it
    unconditionally is safe -- no genuine CSV/TSV/XLSX content would ever
    contain these bytes in its first 512 bytes.
    """
    try:
        head = path.read_bytes()[:512]
    except Exception:                                # noqa: BLE001
        return False
    return any(sig in head for sig in CHALLENGE_SIGNATURES)


# 3. wide matrix: pivot to_long_rows' output into sample_id x taxon, with
#    virus_/prok_-prefixed column names and a `source` column
class MatrixTooLargeError(Exception):
    """Raised by build_table_matrix() instead of letting pivot_table()
    attempt a combinatorially huge reshape. Even a SINGLE table (let alone
    several merged together, which this stage no longer does by default)
    can in principle have a pathological sample x taxon cardinality; this
    guard stays per-table so one huge table cannot OOM the whole paper's
    processing."""


def col_name_for(taxon_string) -> str:
    """Sanitize (filesystem/TSV-safe characters only) the RAW taxon string
    exactly as the source deposited it -- e.g. "d__Bacteria;p__
    Cloacimonadota" stays "d__Bacteria;p__Cloacimonadota" (semicolons ->
    underscores), it is never shortened to just "Cloacimonadota" or any
    other derived label. If a user wants to know what a column means, the
    sibling *.taxonomy.tsv (build_taxonomy_table) is where that lookup
    happens -- the matrix header itself is never "cleaned up" on their
    behalf. Shared by build_table_matrix and build_taxonomy_table so the
    matrix column name and the taxonomy sidecar's col_name always agree
    (a genuine source duplicate -- identical raw string after to_long_rows'
    pandas-dedup-suffix normalization -- still collapses to one col_name
    and is still caught by safe_pivot's dup_report, same as before)."""
    import re
    return re.sub(r"[^\w]+", "_", str(taxon_string).strip())


def build_taxonomy_table(table_df):
    """One row per DISTINCT col_name in this table (the same col_name that
    ends up as a matrix column header), recording what that name actually
    MEANS -- the full raw taxon string as the source deposited it, the
    per-rank breakdown split_taxonomy() already parsed out of it, and the
    domain call + its evidence. This is the reference a user consults to
    answer "what is this column", NOT the matrix header itself -- the
    header is deliberately left as the bare original name (see
    build_table_matrix), so this file is the only place the domain
    classification and its evidence are recorded at all.

    lineage_source is 'explicit_in_source' when the raw taxon string
    itself yielded at least one non-empty rank (i.e. it used a recognized
    rank-prefix convention like d__/p__/c__ or an unprefixed positional
    lineage) -- in that case the source already told us the lineage and
    no external reference lookup is needed. It is
    'no_lineage_info_in_source' when NONE of the ranks could be parsed
    (e.g. a bare OTU/ASV id with no attached taxonomy string at all) --
    these are the only rows where cross-checking an external authority
    (e.g. SILVA) would add information the source itself doesn't give.
    """
    import pandas as pd
    if table_df.empty:
        return pd.DataFrame()
    df = table_df.copy()
    df["col_name"] = df["taxon"].map(col_name_for)
    rank_fields = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    for r in rank_fields:
        if r not in df.columns:
            df[r] = ""
    df["lineage_source"] = df[rank_fields].apply(
        lambda r: "explicit_in_source" if any(str(v).strip() for v in r) else "no_lineage_info_in_source",
        axis=1)
    cols = ["col_name", "taxon"] + rank_fields + ["domain", "domain_evidence", "lineage_source"]
    tax = df[cols].drop_duplicates(subset=["col_name", "taxon"]).rename(columns={"taxon": "raw_taxon_string"})
    return tax.sort_values(["lineage_source", "col_name"]).reset_index(drop=True)


def build_table_matrix(table_df, max_cells: int = 20_000_000):
    """table_df: to_long_rows() output for ONE (dataset_id, source_file,
    sub_table) -- NEVER multiple tables concatenated together; this
    function does not merge anything across files, sheets, or sequencing
    types. Returns (wide_df, dup_report):
      wide_df: sample_id, source, then one column per taxon, using the
        RAW taxon string exactly as the source deposited it (e.g.
        "d__Bacteria;p__Cloacimonadota" stays exactly that, sanitized only
        for filesystem/TSV-safe characters -- never shortened to just
        "Cloacimonadota" or any other derived label, and never given a
        virus_/prok_/unk_ domain prefix). What a column actually MEANS
        (full per-rank breakdown, domain call + evidence) lives in the
        sibling *.taxonomy.tsv from build_taxonomy_table(), not baked into
        or inferred from the header -- a wrong classification no longer
        means renaming a column retroactively invalidates anything that
        already keyed off its name, and a user who wants to know what a
        column represents always looks it up there rather than relying on
        an abbreviated name we chose on their behalf. A (sample, taxon)
        cell with NO underlying row is left BLANK (NaN), never silently
        zero-filled.
      dup_report: rows where the SAME (sample, taxon) pair had more than
        one raw value within this one table -- these are surfaced, not
        summed; the corresponding wide_df cell is left blank pending
        manual resolution.
    Raises MatrixTooLargeError before ever calling pivot_table if
    n_unique_samples * n_unique_taxa would exceed max_cells (20M cells by
    default, ~160MB of float64 -- generous for any single real table)."""
    import pandas as pd
    if table_df.empty:
        return pd.DataFrame(), pd.DataFrame()
    df = table_df.copy()
    df["col_name"] = df["taxon"].map(col_name_for)

    n_samples = df["sample_id"].nunique()
    n_cols = df["col_name"].nunique()
    est_cells = n_samples * n_cols
    if est_cells > max_cells:
        raise MatrixTooLargeError(
            f"{n_samples} unique samples x {n_cols} unique taxa = {est_cells:,} estimated cells "
            f"(cap {max_cells:,}) -- refusing to pivot this single table")

    piv, dup = safe_pivot(df, "sample_id", "col_name", "value")
    piv = piv.reset_index()
    src = df.groupby("sample_id")["dataset_id"].agg(lambda s: ";".join(sorted(set(s)))).rename("source")
    piv = piv.merge(src, on="sample_id", how="left")
    cols = ["sample_id", "source"] + [c for c in piv.columns if c not in ("sample_id", "source")]
    return piv[cols], dup


# 4. extract phase: one paper at a time -- download every abundance_ready
#    dataset it has, accumulate long rows, write the wide matrix, record
#    the outcome (including which files were anti-bot-blocked, if any)
def process_paper(paper: dict, datasets: list[dict], scratch: Path,
                   max_file_mb: float, max_study_mb: float, data_dir: Path,
                   taken_study_ids: set[str], lock: threading.Lock,
                   raw_store: Path, raw_manifest_path: Path, raw_manifest_lock: threading.Lock) -> dict:
    """Download every abundance_ready dataset/file this paper has, and for
    EACH (dataset, file, sub_table) that passes sniff_frame as a real
    abundance table, write its OWN matrix + long-format output -- never
    merged with any other file, sheet, or sub_table, even from the same
    paper or dataset. A paper that ends up with several tables ends up
    with several output files side by side under data_dir/<study_id>/; it
    is the user's call, made with the per-table manifest in hand, whether
    any of them actually belong combined -- this function never guesses.
    """
    paper_id = paper["paper_id"]
    study_dir = scratch / re.sub(r"[^A-Za-z0-9._-]", "_", paper_id)[:100]
    blocked_files = []     # [(dataset_id, filename, reason)]
    needs_review = []      # [{dataset_id, source_file, reason, score}]
    tables_out = []        # one entry per (dataset, file, sub_table) actually written
    study_id_box = [None]  # assigned lazily -- only once there is something to write

    # Per-table row cap (NOT a whole-paper accumulation cap any more -- since
    # tables are no longer concatenated across files, the OOM mechanism that
    # motivated the old whole-paper cap no longer applies the same way; a
    # single pathological table can still be huge, so it still gets capped).
    MAX_TABLE_ROWS = 2_000_000

    def get_study_id() -> str:
        if study_id_box[0] is not None:
            return study_id_box[0]
        surname, year, err = fetch_first_author_year(paper.get("pmid", ""), paper.get("doi", ""))
        with lock:
            if surname and year:
                sid = assign_study_id(surname, year, taken_study_ids)
            else:
                sid = f"paper{paper_id.replace('/', '_').replace(':', '_')[:40]}"
            taken_study_ids.add(sid)
        study_id_box[0] = sid
        study_id_box.append(err)  # stash the author-lookup error alongside, read back below
        return sid

    for ds in datasets:
        repo, acc, dsid = ds["repository"], ds["accession"], ds["dataset_id"]
        try:
            files = list_deposit_files(repo, acc)
        except Exception as e:                      # noqa: BLE001
            blocked_files.append((dsid, "<listing>", f"list_deposit_files error: {e}"))
            continue
        for name, url, _size in files:
            if not url:
                continue
            local = study_dir / re.sub(r"[^A-Za-z0-9._-]", "_", name)[:100]
            ok, status = stream_download(url, local, max_file_mb)
            if not ok:
                blocked_files.append((dsid, name, f"download_failed:{status}"))
                continue
            # Raw-file retention: EVERY downloaded file is kept byte-for-byte,
            # including ones that turn out to be anti-bot challenge pages --
            # the retained copy is the proof of what was actually served.
            manifest_row = retain_raw_file(local, raw_store, paper_id, dsid, repo, acc, url, name)
            append_raw_manifest(raw_manifest_path, manifest_row, raw_manifest_lock)
            local = Path(manifest_row["retained_path"])

            if looks_like_challenge_page(local):
                blocked_files.append((dsid, name, "anti_bot_challenge_page"))
                continue
            # iter_tables (stage6, shared code) SWALLOWS its own top-level
            # exceptions internally (prints to stderr, yields nothing) --
            # capture that stderr so a nested-member failure (corrupt zip
            # member, mislabeled ".gz", ...) gets a real, traceable reason
            # instead of silently vanishing. preview=False means this is a
            # FULL read, not the row-capped structural preview -- the whole
            # point of splitting the two apart.
            stderr_buf = io.StringIO()
            bad_lines_log: list[str] = []
            try:
                with contextlib.redirect_stderr(stderr_buf):
                    for sub_id, _hr, d in iter_tables(local, preview=False, bad_lines_log=bad_lines_log):
                        v = sniff_frame(d)
                        if not v["is_abundance"]:
                            continue
                        source_file = f"{name}!{sub_id}"
                        if v["needs_review"]:
                            # Low-confidence tables are NEVER auto-included
                            # in the formal output -- they go to a human
                            # review queue instead.
                            needs_review.append({"paper_id": paper_id, "dataset_id": dsid,
                                                "source_file": source_file, "reason": v.get("reason", ""),
                                                "score": round(v.get("score", 0), 3)})
                            continue
                        rows = to_long_rows(d, v, paper_id, dsid, source_file)
                        if rows.empty:
                            continue
                        if len(rows) > MAX_TABLE_ROWS:
                            blocked_files.append((dsid, source_file,
                                f"table_too_large:{len(rows)}_rows_exceeds_{MAX_TABLE_ROWS}"))
                            continue
                        try:
                            wide, dup = build_table_matrix(rows)
                        except MatrixTooLargeError as e:
                            blocked_files.append((dsid, source_file, f"matrix_too_large:{e}"))
                            continue
                        if wide.empty:
                            blocked_files.append((dsid, source_file, "pivot_produced_empty_matrix"))
                            continue

                        study_id = get_study_id()
                        study_out_dir = data_dir / study_id
                        study_out_dir.mkdir(parents=True, exist_ok=True)
                        safe_name = re.sub(r"[^A-Za-z0-9._-]", "_", f"{dsid}__{source_file}")[:150]
                        matrix_path = study_out_dir / f"{safe_name}.matrix.tsv"
                        long_path = study_out_dir / f"{safe_name}.long.tsv.gz"
                        taxonomy_path = study_out_dir / f"{safe_name}.taxonomy.tsv"
                        wide.to_csv(matrix_path, sep="\t", index=False)
                        rows.to_csv(long_path, sep="\t", index=False, compression="gzip")
                        if not dup.empty:
                            dup.to_csv(study_out_dir / f"{safe_name}.duplicates.tsv", sep="\t", index=False)
                        taxonomy = build_taxonomy_table(rows)
                        if not taxonomy.empty:
                            taxonomy.to_csv(taxonomy_path, sep="\t", index=False)

                        blank_ids = sorted(rows.loc[rows["sample_flag"] == "blank_or_control", "sample_id"].unique().tolist())
                        seq_hints = sorted({h for h in rows["sequencing_type_hint"] if h})
                        tables_out.append({
                            "dataset_id": dsid, "source_file": source_file,
                            "matrix_path": str(matrix_path), "long_path": str(long_path),
                            "n_sample_ids_raw": int(rows["sample_id"].nunique()),
                            "n_verified_biological_samples": "",  # never auto-filled -- needs metadata to confirm
                            "n_taxa": int(wide.shape[1] - 2),
                            "domains": ";".join(sorted(rows["domain"].unique())),
                            "sequencing_type_hints": ";".join(seq_hints),
                            "blank_or_control_sample_ids": ";".join(blank_ids),
                            "n_duplicate_cells": int(dup.shape[0]) if not dup.empty else 0,
                            "raw_file_sha256": manifest_row["sha256"],
                        })
            except Exception as e:                   # noqa: BLE001
                blocked_files.append((dsid, name, f"parse_error:{type(e).__name__}"))
            captured = stderr_buf.getvalue().strip()
            if captured:
                for eline in captured.splitlines():
                    if eline.strip():
                        blocked_files.append((dsid, name, f"iter_tables_internal_error:{eline.strip()}"))
            if bad_lines_log:
                for bl in bad_lines_log:
                    blocked_files.append((dsid, name, f"malformed_lines_recorded:{bl}"))
        gc.collect()

    shutil.rmtree(study_dir, ignore_errors=True)  # scratch cleanup only -- raw_store is untouched

    if not tables_out:
        return {"paper_id": paper_id, "status": "blocked", "study_id": study_id_box[0] or "",
                "tables": [], "needs_review": needs_review, "blocked_files": blocked_files, "error": ""}

    return {"paper_id": paper_id, "status": "ok", "study_id": study_id_box[0],
            "tables": tables_out, "needs_review": needs_review, "blocked_files": blocked_files,
            "author_lookup_error": study_id_box[1] if len(study_id_box) > 1 else ""}

def parse_range(spec: str | None, n_pending: int) -> tuple[int, int]:
    """'10' -> (0, 10) i.e. the first 10 pending papers (1-indexed display,
    0-indexed slice). '20-30' -> (19, 30) i.e. pending papers 20 through 30
    inclusive, by POSITION in the sorted pending list (not by paper_id).
    None -> the whole pending list."""
    if not spec:
        return 0, n_pending
    spec = spec.strip()
    if "-" in spec:
        a, b = spec.split("-", 1)
        start, end = int(a), int(b)
        if start < 1 or end < start:
            raise ValueError(f"invalid --range {spec!r}: expected 'A-B' with 1 <= A <= B")
        return start - 1, min(end, n_pending)
    n = int(spec)
    if n < 1:
        raise ValueError(f"invalid --range {spec!r}: expected a positive integer or 'A-B'")
    return 0, min(n, n_pending)


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

    if args.paper_id:
        papers = [p for p in papers if p["paper_id"] == args.paper_id]
        print(f"--paper-id given: restricting to {len(papers)} matching paper(s)")
    papers = sorted(papers, key=lambda p: p["paper_id"])  # deterministic order across resumed runs
    print(f"{len(papers)} abundance_ready papers")

    progress_path = Path(args.progress)
    done_ids, taken_study_ids = set(), set()
    if not args.rebuild and progress_path.exists():
        with open(progress_path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                done_ids.add(rec["paper_id"])
                if rec.get("study_id"):
                    taken_study_ids.add(rec["study_id"])
    elif args.rebuild:
        print(f"--rebuild: ignoring {progress_path} entirely -- every paper in --range is reprocessed "
             f"from scratch, written to a NEW progress file; the old one is left untouched.")
    pending_all = [p for p in papers if p["paper_id"] not in done_ids]
    print(f"{len(done_ids)} already processed, {len(pending_all)} pending")

    start, end = parse_range(args.range, len(pending_all))
    pending = pending_all[start:end]
    if args.range:
        print(f"--range {args.range!r}: processing pending papers {start + 1}-{end} of {len(pending_all)} "
             f"({len(pending)} papers this run)")
    if not pending:
        print("nothing to do -- run the build phase.")
        return

    scratch = Path(args.scratch); scratch.mkdir(parents=True, exist_ok=True)
    data_dir = Path(args.data_dir)
    raw_store = Path(args.raw_store); raw_store.mkdir(parents=True, exist_ok=True)
    raw_manifest_path = Path(args.raw_manifest) if args.raw_manifest else raw_store / "raw_files_manifest.csv"
    progress_path.parent.mkdir(parents=True, exist_ok=True)
    out_fh = open(progress_path, "a", encoding="utf-8")
    write_lock = threading.Lock()
    study_id_lock = threading.Lock()
    raw_manifest_lock = threading.Lock()
    n_done = [0]
    t0 = time.time()

    def one(p: dict) -> None:
        ds = datasets_for(p)
        if not ds:
            rec = {"paper_id": p["paper_id"], "status": "blocked", "study_id": "",
                   "tables": [], "needs_review": [],
                   "blocked_files": [], "error": "no abundance_ready dataset resolved for this paper"}
        else:
            # process_paper can, in principle, hit an exception its own
            # internal try/excepts don't cover (a pandas bug, a disk-full
            # write, an unanticipated data shape, ...). Without this catch,
            # that exception propagates through the worker thread's Future
            # and re-raises in the main thread at fu.result() below --
            # crashing the ENTIRE extract run (and, in submit_all_steps.sh,
            # likely preventing the downstream build phase from running at
            # all this submission) over ONE paper, while that paper's
            # record never reaches progress.jsonl either. Catch it here so
            # one paper's crash becomes a normal, traceable "error" record
            # -- exactly like every other failure mode fixed today -- and
            # every other paper still gets processed and recorded.
            try:
                rec = process_paper(p, ds, scratch, args.max_file_mb, args.max_study_mb,
                                    data_dir, taken_study_ids, study_id_lock,
                                    raw_store, raw_manifest_path, raw_manifest_lock)
            except Exception as e:                   # noqa: BLE001
                rec = {"paper_id": p["paper_id"], "status": "blocked", "study_id": "",
                       "tables": [], "needs_review": [], "blocked_files": [],
                       "error": f"process_paper_crashed:{type(e).__name__}:{e}"}
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
    print(f"done. {n_done[0]} papers processed this run. progress written to {progress_path}")


# 5. build phase: local-only merge into studies.tsv / sample.tsv / blocked
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
        print(f"FATAL: {args.progress} does not exist -- run \'extract\' first", file=sys.stderr)
        sys.exit(1)

    studies_rows, sample_rows, blocked_rows, table_rows, review_rows = [], [], [], [], []
    status_by_paper: dict[str, str] = {}
    with open(args.progress, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            status_by_paper[rec["paper_id"]] = rec["status"]
            paper = papers_by_id.get(rec["paper_id"])
            if paper is None:
                continue

            for rv in rec.get("needs_review", []):
                review_rows.append({**{k: rv.get(k, "") for k in NEEDS_REVIEW_FIELDS}})

            tables = rec.get("tables", [])
            if rec["status"] == "ok" and tables:
                repo_accs = sorted({(d.split(":", 1)[0], d.split(":", 1)[1])
                                    for d in (paper.get("dataset_ids") or "").split(";") if ":" in d})
                n_raw_total = sum(t.get("n_sample_ids_raw", 0) for t in tables)
                all_domains = sorted({dm for t in tables for dm in t.get("domains", "").split(";") if dm})
                studies_rows.append({
                    "study_id": rec["study_id"], "title": paper.get("title", ""),
                    "doi": paper.get("doi", ""), "pmid": paper.get("pmid", ""),
                    "year": paper.get("publication_year", ""), "disease": "",
                    "body_site": "", "sample_type": "",
                    "n_tables": len(tables),
                    "n_sample_ids_raw_total": n_raw_total,  # SUM across tables -- NOT deduplicated; the
                                                             # same sample_id in two different tables is
                                                             # counted twice here on purpose, since this
                                                             # stage never asserts those are the same
                                                             # biological sample across different tables.
                    "n_verified_biological_samples": "",    # never auto-filled -- needs metadata to confirm
                    "data_available": "processed",
                    "sequencing_type": "", "sequencing_platform": "",
                    "target_type": "+".join(all_domains),
                    "available_abundance": "yes", "profiling_possible": "from_per_table_matrices",
                    "metadata_available": "no",
                    "repository": ";".join(sorted({r for r, _ in repo_accs})),
                    "submission_accession": ";".join(a for _, a in repo_accs),
                    "url": f"https://pmc.ncbi.nlm.nih.gov/articles/{paper.get('pmcid','')}/" if paper.get("pmcid") else "",
                    "include": "yes",
                    "notes": (f"{len(tables)} separate table(s) extracted, NOT merged -- see "
                              f"table_manifest.tsv for per-table sample/taxon counts, domains, and flags."
                              + (f" author lookup failed ({rec.get('author_lookup_error')}) "
                                 "-- study_id falls back to paper_id." if rec.get("author_lookup_error") else "")),
                })
                for t in tables:
                    table_rows.append({"study_id": rec["study_id"], "paper_id": rec["paper_id"],
                                       **{k: t.get(k, "") for k in TABLE_MANIFEST_FIELDS
                                          if k not in ("study_id", "paper_id")}})
                    try:
                        with open(t["matrix_path"], newline="", encoding="utf-8") as mf:
                            for row in csv.DictReader(mf, delimiter="\t"):
                                sample_rows.append({
                                    "study_id": rec["study_id"], "dataset_id": t["dataset_id"],
                                    "source_file": t["source_file"],
                                    "submission_accession": ";".join(a for _, a in repo_accs),
                                    "sample_id": row.get("sample_id", ""),
                                    "sample_flag": "blank_or_control" if row.get("sample_id", "") in
                                                  set(t.get("blank_or_control_sample_ids", "").split(";")) else "",
                                    "sequencing_type_hint": t.get("sequencing_type_hints", ""),
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

    print(f"studies: {len(studies_rows)}   tables: {len(table_rows)}   samples: {len(sample_rows)}   "
          f"needs_review: {len(review_rows)}   blocked entries: {len(blocked_rows)}")
    n1 = append_tsv(Path(args.studies_tsv), STUDIES_FIELDS, studies_rows, key_fields=("study_id",))
    n2 = append_tsv(Path(args.sample_tsv), SAMPLE_FIELDS, sample_rows,
                    key_fields=("study_id", "dataset_id", "source_file", "sample_id"))
    n3 = append_tsv(Path(args.blocked_tsv), BLOCKED_FIELDS, blocked_rows,
                    key_fields=("paper_id", "dataset_id", "attempted_files"))
    n4 = append_tsv(Path(args.table_manifest_tsv), TABLE_MANIFEST_FIELDS, table_rows,
                    key_fields=("study_id", "dataset_id", "source_file"))
    n5 = append_tsv(Path(args.needs_review_tsv), NEEDS_REVIEW_FIELDS, review_rows,
                    key_fields=("paper_id", "dataset_id", "source_file"))
    print(f"wrote {n1} new study rows, {n2} new sample rows, {n3} new blocked entries, "
          f"{n4} new table-manifest rows, {n5} new needs-review rows "
          f"(rows already present from an earlier build run were skipped)")
    # blocked_rows mixes two DIFFERENT situations that must not be reported
    # as one number: a paper with status=="blocked" never produced ANY
    # table at all, while a paper with status=="ok" but non-empty
    # blocked_files DID get at least one table -- just possibly missing
    # whatever those specific blocked files would have contributed.
    # Conflating them as "had every file blocked" is simply false for the
    # second group.
    fully_blocked_papers = {r["paper_id"] for r in blocked_rows
                            if status_by_paper.get(r["paper_id"]) == "blocked"}
    partial_papers = {r["paper_id"] for r in blocked_rows} - fully_blocked_papers
    print(f"\n{len(studies_rows)} studies got >=1 table written to disk ({len(table_rows)} tables total, "
          f"never merged across files/sheets).")
    print(f"{len(fully_blocked_papers)} abundance_ready papers had EVERY file blocked "
          f"(no table produced at all) -- see {args.blocked_tsv} for which ones and why "
          f"(most commonly anti_bot_challenge_page = manual browser download needed).")
    if partial_papers:
        print(f"{len(partial_papers)} more papers DID get at least one table, but had one or more "
              f"individual files blocked (possible missing data, not a total failure) -- "
              f"same {args.blocked_tsv}, look up these paper_ids to see which files/reasons.")
    if review_rows:
        print(f"{len(review_rows)} low-confidence tables need a human to open and check -- "
              f"see {args.needs_review_tsv}.")

def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    ep = sub.add_parser("extract", help="network phase: download + parse, write matrices + progress")
    ep.add_argument("--abundance-ready", required=True)
    ep.add_argument("--datasets", required=True)
    ep.add_argument("--data-dir", required=True)
    ep.add_argument("--progress", required=True)
    ep.add_argument("--scratch", required=True)
    ep.add_argument("--raw-store", required=True,
                    help="permanent, byte-for-byte retained copy of every downloaded file "
                         "(including anti-bot challenge pages), keyed by paper_id/dataset_id/filename. "
                         "Never deleted by this script; this is the audit trail behind every row.")
    ep.add_argument("--raw-manifest", default=None,
                    help="CSV recording source URL/repository/accession/filename/SHA256 for every "
                         "retained raw file. Defaults to <raw-store>/raw_files_manifest.csv.")
    ep.add_argument("--workers", type=int, default=8)
    ep.add_argument("--max-file-mb", type=float, default=200)
    ep.add_argument("--max-study-mb", type=float, default=500)
    ep.add_argument("--range", default=None,
                    help="Scope this run to a slice of the PENDING paper list (after skipping "
                         "already-done ones), by position -- not paper_id. A bare integer N processes "
                         "the first N pending papers (e.g. '10' = pending papers 1-10). 'A-B' processes "
                         "pending papers A through B inclusive (e.g. '20-30'). Omit to process everything "
                         "pending. Combine with --rebuild to force-reprocess a specific slice.")
    ep.add_argument("--paper-id", default=None,
                    help="Restrict this run to exactly one paper_id (e.g. for the Goodall_2026-style "
                         "single-paper verification pass). Applied before --range.")
    ep.add_argument("--rebuild", action="store_true",
                    help="Ignore --progress entirely when deciding what is 'already done' -- every "
                         "paper selected by --paper-id/--range is reprocessed from scratch. Results are "
                         "APPENDED to the same --progress file you pass (so point --progress at a NEW "
                         "path, e.g. progress_rebuild_<date>.jsonl, to keep the old progress.jsonl and "
                         "its old matrices/long tables completely untouched for side-by-side comparison; "
                         "pointing --rebuild at the OLD progress file would duplicate rows for any "
                         "paper_id appearing in both old and new runs instead of replacing them).")

    bp = sub.add_parser("build", help="local phase: assemble studies.tsv/sample.tsv/blocked list")
    bp.add_argument("--abundance-ready", required=True)
    bp.add_argument("--progress", required=True)
    bp.add_argument("--studies-tsv", required=True)
    bp.add_argument("--sample-tsv", required=True)
    bp.add_argument("--blocked-tsv", required=True)
    bp.add_argument("--table-manifest-tsv", required=True,
                    help="one row per (study, dataset, source_file) actually extracted -- the "
                         "authoritative per-table manifest (n_sample_ids_raw, domains, flags, SHA256).")
    bp.add_argument("--needs-review-tsv", required=True,
                    help="low-confidence tables (weak taxonomic signal) that were NOT auto-included "
                         "in any output -- queue for a human to open and judge.")

    args = ap.parse_args()
    if args.cmd == "extract":
        cmd_extract(args)
    else:
        cmd_build(args)


if __name__ == "__main__":
    main()
