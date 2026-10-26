#!/usr/bin/env python3
"""Incrementally widen europepmc_supp rows in datasets_master.csv WITHOUT
re-downloading the supplementary zip for every one of them.

Why this exists: commit 736a9e6 fixed europepmc_supp's file listing to
also scrape the live PMC article page's Associated Data links (a PMCID
can have multiple NIHMS manuscript submissions; Europe PMC's
supplementaryFiles zip mirrors only ONE of them, and can miss files a
later revision added -- see that commit for the full story). But
stage4's resume can't retroactively apply this to rows already written
by the old code -- deleting and letting stage4 re-resolve them (see
clear_europepmc_supp.py) would re-download every one of those zips
again, which is real cost at this scale (~48k papers, observed average
~16.5MB/zip from an earlier 58-paper sample -> roughly 800GB of repeat
transfer for data we already have on disk from the first pass).

This script does only the NEW, cheap part: for each europepmc_supp row,
fetch the live PMC page (one small HTML GET, not a multi-MB zip) and
union any newly-found filenames into the row's existing listed_files,
leaving the zip-derived names (and everything else about the row) alone.

Two-phase, independently resumable:
  1. fetch    -- network phase. Appends one JSON line per PMCID to a
                 progress file as it resolves. Safe to Ctrl-C and rerun:
                 already-recorded PMCIDs are skipped.
  2. merge    -- local, no network. Reads the progress file + the
                 original datasets_master.csv and writes an augmented
                 copy (backing up the original first). Rerun this alone
                 as many times as you like; it's a fast, pure merge.

Usage:
    # phase 1 (network, can take a while -- rerun to resume)
    python3 augment_europepmc_supp.py fetch \\
        --master /path/to/datasets_master.csv \\
        --progress /path/to/augment_progress.jsonl \\
        --workers 8

    # phase 2 (fast, local merge -- run once fetch looks complete, or
    # partway through if you want a partial merge now)
    python3 augment_europepmc_supp.py merge \\
        --master /path/to/datasets_master.csv \\
        --progress /path/to/augment_progress.jsonl

Then rerun stage5 to reclassify with the widened listings (stage5 makes
no network calls and is fast at this scale -- no --workers needed).
"""
import argparse
import csv
import json
import re
import shutil
import sys
import threading
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

csv.field_size_limit(sys.maxsize)  # listed_files can exceed the 131072-byte default

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _netutil import GLOBAL_THROTTLE  # noqa: E402

UA = {"User-Agent": "microbiome-dataset-resolver/1.0 (academic research)"}
PMC_BIN_RE = re.compile(r'/articles/instance/(\d+)/bin/([^"\']+)')

_print_lock = threading.Lock()


def list_pmc_associated_data(pmcid: str) -> tuple[list[str], str]:
    """(filenames, error). Same scrape as stage4/stage6's fix -- kept as
    an independent copy here so this script has no import-time coupling
    to stage4's other network calls (zenodo/figshare/etc.)."""
    url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"
    try:
        GLOBAL_THROTTLE.wait(url)
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=30) as r:
            html_text = r.read().decode("utf-8", "replace")
    except Exception as e:                          # noqa: BLE001
        return [], f"{type(e).__name__}: {e}"[:300]
    seen, names = set(), []
    for _num, fname in PMC_BIN_RE.findall(html_text):
        if fname not in seen:
            seen.add(fname); names.append(fname)
    return names, ""


def cmd_fetch(args: argparse.Namespace) -> None:
    with open(args.master, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    pmcids = sorted({r["accession"] for r in rows if r.get("repository") == "europepmc_supp" and r.get("accession")})
    print(f"{len(pmcids)} europepmc_supp PMCIDs in {args.master}")

    done = set()
    if Path(args.progress).exists():
        with open(args.progress, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    done.add(json.loads(line)["pmcid"])
                except Exception:                    # noqa: BLE001
                    pass
    pending = [p for p in pmcids if p not in done]
    print(f"{len(done)} already fetched, {len(pending)} pending")
    if not pending:
        print("nothing to do -- run the merge phase.")
        return

    out_fh = open(args.progress, "a", encoding="utf-8")
    write_lock = threading.Lock()
    n_done = [0]
    n_new_files = [0]
    t0 = time.time()

    def one(pmcid: str) -> None:
        names, err = list_pmc_associated_data(pmcid)
        rec = {"pmcid": pmcid, "names": names, "error": err, "ts": time.time()}
        with write_lock:
            out_fh.write(json.dumps(rec) + "\n")
            out_fh.flush()
            n_done[0] += 1
            if names:
                n_new_files[0] += len(names)
            if n_done[0] % 500 == 0:
                elapsed = time.time() - t0
                rate = n_done[0] / elapsed if elapsed > 0 else 0
                eta_min = (len(pending) - n_done[0]) / rate / 60 if rate > 0 else float("inf")
                with _print_lock:
                    print(f"  fetched {n_done[0]}/{len(pending)}  "
                          f"({rate:.2f}/s, ETA {eta_min:.0f} min)  "
                          f"files_found_so_far={n_new_files[0]}", flush=True)

    with ThreadPoolExecutor(args.workers) as ex:
        futs = [ex.submit(one, p) for p in pending]
        for fu in as_completed(futs):
            fu.result()  # surface any unexpected exception immediately
    out_fh.close()
    print(f"done. {n_done[0]} fetched this run, {n_new_files[0]} associated-data filenames found "
          f"(includes files already known from the zip -- merge phase dedupes).")


def cmd_merge(args: argparse.Namespace) -> None:
    scraped: dict[str, list[str]] = {}
    errors = 0
    if not Path(args.progress).exists():
        print(f"FATAL: {args.progress} does not exist -- run 'fetch' first", file=sys.stderr)
        sys.exit(1)
    with open(args.progress, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            if rec.get("error"):
                errors += 1
            scraped[rec["pmcid"]] = rec.get("names", [])
    print(f"{len(scraped)} PMCIDs in progress file ({errors} had a fetch error and contribute no new names)")

    with open(args.master, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    n_widened = 0
    n_files_added = 0
    for row in rows:
        if row.get("repository") != "europepmc_supp":
            continue
        extra = scraped.get(row.get("accession", ""))
        if not extra:
            continue
        existing = [x for x in row.get("listed_files", "").split(";") if x]
        merged = list(dict.fromkeys(existing + extra))
        added = len(merged) - len(existing)
        if added > 0:
            row["listed_files"] = ";".join(merged)
            row["n_listed_files"] = str(len(merged))
            page_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{row['accession']}/"
            if page_url not in (row.get("metadata_api_url") or ""):
                row["metadata_api_url"] = ";".join(v for v in (row.get("metadata_api_url"), page_url) if v)
            n_widened += 1
            n_files_added += added

    print(f"rows widened: {n_widened}  (total new filenames added: {n_files_added})")
    if n_widened == 0:
        print("nothing changed -- not writing a new file.")
        return

    master_path = Path(args.master)
    backup = master_path.with_name(f"{master_path.stem}.backup_{time.strftime('%Y%m%d_%H%M%S')}{master_path.suffix}")
    shutil.copy2(master_path, backup)
    print(f"backup written: {backup}")

    with master_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"rewrote {master_path} in place ({len(rows)} rows, {n_widened} widened).")
    print("\nNext: rerun stage5_classify_candidates.py to reclassify with the widened listings.")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    fp = sub.add_parser("fetch", help="network phase: scrape PMC pages, append results to --progress")
    fp.add_argument("--master", required=True, help="path to datasets_master.csv")
    fp.add_argument("--progress", required=True, help="JSONL file to append fetch results to (resumable)")
    fp.add_argument("--workers", type=int, default=8)

    mp = sub.add_parser("merge", help="local phase: merge --progress results into datasets_master.csv")
    mp.add_argument("--master", required=True, help="path to datasets_master.csv")
    mp.add_argument("--progress", required=True, help="JSONL file written by the fetch phase")

    args = ap.parse_args()
    if args.cmd == "fetch":
        cmd_fetch(args)
    else:
        cmd_merge(args)


if __name__ == "__main__":
    main()
