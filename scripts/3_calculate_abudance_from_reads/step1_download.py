#!/usr/bin/env python3
"""Download raw FASTQ reads for every paper in 1_2_status_raw_reads_only.csv.
"""

import argparse
import sys
from pathlib import Path

MYPROJECT = Path(
    "/hickory/proj/didonglab/dataset/virus/yanting/microbiome_dataset/myproject"
)
sys.path.insert(0, str(MYPROJECT))
from download_raw_reads import (  # noqa: E402
    COMPLETED_CSV, FAILED_CSV, Progress, install_signal_handlers,
    process_paper, read_status_table,
)

STATUS_CSV = Path(
    "/hickory/proj/didonglab/dataset/virus/Database/data/abudance_integration/"
    "1_2_status_raw_reads_only.csv"
)
SAVE_DIR = Path(
    "/hickory/proj/didonglab/dataset/virus/Database/data/abundance_matrix_precise_v2"
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", type=Path, default=STATUS_CSV)
    ap.add_argument("--save-dir", type=Path, default=SAVE_DIR)
    ap.add_argument("--paper", action="append", help="only this paper_folder (repeatable)")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--no-checksum", action="store_true",
                    help="verify size only; recorded as such, not as MD5 verified")
    ap.add_argument("--retry-failed", action="store_true",
                    help="also retry papers recorded as failed or partial")
    args = ap.parse_args()

    rows = read_status_table(args.csv)
    papers = rows
    if args.paper:
        papers = [p for p in papers if p["paper_folder"].strip() in set(args.paper)]
    if args.limit:
        papers = papers[: args.limit]

    progress = Progress(args.save_dir, rows)
    install_signal_handlers(progress)
    print(f"{len(papers)} papers -> {args.save_dir}")
    print(f"already complete: {len(progress.done)}, unfinished: {len(progress.bad)}")

    counts: dict[str, int] = {}
    statuses: list[str] = []
    for i, paper in enumerate(papers, 1):
        print(f"[{i}/{len(papers)}] {paper['paper_folder'].strip()}", flush=True)
        status, summary = process_paper(
            paper, args.save_dir, progress, args.csv,
            workers=args.workers, checksum=not args.no_checksum,
            retry_failed=args.retry_failed)
        counts[status] = counts.get(status, 0) + 1
        statuses.append(status)
        print(f"    {status}: {summary}", flush=True)

    print("\n" + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    print(f"completed -> {args.save_dir / COMPLETED_CSV} ({len(progress.done)})")
    print(f"failed    -> {args.save_dir / FAILED_CSV} ({len(progress.bad)})")
    # The exit code is about THIS run's papers only. progress.bad counts every
    # unfinished paper in the whole table, so using it would fail a successful
    # single-paper run just because other papers are still outstanding.
    ok = {"complete", "skipped_complete"}
    bad_here = [p for p, s in zip(papers, statuses) if s not in ok]
    if bad_here:
        print(f"not downloaded in this run: "
              f"{', '.join(p['paper_folder'].strip() for p in bad_here)}")
    return 1 if bad_here else 0


if __name__ == "__main__":
    raise SystemExit(main())
