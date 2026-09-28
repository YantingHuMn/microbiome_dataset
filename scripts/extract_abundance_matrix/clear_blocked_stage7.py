#!/usr/bin/env python3
"""Remove status=="blocked" entries from stage7's progress.jsonl so a
resumed `extract` run re-checks just those papers with the current code
(commit e06f939's anti-bot detection fix -- see that commit for why
"blocked" entries recorded before it can carry an uninformative, blank
reason instead of the real anti_bot_challenge_page cause).

Leaves status=="ok" entries untouched -- those papers already have a
real, downloaded abundance matrix on disk; this script does not force
them to be re-checked (a paper with multiple candidate files where only
SOME hit the bug would already have real data from whichever file(s)
parsed successfully, so re-checking is a completeness nice-to-have, not
a correctness fix, for that subset).

Usage:
    python3 clear_blocked_stage7.py /path/to/stage7_extract/progress.jsonl
"""
import argparse
import json
import shutil
import sys
import time
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("progress_path")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    path = Path(args.progress_path)
    if not path.exists():
        print(f"FATAL: {path} does not exist", file=sys.stderr)
        sys.exit(1)

    kept, removed = [], []
    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            (removed if rec.get("status") == "blocked" else kept).append(line)

    print(f"total records:     {len(kept) + len(removed)}")
    print(f"status=ok kept:    {len(kept)}  (untouched, will be skipped on resume)")
    print(f"status=blocked:    {len(removed)}  (will be re-checked on next extract run)")

    if not removed:
        print("nothing to clear.")
        return
    if args.dry_run:
        print("\n--dry-run: no file was modified.")
        return

    backup = path.with_name(f"{path.stem}.backup_{time.strftime('%Y%m%d_%H%M%S')}{path.suffix}")
    shutil.copy2(path, backup)
    print(f"\nbackup written: {backup}")

    with path.open("w", encoding="utf-8") as f:
        for line in kept:
            f.write(line + "\n")
    print(f"rewrote {path} with {len(kept)} records ({len(removed)} blocked records removed)")
    print("\nNext: resubmit the extract phase against the SAME --progress/--data-dir. "
          "Only the removed (previously blocked) papers will be reprocessed.")


if __name__ == "__main__":
    main()
