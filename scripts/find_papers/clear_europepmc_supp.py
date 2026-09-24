#!/usr/bin/env python3
"""Remove europepmc_supp rows from datasets_master.csv so a resumed stage4
run reprocesses them with the current code (commit 736a9e6 and later --
the fix that widens europepmc_supp's file listing to also pick up files
from a newer PMC manuscript revision, not just whatever Europe PMC's
supplementaryFiles zip happens to mirror).

Everything else in datasets_master.csv (zenodo/figshare/dryad/osf/mgnify/
bioproject/... rows) is left untouched -- stage4's resume logic will see
those dataset_ids are still present and skip them, so only the
europepmc_supp subset gets redone.

Usage:
    python3 clear_europepmc_supp.py /path/to/stage4_datasets/datasets_master.csv

Writes a timestamped backup next to the original before overwriting it.
Safe to run more than once (it's idempotent -- if there's nothing left to
clear, it says so and exits without touching the file).
"""
import argparse
import csv
import shutil
import sys
import time
from pathlib import Path

csv.field_size_limit(sys.maxsize)  # listed_files can exceed the 131072-byte default


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("csv_path", help="path to datasets_master.csv")
    ap.add_argument("--dry-run", action="store_true",
                    help="just report counts, don't modify the file")
    args = ap.parse_args()

    path = Path(args.csv_path)
    if not path.exists():
        print(f"FATAL: {path} does not exist", file=sys.stderr)
        sys.exit(1)

    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    removed = [r for r in rows if r.get("repository") == "europepmc_supp"]
    kept = [r for r in rows if r.get("repository") != "europepmc_supp"]

    print(f"total rows:        {len(rows)}")
    print(f"europepmc_supp:    {len(removed)}  (will be re-resolved on next stage4 run)")
    print(f"everything else:   {len(kept)}  (untouched, stage4 will skip these on resume)")

    if not removed:
        print("nothing to clear -- already clean, or stage4 hasn't run against this file yet.")
        return

    if args.dry_run:
        print("\n--dry-run: no file was modified.")
        return

    backup = path.with_name(f"{path.stem}.backup_{time.strftime('%Y%m%d_%H%M%S')}{path.suffix}")
    shutil.copy2(path, backup)
    print(f"\nbackup written: {backup}")

    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(kept)
    print(f"rewrote {path} with {len(kept)} rows ({len(removed)} europepmc_supp rows removed)")
    print("\nNext: resubmit stage4 against the SAME --outdir. Resume will treat only")
    print("the removed europepmc_supp entries as pending -- everything else is skipped.")


if __name__ == "__main__":
    main()
