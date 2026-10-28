#!/usr/bin/env python3
"""Claim the next N unprocessed papers, so several jobs can run side by side.

Each submitted job takes a disjoint slice of the work. The pipeline itself still
runs one paper at a time inside a job - nothing is parallelised within a job -
but several people can submit at once and never collide.

How a paper is chosen: everything in the source table, minus
    - papers another job has already claimed and not yet finished
    - papers in pipeline_completed.csv
    - papers in pipeline_failed.csv            (fix them deliberately, not by accident)
    - papers in method_missing.csv             (no pipeline exists for their assay)
    - papers in any other csv named with --exclude

Risk files are NOT excluded: a paper in risk_no_grouping.csv finished
successfully and should not be handed out again.

The claim itself is a file under claims/. Creating it is atomic (O_EXCL), so two
jobs racing for the same slice cannot both win. A paper is removed from its
claim file as soon as its own run ends, which is what lets the next submission
see it as finished rather than as still-running.

    claim_papers.py --n 10 --out papers_job123.txt      # take 10, write a list
    claim_papers.py --release Ambat2024 --claim-file claims/job123.csv
    claim_papers.py --status                            # who holds what
    claim_papers.py --release-stale 24                  # free claims older than 24h
"""

import argparse
import csv
import os
import subprocess
import time
from pathlib import Path

OUT_ROOT = Path("/hickory/proj/didonglab/dataset/virus/Database/data/reads_abundance")
SOURCE = Path("/hickory/proj/didonglab/dataset/virus/Database/data/abudance_integration/"
              "1_2_status_raw_reads_only.csv")
CLAIMS = OUT_ROOT / "claims"
#: Finishing or failing takes a paper out of circulation; a risk does not.
# These are the registry's level-1 files. Their names must track registry.py:
# a name that no longer exists would read as "nothing to exclude" and the next
# claim would hand out papers that are already finished.
EXCLUDE_DEFAULT = ["1_completed.csv", "1_pipeline_failed.csv",
                   "1_method_missing.csv"]


def folders(path, column="paper_folder"):
    if not Path(path).exists():
        return set()
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        return {(r.get(column) or "").strip() for r in csv.DictReader(fh)
                if (r.get(column) or "").strip()}


def job_is_gone(job_id):
    """True only when slurm itself says that job id is no longer in the queue.

    Claim files are named after the job that took them, so a job killed by the
    wall clock can be told apart from one still working. Only slurm's own clean
    "Invalid job id" counts: a scheduler hiccup, a non-numeric name, or a local
    run must leave the claim standing, because freeing a live job's papers lets
    a second job write the same output directory.
    """
    if not job_id.isdigit():
        return False
    try:
        r = subprocess.run(["squeue", "-h", "-j", job_id, "-o", "%T"],
                           capture_output=True, text=True, timeout=60)
    except Exception:
        return False
    if r.returncode != 0:
        return "Invalid job id" in (r.stderr or "")
    return not r.stdout.strip()


def claimed_now(auto_free=True):
    """paper -> claim file currently holding it.

    Claims whose job has left the queue are dropped on sight, so one timeout
    does not strand its papers for every later submission.
    """
    held = {}
    if not CLAIMS.exists():
        return held
    for f in sorted(CLAIMS.glob("*.csv")):
        if auto_free and job_is_gone(f.stem.split(".")[0]):
            print(f"  freeing {f.name}: slurm job {f.stem.split('.')[0]} has "
                  f"left the queue, so its claim cannot still be in use")
            f.unlink(missing_ok=True)
            continue
        for p in folders(f):
            held[p] = f
    return held


def write_claim(path, papers, owner, note):
    """Create the claim file atomically; refuse if another job got there first."""
    fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    with os.fdopen(fd, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["paper_folder", "claimed_by",
                                           "claimed_utc", "note"])
        w.writeheader()
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        for p in papers:
            w.writerow({"paper_folder": p, "claimed_by": owner,
                        "claimed_utc": now, "note": note})


def release(paper, claim_file):
    """Drop one paper from a claim file. Removes the file when it empties."""
    claim_file = Path(claim_file)
    if not claim_file.exists():
        return False
    with claim_file.open(encoding="utf-8-sig", newline="") as fh:
        rows = [r for r in csv.DictReader(fh)]
    keep = [r for r in rows if (r.get("paper_folder") or "").strip() != paper]
    if len(keep) == len(rows):
        return False
    if not keep:
        claim_file.unlink()
        return True
    tmp = claim_file.with_suffix(".csv.tmp")
    with tmp.open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(keep)
    os.replace(tmp, claim_file)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, help="how many papers to claim")
    ap.add_argument("--out", type=Path, help="write the claimed list here, one per line")
    ap.add_argument("--source", type=Path, default=SOURCE)
    ap.add_argument("--out-root", type=Path, default=OUT_ROOT)
    ap.add_argument("--owner", default=os.environ.get("USER", "unknown"))
    ap.add_argument("--job", default=os.environ.get("SLURM_JOB_ID", f"local{os.getpid()}"))
    ap.add_argument("--note", default="")
    ap.add_argument("--exclude", action="append", default=[],
                    help="extra csv of papers to skip (repeatable)")
    ap.add_argument("--release", help="mark this paper as no longer being worked on")
    ap.add_argument("--claim-file", type=Path, help="which claim file --release touches")
    ap.add_argument("--release-stale", type=float, metavar="HOURS",
                    help="free claims older than this, e.g. after a job was killed")
    ap.add_argument("--status", action="store_true")
    args = ap.parse_args()

    global CLAIMS
    CLAIMS = args.out_root / "claims"
    CLAIMS.mkdir(parents=True, exist_ok=True)

    if args.release:
        if not args.claim_file:
            cf = claimed_now(auto_free=False).get(args.release)
            if not cf:
                print(f"{args.release} is not claimed by anyone")
                return 0
            args.claim_file = cf
        print(("released " if release(args.release, args.claim_file) else
               "not in ") + f"{args.claim_file}: {args.release}")
        return 0

    if args.release_stale is not None:
        cutoff = time.time() - args.release_stale * 3600
        freed = 0
        for f in sorted(CLAIMS.glob("*.csv")):
            if f.stat().st_mtime < cutoff:
                print(f"freeing stale claim {f.name} "
                      f"({(time.time()-f.stat().st_mtime)/3600:.1f}h old): "
                      f"{sorted(folders(f))}")
                f.unlink(); freed += 1
        print(f"{freed} stale claim file(s) freed")
        return 0

    held = claimed_now()
    if args.status:
        print(f"claims dir: {CLAIMS}")
        for f in sorted(CLAIMS.glob("*.csv")):
            age = (time.time() - f.stat().st_mtime) / 3600
            print(f"  {f.name}  {age:5.1f}h  {sorted(folders(f))}")
        print(f"  {len(held)} paper(s) claimed in total")
        return 0

    if not args.n:
        ap.error("give --n, or use --status / --release / --release-stale")

    everything = [r for r in csv.DictReader(
        args.source.open(encoding="utf-8-sig", newline=""))]
    all_papers = [(r.get("paper_folder") or "").strip() for r in everything]

    skip, why = set(), {}
    missing = [n for n in EXCLUDE_DEFAULT
               if not (Path(n) if Path(n).is_absolute()
                       else args.out_root / n).exists()]
    if missing:
        raise SystemExit(
            "cannot claim safely: these status files are missing from "
            f"{args.out_root}: {', '.join(missing)}.\nWithout them a claim "
            "cannot tell finished papers from unstarted ones and would hand "
            "out work that is already done. Run `python registry.py` to "
            "rebuild them from pipeline_state.csv, then claim again.")
    for name in EXCLUDE_DEFAULT + args.exclude:
        p = Path(name) if Path(name).is_absolute() else args.out_root / name
        got = folders(p)
        for g in got:
            why.setdefault(g, p.name)
        skip |= got
    for p in held:
        why.setdefault(p, f"claimed by {held[p].name}")
    skip |= set(held)

    free = [p for p in all_papers if p and p not in skip]
    take = free[:args.n]
    print(f"source {args.source.name}: {len(all_papers)} papers")
    for name in EXCLUDE_DEFAULT + args.exclude:
        p = Path(name) if Path(name).is_absolute() else args.out_root / name
        print(f"  minus {len(folders(p)):4d} in {p.name}")
    print(f"  minus {len(held):4d} claimed by a running job")
    print(f"  -> {len(free)} free, claiming {len(take)}")
    if not take:
        print("nothing left to claim")
        if args.out:
            args.out.write_text("")
        return 0

    claim = CLAIMS / f"{args.job}.csv"
    try:
        write_claim(claim, take, args.owner, args.note)
    except FileExistsError:
        raise SystemExit(f"{claim} already exists - this job id has already "
                         "claimed work. Use a different --job or release it first.")
    print(f"claimed in {claim}:")
    for p in take:
        print(f"  {p}")
    if args.out:
        args.out.write_text("\n".join(take) + "\n", encoding="utf-8")
        print(f"list -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
