#!/usr/bin/env python3
"""Per-paper pipeline: download -> metadata -> abundance. One paper at a time.

Three steps, each runnable on its own so no block depends on the one before it
being run in the same invocation:

    1 download   reuses myproject/download_raw_reads.py
    2 metadata   run -> BioSample from ENA, then one unattended model session
                 per paper to finish the chain to the paper's own sample ids,
                 groups and tissue (step2_metadata.py + llm_extract.py)
    3 abundance  branch by assay type (step3_abundance.py)

State lives in pipeline_state.csv (one row per paper x step, atomic writes) and
is derived into pipeline_completed.csv / pipeline_failed.csv. A step that fails
is never recorded as done and never blocks other papers.

    python pipeline.py --paper Ambat2024                 # all three steps
    python pipeline.py --paper Ambat2024 --steps 2       # just metadata
    python pipeline.py --paper Ambat2024 --steps 3 --force   # redo abundance
    python pipeline.py --papers-file list.txt            # a list, in order
    python pipeline.py --paper X --steps 2 --dry-run     # show, do nothing
"""

import argparse
import csv
import json
import os
import subprocess
import socket
import sys
import time
from pathlib import Path

import registry

HERE = Path(__file__).resolve().parent
MYPROJECT = Path("/hickory/proj/didonglab/dataset/virus/yanting/microbiome_dataset/myproject")
STATUS_CSV = Path("/hickory/proj/didonglab/dataset/virus/Database/data/abudance_integration/"
                  "1_2_status_raw_reads_only.csv")
READS_ROOT = Path("/hickory/proj/didonglab/dataset/virus/Database/data/abundance_matrix_precise_v2")
OUT_ROOT = Path("/hickory/proj/didonglab/dataset/virus/Database/data/reads_abundance")
STATE = OUT_ROOT / "pipeline_state.csv"
COLUMNS = ["paper_folder", "doi", "title", "step", "status", "reason",
           "updated_utc", "result_path"]
STEPS = {1: "download", 2: "metadata", 3: "abundance"}
# Statuses that let the next step start. `done_quantification_only` does: the
# samples are sound even though their group labels are not confirmed.
DONE = {"done", "done_quantification_only", "done_partial"}


def now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def load_state():
    if not STATE.exists():
        return {}
    with STATE.open(encoding="utf-8-sig", newline="") as fh:
        return {(r["paper_folder"], r["step"]): r for r in csv.DictReader(fh)}


def save_state(state):
    """Write the audit trail, then rebuild every per-paper view from it.

    pipeline_state.csv keeps one row per (paper, step) because that is what the
    resume logic reads. Everything a human reads is one row per PAPER and is
    derived here, so the two can never disagree.
    """
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    rows = sorted(state.values(), key=lambda r: (r["paper_folder"].lower(), r["step"]))
    tmp = STATE.with_suffix(".csv.tmp")
    with tmp.open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, STATE)
    registry.rebuild(OUT_ROOT, rows)


def record(state, paper, step, status, reason="", result_path=""):
    meta = paper_row(paper)
    state[(paper, str(step))] = {
        "paper_folder": paper, "doi": meta.get("doi", ""),
        "title": (meta.get("title") or "")[:150],
        "step": str(step), "status": status, "reason": str(reason)[:500],
        "updated_utc": now(), "result_path": str(result_path),
    }
    save_state(state)


_TABLE = None


def paper_row(paper):
    global _TABLE
    if _TABLE is None:
        with STATUS_CSV.open(encoding="utf-8-sig", newline="") as fh:
            _TABLE = {r["paper_folder"].strip(): r for r in csv.DictReader(fh)}
    return _TABLE.get(paper, {})


def run(cmd, log):
    """Run a child step. Returns (ok, tail-of-output)."""
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as fh:
        p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, text=True)
    tail = log.read_text(encoding="utf-8", errors="replace").strip().splitlines()[-3:]
    return p.returncode == 0, " | ".join(tail)[:400]


def step_download(paper, args):
    cmd = [sys.executable, str(HERE / "step1_download.py"),
           "--paper", paper, "--csv", str(STATUS_CSV), "--save-dir", str(READS_ROOT)]
    ok, tail = run(cmd, OUT_ROOT / paper / "logs" / "step1_download.log")
    return ("done" if ok else "failed"), tail, READS_ROOT / paper / "raw_reads"


def step_metadata(paper, args):
    out = OUT_ROOT / paper
    cmd = [sys.executable, str(HERE / "step2_metadata.py"), "--paper", paper,
           "--reads-root", str(READS_ROOT), "--out-dir", str(out),
           "--status-csv", str(STATUS_CSV), "--agent", args.agent,
           "--llm-timeout", str(args.llm_timeout),
           "--llm-retries", str(args.llm_retries),
           "--max-turns", str(args.max_turns)]
    if args.model:
        cmd += ["--model", args.model]
    if args.refresh_metadata:
        cmd.append("--refresh-metadata")
    if args.no_llm:
        cmd.append("--no-llm")
    ok, tail = run(cmd, out / "logs" / "step2_metadata.log")
    if not ok:
        return "failed", tail, out
    gate = out / "metadata_gate.txt"
    verdict = gate.read_text(encoding="utf-8").strip() if gate.exists() else "blocked:no gate file"
    # Three outcomes, kept distinct: a paper that can only be quantified is not
    # the same as one that is fully resolved, and must not be filed as such.
    if verdict.startswith("pass_quantification_only"):
        status = "done_quantification_only"
    elif verdict.startswith("pass"):
        status = "done"
    else:
        status = "metadata_blocked"
    return status, verdict[:400], out


def step_abundance(paper, args):
    out = OUT_ROOT / paper
    cmd = [sys.executable, str(HERE / "step3_abundance.py"), "--paper", paper,
           "--reads-root", str(READS_ROOT), "--out-dir", str(out)]
    if args.silva:
        cmd += ["--silva", str(args.silva)]
    if args.silva_species:
        cmd += ["--silva-species", str(args.silva_species)]
    cmd += ["--quantify", args.quantify]
    for flag, val in (("--sylph-db", args.sylph_db), ("--sylph-tax", args.sylph_tax),
                      ("--viral-index", args.viral_index),
                      ("--viral-catalogue-name", args.viral_catalogue_name),
                      ("--viral-catalogue-scope", args.viral_catalogue_scope),
                      ("--host-index", args.host_index)):
        if val:
            cmd += [flag, str(val)]
    if args.archive_old:
        cmd.append("--archive-old")
    if args.preflight:
        cmd.append("--preflight")
    ok, tail = run(cmd, out / "logs" / "step3_abundance.log")
    if not ok:
        return "failed", tail, out
    # A preflight proves nothing about abundance: record it as its own state so
    # it can never satisfy "step 3 done".
    if args.preflight:
        return "preflight_ok", tail, out
    # Exit code 0 is not evidence that a matrix exists. A --probe-only run also
    # exits 0, and so does a run that stopped after writing only some outputs.
    # Which branches ran, and did each leave a real matrix behind? A paper with
    # a 16S branch and a viral branch is only done when both are.
    try:
        plan = json.loads((out / "analysis_plan.json").read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return "failed", f"no readable analysis_plan.json: {exc}", out
    expected = {
        "amplicon_16s": ("dada2", "asv_counts.csv",
                         ["asv_relative_abundance.csv", "asv_sequences.fasta",
                          "asv_taxonomy.csv", "sample_metadata.csv",
                          "read_counts_per_step.csv"]),
        "dna_shotgun_bacteria": ("shotgun", "bacteria_taxonomic_abundance_pct.csv",
                                 ["feature_metadata.csv", "qc_stats.csv", "README.txt"]),
        "viral_mapping": ("virome", "viral_mean_depth.csv",
                          ["viral_relative_abundance_within_virome.csv",
                           "feature_metadata.csv", "mapping_stats.csv", "README.txt"]),
    }
    ran = [b for b, st in plan.get("branch_status", {}).items() if st in ("done", "failed")]
    if not ran:
        return "failed", "exit 0 but no branch reported a result", out
    notes, bad = [], []
    for b in ran:
        if plan["branch_status"][b] == "failed":
            bad.append(f"{b}: reported failed")
            continue
        sub, matrix, others = expected[b]
        d = out / sub
        miss = [f for f in [matrix] + others if not (d / f).exists()]
        if miss:
            bad.append(f"{b}: missing {miss}")
            continue
        try:
            with (d / matrix).open(encoding="utf-8-sig", newline="") as fh:
                header = next(csv.reader(fh))
                n = sum(1 for _ in fh)
        except (OSError, StopIteration) as exc:
            bad.append(f"{b}: {matrix} unreadable ({exc})")
            continue
        if n == 0 or len(header) < 2:
            bad.append(f"{b}: {matrix} has {n} rows x {max(len(header) - 1, 0)} features")
            continue
        notes.append(f"{b}={n}x{len(header) - 1}")
    if bad:
        return "failed", "; ".join(bad)[:400], out
    # A branch that never ran is not a failure, but the paper is not finished
    # either. Give it its own status so "done" keeps meaning "everything this
    # paper has was produced".
    deferred = {b: st for b, st in plan.get("branch_status", {}).items()
                if st in ("missing_dependencies", "unsupported")}
    if deferred:
        return ("done_partial",
                " | ".join(notes) + "; NOT RUN: " +
                ", ".join(f"{b}({st})" for b, st in sorted(deferred.items())), out)
    return "done", " | ".join(notes) + f" | {tail}", out


RUNNERS = {1: step_download, 2: step_metadata, 3: step_abundance}


def lock_id():
    """Who this process is, in terms another job can check from another node."""
    return json.dumps({"host": socket.gethostname(), "pid": os.getpid(),
                       "job": os.environ.get("SLURM_JOB_ID", ""), "at": now()})


def holder_is_gone(raw):
    """True when the process named in a lock file is provably no longer running.

    A lock is only ever stolen on PROOF of death, never on age: a 16S paper can
    legitimately sit in dada2 for many hours, and stealing from a live run would
    have two jobs writing the same output directory. Unparseable or unprovable
    means held.
    """
    try:
        h = json.loads(raw)
    except Exception:
        return False                      # old plain-text lock: leave it alone
    job = h.get("job") or ""
    if job:
        # squeue is authoritative and works from any node. Only its own clean
        # "not in queue" answer counts; a scheduler hiccup must not free a lock.
        try:
            r = subprocess.run(["squeue", "-h", "-j", job, "-o", "%T"],
                               capture_output=True, text=True, timeout=60)
        except Exception:
            return False
        if r.returncode != 0:
            # slurm says the job id is unknown -> it has left the queue.
            return "Invalid job id" in (r.stderr or "")
        return not r.stdout.strip()
    if h.get("host") != socket.gethostname():
        return False                      # another node's pid means nothing here
    try:
        os.kill(int(h["pid"]), 0)
    except ProcessLookupError:
        return True
    except Exception:
        return False
    return False


def process(paper, args, state):
    lock = OUT_ROOT / paper / ".lock"
    lock.parent.mkdir(parents=True, exist_ok=True)
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.write(fd, (lock_id() + "\n").encode())
        os.close(fd)
    except FileExistsError:
        try:
            held = lock.read_text().strip()
            age = int(time.time() - lock.stat().st_mtime)
        except OSError:
            held, age = "", -1
        if holder_is_gone(held):
            # A job killed by the wall clock leaves its lock behind. Taking it
            # over is the whole point: otherwise one timeout strands every
            # paper that job held, silently, for every later submission.
            print(f"  {paper}: taking over a lock left by a job that is gone ({held})")
            try:
                lock.write_text(lock_id() + "\n", encoding="utf-8")
            except OSError as e:
                print(f"  {paper}: could not take the lock over ({e}); skipped")
                return
        else:
            print(f"  {paper}: locked by '{held}' for {age}s ({lock}). If that "
                  f"process is gone, delete the file to let this paper run.")
            return
    try:
        for step in args.steps:
            prev = state.get((paper, str(step)), {})
            redo = args.force or (step == 2 and args.refresh_metadata)
            if prev.get("status") in DONE and not redo:
                print(f"  [{step} {STEPS[step]}] already {prev['status']}, skipped")
                continue
            if args.dry_run:
                print(f"  [{step} {STEPS[step]}] would run")
                continue
            record(state, paper, step, "running")
            status, reason, path = RUNNERS[step](paper, args)
            record(state, paper, step, status, reason, path)
            print(f"  [{step} {STEPS[step]}] {status}: {reason[:160]}")
            if status not in DONE:
                break          # this paper stops; others continue
    finally:
        lock.unlink(missing_ok=True)
        # Hand the paper back as soon as ITS run ends, not when the whole job
        # ends: the next submission must see it as finished rather than as
        # still being worked on.
        if args.release_claims:
            subprocess.run([sys.executable, str(HERE / "claim_papers.py"),
                            "--release", paper, "--out-root", str(OUT_ROOT)],
                           capture_output=True, text=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--paper", action="append")
    g.add_argument("--papers-file", type=Path)
    ap.add_argument("--steps", default="1,2,3",
                    help="comma list of 1=download 2=metadata 3=abundance")
    ap.add_argument("--force", action="store_true", help="redo steps already marked done")
    ap.add_argument("--preflight", action="store_true",
                    help="step 3 only checks tools/inputs and estimates resources")
    ap.add_argument("--refresh-metadata", action="store_true",
                    help="redo step 2's model session even if a validated answer is cached")
    ap.add_argument("--no-llm", action="store_true",
                    help="step 2 uses ENA only (cannot reach a full pass)")
    ap.add_argument("--agent", default="claude", choices=("claude", "codex"))
    ap.add_argument("--model", default="")
    ap.add_argument("--llm-timeout", type=int, default=2400)
    ap.add_argument("--llm-retries", type=int, default=1)
    ap.add_argument("--max-turns", type=int, default=120)
    ap.add_argument("--silva", type=Path, default=None)
    ap.add_argument("--silva-species", type=Path, default=None)
    ap.add_argument("--quantify", default="bacteria",
                    choices=("bacteria", "viruses", "both"))
    ap.add_argument("--sylph-db", default="")
    ap.add_argument("--sylph-tax", default="")
    ap.add_argument("--viral-index", default="")
    ap.add_argument("--viral-catalogue-name", default="unnamed")
    ap.add_argument("--viral-catalogue-scope", default="")
    ap.add_argument("--host-index", default="")
    ap.add_argument("--archive-old", action="store_true",
                    help="move previous step-3 outputs into _archive/<utc>/")
    ap.add_argument("--release-claims", action="store_true",
                    help="remove each paper from its claims/ file when its run ends")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    args.steps = [int(s) for s in str(args.steps).split(",") if s.strip()]

    papers = args.paper or [l.strip() for l in args.papers_file.read_text().splitlines()
                            if l.strip() and not l.startswith("#")]
    state = load_state()
    for i, paper in enumerate(papers, 1):
        print(f"[{i}/{len(papers)}] {paper}", flush=True)
        process(paper, args, state)
    print(f"\nstate -> {STATE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
