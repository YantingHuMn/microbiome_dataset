#!/usr/bin/env python3
"""Measure a paper's amplicon facts from its reads and write paper_facts/<paper>.json.

The first few papers had this file written by hand: read the Methods, find the
primer pair, work out the amplicon length, check it against the reads. That does
not scale to a hundred papers, and most of it was measurement anyway. This does
the measuring.

What it MEASURES from the reads (strong, direct evidence):
    which primer pair is on the reads, and in what fraction of them
    whether the primers are still there at all (submitters often trim them)
    the primer-free insert length distribution, from R1 vs revcomp(R2) overlap
    3' read-through, which caps truncLen

What it does NOT do: read the paper. A measured primer match of 98% is direct
evidence of which primer was used, but it is not the same as the Methods saying
so. Each section records its provenance, and `amplicon.status` is
`measured_from_reads`, never `verified`, unless a human confirms it against the
paper and says so in the file.

    python make_paper_facts.py --paper X --reads-root ... --out-dir ...
    python make_paper_facts.py --paper X ... --overwrite      # redo
"""

import argparse
import collections
import csv
import gzip
import json
import statistics
import sys
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from step3_abundance import IUPAC, PRIMERS, detect_region, matches  # noqa: E402

COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N",
        "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
        "B": "V", "V": "B", "D": "H", "H": "D"}


def rc(s):
    return "".join(COMP.get(c, "N") for c in reversed(s))


def read_pairs(r1, r2, n, cut1, cut2, skip1, skip2):
    """Yield (R1 prefix, revcomp(R2) prefix) for up to n pairs."""
    with gzip.open(r1, "rt") as a, gzip.open(r2, "rt") as b:
        for i, (x, y) in enumerate(zip(a, b)):
            if i % 4 == 1:
                yield x.strip()[skip1:cut1], rc(y.strip()[skip2:cut2])
                n -= 1
                if n <= 0:
                    return


def overlap_length(a, b, min_ov=20, min_id=0.80):
    """Insert length implied by the best overlap of a and rc(b), else None."""
    best = (0.0, None)
    for ov in range(min(len(a), len(b)), min_ov - 1, -1):
        x, y = a[-ov:], b[:ov]
        ident = sum(p == q for p, q in zip(x, y)) / ov
        if ident > best[0]:
            best = (ident, len(a) + len(b) - ov)
        if ident > 0.97 and ov > 30:
            break
    return best[1] if best[0] >= min_id else None


def find(seq, pat, allow=1):
    for i in range(len(seq) - len(pat) + 1):
        if matches(seq[i:], pat, allow):
            return i
    return -1


def probe_runs(out_dir, n):
    """Representative runs spread across the sample list, from step 2's map."""
    with (out_dir / "run_sample_mapping.csv").open(encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))
    by_run = collections.defaultdict(dict)
    for r in rows:
        by_run[r["run_accession"]][r["mate"]] = r["fastq_path"]
    runs = sorted(by_run)
    if not runs:
        raise SystemExit("run_sample_mapping.csv has no runs; run step 2 first")
    idx = sorted({round(i * (len(runs) - 1) / max(n - 1, 1)) for i in range(min(n, len(runs)))})
    return [(runs[i], by_run[runs[i]]) for i in idx], len(runs)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paper", required=True)
    ap.add_argument("--reads-root", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--probe-runs", type=int, default=6)
    ap.add_argument("--probe-pairs", type=int, default=800,
                    help="read pairs per run used for the length measurement")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    dest = HERE / "paper_facts" / f"{args.paper}.json"
    if dest.exists() and not args.overwrite:
        raise SystemExit(f"{dest} exists; pass --overwrite to replace it. A file "
                         "written or corrected by hand is not overwritten by "
                         "accident.")

    probes, n_runs = probe_runs(args.out_dir, args.probe_runs)
    print(f"{args.paper}: probing {len(probes)} of {n_runs} runs")

    # ---- 1. which primers, and are they still on the reads -----------------
    calls = {}
    for run, mates in probes:
        r1, r2 = mates.get("R1") or mates.get("single"), mates.get("R2")
        calls[run] = detect_region(r1, r2) if r1 else ("unknown", False, "no R1")
    regions = {v[0] for v in calls.values()}
    if len(regions) > 1:
        raise SystemExit(f"the probed runs disagree on the region: {sorted(regions)}. "
                         "This paper may mix primer sets; resolve it by hand.")
    region, primers_present, evidence = calls[probes[0][0]]
    print(f"  region {region}, primers on reads: {primers_present}")

    base = region[len("primer_like_"):] if region.startswith("primer_like_") else region
    fwd, rev = PRIMERS.get(base, ("", ""))
    trim_left = [len(fwd), len(rev)] if primers_present and fwd else [0, 0]

    # ---- 2. insert length, measured ----------------------------------------
    # Read length is measured on UNTRUNCATED reads. Taking it from the probe
    # strings below would report the probe's own cut (290/280) as the read
    # length and quietly throw away the last 20-30 usable bases of every read.
    read_len = [0, 0]
    for run, mates in probes:
        for i, key in ((0, "R1"), (1, "R2")):
            path = mates.get(key)
            if not path:
                continue
            with gzip.open(path, "rt") as fh:
                for j, line in enumerate(fh):
                    if j % 4 == 1:
                        read_len[i] = max(read_len[i], len(line.strip()))
                    if j >= 4000:
                        break

    lens, through = [], 0
    tried = 0
    for run, mates in probes:
        r1, r2 = mates.get("R1"), mates.get("R2")
        if not (r1 and r2):
            continue
        for a, b in read_pairs(r1, r2, args.probe_pairs, 290, 280,
                               trim_left[0], trim_left[1]):
            tried += 1
            L = overlap_length(a, b)
            if L:
                lens.append(L)
            if rev and find(a, rc(rev)) >= 0:
                through += 1
    if not lens:
        raise SystemExit("could not measure an insert length from any read pair; "
                         "this paper needs a hand-written facts file")
    q = statistics.quantiles(lens, n=100)
    p99 = round(q[98])
    frac_through = through / max(tried, 1)
    print(f"  insert length: median {statistics.median(lens):.0f}, p99 {p99}, "
          f"max {max(lens)}  ({len(lens)}/{tried} pairs overlapped)")
    print(f"  3' read-through into the reverse primer: {frac_through:.1%}")

    # Read-through means the read runs off the far end into adaptor, so truncLen
    # must stop at the insert; otherwise the read length is the only cap.
    max_trunc = ([p99 + trim_left[0], p99 + trim_left[1]] if frac_through > 0.05
                 else list(read_len))
    max_trunc = [min(max_trunc[i], read_len[i]) for i in (0, 1)]
    need = p99 + 12 + sum(trim_left)
    if sum(max_trunc) < need:
        print(f"  WARNING reads cannot span this amplicon: {sum(max_trunc)} < {need}")

    # ---- 3. sample identity, from step 2 -----------------------------------
    gate = (args.out_dir / "metadata_gate.txt")
    verdict = gate.read_text(encoding="utf-8").strip() if gate.exists() else "unknown"
    with (args.out_dir / "sample_metadata.csv").open(encoding="utf-8-sig", newline="") as fh:
        samples = list(csv.DictReader(fh))
    ctrl = [s for s in samples if s.get("is_technical_control") == "True"]
    bio = [s for s in samples if s.get("is_technical_control") != "True"]
    conf = lambda f: sum(1 for s in bio if s.get(f) not in ("unresolved", "not_reported", ""))

    facts = {
        "paper_folder": args.paper,
        "doi": samples[0].get("doi", "") if samples else "",
        "generated_by": "make_paper_facts.py (measured from the reads)",
        "generated_on": date.today().isoformat(),
        "sources": [
            "the reads themselves, via step 2's run_sample_mapping.csv",
            "step 2's metadata extraction (sample identity section)",
        ],
        "amplicon": {
            # NOT "verified": nobody has read the Methods. Measured is strong
            # evidence of which primer was used, but it is a different claim.
            "status": "measured_from_reads",
            "region": base,
            "paper_statement": "NOT CHECKED - no one has read this paper's Methods. "
                               "Confirm the primer pair and set status to 'verified' "
                               "before relying on this for anything but denoising.",
            "forward_primer": {"name": base + "_fwd", "seq": fwd, "length": len(fwd)},
            "reverse_primer": {"name": base + "_rev", "seq": rev, "length": len(rev)},
            "primers_in_reads": bool(primers_present),
            "primer_evidence": f"probed {len(probes)} runs: {evidence}",
            "read_length": read_len,
            "trim_left": trim_left,
            "insert_length_primer_free": {
                "method": f"R1 vs revcomp(R2) overlap over {len(lens)} pairs from "
                          f"{len(probes)} runs (>=80% identity in the overlap)",
                "min": min(lens), "p1": round(q[0]), "p5": round(q[4]),
                # p95/p99/max drive the merge budget and are stable across
                # probe sets; the median moves with which runs are sampled on a
                # multi-modal distribution, so treat it as descriptive only.
                "median": round(statistics.median(lens)), "p95": round(q[94]),
                "p99": p99, "max": max(lens),
                "modes": [L for L, _ in collections.Counter(lens).most_common(5)],
            },
            "max_trunc_len": max_trunc,
            "max_trunc_len_reason": (
                f"{frac_through:.0%} of reads run through into the reverse primer, so "
                "truncLen is capped at the insert length" if frac_through > 0.05
                else "read length; no 3' read-through to avoid"),
            "trunc_len_basis": (
                "RAW read bases. dada2::fastqPairedFilter removes trimLeft first and "
                "then truncates to truncLen - trimLeft, keeping only reads at least "
                "truncLen long, so a read contributes bases trimLeft+1 .. truncLen."),
            "merge_budget": {
                "min_overlap": 12,
                "required_raw_sum": need,
                "formula": f"truncLen[1]+truncLen[2] >= {p99} + 12 + "
                           f"{trim_left[0]}+{trim_left[1]}",
                "available_raw_sum": sum(max_trunc),
                "slack": sum(max_trunc) - need,
            },
            "read_through": f"{frac_through:.1%} of probed reads reached revcomp of the "
                            f"reverse primer",
        },
        "sample_identity": {
            "status": ("verified" if verdict == "pass"
                       else "partial" if verdict.startswith("pass") else "unresolved"),
            "from_step2_gate": verdict[:300],
            "n_biological_samples": len(bio),
            "n_technical_controls": len(ctrl),
            "exclude_alias_prefixes": [],
            "paper_sample_id_confirmed": f"{conf('paper_sample_id')}/{len(bio)}",
            "experimental_group_confirmed": f"{conf('experimental_group')}/{len(bio)}",
            "usable_for": (["ASV / relative-abundance tables",
                            "grouped comparisons"] if verdict == "pass"
                           else ["ASV / relative-abundance tables"]),
            "not_usable_for": ([] if verdict == "pass"
                               else ["grouped comparisons: see from_step2_gate"]),
            "note": "Technical controls come from step 2's per-sample "
                    "is_technical_control. Add alias prefixes here only if step 2 "
                    "missed a control that the paper names as one.",
        },
    }
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(facts, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")
    print(f"  wrote {dest}")
    print(f"  merge budget: {sum(max_trunc)} available vs {need} needed "
          f"(slack {sum(max_trunc) - need})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
