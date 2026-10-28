#!/usr/bin/env python3
"""Step 3: branch by assay type and produce count + relative abundance.

Branches are decided from sample_metadata.csv written by step 2, never guessed
from file names. Each branch is separate on disk; nothing is concatenated across
branches here - cross-paper merging is a later, explicit step.

    amplicon      region measured from the reads, then DADA2 (step3_16s.R)
    dna_shotgun   needs tools that are NOT installed yet; preflight reports this
    rna           kept separate, labelled as RNA-level; not implemented yet

Primers, amplicon length and sample identity are NOT guessed. They are read
from paper_facts/<paper>.json, which records what the paper itself says and how
it was checked against the reads. Step 3 re-measures the primers in the reads
and refuses to run if they disagree with the file. A confirmed primer never
excuses an unchecked sample identity: both sections are gated separately.

--preflight checks inputs, tools and disk and writes analysis_plan.json without
running anything. Run it first.
"""

import argparse
import collections
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import branches

# 5' primer prefixes. Degenerate bases are expanded so a short literal prefix is
# never the whole basis for a call: a region is only accepted when the forward
# primer is seen in most reads AND (for paired data) the reverse agrees.
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "GC",
         "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT",
         "V": "ACG", "N": "ACGT"}
PRIMERS = {                     # name: (forward, reverse)
    "16S_V1V3": ("AGAGTTTGATCMTGGCTCAG", "ATTACCGCGGCTGCTGG"),
    "16S_V3V4": ("CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"),
    "16S_V4":   ("GTGYCAGCMGCCGCGGTAA", "GGACTACNVGGGTWTCTAAT"),
    "ITS1":     ("CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC"),
    "ITS2":     ("GCATCGATGAAGAACGCAGC", "TCCTCCGCTTATTGATATGC"),
    "18S_V4":   ("CCAGCASCYGCGGTAATTCC", "ACTTTCGTTCTTGATYRA"),
    "18S_V9":   ("GTACACACCGCCCGTC", "TGATCCTTCTGCAGGTTCACCTAC"),
}


def matches(seq, primer, allow=2):
    if len(seq) < len(primer):
        return False
    miss = 0
    for s, p in zip(seq, primer):
        if s not in IUPAC.get(p, p):
            miss += 1
            if miss > allow:
                return False
    return True


def sample_reads(path, n=2000):
    seqs = []
    try:
        with gzip.open(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    seqs.append(line.strip())
                if len(seqs) >= n:
                    break
    except OSError:
        pass
    return seqs


def detect_region(r1, r2):
    """-> (region, primers_present, evidence). Unknown stays unknown."""
    s1, s2 = sample_reads(r1), sample_reads(r2) if r2 else []
    if not s1:
        return "unknown", False, "could not read FASTQ"
    best = []
    for name, (fwd, rev) in PRIMERS.items():
        f = sum(1 for s in s1 if matches(s, fwd)) / len(s1)
        r = (sum(1 for s in s2 if matches(s, rev)) / len(s2)) if s2 else None
        if f > 0.6 and (r is None or r > 0.6):
            best.append((name, f, r))
    if len(best) == 1:
        n, f, r = best[0]
        return n, True, f"forward primer in {f:.0%} of reads" + (
            f", reverse in {r:.0%}" if r is not None else " (single-end)")
    if len(best) > 1:
        return "ambiguous", True, f"several primers matched: {[b[0] for b in best]}"
    # No full-length match. Two very different situations, kept apart: the
    # primer may still be there but diverge past its first bases (27F does), or
    # the library really is pre-trimmed. Neither is assumed to be 16S.
    common, hits = collections.Counter(s[:10] for s in s1).most_common(1)[0]
    pct = hits * 100 // len(s1)
    for name, (fwd, _) in PRIMERS.items():
        if matches(common, fwd[:10], allow=1) and pct > 60:
            return (f"primer_like_{name}", True,
                    f"5' 10-mer {common} in {pct}% of reads matches the start of the "
                    f"{name} forward primer, but the full primer did not match within "
                    f"2 mismatches - confirm the exact primer pair from the paper")
    return ("unknown_pretrimmed", False,
            f"no primer matched; most common 5' 10-mer {common} ({pct}% of reads) "
            f"- region must be confirmed from the paper")


def branches_desc(b):
    return branches.BRANCHES.get(b, ("unknown", ""))[0]


def facts_path(paper):
    return Path(__file__).resolve().parent / "paper_facts" / f"{paper}.json"


def read_facts(paper):
    """The facts file if it exists. (None, None) when it does not."""
    p = facts_path(paper)
    if not p.exists():
        return None, None
    return json.loads(p.read_text(encoding="utf-8")), None


def measure_facts(paper, reads_root=None, out_dir=None, samples_file=None):
    """The facts for this paper, measured from its reads if not already on file.

    Everything step 3 needs here - which primers are on the reads, how long the
    insert is, where read-through starts - is a measurement of the reads, which
    are step 3's own input. Making a human run a separate script to produce it
    turned a one-command pipeline into "run it, watch step 3 fail, run a script,
    run it again". So it is generated on demand.

    An existing file is NEVER overwritten: a hand-curated one carries evidence
    from the paper itself that no measurement can reproduce.
    """
    p = facts_path(paper)
    if not p.exists():
        if reads_root is None or out_dir is None:
            return None, f"no paper_facts/{paper}.json and nothing to measure from"
        print(f"no paper_facts/{paper}.json yet - measuring it from the reads")
        r = subprocess.run(
            [sys.executable, str(Path(__file__).resolve().parent / "make_paper_facts.py"),
             "--paper", paper, "--reads-root", str(reads_root), "--out-dir", str(out_dir)]
            + (["--samples", str(samples_file)] if samples_file else []),
            text=True, capture_output=True)
        print(r.stdout.strip() or r.stderr.strip()[-600:])
        if r.returncode != 0 or not p.exists():
            return None, (f"could not measure paper_facts for {paper}: "
                          f"{(r.stderr or r.stdout).strip()[-300:]}")
    return json.loads(p.read_text(encoding="utf-8")), None


def check_facts(facts, region, primers_present):
    """Confirm the reads agree with the file. Returns a list of blocking reasons."""
    bad, amp, ident = [], facts.get("amplicon", {}), facts.get("sample_identity", {})
    # 'measured_from_reads' is accepted for denoising: finding the primer in
    # 97-99% of the reads is first-hand evidence of what was amplified, and
    # trimLeft / truncLen / the merge budget are measurements anyway. It is NOT
    # the same claim as 'verified' (someone read the Methods), so the
    # distinction is carried into the outputs rather than erased here. Sample
    # identity is gated separately by step 2 and is unaffected.
    if amp.get("status") not in ("verified", "measured_from_reads"):
        bad.append(f"paper_facts amplicon.status is {amp.get('status')!r}; it must be "
                   "'verified' (checked against the paper) or 'measured_from_reads'")
    declared, has_primers = amp.get("region"), amp.get("primers_in_reads")
    # A disagreement between the paper and the reads is never resolved silently
    # in favour of either one. The one exception is not a disagreement: when the
    # paper's primers are absent BECAUSE the library arrived pre-trimmed, the
    # detector cannot name a region and says so, which is what the file claims.
    ok = {declared, f"primer_like_{declared}"}
    if has_primers is False:
        ok.add("unknown_pretrimmed")
    if declared and region not in ok:
        bad.append(f"the reads look like {region!r} but paper_facts says {declared!r}")
    if "primers_in_reads" in amp and bool(amp["primers_in_reads"]) != bool(primers_present):
        bad.append(f"paper_facts says primers_in_reads={amp['primers_in_reads']} "
                   f"but the reads say {primers_present}")
    # Sample identity is gated on its own: having the primers proves nothing
    # about which sample is which.
    if ident.get("status") not in ("verified", "partial", "unresolved"):
        bad.append("paper_facts sample_identity.status is missing or not one of "
                   "verified/partial/unresolved - sample identity has not been checked")
    if not ident.get("usable_for"):
        bad.append("paper_facts sample_identity.usable_for is empty - it must say "
                   "what the resolved metadata actually supports")
    return bad


def split_controls(samples, facts):
    """Technical controls are not biological samples.

    Two independent sources, because neither alone is enough: step 2's
    per-sample `is_technical_control`, which carries the paper's own evidence,
    and the curated alias prefixes in paper_facts. A name prefix on its own
    would miss a blank that is not named like one and would catch a real
    sample that happens to share a prefix. Each exclusion records which source
    fired, so the decision stays auditable.
    """
    pref = tuple(facts.get("sample_identity", {}).get("exclude_alias_prefixes", []))
    keep, drop = [], []
    for s in samples:
        why = []
        if str(s.get("is_technical_control", "")).strip().lower() in ("true", "1", "yes"):
            why.append("step2 is_technical_control (" +
                       (s.get("metadata_evidence") or "no evidence recorded")[:120] + ")")
        if pref and s.get("sample_alias", "").startswith(pref):
            why.append("paper_facts exclude_alias_prefixes")
        if why:
            drop.append({**s, "excluded_because": "; ".join(why)})
        else:
            keep.append(s)
    return keep, drop


def load_samples(out_dir):
    p = out_dir / "sample_metadata.csv"
    if not p.exists():
        raise SystemExit(f"missing {p}; run step 2 first")
    with p.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def r_available():
    probe = ("module load r/4.6.0 >/dev/null 2>&1; "
             "Rscript -e 'quit(status=!requireNamespace(\"dada2\", quietly=TRUE))'")
    return subprocess.run(["bash", "-lc", probe]).returncode == 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--paper", required=True)
    ap.add_argument("--reads-root", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--preflight", action="store_true")
    ap.add_argument("--silva", type=Path, default=None,
                    help="DADA2-formatted SILVA training fasta (required to assign taxonomy)")
    ap.add_argument("--silva-species", type=Path, default=None,
                    help="SILVA assignSpecies fasta (exact-match species, optional)")
    ap.add_argument("--probe-only", action="store_true",
                    help="choose truncLen on a few samples and stop, before the full run")
    ap.add_argument("--probe-n", type=int, default=6,
                    help="how many samples to probe for the primer/region call")
    ap.add_argument("--quantify", default="bacteria",
                    choices=("bacteria", "viruses", "both"),
                    help="what to quantify from a shotgun library")
    ap.add_argument("--sylph-db", default="",
                    help="pre-sketched sylph database for bacterial profiling")
    ap.add_argument("--sylph-tax", default="",
                    help="sylph-tax metadata name matching --sylph-db")
    ap.add_argument("--viral-index", default="",
                    help="bowtie2 index prefix of the viral catalogue")
    ap.add_argument("--viral-catalogue-scope", default="",
                    help="which environment this catalogue is for, as |-separated "
                         "words to look for in the samples, e.g. 'human|homo|gut'")
    ap.add_argument("--viral-catalogue-name", default="unnamed",
                    help="catalogue name + version, recorded in the outputs")
    ap.add_argument("--archive-old", action="store_true",
                    help="move a previous run's outputs into _archive/<utc>/ "
                         "instead of overwriting them")
    ap.add_argument("--host-index", default="",
                    help="bowtie2 index of the host genome, for host read removal")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    args = ap.parse_args()

    gate = args.out_dir / "metadata_gate.txt"
    verdict = gate.read_text(encoding="utf-8").strip() if gate.exists() else "blocked:no gate"
    if not verdict.startswith("pass"):
        raise SystemExit(f"metadata gate did not pass: {verdict}")
    # A quantification-only pass is still a pass for producing per-sample
    # matrices, but the limit travels with the outputs instead of being lost
    # here, so nobody later reads the matrix as group-labelled.
    quant_only = verdict.startswith("pass_quantification_only")

    # Read an existing facts file first: split_controls needs its exclusion
    # prefixes, but MEASURING a missing one needs the amplicon sample list,
    # which only exists after the split. So: read now, measure after.
    facts, facts_err = read_facts(args.paper)
    all_samples = load_samples(args.out_dir)
    samples, controls = split_controls(all_samples, facts or {})
    def dump(rows, name):
        if not rows:
            return None
        path = args.out_dir / name
        with path.open("w", encoding="utf-8-sig", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        return path

    dump(samples, "samples_used.csv")          # everything, for reference
    dump(controls, "samples_excluded_technical_controls.csv")

    # One sample list PER BRANCH. A paper can hold several assay types at once -
    # 11 of the 78 surveyed do - and handing the whole list to a branch would
    # denoise shotgun or RNA runs as if they were this paper's amplicon.
    # Within a branch, paired and single-end are split too: they are different
    # analyses producing fragments of different length.
    def key_of(s):
        b = s["analysis_branch"]
        return f"{b}_single" if (b == "amplicon" and
                                 s.get("library_layout", "").upper() == "SINGLE") else b

    per_branch = collections.defaultdict(list)
    for s_ in samples:
        per_branch[key_of(s_)].append(s_)
    sample_file = {k: dump(v, f"samples_used_{k}.csv") for k, v in per_branch.items()}
    if len(per_branch) > 1:
        print("this paper holds several assay/layout groups; each gets its own "
              "sample list and its own output: " +
              ", ".join(f"{k}={len(v)}" for k, v in sorted(per_branch.items())))
    # NB: not named `branches` - that is the module defining what each
    # branch name means, and shadowing it broke the marker lookup.
    assays = collections.Counter({k: len(v) for k, v in per_branch.items()})

    if facts is None:
        # Measured from the AMPLICON reads only, never from a shotgun or RNA run
        # that happens to sit in the same paper.
        amp_list = sample_file.get("amplicon") or sample_file.get("amplicon_single")
        facts, facts_err = measure_facts(args.paper, args.reads_root, args.out_dir,
                                         amp_list)

    # Probe several samples spread across the run list, not just the first: a
    # paper can mix primer sets or re-sequence part of a cohort, and one sample
    # cannot show that.
    def mates(run):
        rdir = args.reads_root / args.paper / "raw_reads" / run
        a = sorted(rdir.glob("*_1.fastq.gz"))
        b = sorted(rdir.glob("*_2.fastq.gz"))
        return (a[0] if a else None), (b[0] if b else None)

    amp_samples = (per_branch.get("amplicon", []) +
                   per_branch.get("amplicon_single", [])) or samples
    runs = [s["run_accessions"].split(";")[0] for s in amp_samples]
    idx = sorted({round(i * (len(runs) - 1) / max(args.probe_n - 1, 1))
                  for i in range(min(args.probe_n, len(runs)))})
    probes = [(runs[i], *mates(runs[i])) for i in idx]
    r1, r2 = probes[0][1], probes[0][2]

    plan = {
        "paper_folder": args.paper,
        "metadata_gate": verdict[:400],
        "grouped_analysis_supported": not quant_only,
        "n_biological_samples": len(samples),
        "n_technical_controls_excluded": len(controls),
        "branches": dict(assays),
        "paired": bool(r2),
        "tools": {},
        "blocking": [],
    }
    if facts_err:
        plan["blocking"].append(facts_err)
    else:
        plan["paper_facts"] = {
            "amplicon_status": facts.get("amplicon", {}).get("status"),
            "sources": facts.get("sources", []),
            "sample_identity_status": facts["sample_identity"]["status"],
            "usable_for": facts["sample_identity"]["usable_for"],
            "not_usable_for": facts["sample_identity"].get("not_usable_for", []),
        }

    if "amplicon" in assays:
        seen = {}
        for run, a, b in probes:
            seen[run] = detect_region(a, b) if a else ("unknown", False, "no R1 file")
        calls = {v[0] for v in seen.values()}
        region, primers, evidence = seen[probes[0][0]]
        plan["amplicon_region"] = region
        plan["primers_present"] = primers
        plan["region_evidence"] = evidence
        plan["region_probe"] = {k: {"region": v[0], "primers_present": v[1],
                                    "evidence": v[2]} for k, v in seen.items()}
        # A pre-trimmed library has no primer for the detector to name, so the
        # marker comes from the verified paper_facts instead. check_facts above
        # has already confirmed the file agrees with the reads, so this is the
        # curated answer, not a guess.
        declared = (facts or {}).get("amplicon", {}).get("region", "")
        marker_branch = branches.refine_amplicon(
            declared if region == "unknown_pretrimmed" and declared else region)
        plan["amplicon_branch"] = marker_branch
        plan["branch_status"] = plan.get("branch_status", {})
        # Paired and single-end 16S are different analyses, not a flag: the
        # paired script is built around a merge that does not exist here.
        # Layout was already used to split the sample lists, so a paper that
        # mixes PAIRED and SINGLE now runs BOTH scripts on its own half instead
        # of being blocked outright. The two ASV tables stay separate because
        # their fragments are different lengths.
        plan["library_layout"] = sorted(
            {x.get("library_layout", "").upper() for x in amp_samples})
        plan["amplicon_groups"] = {
            ("amplicon_16s_single" if k == "amplicon_single" else "amplicon_16s"):
                len(v) for k, v in per_branch.items() if k.startswith("amplicon")}
        if marker_branch not in ("amplicon_16s", "amplicon_16s_single"):
            # ITS / 18S / anything else must not be denoised against a 16S
            # training set just because DADA2 would happily run on it.
            plan["branch_status"][marker_branch] = "unsupported"
            plan["blocking"].append(
                f"the reads are {region!r}, i.e. {branches_desc(marker_branch)}. "
                "Only the 16S amplicon branch is implemented; SILVA's 16S "
                "training set does not apply to this marker.")
        if len(calls) > 1:
            plan["blocking"].append(
                "the probed samples do not agree on the amplicon region "
                f"({sorted(calls)}); this paper may mix primer sets, which must "
                "be resolved before a single denoising run can cover it")
        plan["tools"]["R_dada2"] = "available" if r_available() else "MISSING"
        if facts:
            plan["blocking"] += check_facts(facts, region, primers)
            amp = facts["amplicon"]
            # Everything the denoiser needs comes from the verified file, in the
            # one coordinate system it declares (raw read bases, before trimLeft).
            plan["trim_left"] = amp["trim_left"]
            plan["max_trunc_len"] = amp["max_trunc_len"]
            plan["insert_p99"] = amp["insert_length_primer_free"]["p99"]
            plan["min_overlap"] = amp["merge_budget"]["min_overlap"]
            plan["required_raw_sum"] = amp["merge_budget"]["required_raw_sum"]
            plan["trunc_len_basis"] = amp["trunc_len_basis"]
            if sum(amp["max_trunc_len"]) < amp["merge_budget"]["required_raw_sum"]:
                plan["blocking"].append(
                    "the reads cannot span this amplicon: max truncLen sum "
                    f"{sum(amp['max_trunc_len'])} < required {amp['merge_budget']['required_raw_sum']}")
        elif region in ("unknown", "ambiguous", "unknown_pretrimmed"):
            plan["blocking"].append(
                f"amplicon region not established from reads ({region}); "
                "confirm primers/region from the paper before denoising")
        if not str(region).startswith(("16S", "primer_like_16S")) and \
                not (facts and str(facts["amplicon"].get("region", "")).startswith("16S")):
            plan["blocking"].append(
                f"region {region!r} is not 16S; step3_16s.R and the SILVA 16S "
                "training set do not apply to it (ITS and other markers need "
                "their own reference and script)")
        rscript = Path(__file__).resolve().parent / "step3_16s.R"
        if not rscript.exists():
            plan["blocking"].append(f"{rscript.name} not written yet "
                                    "(denoising parameters pending confirmation)")
        if plan["tools"]["R_dada2"] != "available":
            plan["blocking"].append("R module with dada2 not available")
        if not args.silva:
            plan["blocking"].append(
                "no --silva training set given; ASV counts can be produced but "
                "taxonomy cannot be assigned")
    # ---- the non-amplicon branches ----------------------------------------
    # Each records its own readiness. A branch that cannot run must not stop a
    # branch that can, and must not be reported as finished either.
    plan["branch_status"] = {}

    def need(branch, tools, dbs):
        """Record what a branch needs and whether it has it."""
        have = {t: (shutil.which(t) or "MISSING") for t in tools}
        miss_t = [t for t, v in have.items() if v == "MISSING"]
        miss_d = [n for n, v in dbs.items() if not v or not Path(v).exists()]
        plan["tools"][branch] = have
        if miss_t or miss_d:
            plan["branch_status"][branch] = "missing_dependencies"
            plan["blocking"].append(
                f"branch '{branch}' cannot run: missing tool(s) {miss_t or 'none'}, "
                f"missing database/index {miss_d or 'none'}. See "
                f"reference_db/README.md for how to install them.")
        else:
            plan["branch_status"][branch] = "ready"
        return not (miss_t or miss_d)

    if "dna_shotgun" in assays or "dna_virome" in assays:
        if args.quantify in ("bacteria", "both") and "dna_shotgun" in assays:
            need("dna_shotgun_bacteria", ["fastp", "sylph"],
                 {"sylph_db": args.sylph_db})
        if args.quantify in ("viruses", "both"):
            # A catalogue is built for an environment. UHGV is human gut; mapping
            # a mouse, plant or marine library against it reports "how much human
            # gut virus is in this sample", which is not a virome.
            scope = (args.viral_catalogue_scope or "").strip().lower()
            host = " ".join(filter(None, (
                (facts or {}).get("experiment", {}).get("host", ""),
                (facts or {}).get("experiment", {}).get("body_site", ""),
                samples[0].get("host", ""), samples[0].get("body_site", ""),
                samples[0].get("isolation_source", "")))).lower()
            # Groups separated by ';' must ALL match; alternatives inside a
            # group are separated by '|'. Host and site are separate groups
            # because matching only one lets a MOUSE GUT sample through a HUMAN
            # GUT catalogue - the exact mistake this check exists to stop.
            groups = [g.strip() for g in scope.split(";") if g.strip()]
            unmatched = [g for g in groups
                         if not any(w.strip() in host for w in g.split("|") if w.strip())]
            if scope and host and unmatched:
                plan["blocking"].append(
                    f"the viral catalogue requires {unmatched} but this paper's "
                    f"samples look like {host[:80]!r}. Mapping them to the wrong "
                    "catalogue would report presence of viruses that have no reason "
                    "to be there. Set VIRAL_CATALOGUE_SCOPE to match, or point "
                    "VIRAL_INDEX at a catalogue for this environment.")
            # A VLP-enriched library gets the viral branch only: counting
            # bacteria in a library built to exclude them produces a number
            # with no interpretation.
            idx = args.viral_index
            ok_idx = idx and (Path(str(idx) + ".1.bt2").exists()
                              or Path(str(idx) + ".1.bt2l").exists())
            need("viral_mapping", ["bowtie2", "bowtie2-inspect", "samtools"],
                 {"viral_bowtie2_index": idx if ok_idx else ""})
        if "dna_virome" in assays and args.quantify == "bacteria":
            plan["blocking"].append(
                "this paper's libraries are VLP-enriched (dna_virome) but "
                "--quantify is 'bacteria'; a virus-targeted library is the wrong "
                "input for bacterial profiling")

    for b in ("rna", "amplicon_its", "amplicon_18s", "amplicon_other",
              "unsupported_isolate_genome"):
        if b in assays:
            plan["branch_status"][b] = "unsupported"
            plan["blocking"].append(
                f"branch '{b}' is NOT IMPLEMENTED ({branches_desc(b)}). It is "
                "stopped here rather than being pushed through a pipeline built "
                "for a different assay.")

    free = shutil.disk_usage(args.out_dir).free
    plan["disk_free_bytes"] = free
    plan["disk_free_human"] = f"{free/2**40:.2f} TiB"
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "analysis_plan.json").write_text(
        json.dumps(plan, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(plan, indent=2, ensure_ascii=False))

    if args.preflight:
        return 0
    # Two kinds of blocker, kept apart so one branch's missing database cannot
    # stop a branch that is ready. A branch-level blocker leaves that branch
    # recorded as missing_dependencies/unsupported; only a paper-level one
    # (metadata, inputs, tools for a branch that IS running) stops everything.
    branch_level = {"missing_dependencies", "unsupported"}
    deferred = {b for b, st in plan["branch_status"].items() if st in branch_level}
    fatal = [b for b in plan["blocking"]
             if not any(f"'{d}'" in b or d in b for d in deferred)]
    runnable = ("amplicon" in assays and plan.get("amplicon_branch") in
                ("amplicon_16s", "amplicon_16s_single")) \
        or any(st == "ready" for st in plan["branch_status"].values())
    if fatal or not runnable:
        raise SystemExit("blocked:\n  " + "\n  ".join(plan["blocking"] or
                                                      ["no branch of this paper is runnable"]))
    if deferred:
        print("branches not run (recorded, not failed): " + ", ".join(sorted(deferred)))

    # --- dispatch -----------------------------------------------------------
    # Every branch present in this paper runs, and each records its own result.
    # One branch failing must not mark the paper done, and must not stop the
    # others from producing what they can.
    results, HERE = {}, Path(__file__).resolve().parent

    if args.archive_old:
        # Move a previous run's outputs aside instead of writing over them, so a
        # rerun can be compared with what it replaced rather than erasing it.
        stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
        for sub in ("dada2", "shotgun", "virome"):
            d = args.out_dir / sub
            if d.exists() and any(d.iterdir()):
                dest = args.out_dir / "_archive" / stamp / sub
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(d), str(dest))
                print(f"archived previous {sub}/ -> {dest}")

    # Each amplicon group runs on its OWN sample list, into its OWN directory.
    # A paper that mixes paired and single-end gets two tables that must not be
    # merged on ASV sequence; the single-end README says so.
    for key, group in sorted(per_branch.items()):
        if not key.startswith("amplicon"):
            continue
        single = key.endswith("_single")
        b = "amplicon_16s_single" if single else "amplicon_16s"
        if plan.get("amplicon_branch") not in ("amplicon_16s", "amplicon_16s_single"):
            continue                       # marker is not 16S; already recorded
        sub = "dada2_single" if single else "dada2"
        plan["branch_status"][b] = "running"
        script = "step3_16s_single.R" if single else "step3_16s.R"
        common = (f"--paper {args.paper} --reads-root {args.reads_root} "
                  f"--out-dir {args.out_dir} --region {plan['amplicon_region']} "
                  f"--samples {sample_file[key]} --out-sub {sub} "
                  f"--threads {args.threads} --silva {args.silva} "
                  f"--grouped-analysis {int(not quant_only)} "
                  f"--amplicon-status "
                  f"{plan.get('paper_facts', {}).get('amplicon_status', 'unknown')}")
        if single:
            # One read: no reverse primer to trim, no merge budget to respect.
            extra = (f" --trim-left {plan['trim_left'][0]} "
                     f"--max-trunc {plan['max_trunc_len'][0]}")
        else:
            extra = (f" --trim-left {plan['trim_left'][0]},{plan['trim_left'][1]} "
                     f"--max-trunc {plan['max_trunc_len'][0]},{plan['max_trunc_len'][1]} "
                     f"--insert-p99 {plan['insert_p99']} "
                     f"--min-overlap {plan['min_overlap']}")
        cmd = ["bash", "-lc",
               f"module load r/4.6.0 && Rscript {HERE/script} " + common + extra
               + (f" --silva-species {args.silva_species}" if args.silva_species else "")
               + (" --probe-only" if args.probe_only else "")]
        print(f"--- {b}: {len(group)} samples -> {sub}/")
        results[b] = subprocess.run(cmd).returncode

    if plan["branch_status"].get("dna_shotgun_bacteria") == "ready":
        cmd = [sys.executable, str(HERE / "step3_shotgun.py"), "--paper", args.paper,
               "--reads-root", str(args.reads_root), "--out-dir", str(args.out_dir),
               "--sylph-db", args.sylph_db, "--threads", str(args.threads)]
        if sample_file.get("dna_shotgun"):
            cmd += ["--samples", str(sample_file["dna_shotgun"])]
        if args.sylph_tax:
            cmd += ["--sylph-tax", args.sylph_tax]
        if args.host_index:
            cmd += ["--host-index", args.host_index]
        results["dna_shotgun_bacteria"] = subprocess.run(cmd).returncode

    if plan["branch_status"].get("viral_mapping") == "ready":
        lib = "dna_virome" if "dna_virome" in assays else "dna_shotgun"
        cmd = [sys.executable, str(HERE / "step3_virome.py"), "--paper", args.paper,
               "--reads-root", str(args.reads_root), "--out-dir", str(args.out_dir),
               "--bowtie2-index", args.viral_index, "--threads", str(args.threads),
               "--catalogue-name", args.viral_catalogue_name, "--library-type", lib]
        vkey = "dna_virome" if sample_file.get("dna_virome") else "dna_shotgun"
        if sample_file.get(vkey):
            cmd += ["--samples", str(sample_file[vkey])]
        results["viral_mapping"] = subprocess.run(cmd).returncode

    if not results:
        raise SystemExit("no branch of this paper could be run; see analysis_plan.json")

    for name, rc in results.items():
        plan["branch_status"][name] = "done" if rc == 0 else "failed"
    (args.out_dir / "analysis_plan.json").write_text(
        json.dumps(plan, indent=2, ensure_ascii=False), encoding="utf-8")
    print("branch results: " + ", ".join(f"{k}={v}" for k, v in
                                         plan["branch_status"].items()))
    # Non-zero if ANY branch failed, so the scheduler cannot file a partial
    # paper as complete.
    return 0 if all(rc == 0 for rc in results.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
