#!/usr/bin/env python3
"""Quantify viruses that are ALREADY IN A REFERENCE CATALOGUE, by read mapping.

What this does and does not do
------------------------------
It maps each sample's reads to a catalogue of viral genomes and reports, per
genome, how much of it is covered and how deeply. That measures the abundance
of viruses the catalogue already contains.

It does NOT discover viruses. A virus absent from the catalogue is invisible
here, and no amount of mapping will reveal it. Finding new ones needs assembly,
viral identification (geNomad/VirSorter2/CheckV) and clustering into new vOTUs
- a separate pipeline that this script deliberately does not pretend to be.

Detection thresholds
--------------------
Defaults follow Roux et al. 2017 (PeerJ 5:e3817), the benchmarking study that
established the thresholds most viromics work now cites: reads mapped at >=90%
identity, and a genome called present only when >=75% of its length is covered
at >=1x. Mapping a read somewhere is not detection; a few conserved genes will
recruit reads from relatives that are not there.

    python step3_virome.py --paper X --reads-root ... --out-dir ... \\
        --catalogue /path/votus.fna --bowtie2-index /path/votus \\
        [--min-identity 0.90] [--min-breadth 0.75] [--min-depth 1]
"""

import argparse
import collections
import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=isinstance(cmd, str), text=True,
                          capture_output=True, **kw)


def load_samples(out_dir, path=None):
    # The caller passes the list for THIS branch; a paper can hold several assay
    # types and the others' runs must not come along.
    p = Path(path) if path else out_dir / "samples_used.csv"
    if not p.exists():
        raise SystemExit(f"missing {p}; run step 3 preflight first")
    with p.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def fastqs_for(out_dir, sample_id):
    """Every FASTQ of this sample, from run_sample_mapping.csv - not a glob."""
    with (out_dir / "run_sample_mapping.csv").open(encoding="utf-8-sig", newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r["biosample_accession"] == sample_id]
    r1 = sorted(r["fastq_path"] for r in rows if r["mate"] == "R1")
    r2 = sorted(r["fastq_path"] for r in rows if r["mate"] == "R2")
    single = sorted(r["fastq_path"] for r in rows if r["mate"] == "single")
    return r1, r2, single


def genome_lengths(index_prefix):
    """Reference name -> length, straight from the bowtie2 index."""
    out = sh(["bowtie2-inspect", "-s", str(index_prefix)])
    if out.returncode != 0:
        raise SystemExit(f"bowtie2-inspect failed: {out.stderr[:300]}")
    lens = {}
    for line in out.stdout.splitlines():
        if line.startswith("Sequence-"):
            _, name, length = line.split("\t")[:3]
            lens[name.split()[0]] = int(length)
    return lens


def coverage_from_bam(bam, threads):
    """Per reference: covered bases at >=1x, summed depth, mapped reads.

    Computed with `samtools coverage`, which reports exactly these without
    materialising a per-base table.
    """
    # No -@ here: `samtools coverage` has no threads option, prints its usage
    # and still EXITS 0, so passing it silently yields an empty table that
    # looks like "nothing detected".
    out = sh(["samtools", "coverage", str(bam)])
    if out.returncode != 0 or "invalid option" in out.stderr:
        raise SystemExit(f"samtools coverage failed: {out.stderr[:300]}")
    cov, parsed = {}, 0
    for line in out.stdout.splitlines():
        f = line.split("\t")
        # Skip the header and anything that is not a full coverage row, rather
        # than indexing blindly into it.
        if line.startswith("#") or len(f) < 7:
            continue
        parsed += 1
        # rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq
        cov[f[0]] = {"reads": int(f[3]), "covered_bases": int(f[4]),
                     "breadth": float(f[5]) / 100.0, "mean_depth": float(f[6])}
    if parsed == 0:
        # An unparsable table is a tooling failure, not a biological result.
        raise SystemExit(
            "samtools coverage produced no parsable rows for "
            f"{bam}; stdout began {out.stdout[:200]!r}, stderr {out.stderr[:200]!r}")
    return cov


def map_sample(sid, r1, r2, single, index, out_dir, threads, min_identity, keep_bam):
    """bowtie2 -> sorted BAM -> coverage. Returns (coverage dict, stats)."""
    work = out_dir / "virome" / "bam"
    work.mkdir(parents=True, exist_ok=True)
    bam = work / f"{sid}.bam"
    # --no-unal keeps the BAM small; end-to-end (not --local) so that an
    # alignment must span the read, which is what the identity filter assumes.
    args = ["bowtie2", "-x", str(index), "--threads", str(threads),
            "--no-unal", "--end-to-end", "--sensitive"]
    if r1 and r2:
        args += ["-1", ",".join(r1), "-2", ",".join(r2)]
    else:
        args += ["-U", ",".join(single or r1)]
    # NM/MD let us filter on identity; samtools view -e needs the NM tag.
    cmd = (" ".join(args) + f" 2> {work/(sid + '.bowtie2.log')}"
           f" | samtools sort -@ {threads} -o {bam} -")
    res = sh(cmd)
    if res.returncode != 0:
        raise SystemExit(f"{sid}: bowtie2/samtools failed: {res.stderr[:300]}")
    sh(["samtools", "index", "-@", str(threads), str(bam)])

    # Identity filter: drop alignments with more than (1 - min_identity) of the
    # aligned length mismatching. Applied before coverage so a 70%-identity hit
    # cannot contribute breadth.
    filt = work / f"{sid}.filt.bam"
    # samtools filter expressions have no int(); the float comparison is what
    # works and is what the thresholds mean anyway.
    max_div = 1.0 - min_identity
    expr = f"[NM] <= {max_div:.4f} * qlen"
    res = sh(["samtools", "view", "-@", str(threads), "-b", "-e", expr,
              "-o", str(filt), str(bam)])
    if res.returncode != 0:
        raise SystemExit(f"{sid}: identity filter failed: {res.stderr[:300]}")
    sh(["samtools", "index", "-@", str(threads), str(filt)])

    cov = coverage_from_bam(filt, threads)
    total = sh(["samtools", "view", "-c", "-@", str(threads), str(bam)]).stdout.strip()
    kept = sh(["samtools", "view", "-c", "-@", str(threads), str(filt)]).stdout.strip()
    log = (work / f"{sid}.bowtie2.log").read_text(errors="replace")
    if not keep_bam:
        for f in (bam, filt, Path(str(bam) + ".bai"), Path(str(filt) + ".bai")):
            f.unlink(missing_ok=True)
    return cov, {"alignments_total": int(total or 0),
                 "alignments_passing_identity": int(kept or 0),
                 "bowtie2_summary": log.strip().splitlines()[-1] if log.strip() else ""}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paper", required=True)
    ap.add_argument("--reads-root", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--samples", default="",
                    help="sample list for this branch; defaults to samples_used.csv")
    ap.add_argument("--bowtie2-index", required=True,
                    help="prefix of a prebuilt bowtie2 index of the viral catalogue")
    ap.add_argument("--catalogue-name", default="unnamed",
                    help="catalogue name and version, recorded in the README")
    ap.add_argument("--library-type", default="dna_shotgun",
                    help="dna_shotgun or dna_virome; recorded, never compared across")
    ap.add_argument("--min-identity", type=float, default=0.90)
    ap.add_argument("--min-breadth", type=float, default=0.75)
    ap.add_argument("--min-depth", type=float, default=1.0)
    ap.add_argument("--threads", type=int,
                    default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    ap.add_argument("--keep-bam", action="store_true")
    args = ap.parse_args()

    for tool in ("bowtie2", "bowtie2-inspect", "samtools"):
        if not shutil.which(tool):
            raise SystemExit(f"{tool} not on PATH")
    if not Path(str(args.bowtie2_index) + ".1.bt2").exists() and \
            not Path(str(args.bowtie2_index) + ".1.bt2l").exists():
        raise SystemExit(f"no bowtie2 index at {args.bowtie2_index}")

    out = args.out_dir / "virome"
    out.mkdir(parents=True, exist_ok=True)
    samples = load_samples(args.out_dir, args.samples)
    lens = genome_lengths(args.bowtie2_index)
    print(f"{args.paper}: {len(samples)} samples vs {len(lens)} reference genomes")

    detected = {}        # sample -> {genome: depth}
    per_sample_stats = []
    for i, s in enumerate(samples, 1):
        sid = s["sample_id"]
        r1, r2, single = fastqs_for(args.out_dir, sid)
        if not (r1 or single):
            raise SystemExit(f"{sid}: no FASTQ in run_sample_mapping.csv")
        cov, st = map_sample(sid, r1, r2, single, args.bowtie2_index, args.out_dir,
                             args.threads, args.min_identity, args.keep_bam)
        hits = {g: v for g, v in cov.items()
                if v["breadth"] >= args.min_breadth and v["mean_depth"] >= args.min_depth}
        detected[sid] = hits
        st.update({"sample_id": sid, "genomes_with_any_read": len(cov),
                   "genomes_detected": len(hits)})
        per_sample_stats.append(st)
        print(f"  [{i}/{len(samples)}] {sid}: {len(hits)} genomes pass "
              f"breadth>={args.min_breadth:.0%} depth>={args.min_depth}", flush=True)

    genomes = sorted({g for h in detected.values() for g in h})
    if not genomes:
        raise SystemExit("no genome in the catalogue passed the detection "
                         "thresholds in any sample; the catalogue may be wrong "
                         "for this environment, or the library may hold no "
                         "catalogued virus. Not writing an empty matrix.")

    # ---- matrices ----------------------------------------------------------
    # Length-corrected depth, not raw reads: a 200 kb genome recruits more reads
    # than a 5 kb one at the same true abundance. mean_depth is already
    # per-base, so it is the length-corrected quantity.
    sids = [s["sample_id"] for s in samples]

    def write_matrix(name, value_of):
        with (out / name).open("w", encoding="utf-8-sig", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["sample_id"] + genomes)
            for sid in sids:
                w.writerow([sid] + [value_of(sid, g) for g in genomes])

    write_matrix("viral_mean_depth.csv",
                 lambda s, g: round(detected[s].get(g, {}).get("mean_depth", 0.0), 6))
    write_matrix("viral_mapped_reads.csv",
                 lambda s, g: detected[s].get(g, {}).get("reads", 0))

    # Relative abundance WITHIN the detected viral fraction. The denominator is
    # named in the filename and in the README because a viral-internal
    # proportion is not comparable with a whole-community proportion.
    def rel(sid, g):
        tot = sum(v["mean_depth"] for v in detected[sid].values())
        return round(detected[sid].get(g, {}).get("mean_depth", 0.0) / tot, 8) if tot else 0.0
    write_matrix("viral_relative_abundance_within_virome.csv", rel)

    with (out / "feature_metadata.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["genome_id", "length_bp", "catalogue", "n_samples_detected"])
        for g in genomes:
            w.writerow([g, lens.get(g, ""), args.catalogue_name,
                        sum(1 for s in sids if g in detected[s])])

    with (out / "mapping_stats.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_sample_stats[0].keys()))
        w.writeheader(); w.writerows(per_sample_stats)

    empty = [s for s in sids if not detected[s]]
    (out / "README.txt").write_text("\n".join([
        f"paper: {args.paper}",
        f"branch: viral quantification by read mapping ({args.library_type})",
        f"catalogue: {args.catalogue_name}  ({len(lens)} genomes)",
        f"samples: {len(sids)}, of which {len(empty)} detected nothing",
        f"genomes detected in at least one sample: {len(genomes)}",
        "",
        "method:",
        "  bowtie2 --end-to-end --sensitive, then samtools sort/index",
        f"  alignments kept when NM <= {1 - args.min_identity:.2f} x read length "
        f"(>= {args.min_identity:.0%} identity)",
        f"  a genome is PRESENT only when >= {args.min_breadth:.0%} of its length is "
        f"covered at >= {args.min_depth}x",
        "  thresholds follow Roux et al. 2017, PeerJ 5:e3817 (benchmarking viromics)",
        "",
        "multi-mapping:",
        "  bowtie2 reports one alignment per read (its best; ties broken at random).",
        "  Reads shared between near-identical catalogue genomes are therefore",
        "  assigned to one of them, not split. The breadth threshold is what stops",
        "  a relative from being called present on a handful of conserved genes.",
        "",
        "files and what the numbers mean:",
        "  viral_mean_depth.csv    mean per-base depth over the genome. This is the",
        "                          length-corrected abundance; a long genome does not",
        "                          score higher merely for being long.",
        "  viral_mapped_reads.csv  raw mapped read counts, NOT length corrected.",
        "  viral_relative_abundance_within_virome.csv",
        "                          each sample's depths divided by that SAMPLE'S OWN",
        "                          total over detected catalogue genomes. Rows sum to 1.",
        "                          DENOMINATOR = the detected viral fraction only.",
        "                          It is NOT a fraction of all reads and NOT comparable",
        "                          with a bacterial relative abundance. Do not",
        "                          concatenate it with another branch's table.",
        "  mapping_stats.csv       per sample: alignments before/after the identity",
        "                          filter, genomes with any read, genomes detected.",
        "",
        "limits:",
        "  Only viruses already in the catalogue can be seen. This is not virus",
        "  discovery; novel viruses need assembly + geNomad/VirSorter2 + CheckV +",
        "  vOTU clustering.",
        f"  The library type is recorded as {args.library_type}. VLP-enriched and bulk",
        "  shotgun libraries are not directly comparable even when both are",
        "  quantified this way.",
    ]) + "\n", encoding="utf-8")

    print(f"done: {len(genomes)} genomes x {len(sids)} samples -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
