#!/usr/bin/env python3
"""Bacterial/archaeal profiling of bulk DNA shotgun metagenomes, with sylph.

Why sylph
---------
sylph estimates containment ANI between reference genomes and a sample by
k-mer sketching, and was the most accurate of seven profilers on the CAMI2
marine benchmark while using far less CPU and memory than Kraken2 for
multi-sample profiling (Shaw & Yu, Nature Biotechnology 2024). Memory matters
here: a Kraken2 standard database wants tens of GB resident, and this cluster
has no bioinformatics database storage set up yet.

What the numbers ARE
--------------------
sylph reports two percentages per genome and they are NOT read counts:

    Taxonomic_abundance   genome-size-normalised - the proportion of the
                          community's CELLS, which is what "relative abundance"
                          normally means in a microbiome paper
    Sequence_abundance    the proportion of READS, which a large genome inflates

Neither is a count. This script never writes a `count` column for sylph output,
because rounding a percentage into integers would invent data that does not
exist. If you need counts, that is a different tool (Kraken2 + Bracken produces
*estimated* counts, still not observations).

Steps
-----
    1. fastp      adapter and quality trimming, per sample
    2. bowtie2    optional host-read removal against a host index
    3. sylph      profile against a pre-sketched GTDB database
    4. matrices   one row per sample

    python step3_shotgun.py --paper X --reads-root ... --out-dir ... \\
        --sylph-db /path/gtdb-r220-c200.syldb [--host-index /path/host] \\
        [--sylph-tax GTDB_r220]
"""

import argparse
import csv
import os
import shutil
import subprocess
from pathlib import Path


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=isinstance(cmd, str), text=True,
                          capture_output=True, **kw)


def load_samples(out_dir):
    p = out_dir / "samples_used.csv"
    if not p.exists():
        raise SystemExit(f"missing {p}; run step 3 preflight first")
    with p.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def fastqs_for(out_dir, sample_id):
    with (out_dir / "run_sample_mapping.csv").open(encoding="utf-8-sig", newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r["biosample_accession"] == sample_id]
    return (sorted(r["fastq_path"] for r in rows if r["mate"] == "R1"),
            sorted(r["fastq_path"] for r in rows if r["mate"] == "R2"),
            sorted(r["fastq_path"] for r in rows if r["mate"] == "single"))


def qc_sample(sid, r1, r2, single, work, threads):
    """fastp, then optional host removal. Returns (r1, r2, stats)."""
    work.mkdir(parents=True, exist_ok=True)
    o1, o2 = work / f"{sid}_1.fq.gz", work / f"{sid}_2.fq.gz"
    rep = work / f"{sid}.fastp.json"
    if r1 and r2:
        cmd = ["fastp", "-i", r1[0], "-I", r2[0], "-o", str(o1), "-O", str(o2)]
        if len(r1) > 1:          # several runs: fastp takes one pair, so pool first
            cat1, cat2 = work / f"{sid}_pool_1.fq.gz", work / f"{sid}_pool_2.fq.gz"
            for dest, parts in ((cat1, r1), (cat2, r2)):
                with dest.open("wb") as out:
                    for p in parts:
                        out.write(Path(p).read_bytes())
            cmd = ["fastp", "-i", str(cat1), "-I", str(cat2), "-o", str(o1), "-O", str(o2)]
    else:
        cmd = ["fastp", "-i", (single or r1)[0], "-o", str(o1)]
        o2 = None
    cmd += ["--thread", str(min(threads, 16)), "--json", str(rep), "--html", "/dev/null",
            "--detect_adapter_for_pe" if o2 else "--dont_eval_duplication"]
    res = sh(cmd)
    if res.returncode != 0:
        raise SystemExit(f"{sid}: fastp failed: {res.stderr[-400:]}")
    return o1, o2, rep


def remove_host(sid, r1, r2, host_index, work, threads):
    """Keep the pairs that do NOT align to the host."""
    k1, k2 = work / f"{sid}_nohost_1.fq.gz", work / f"{sid}_nohost_2.fq.gz"
    log = work / f"{sid}.host.log"
    if r2:
        cmd = (f"bowtie2 -x {host_index} -1 {r1} -2 {r2} --threads {threads} "
               f"--very-sensitive --un-conc-gz {work}/{sid}_nohost_%.fq.gz "
               f"-S /dev/null 2> {log}")
    else:
        k2 = None
        cmd = (f"bowtie2 -x {host_index} -U {r1} --threads {threads} "
               f"--very-sensitive --un-gz {k1} -S /dev/null 2> {log}")
    res = sh(cmd)
    if res.returncode != 0:
        raise SystemExit(f"{sid}: host removal failed: {res.stderr[-400:]}")
    return k1, k2, log


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paper", required=True)
    ap.add_argument("--reads-root", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--sylph-db", required=True, help="a pre-sketched .syldb")
    ap.add_argument("--sylph-tax", default="",
                    help="sylph-tax metadata name, e.g. GTDB_r220; must match --sylph-db")
    ap.add_argument("--host-index", default="",
                    help="bowtie2 index of the host genome; host reads are removed if given")
    ap.add_argument("--threads", type=int,
                    default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    ap.add_argument("--keep-intermediate", action="store_true")
    args = ap.parse_args()

    for tool in ("fastp", "sylph"):
        if not shutil.which(tool):
            raise SystemExit(f"{tool} not on PATH")
    if args.host_index and not shutil.which("bowtie2"):
        raise SystemExit("bowtie2 not on PATH but --host-index was given")
    if not Path(args.sylph_db).exists():
        raise SystemExit(f"no sylph database at {args.sylph_db}")

    out = args.out_dir / "shotgun"
    work = out / "work"
    out.mkdir(parents=True, exist_ok=True)
    samples = load_samples(args.out_dir)
    print(f"{args.paper}: {len(samples)} shotgun samples")

    qc_rows, prof_inputs = [], {}
    for i, s in enumerate(samples, 1):
        sid = s["sample_id"]
        r1, r2, single = fastqs_for(args.out_dir, sid)
        if not (r1 or single):
            raise SystemExit(f"{sid}: no FASTQ in run_sample_mapping.csv")
        o1, o2, rep = qc_sample(sid, r1, r2, single, work, args.threads)
        host_log = ""
        if args.host_index:
            o1, o2, hl = remove_host(sid, o1, o2, args.host_index, work, args.threads)
            host_log = str(hl)
        prof_inputs[sid] = (o1, o2)
        qc_rows.append({"sample_id": sid, "fastp_report": str(rep),
                        "host_removal_log": host_log,
                        "host_index": args.host_index or "none"})
        print(f"  [{i}/{len(samples)}] {sid} QC done", flush=True)

    # ---- sylph -------------------------------------------------------------
    tsv = out / "sylph_profile.tsv"
    paired = all(v[1] for v in prof_inputs.values())
    cmd = ["sylph", "profile", args.sylph_db, "-t", str(args.threads), "-o", str(tsv)]
    if paired:
        cmd += ["-1"] + [str(v[0]) for v in prof_inputs.values()]
        cmd += ["-2"] + [str(v[1]) for v in prof_inputs.values()]
    else:
        cmd += [str(v[0]) for v in prof_inputs.values()]
    res = sh(cmd)
    if res.returncode != 0:
        raise SystemExit(f"sylph profile failed: {res.stderr[-600:]}")

    with tsv.open(encoding="utf-8", newline="") as fh:
        prof = list(csv.DictReader(fh, delimiter="\t"))
    if not prof:
        raise SystemExit("sylph detected no genome in any sample. Not writing an "
                         "empty matrix - check that the database matches these "
                         "samples' environment and that QC left enough reads.")

    # sylph names each sample by its input file; map back to sample_id.
    byfile = {str(v[0]): sid for sid, v in prof_inputs.items()}
    def sid_of(row):
        f = row.get("Sample_file", "")
        return byfile.get(f) or byfile.get(str(Path(f))) or f

    genomes = sorted({r["Genome_file"] for r in prof})
    sids = [s["sample_id"] for s in samples]
    cell = {(sid_of(r), r["Genome_file"]): r for r in prof}

    def write(name, col):
        with (out / name).open("w", encoding="utf-8-sig", newline="") as fh:
            w = csv.writer(fh); w.writerow(["sample_id"] + genomes)
            for sid in sids:
                w.writerow([sid] + [cell.get((sid, g), {}).get(col, 0) or 0
                                    for g in genomes])

    # Both are percentages from sylph. Neither is a count, and they are written
    # under names that say which is which.
    write("bacteria_taxonomic_abundance_pct.csv", "Taxonomic_abundance")
    write("bacteria_sequence_abundance_pct.csv", "Sequence_abundance")

    with (out / "feature_metadata.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        keys = ["Genome_file", "Contig_name", "Adjusted_ANI", "Eff_cov",
                "ANI_5-95_percentile", "Median_cov", "Mean_cov_geq1"]
        w.writerow(keys)
        seen = set()
        for r in prof:
            if r["Genome_file"] in seen:
                continue
            seen.add(r["Genome_file"])
            w.writerow([r.get(k, "") for k in keys])

    with (out / "qc_stats.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(qc_rows[0].keys()))
        w.writeheader(); w.writerows(qc_rows)

    if args.sylph_tax:
        r = sh(["sylph-tax", "taxprof", str(tsv), "-t", args.sylph_tax,
                "-o", str(out / "taxonomy_")])
        (out / "sylph_tax.log").write_text(r.stdout + r.stderr, encoding="utf-8")

    if not args.keep_intermediate:
        shutil.rmtree(work, ignore_errors=True)

    (out / "README.txt").write_text("\n".join([
        f"paper: {args.paper}",
        "branch: bacteria/archaea from bulk DNA shotgun",
        f"samples: {len(sids)}; genomes reported: {len(genomes)}",
        "",
        "method:",
        "  1. fastp (adapter + quality trimming)",
        f"  2. host removal: {args.host_index or 'NOT PERFORMED (no --host-index given)'}",
        f"  3. sylph profile against {args.sylph_db}",
        f"  4. taxonomy: {args.sylph_tax or 'not annotated (no --sylph-tax given)'}",
        "",
        "what the numbers are - read this before using them:",
        "  bacteria_taxonomic_abundance_pct.csv",
        "      sylph Taxonomic_abundance: genome-size-normalised percentage, i.e.",
        "      the proportion of CELLS. This is the column that corresponds to what",
        "      a microbiome paper usually calls relative abundance.",
        "  bacteria_sequence_abundance_pct.csv",
        "      sylph Sequence_abundance: percentage of READS. Large genomes inflate",
        "      this. Do not mix it with the column above.",
        "",
        "  NEITHER IS A READ COUNT. sylph is a containment-ANI profiler; it does not",
        "  assign individual reads to genomes, so there is no count to report and",
        "  none is written. Rounding these percentages to integers would fabricate",
        "  observations.",
        "",
        "  Denominator: sylph's percentages are over the DETECTED portion of the",
        "  community in that sample. Reads from organisms absent from the database",
        "  are not in the denominator, so this is not a fraction of all reads.",
        "",
        "limits:",
        "  Only genomes in the database can be detected; database choice decides",
        "  what you can see. Bacterial and viral tables come from different methods",
        "  and different denominators and must not be concatenated.",
    ]) + "\n", encoding="utf-8")

    print(f"done: {len(genomes)} genomes x {len(sids)} samples -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
