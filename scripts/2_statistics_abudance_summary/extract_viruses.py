#!/usr/bin/env python3
"""Pull every viral feature out of the integrated matrices into one viral dataset.

Scope: the 17 papers that carry viral features - 6 whose main matrix IS viral
(vOTU / viral genome) and 11 whose bacterial matrix happens to include named
viruses. Viral features are identified with the same per-feature logic the
summary statistics use, imported from summary_statistics.py so the two can never
drift apart.

Why long format
---------------
Feature IDs are paper-local: Fan2023's vOTU_1 and Wallace2025's vOTU_1 are
unrelated sequences. The combined matrix is therefore block-diagonal -
5,542 samples x 106,434 features is 590M cells of which ~95% are structurally
empty. The data is written long (one row per observation) and only non-zero
values are stored; `viral_blocks.csv` records each paper's block shape and
whether its absent cells are explicit zeros or unreported blanks, so nothing is
lost.

Units are NOT harmonised. The blocks are read counts, coverage, RPKM and
relative abundance; `abundance_unit` travels with every row and the blocks must
not be compared numerically without a deliberate normalisation step.

Usage:  python3 extract_viruses.py [--out-dir DIR]
"""

from __future__ import annotations

import argparse
import collections
import csv
import importlib.util
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = Path("/hickory/proj/didonglab/dataset/virus/Database/data/abudance_integration")
STATS = Path("/hickory/proj/didonglab/dataset/virus/Database/results/2_statistics_abudance_summary")
OUT_DIR = Path("/hickory/proj/didonglab/dataset/virus/Database/data/viruses")

SEP = "::"          # namespace separator; feature ids already use "|" and "||"

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2 ** 31 - 1)

# reuse the exact classification the statistics use
_spec = importlib.util.spec_from_file_location("stats", HERE / "summary_statistics.py")
stats = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(stats)

# Established from each paper's Methods: only Feng2024 enriched virus-like
# particles. Everything else recovers viruses from bulk shotgun data, which is
# prophage-biased and DNA-only. This travels with the data because it decides
# what the numbers can be used for.
RECOVERY = {
    "Feng2024": "VLP-enriched virome",
    "Garmaeva2024": "VLP-enriched virome",
    "Fan2023": "bulk shotgun (viral contigs + CheckV)",
    "Wallace2025": "bulk shotgun (VIBRANT/geNomad + CheckV)",
    "Merenstein2026": "bulk shotgun (Kraken2 vs vOTU DB, prophage regions masked)",
    "Lin2025": "bulk shotgun (viral genomes, read mapping)",
    "Hedžet2021": "bulk public metagenomes mapped to one phage genome",
}
DEFAULT_RECOVERY = "incidental: viruses reported by a bacterial profiler"

TAX_COLS = ["realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
SAMPLE_COLS = ["sample_type", "host", "subject_id", "disease", "disease_name",
               "disease_status", "experimental_group", "timepoint",
               "measurement_method", "batch", "is_pooled_sample", "abundance_unit"]

# Papers name the disease in whichever column their integration chose, and carry
# extra grouping columns of their own (Grazioli2022 holds all 8 of its cohort
# diseases in `disease_name`, which an earlier version of this script dropped,
# leaving 1,793 samples looking unlabelled). Rather than enumerate every variant,
# the fixed columns above are filled from the first that has a value, and every
# other label-bearing column is preserved verbatim in `source_label_fields`.
DISEASE_SOURCES = ["disease", "disease_name", "disease_label_original",
                   "paper_cohort_name", "group_description"]
LABEL_COL = re.compile(r"disease|cohort|diagnos|condition|group|phenotype|clinical|"
                       r"salinity|treatment", re.I)


def resolve_labels(src: dict) -> tuple[str, str]:
    """-> (disease name resolved across the known column variants, overflow string)."""
    disease = ""
    for c in DISEASE_SOURCES:
        v = (src.get(c) or "").strip()
        if v and not re.match(r"^(not |none$|n/?a$|unknown)", v, re.I):
            disease = v
            break
    extra = []
    for c, v in src.items():
        if not c or c in SAMPLE_COLS or not LABEL_COL.search(c):
            continue
        v = (v or "").strip()
        if v:
            extra.append(f"{c}={v}")
    return disease, "; ".join(extra)[:600]


def write_figures(dest: Path) -> None:
    """Summary tables + pies, built from the files just written.

    These live in the script rather than in a one-off notebook: the first
    version of them was generated ad hoc and went stale the moment the dataset
    grew from 17 to 20 papers.
    """
    fig = dest / "figures"
    fig.mkdir(exist_ok=True)
    feats = stats.read_csv(dest / "viral_feature_metadata.csv")
    samps = stats.read_csv(dest / "viral_sample_metadata.csv")

    deepest = collections.Counter()
    for r in feats:
        got = [c for c in TAX_COLS if (r.get(c) or "").strip()]
        deepest[got[-1].capitalize() if got else "No taxonomy"] += 1
    rank_items = [("No taxonomy", deepest["No taxonomy"]), ("Family", deepest["Family"]),
                  ("Kingdom", deepest["Kingdom"]), ("Species", deepest["Species"]),
                  ("Class", deepest["Class"]),
                  ("Other ranks", sum(deepest[k] for k in
                                      ("Realm", "Phylum", "Order", "Genus")))]
    rank_items = [(k, v) for k, v in rank_items if v]

    dis = collections.Counter(stats.disease_bucket(r.get("disease_status") or "")
                              for r in samps)
    rec = collections.Counter()
    for r in feats:
        rec["VLP-enriched virome" if "VLP" in r["viral_recovery_method"]
            else "bulk shotgun / incidental"] += 1
    comp = collections.Counter(r["metadata_completeness"] for r in samps)

    sub = f"{len({r['paper_folder'] for r in feats})} papers, {len(feats):,} viral features"
    ssub = f"{len({r['paper_folder'] for r in samps})} papers, {len(samps):,} samples"
    for fname, title, items, st_ in [
            ("vfig01_rank_depth_pie", "Viral features by deepest taxonomic rank",
             rank_items, sub),
            ("vfig02_disease_pie", "Viral samples by disease status",
             dis.most_common(), ssub),
            ("vfig03_recovery_pie", "Viral features by recovery method",
             rec.most_common(), sub),
            ("vfig04_metadata_completeness_pie", "Viral samples by metadata completeness",
             comp.most_common(), ssub)]:
        stats.svg_pie(fig / f"{fname}.svg", title, items, subtitle=st_)

    for fname, counter, key, total in [
            ("viral_rank_depth.csv", dict(rank_items), "deepest_rank", len(feats)),
            ("viral_disease_status.csv", dis, "disease_status", len(samps)),
            ("viral_recovery_method.csv", rec, "recovery_method", len(feats)),
            ("viral_metadata_completeness.csv", comp, "metadata_completeness", len(samps))]:
        stats.write_table(dest / fname, [key, "n", "share"],
                          [[k, v, stats.pct(v, total)]
                           for k, v in sorted(counter.items(), key=lambda x: -x[1])])
    print(f"  wrote 4 figures + 4 summary tables")


def usability(row: dict) -> tuple[str, str, str]:
    """-> (completeness level, what the sample can answer, what it cannot).

    A sample whose attributes are missing is still a distinct observation - its
    identity is the abundance vector, not the metadata row. Garmaeva2024's 205
    samples carry only an id, yet detect between 19 and 19,961 vOTUs and share a
    median 1.8% of them, so they must not be merged. What they lose is the
    ability to be STRATIFIED, and that is what these columns state explicitly,
    rather than leaving a reader to infer it from blank cells.
    """
    def has(*names: str) -> bool:
        for n in names:
            v = (row.get(n) or "").strip()
            if v and not re.match(r"^(not |unknown|n/?a$|none$)", v, re.I):
                return True
        return False

    subject = has("subject_id")
    timepoint = has("timepoint")
    group = has("experimental_group")
    if not group:
        bucket = stats.disease_bucket(row.get("disease_status") or "")
        group = bucket in ("Case / disease", "Control / healthy") or has(*DISEASE_SOURCES)

    can = ["prevalence", "alpha/beta diversity", "co-occurrence", "vOTU catalogue"]
    cannot = []
    (can if group else cannot).append("group comparison")
    (can if subject and timepoint else cannot).append("longitudinal / repeated measures")
    (can if subject else cannot).append("within-subject pairing")

    if subject and timepoint and group:
        level = "full"
    elif group and timepoint:
        level = "grouped_longitudinal"
    elif group:
        level = "grouped"
    elif subject or timepoint:
        level = "unlabelled_but_structured"
    else:
        level = "ids_only"
    return level, "; ".join(can), "; ".join(cannot)


def as_int(v, default=0):
    try:
        return int(float(str(v).strip()))
    except (TypeError, ValueError):
        return default


def viral_features(folder: Path, paper_domain: str) -> tuple[dict, list]:
    """-> ({original feature id: {taxonomy...}}, header of feature_metadata)."""
    path = folder / "feature_metadata.csv"
    if not path.exists():
        return {}, []
    with path.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header:
            return {}, []

        def find(pattern: str) -> list[int]:
            return [i for i, c in enumerate(header)
                    if c and re.search(pattern, c.strip(), re.I)
                    and not re.search(r"(^|_)host(_|$)", c.strip(), re.I)]

        king = find(r"(^|_)(realm|kingdom|domain|superkingdom)(_|$)")
        lineage = find(r"taxonomy|lineage|taxon_path|classification|full_lineage")
        name_i = find(r"^(cleaned_name|original_name|feature_id)$")
        rank_cols = {r: find(rf"(^|_)({r}|{r[0]})$" if r not in ("kingdom", "realm")
                             else r"(^|_)(realm|kingdom|domain|superkingdom)(_|$)")
                     for r in TAX_COLS}
        id_i = header.index("feature_id") if "feature_id" in header else 0
        clean_i = header.index("cleaned_name") if "cleaned_name" in header else id_i

        out = {}
        for row in reader:
            domain = None
            for group in (king, lineage):
                for i in group:
                    if i < len(row):
                        domain = stats.classify_domain(row[i])
                        if domain:
                            break
                if domain:
                    break
            if not domain:
                for i in name_i:
                    if i < len(row) and (";" in row[i] or "__" in row[i]):
                        domain = stats.classify_domain(row[i])
                        if domain:
                            break
            if not domain:
                domain = paper_domain
            if domain != "Viruses":
                continue

            tax = {}
            for rank, idx in rank_cols.items():
                for i in idx:
                    if i < len(row):
                        name = stats.clean_taxon(row[i])
                        if name and not re.match(r"^(unassigned|unclassified|unknown|"
                                                 r"norank|uncultured|na|none)$", name, re.I):
                            tax[stats.ictv_rank(name) or rank] = name
                            break
            fid = row[id_i] if id_i < len(row) else ""
            tax["cleaned_name"] = row[clean_i] if clean_i < len(row) else fid
            out[fid] = tax
    return out, header


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--figures-only", action="store_true",
                    help="rebuild figures and summary tables from the files already "
                         "written, without re-scanning the matrices")
    args = ap.parse_args()
    dest = args.out_dir
    dest.mkdir(parents=True, exist_ok=True)
    if args.figures_only:
        write_figures(dest)
        print(f"done -> {dest}")
        return 0

    # Inclusion is decided HERE, not inherited from the statistics cohort.
    # That cohort drops papers whose SAMPLES lack grouping labels - a rule made
    # for disease-association work. Applied to a viral sequence catalogue it
    # silently removed Garmaeva2024, the only large VLP virome in the database
    # (205 samples x 102,170 vOTUs). A viral matrix is usable here whenever it
    # contains viral features, whatever its sample labelling.
    # Derived here, not looked up: 01_papers_included.csv only covers the
    # statistics cohort, so the three papers it excludes came back blank.
    def row_rank_of(folder: str, flev: str) -> str:
        fp = ROOT / folder / "feature_metadata.csv"
        if fp.exists():
            with fp.open(encoding="utf-8-sig", newline="") as fh:
                rd = csv.reader(fh)
                hdr = next(rd, None) or []
                idx = [i for i, c in enumerate(hdr)
                       if c and c.strip().lower() in ("taxonomic_rank", "feature_level")]
                counts: collections.Counter = collections.Counter()
                for row in rd:
                    for i in idx:
                        if i < len(row) and row[i].strip():
                            counts[stats.bucket(row[i], stats.RANK_RULES,
                                                "Other / mixed")] += 1
                            break
                if counts:
                    return counts.most_common(1)[0][0]
        return stats.bucket(flev, stats.RANK_RULES, "Other / mixed")
    papers = []
    allp = stats.read_csv(stats.status_integrated)
    print(f"scanning {len(allp)} integrated papers for viral features ...")
    for i, r in enumerate(allp, 1):
        r["paper_folder"] = r["paper_folder"].strip()
        r["matrix_row_rank"] = row_rank_of(r["paper_folder"], r["feature_level"])
        pd = stats.paper_domain_of(f"{r['feature_level']} {r['main_matrix_unit']}")
        feats, _ = viral_features(ROOT / r["paper_folder"], pd)
        if feats:
            r["_feats"] = feats
            papers.append(r)
        if i % 50 == 0:
            print(f"  ...{i}/{len(allp)}", flush=True)
    papers.sort(key=lambda r: -len(r["_feats"]))
    print(f"{len(papers)} papers carry viral features")

    ab_fh = (dest / "viral_abundance_long.csv").open("w", encoding="utf-8-sig", newline="")
    ab = csv.writer(ab_fh)
    ab.writerow(["paper_folder", "sample_uid", "feature_uid", "value", "abundance_unit"])

    feat_rows, samp_rows, block_rows = [], [], []
    total_obs = 0
    collisions = 0

    for n, p in enumerate(papers, 1):
        name = p["paper_folder"]
        folder = ROOT / name
        feats = p["_feats"]

        unit = p["main_matrix_unit"]
        # "Was this matrix built to measure viruses?" A matrix whose features are
        # ALL viral, or whose row unit is a viral one, is viral-dedicated; viral
        # rows inside a bacterial profile are incidental. Deciding this from the
        # data beats a hand-written list, which had filed Garmaeva2024 (102,170
        # of 102,170 features viral) as incidental.
        n_total = as_int(p.get("n_features"), len(feats))
        viral_rank = bool(re.search(r"votu|viral|phage", p["matrix_row_rank"], re.I))
        viral_matrix = viral_rank or (n_total and len(feats) == n_total)
        p["viral_matrix_type"] = "viral-dedicated" if viral_matrix else "incidental"
        recovery = RECOVERY.get(name, "VLP/virome pipeline (see paper)"
                                if viral_matrix else DEFAULT_RECOVERY)

        for fid, tax in feats.items():
            if SEP in fid:
                collisions += 1
            feat_rows.append({
                "feature_uid": f"{name}{SEP}{fid}", "paper_folder": name,
                "feature_id_original": fid, "cleaned_name": tax.get("cleaned_name", ""),
                "matrix_row_rank": p["matrix_row_rank"],
                "viral_recovery_method": recovery,
                **{c: tax.get(c, "") for c in TAX_COLS},
            })

        # sample metadata
        smeta = {}
        spath = folder / "sample_metadata.csv"
        if spath.exists():
            for row in stats.read_csv(spath):
                sid = (row.get("sample_id") or "").strip()
                if sid:
                    smeta[sid] = row

        # stream the matrix, keeping only viral columns and non-zero values
        mpath = folder / "abundance_matrix.csv"
        if not mpath.exists():
            print(f"  ! {name}: no abundance_matrix.csv, skipped", file=sys.stderr)
            continue
        n_obs = n_blank = n_zero = 0
        with mpath.open(encoding="utf-8-sig", newline="") as fh:
            rd = csv.reader(fh)
            header = next(rd)
            keep = [(i, header[i]) for i in range(1, len(header)) if header[i] in feats]
            seen_samples = set()
            for row in rd:
                sid = row[0].strip()
                seen_samples.add(sid)
                suid = f"{name}{SEP}{sid}"
                for i, fname in keep:
                    if i >= len(row):
                        continue
                    v = (row[i] or "").strip()
                    if not v:
                        n_blank += 1
                        continue
                    try:
                        if float(v) == 0:
                            n_zero += 1
                            continue
                    except ValueError:
                        pass
                    ab.writerow([name, suid, f"{name}{SEP}{fname}", v, unit])
                    n_obs += 1

            for sid in sorted(seen_samples):
                src = smeta.get(sid, {})
                level, can, cannot = usability(src)
                disease_resolved, extra_labels = resolve_labels(src)
                samp_rows.append({
                    "sample_uid": f"{name}{SEP}{sid}", "paper_folder": name,
                    "sample_id_original": sid, "viral_recovery_method": recovery,
                    "abundance_unit": unit,
                    "metadata_completeness": level,
                    "usable_for": can, "not_usable_for": cannot,
                    "disease_resolved": disease_resolved,
                    "source_label_fields": extra_labels,
                    **{c: (src.get(c) or "") for c in SAMPLE_COLS if c != "abundance_unit"},
                })

            levels = collections.Counter(
                usability(smeta.get(sid, {}))[0] for sid in seen_samples)
            block_rows.append({
                "paper_folder": name, "n_samples": len(seen_samples),
                "metadata_completeness": "; ".join(
                    f"{k}:{v}" for k, v in levels.most_common()),
                "n_viral_features": len(feats),
                "n_cells_in_block": len(seen_samples) * len(feats),
                "n_nonzero_written": n_obs, "n_explicit_zero": n_zero,
                "n_blank_unreported": n_blank,
                "absent_cell_means": ("explicit zero" if n_blank == 0
                                      else "zero OR unreported - see the paper's "
                                           "missing_value_reasons.csv"),
                "abundance_unit": unit, "matrix_row_rank": p["matrix_row_rank"],
                "viral_matrix_type": p["viral_matrix_type"],
                "n_features_total": n_total,
                "viral_share_of_matrix": stats.pct(len(feats), n_total),
                "viral_recovery_method": recovery, "doi": p["doi"],
            })
        total_obs += n_obs
        print(f"  [{n:>2}/{len(papers)}] {name:<20} {len(seen_samples):>5} samples x "
              f"{len(feats):>6} viral features -> {n_obs:>9,} non-zero", flush=True)

    ab_fh.close()
    if collisions:
        print(f"  ! {collisions} feature ids already contain {SEP!r}", file=sys.stderr)

    for fname, rows in [("viral_feature_metadata.csv", feat_rows),
                        ("viral_sample_metadata.csv", samp_rows),
                        ("viral_blocks.csv", block_rows)]:
        with (dest / fname).open("w", encoding="utf-8-sig", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"  wrote {fname}: {len(rows):,} rows")

    write_figures(dest)
    print(f"\ndone -> {dest}")
    print(f"  {len(block_rows)} papers | {len(samp_rows):,} samples | "
          f"{len(feat_rows):,} viral features | {total_obs:,} non-zero observations")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
