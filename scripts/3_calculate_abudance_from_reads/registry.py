#!/usr/bin/env python3
"""The per-paper registry: what finished, what failed, and what the risks are.

pipeline_state.csv stays as it is - one row per paper x step, the audit trail.
Everything here is derived from it, one row per PAPER, so a reader does not have
to reassemble three rows to learn what happened to one paper.

The point of the risk files is that a paper can be perfectly usable for one
thing and unusable for another. A matrix whose group labels were never
confirmed is still a correct matrix; it just must not be fed into a case/control
test. Rather than silently dropping such a paper or silently including it, each
limitation gets its own file, named for the limitation, so whoever uses the
database can see exactly what they are taking on.

Files written into <out_root>:

    0_counts.csv                  the row counts, and whether they balance

    1_completed.csv               produced at least one verified matrix
    1_pipeline_failed.csv         something broke; which step and why
    1_method_missing.csv          no implemented pipeline for this assay
      these three are disjoint and add up to the paper count

    methods/2_method_<name>.csv   one file per method; every completed paper is
                                  in exactly one, so they add up to 1_completed

    3_risk_no_grouping.csv        quantification is sound, group labels are NOT
    3_risk_noisy_rare_asvs.csv    chimera removal took an abnormal share of reads
    3_risk_partial_branches.csv   one branch produced output, another never ran
    3_risk_none.csv               no limitation fired
      a paper lands in exactly one, by severity, so these add up to 1_completed.
      `all_risks` lists every risk it has, for the papers that have several.

A paper can appear in `pipeline_completed.csv` AND in a risk file. That is the
intended behaviour: it finished, and here is what is wrong with it.
"""

import csv
import json
import os
from pathlib import Path

STEP_NAMES = {"1": "download", "2": "metadata", "3": "abundance"}

#: status -> (which bucket, human wording)
DONE_STATUSES = {
    "done": "done",
    "done_quantification_only": "done (quantification only, groups unconfirmed)",
    "done_partial": "done (some branches did not run)",
    "preflight_ok": "preflight only - NOT a finished analysis",
}
FAIL_STATUSES = {"failed", "metadata_blocked", "no_space", "skipped_unsuitable"}

#: method name -> (script, what it produces, which reference)
METHODS = {
    "16s_paired": ("step3_16s.R", "DADA2 ASVs from paired reads", "SILVA 138.2"),
    "16s_single": ("step3_16s_single.R", "DADA2 ASVs from single-end reads",
                   "SILVA 138.2"),
    "shotgun_bacteria": ("step3_shotgun.py", "sylph containment-ANI profile",
                         "GTDB sketch"),
    "viral_mapping": ("step3_virome.py", "read mapping + coverage",
                      "viral catalogue"),
}


def write_csv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


def plan_of(out_root, paper):
    try:
        return json.loads((out_root / paper / "analysis_plan.json")
                          .read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def gate_of(out_root, paper):
    try:
        return (out_root / paper / "metadata_gate.txt").read_text(
            encoding="utf-8").strip()
    except OSError:
        return ""


def chimera_retention(out_root, paper):
    """Fraction of merged reads surviving chimera removal, or None."""
    p = out_root / paper / "dada2" / "read_counts_per_step.csv"
    if not p.exists():
        return None
    try:
        with p.open(encoding="utf-8-sig", newline="") as fh:
            rows = list(csv.DictReader(fh))
        merged = sum(float(r["merged"]) for r in rows)
        nochim = sum(float(r["nonchimeric"]) for r in rows)
        return nochim / merged if merged else None
    except (OSError, ValueError, KeyError):
        return None


def matrix_shape(out_root, paper, sub, name):
    p = out_root / paper / sub / name
    if not p.exists():
        return None
    try:
        with p.open(encoding="utf-8-sig", newline="") as fh:
            header = next(csv.reader(fh))
            return sum(1 for _ in fh), len(header) - 1
    except (OSError, StopIteration):
        return None


def rebuild(out_root, state_rows):
    """Rewrite every derived file from pipeline_state.csv.

    `state_rows` is the full list of per-(paper, step) rows.
    """
    out_root = Path(out_root)
    by_paper = {}
    for r in state_rows:
        by_paper.setdefault(r["paper_folder"], {})[r["step"]] = r

    completed, failed, in_progress = [], [], []
    no_group, noisy, partial, missing = [], [], [], []
    per_method = {m: [] for m in METHODS}

    for paper, steps in sorted(by_paper.items(), key=lambda kv: kv[0].lower()):
        latest = max(steps.values(), key=lambda r: r.get("updated_utc", ""))
        meta = {"paper_folder": paper, "doi": latest.get("doi", ""),
                "title": latest.get("title", ""),
                "updated_utc": latest.get("updated_utc", "")}
        plan = plan_of(out_root, paper)
        gate = gate_of(out_root, paper)

        bad = {s: r for s, r in steps.items() if r["status"] in FAIL_STATUSES}
        ok = {s: r for s, r in steps.items() if r["status"] in DONE_STATUSES}

        if bad:
            step = sorted(bad)[0]
            failed.append({**meta, "failed_at_step": f"{step} ({STEP_NAMES.get(step, step)})",
                           "status": bad[step]["status"],
                           "steps_done": ",".join(sorted(ok)) or "none",
                           "reason": bad[step].get("reason", "")[:500]})
        step3 = steps.get("3", {}).get("status")
        if ok and step3 not in DONE_STATUSES:
            # Steps finished, but no matrix yet. Not completed, not failed.
            in_progress.append({**meta, "steps_done": ",".join(sorted(ok)),
                                "next_step": "3 (abundance)" if "2" in ok else
                                             "2 (metadata)",
                                "why": "the pipeline has not reached a finished "
                                       "abundance matrix for this paper yet"})
        elif ok:
            shape = (matrix_shape(out_root, paper, "dada2", "asv_counts.csv")
                     or matrix_shape(out_root, paper, "shotgun",
                                     "bacteria_taxonomic_abundance_pct.csv")
                     or matrix_shape(out_root, paper, "virome", "viral_mean_depth.csv"))
            completed.append({
                **meta,
                "steps_done": ",".join(sorted(ok)),
                "all_three_done": "yes" if set(ok) >= {"1", "2", "3"} else "no",
                "step3_status": steps.get("3", {}).get("status", "not run"),
                "n_samples": shape[0] if shape else "",
                "n_features": shape[1] if shape else "",
                "branches": ";".join(f"{k}={v}" for k, v in
                                     (plan.get("branch_status") or {}).items()),
                "result_path": steps.get("3", {}).get("result_path", ""),
            })

        # ---- risks, independent of success ---------------------------------
        if steps.get("3", {}).get("status") in DONE_STATUSES:
            if gate.startswith("pass_quantification_only") or \
                    plan.get("grouped_analysis_supported") is False:
                no_group.append({**meta,
                    "usable_for": "per-sample abundance; within-paper comparisons by "
                                  "any field that IS confirmed (timepoint, body site)",
                    "not_usable_for": "any grouped comparison: case/control, treatment arm",
                    "why": gate[:400] or "analysis_plan.json: grouped_analysis_supported=false"})
            chim = chimera_retention(out_root, paper)
            if chim is not None and chim < 0.80:
                shape = matrix_shape(out_root, paper, "dada2", "asv_counts.csv")
                noisy.append({**meta,
                    "chimera_read_retention": f"{chim:.1%}",
                    "n_asvs": shape[1] if shape else "",
                    "usable_for": "the abundant fraction of the community",
                    "not_usable_for": "rare/low-prevalence ASVs taken at face value",
                    "why": f"chimera removal discarded {1-chim:.1%} of merged reads "
                           "(a healthy run keeps >90%). Filter by prevalence before "
                           "using the rare tail."})
            deferred = {b: s for b, s in (plan.get("branch_status") or {}).items()
                        if s in ("missing_dependencies", "unsupported")}
            if deferred:
                partial.append({**meta,
                    "branches_not_run": ";".join(f"{b}={s}" for b, s in deferred.items()),
                    "why": "this paper has more than one analysable branch and at "
                           "least one of them did not run"})

        # ---- which method ran, and which assay has none ---------------------
        bs = dict(plan.get("branch_status") or {})
        if not bs and step3 in DONE_STATUSES:
            # An older run wrote no branch_status. Infer it from the outputs that
            # exist, so the method files still account for every completed paper.
            for sub, name in (("dada2", "amplicon_16s"),
                              ("dada2_single", "amplicon_16s_single"),
                              ("shotgun", "dna_shotgun_bacteria"),
                              ("virome", "viral_mapping")):
                if (out_root / paper / sub / "README.txt").exists():
                    bs[name] = "done"
        for b, st in bs.items():
            name = {"amplicon_16s": "16s_paired", "amplicon_16s_single": "16s_single",
                    "dna_shotgun_bacteria": "shotgun_bacteria",
                    "viral_mapping": "viral_mapping"}.get(b)
            if name and st == "done":
                script, produces, ref = METHODS[name]
                sub = {"16s_paired": ("dada2", "asv_counts.csv"),
                       "16s_single": ("dada2_single", "asv_counts.csv"),
                       "shotgun_bacteria": ("shotgun", "bacteria_taxonomic_abundance_pct.csv"),
                       "viral_mapping": ("virome", "viral_mean_depth.csv")}[name]
                shape = matrix_shape(out_root, paper, *sub)
                per_method[name].append({
                    **meta, "method": name, "script": script,
                    "produces": produces, "reference": ref,
                    "region_or_catalogue": plan.get("amplicon_region")
                                           or plan.get("viral_catalogue") or "",
                    "primer_evidence": (plan.get("paper_facts") or {}).get(
                        "amplicon_status", ""),
                    "n_samples": shape[0] if shape else "",
                    "n_features": shape[1] if shape else "",
                    "output_dir": str(out_root / paper / sub[0]),
                })
            elif st == "unsupported":
                missing.append({**meta, "assay_branch": b,
                    "why": "no implemented pipeline for this assay",
                    "what_is_needed": NEEDED.get(b, "a pipeline for this assay type")})

    # ---- level 1: every paper is in exactly one of these ------------------
    done_set = {r["paper_folder"] for r in completed}
    fail_set = {r["paper_folder"] for r in failed}
    miss_set = {r["paper_folder"] for r in missing}
    # A paper that produced something is 'completed' even if another branch had
    # no method; a paper that failed is not 'missing a method'. Resolving the
    # overlaps here is what makes the three files add up to the paper count.
    prog_set = {r["paper_folder"] for r in in_progress}
    # "No implemented method" is not a failure. A paper whose step 3 stopped
    # solely because its assay has no pipeline belongs in 1_method_missing;
    # filing it under 1_pipeline_failed would imply something broke and would
    # hide it from anyone counting what the pipeline still cannot do.
    only_unsupported = {
        r["paper_folder"] for r in missing
        if all(st == "unsupported" for st in
               (plan_of(out_root, r["paper_folder"]).get("branch_status") or
                {"x": "unsupported"}).values())}
    fail_set -= done_set | only_unsupported
    prog_set -= done_set | fail_set | only_unsupported
    miss_set -= done_set | prog_set
    miss_set &= (miss_set - fail_set) | only_unsupported
    in_progress = [r for r in in_progress if r["paper_folder"] in prog_set]
    missing = [r for r in missing if r["paper_folder"] in miss_set]
    failed = [r for r in failed if r["paper_folder"] in fail_set]
    seen_miss = set()
    missing = [r for r in missing
               if not (r["paper_folder"] in seen_miss or seen_miss.add(r["paper_folder"]))]

    write_csv(out_root / "1_completed.csv", completed,
              ["paper_folder", "doi", "title", "steps_done", "all_three_done",
               "step3_status", "n_samples", "n_features", "branches",
               "result_path", "updated_utc"])
    write_csv(out_root / "1_pipeline_failed.csv", failed,
              ["paper_folder", "doi", "title", "failed_at_step", "status",
               "steps_done", "reason", "updated_utc"])
    write_csv(out_root / "1_in_progress.csv", in_progress,
              ["paper_folder", "doi", "title", "steps_done", "next_step", "why",
               "updated_utc"])
    write_csv(out_root / "1_method_missing.csv", missing,
              ["paper_folder", "doi", "title", "assay_branch", "why",
               "what_is_needed", "updated_utc"])

    # ---- level 2: which method ran. One file per method, and every completed
    # paper appears in exactly one of them. -------------------------------
    for name, rows in per_method.items():
        write_csv(out_root / "methods" / f"2_method_{name}.csv", rows,
                  ["paper_folder", "doi", "title", "method", "script", "produces",
                   "reference", "region_or_catalogue", "primer_evidence",
                   "n_samples", "n_features", "output_dir", "updated_utc"])

    # ---- level 3: risks. Each completed paper lands in exactly ONE file, by
    # severity, so the files still add up; `all_risks` keeps the full picture
    # for a paper that has more than one. ----------------------------------
    risk_of = {}
    for bucket, rows in (("no_grouping", no_group), ("noisy_rare_asvs", noisy),
                         ("partial_branches", partial)):
        for r in rows:
            risk_of.setdefault(r["paper_folder"], []).append((bucket, r))
    ORDER = ["no_grouping", "noisy_rare_asvs", "partial_branches"]
    files = {b: [] for b in ORDER}
    clean = []
    for r in completed:
        got = risk_of.get(r["paper_folder"], [])
        if not got:
            clean.append({**{k: r[k] for k in ("paper_folder", "doi", "title",
                                               "updated_utc")},
                          "all_risks": "none",
                          "note": "no limitation recorded by the pipeline. That is not "
                                  "a guarantee the paper is free of problems, only "
                                  "that none of the checks fired."})
            continue
        names = [b for b, _ in got]
        primary = min(names, key=ORDER.index)
        row = dict(next(x for b, x in got if b == primary))
        row["all_risks"] = ";".join(sorted(names, key=ORDER.index))
        files[primary].append(row)

    write_csv(out_root / "3_risk_no_grouping.csv", files["no_grouping"],
              ["paper_folder", "doi", "title", "all_risks", "usable_for",
               "not_usable_for", "why", "updated_utc"])
    write_csv(out_root / "3_risk_noisy_rare_asvs.csv", files["noisy_rare_asvs"],
              ["paper_folder", "doi", "title", "all_risks", "chimera_read_retention",
               "n_asvs", "usable_for", "not_usable_for", "why", "updated_utc"])
    write_csv(out_root / "3_risk_partial_branches.csv", files["partial_branches"],
              ["paper_folder", "doi", "title", "all_risks", "branches_not_run",
               "why", "updated_utc"])
    write_csv(out_root / "3_risk_none.csv", clean,
              ["paper_folder", "doi", "title", "all_risks", "note", "updated_utc"])

    # ---- the identities, checked rather than asserted in a comment ---------
    n_method = sum(len(v) for v in per_method.values())
    n_risk = sum(len(v) for v in files.values()) + len(clean)
    counts = {
        "1_completed": len(completed), "1_pipeline_failed": len(failed),
        "1_method_missing": len(missing), "1_in_progress": len(in_progress),
        "level1_total": len(completed) + len(failed) + len(missing) + len(in_progress),
        "papers_in_state": len(by_paper),
        "2_method_sum": n_method,
        "3_risk_sum": n_risk,
    }
    counts["level1_balances"] = counts["level1_total"] == counts["papers_in_state"]
    counts["method_balances"] = n_method == len(completed)
    counts["risk_balances"] = n_risk == len(completed)
    write_csv(out_root / "0_counts.csv",
              [{"name": k, "value": v} for k, v in counts.items()],
              ["name", "value"])
    return counts

NEEDED = {
    "rna": "an RNA/metatranscriptome pipeline: rRNA depletion handling plus a "
           "transcript reference or assembly. RNA abundance is a different "
           "measurement from DNA abundance and must be stored separately.",
    "amplicon_its": "an ITS pipeline with the UNITE reference. ITS length varies "
                    "far more than 16S, so the merge-budget logic here does not "
                    "transfer.",
    "amplicon_18s": "an 18S pipeline with PR2 or SILVA-18S.",
    "amplicon_other": "identification of which marker this is, then a pipeline and "
                      "reference for it.",
    "unsupported_isolate_genome": "nothing: an isolate or host genome is not a "
                                  "community and has no abundance matrix.",
}


if __name__ == "__main__":
    # Rebuilding the derived views should not require running the pipeline:
    # after an edit to pipeline_state.csv the per-paper files must be able to
    # catch up on their own, or they quietly disagree with the audit trail.
    import argparse
    ap = argparse.ArgumentParser(description="Rebuild the per-paper views from "
                                             "pipeline_state.csv.")
    ap.add_argument("--out-root", default="/hickory/proj/didonglab/dataset/"
                                          "virus/Database/data/reads_abundance")
    a = ap.parse_args()
    root = Path(a.out_root)
    with (root / "pipeline_state.csv").open(encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))
    rebuild(root, rows)
    print(f"rebuilt from {len(rows)} state row(s) -> {root}")
