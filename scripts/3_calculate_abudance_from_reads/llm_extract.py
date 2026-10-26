#!/usr/bin/env python3
"""Run one headless model session to resolve one paper's sample metadata.

Same connection the paper-integration batch already uses: the Claude Code CLI in
`-p` (headless) mode with an explicit tool allow-list and --add-dir, or the Codex
CLI as an alternative. Nothing here is a placeholder - step 2 really calls it.

Lifecycle, which is the part that has to be right:

    start    the CLI is launched in its OWN process group (start_new_session),
             so nothing we do can reach another Claude or Codex session that
             happens to be running on the same node
    wait     a wall-clock timeout, from the caller
    verify   the session's own exit code is not trusted. The structured file it
             was asked to write is parsed and validated; a session that exits 0
             without a usable file is a failure
    end      on timeout, SIGINT/SIGTERM or a crash, the whole process GROUP is
             signalled (SIGTERM, then SIGKILL) so no child survives the parent

A failed or unparsable call leaves `llm_status.json` behind saying what happened
and how to resume. It never writes a result that implies success.

Credentials come from the CLI's own login/environment. Nothing in here reads,
prints or stores a key.
"""

import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROMPT = HERE / "prompt_step2.md"
CONFIDENCE = {"confirmed", "inferred", "unknown"}
# ---------------------------------------------------------------- field spec
# The ONE definition of the per-sample fields. The prompt's tables and JSON
# skeleton are generated from this tuple, and step2_metadata.py imports it for
# its CSV columns, so the three cannot drift apart.
#
# Confidence is per field, not per sample: a sample id can be certain while the
# group beside it is only a guess, and the gate has to see the difference.
#
# required=True means the entry must be present in the answer (its confidence
# may still be "unknown"), because the gate reads it. Everything else is
# extracted when the paper supports it and omitted when it does not: a
# covariate the paper never reports is a fact worth recording, an invented one
# silently corrupts every analysis built on it.
FIELD_SPEC = (
    ("paper_sample_id", True,
     "the identifier the PAPER uses for this sample"),
    ("experimental_group", True,
     "the arm/condition this sample belongs to, named as the paper names it"),
    ("disease", False,
     'the named condition, e.g. "ulcerative colitis", "type 2 diabetes". '
     "`none` for a healthy subject, `not_applicable` for a blank or mock"),
    ("disease_state", False,
     "one of `case`, `control`, `healthy`, `not_applicable`"),
    ("treatment", False,
     "the intervention applied to THIS sample: drug, diet, antibiotic, "
     "probiotic, FMT, with dose/duration if stated"),
    ("timepoint", False,
     "day/week/visit, as the paper labels it"),
    ("subject_id", False,
     "the individual (patient, animal, donor) this sample came from. Two "
     "samples from one subject are not independent, so this is what makes "
     "repeated measures visible"),
    ("replicate_of", False,
     "if this is a technical replicate, the sample it replicates"),
    ("sample_type", False,
     "WHAT THE MATERIAL IS: stool, biopsy, swab, aspirate, lavage, isolate, "
     "culture. Not where it came from"),
    ("body_site", False,
     "WHERE IT CAME FROM, anatomically: colon, cecum, tongue dorsum, nasal "
     "cavity, skin. Not what the material is. (The older name `tissue` is "
     "accepted as a synonym and merged into this field.)"),
    ("host_species", False,
     "the host organism, if there is one"),
    ("sex", False, "as reported"),
    ("age", False, "with its unit, as reported"),
    ("geographic_location", False, "country/region of collection"),
    ("collection_date", False, "as reported"),
)
REQUIRED_SAMPLE_FIELDS = tuple(n for n, req, _ in FIELD_SPEC if req)
ANALYSIS_SAMPLE_FIELDS = tuple(n for n, req, _ in FIELD_SPEC if not req)
SAMPLE_FIELDS = tuple(n for n, _, _ in FIELD_SPEC)
# Names an answer may use for a field we call something else.
FIELD_ALIASES = {"tissue": "body_site", "bodysite": "body_site",
                 "host": "host_species", "group": "experimental_group"}


def field_table(required):
    """The prompt's markdown table for one half of the spec."""
    rows = ["| field | meaning |", "| --- | --- |"]
    rows += [f"| `{n}` | {d} |" for n, req, d in FIELD_SPEC if req is required]
    return "\n".join(rows)


def fields_skeleton(names, indent=8):
    """The prompt's JSON skeleton, so the example can never omit a field."""
    pad, width = " " * indent, max(len(n) for n in names) + 3
    return "\n".join(
        f'{pad}{(chr(34) + n + chr(34) + ":").ljust(width)} '
        '{"value": "", "confidence": "confirmed|inferred|unknown", "evidence": ""}'
        + ("," if i < len(names) - 1 else "")
        for i, n in enumerate(names))


def canonical(fields):
    """Rename any accepted synonym to the canonical field name."""
    out = {}
    for key, value in (fields or {}).items():
        out.setdefault(FIELD_ALIASES.get(str(key).strip().lower(), key), value)
    return out


def resolve_fields(sample, defaults):
    """Merge a sample's own fields over the paper-wide defaults.

    Most covariates are constant across a paper (host, country, often disease),
    so the answer may state them once. A per-sample entry always wins.
    """
    own = canonical(sample.get("fields"))
    defaults = canonical(defaults)
    out = {}
    for name in SAMPLE_FIELDS:
        entry = own.get(name)
        if not isinstance(entry, dict) or not str(entry.get("value") or "").strip():
            fallback = defaults.get(name)
            if isinstance(fallback, dict) and str(fallback.get("value") or "").strip():
                entry = {**fallback,
                         "evidence": (fallback.get("evidence", "") +
                                      " [paper-wide default]").strip()}
        out[name] = entry if isinstance(entry, dict) else {
            "value": "", "confidence": "unknown", "evidence": ""}
    return out


# --------------------------------------------------------------------------- CLI
def find_cli(agent):
    """Locate the CLI, preferring a real install over the VS Code extension copy."""
    if agent == "claude":
        found = shutil.which("claude")
        if found:
            return found
        cands = sorted(Path.home().glob(
            ".vscode-server/extensions/anthropic.claude-code-*/resources/native-binary/claude"))
    elif agent == "codex":
        found = shutil.which("codex")
        if found:
            return found
        cands = sorted(Path.home().glob(
            ".vscode-server/extensions/*/bin/linux-x86_64/codex"))
    else:
        raise ValueError(f"unknown agent {agent!r}")
    return str(cands[-1]) if cands else None


def build_command(agent, cli, prompt, work_dir, extra_dirs, model, max_turns):
    """The argv for one unattended session. Tools are allow-listed, not open."""
    if agent == "claude":
        cmd = [cli, "-p", prompt, "--output-format", "stream-json", "--verbose",
               "--max-turns", str(max_turns),
               "--allowedTools", "Bash Read Write Edit Glob Grep WebFetch WebSearch",
               "--disallowedTools", "Agent"]
        if model:
            cmd += ["--model", model]
        cmd += ["--add-dir", str(work_dir)] + [a for d in extra_dirs
                                               for a in ("--add-dir", str(d))]
        return cmd, None
    cmd = [cli, "--ask-for-approval", "never", "--search", "exec",
           "--sandbox", "workspace-write",
           "-c", "sandbox_workspace_write.network_access=true",
           "--skip-git-repo-check", "--cd", str(work_dir)]
    if model:
        cmd += ["--model", model]
    cmd += [a for d in extra_dirs for a in ("--add-dir", str(d))]
    cmd += ["--json", "-"]
    return cmd, prompt          # Codex takes the prompt on stdin


# ------------------------------------------------------------------- run / stop
def _stop(proc, log):
    """End this session's whole process group, and only this one."""
    if proc.poll() is not None:
        return
    try:
        pgid = os.getpgid(proc.pid)
    except ProcessLookupError:
        return
    for sig, wait in ((signal.SIGTERM, 10), (signal.SIGKILL, 5)):
        try:
            os.killpg(pgid, sig)        # the group we created, never a shared one
        except ProcessLookupError:
            return
        try:
            proc.wait(timeout=wait)
            log(f"session ended after {sig.name}")
            return
        except subprocess.TimeoutExpired:
            continue
    log("session did not die even after SIGKILL")


def run_session(cmd, stdin_text, timeout, jsonl, errlog, log):
    """Launch, wait, and always leave the group dead. Returns (rc, reason)."""
    interrupted = {}

    def on_signal(signum, _frame):
        interrupted["sig"] = signal.Signals(signum).name

    previous = {s: signal.getsignal(s) for s in (signal.SIGINT, signal.SIGTERM)}
    for s in previous:
        signal.signal(s, on_signal)
    proc = None
    try:
        with jsonl.open("w") as out, errlog.open("w") as err:
            proc = subprocess.Popen(
                cmd, stdin=subprocess.PIPE if stdin_text else subprocess.DEVNULL,
                stdout=out, stderr=err, text=True,
                start_new_session=True)      # its own process group
            if stdin_text:
                try:
                    proc.stdin.write(stdin_text)
                finally:
                    proc.stdin.close()
            try:
                rc = proc.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                log(f"timeout after {timeout}s")
                _stop(proc, log)
                return None, f"timeout after {timeout}s"
        if interrupted:
            _stop(proc, log)
            return None, f"interrupted by {interrupted['sig']}"
        return rc, ("" if rc == 0 else f"{cmd[0].rsplit('/', 1)[-1]} exited {rc}")
    except Exception as exc:                 # noqa: BLE001 - reported, not hidden
        if proc:
            _stop(proc, log)
        return None, f"could not run the session: {exc}"
    finally:
        for s, h in previous.items():
            signal.signal(s, h)
        if proc:
            _stop(proc, log)


# ------------------------------------------------------------------- validation
def validate(path, expected_biosamples):
    """Parse and sanity-check the model's file. Returns (data, problems)."""
    if not path.exists():
        return None, [f"the session wrote no {path.name}"]
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:                  # noqa: BLE001
        return None, [f"{path.name} is not valid JSON: {exc}"]

    bad = []
    for key in ("sample_id_scheme", "key_experiment", "samples", "unresolved",
                "sources_consulted"):
        if key not in data:
            bad.append(f"missing required key {key!r}")
    if bad:
        return data, bad

    samples = data.get("samples") or []
    if not isinstance(samples, list) or not samples:
        return data, ["'samples' is empty"]

    seen = set()
    for i, s in enumerate(samples):
        acc = (s.get("biosample_accession") or "").strip()
        if not acc:
            bad.append(f"sample {i} has no biosample_accession")
            continue
        seen.add(acc)
        if not isinstance(s.get("fields"), dict):
            bad.append(f"{acc}: no 'fields' object")
            continue
        fields = resolve_fields(s, data.get("sample_defaults"))
        for name in SAMPLE_FIELDS:
            f = fields[name]
            if name in REQUIRED_SAMPLE_FIELDS and not isinstance(
                    (s.get("fields") or {}).get(name), dict) and not isinstance(
                    (data.get("sample_defaults") or {}).get(name), dict):
                bad.append(f"{acc}: required field {name!r} is missing")
                continue
            conf = (f.get("confidence") or "").strip()
            if conf not in CONFIDENCE:
                bad.append(f"{acc}.{name}: confidence {conf!r} is not one of "
                           f"{sorted(CONFIDENCE)}")
            # A claim of fact without a source is not a fact.
            elif conf == "confirmed" and not (f.get("evidence") or "").strip():
                bad.append(f"{acc}.{name}: marked 'confirmed' with no evidence")
            elif conf == "confirmed" and not str(f.get("value") or "").strip():
                bad.append(f"{acc}.{name}: marked 'confirmed' but the value is empty")
    missing = expected_biosamples - seen
    extra = seen - expected_biosamples
    if missing:
        bad.append(f"{len(missing)} downloaded BioSample(s) absent from the answer: "
                   + ", ".join(sorted(missing)[:6]))
    if extra:
        bad.append(f"{len(extra)} BioSample(s) that were not downloaded: "
                   + ", ".join(sorted(extra)[:6]))
    for key in ("sample_id_scheme", "key_experiment"):
        if (data[key].get("confidence") or "") not in CONFIDENCE:
            bad.append(f"{key}.confidence is not one of {sorted(CONFIDENCE)}")
    return data, bad


def final_text(jsonl, agent):
    """The session's closing message, for the log. Never parsed as the answer."""
    try:
        if agent == "codex":
            return jsonl.read_text(encoding="utf-8", errors="replace")[-800:]
        for line in reversed(jsonl.read_text(encoding="utf-8",
                                             errors="replace").splitlines()):
            try:
                msg = json.loads(line)
            except ValueError:
                continue
            if msg.get("type") == "result":
                return (f"[turns={msg.get('num_turns')} "
                        f"cost_usd={msg.get('total_cost_usd')}] "
                        f"{str(msg.get('result', ''))[:600]}")
    except OSError:
        pass
    return ""


# ------------------------------------------------------------------------ entry
def extract(paper, doi, title, work_dir, reads_dir, orig_dir, paper_facts,
            expected, *, agent="claude", model="", timeout=2400, retries=1,
            max_turns=120, log=print):
    """Resolve one paper's sample metadata. Returns (data_or_None, status dict)."""
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    out_json = work_dir / "metadata_extraction.json"
    status_path = work_dir / "llm_status.json"

    cli = find_cli(agent)
    if not cli:
        st = {"ok": False, "reason": f"{agent} CLI not found on this node",
              "resume": f"install {agent} or set the agent in submit_specific_steps.sh"}
        status_path.write_text(json.dumps(st, indent=2), encoding="utf-8")
        return None, st

    expected = set(expected)
    text = PROMPT.read_text(encoding="utf-8")
    # The prompt's field tables and JSON skeleton are rendered from FIELD_SPEC,
    # so the instructions and the validator are the same list by construction.
    for key, val in (("{PAPER}", paper), ("{DOI}", doi), ("{TITLE}", title),
                     ("{READS_DIR}", str(reads_dir)), ("{ORIG_DIR}", str(orig_dir)),
                     ("{WORK_DIR}", str(work_dir)), ("{OUT_JSON}", str(out_json)),
                     ("{N_SAMPLES}", str(len(expected))),
                     ("{FIELD_TABLE_REQUIRED}", field_table(True)),
                     ("{FIELD_TABLE_ANALYSIS}", field_table(False)),
                     ("{FIELDS_JSON}", fields_skeleton(SAMPLE_FIELDS)),
                     ("{DEFAULTS_JSON}", fields_skeleton(ANALYSIS_SAMPLE_FIELDS, 4)),
                     ("{PAPER_FACTS}", str(paper_facts) if paper_facts
                      else "(none on file for this paper)")):
        text = text.replace(key, val)
    (work_dir / "prompt_sent.md").write_text(text, encoding="utf-8")

    extra = [d for d in (reads_dir, orig_dir,
                         paper_facts.parent if paper_facts else None)
             if d and Path(d).exists()]
    cmd, stdin_text = build_command(agent, cli, text, work_dir, extra, model, max_turns)

    attempts = []
    for attempt in range(1, retries + 2):
        jsonl = work_dir / f"llm_session_{attempt}.jsonl"
        errlog = work_dir / f"llm_session_{attempt}.err"
        if out_json.exists():
            out_json.unlink()            # never validate a previous attempt's file
        started = time.time()
        log(f"[{paper}] {agent} session {attempt}/{retries + 1} "
            f"(model={model or 'CLI default'}, timeout={timeout}s)")
        rc, reason = run_session(cmd, stdin_text, timeout, jsonl, errlog, log)
        took = round(time.time() - started, 1)
        data, problems = validate(out_json, expected)
        note = final_text(jsonl, agent)
        attempts.append({"attempt": attempt, "seconds": took, "exit_code": rc,
                         "run_problem": reason, "output_problems": problems,
                         "session_log": str(jsonl), "final_message": note})
        if data is not None and not problems:
            st = {"ok": True, "agent": agent, "model": model or "CLI default",
                  "cli": cli, "seconds": took, "attempts": attempts,
                  "result": str(out_json),
                  "note": "the session's exit code alone was not treated as success; "
                          "this file was parsed and validated"}
            status_path.write_text(json.dumps(st, indent=2, ensure_ascii=False),
                                   encoding="utf-8")
            log(f"[{paper}] extraction validated: {len(data['samples'])} samples, "
                f"{took}s")
            return data, st
        log(f"[{paper}] attempt {attempt} unusable: "
            f"{reason or '; '.join(problems)[:200]}")
        if "interrupted by" in (reason or ""):
            break                        # the operator stopped us; do not retry

    st = {"ok": False, "agent": agent, "model": model or "CLI default", "cli": cli,
          "reason": attempts[-1]["run_problem"] or "; ".join(attempts[-1]["output_problems"]),
          "attempts": attempts,
          "resume": f"inspect {work_dir}/llm_session_*.jsonl and .err, then rerun "
                    f"step 2 for {paper} with --refresh-metadata"}
    status_path.write_text(json.dumps(st, indent=2, ensure_ascii=False),
                           encoding="utf-8")
    return None, st


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="One metadata extraction session.")
    ap.add_argument("--paper", required=True)
    ap.add_argument("--doi", default="")
    ap.add_argument("--title", default="")
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--reads-dir", type=Path, required=True)
    ap.add_argument("--orig-dir", type=Path, required=True)
    ap.add_argument("--paper-facts", type=Path, default=None)
    ap.add_argument("--agent", default="claude", choices=("claude", "codex"))
    ap.add_argument("--model", default="")
    ap.add_argument("--timeout", type=int, default=2400)
    ap.add_argument("--retries", type=int, default=1)
    ap.add_argument("--max-turns", type=int, default=120)
    a = ap.parse_args()
    import csv as _csv
    tsv = a.work_dir / "ena_runs.tsv"
    with tsv.open(encoding="utf-8", newline="") as fh:
        want = {r["sample_accession"].strip() for r in _csv.DictReader(fh, delimiter="\t")
                if (r.get("in_manifest") or "yes") == "yes"}
    _, st = extract(a.paper, a.doi, a.title, a.work_dir, a.reads_dir, a.orig_dir,
                    a.paper_facts, want, agent=a.agent, model=a.model,
                    timeout=a.timeout, retries=a.retries, max_turns=a.max_turns)
    sys.exit(0 if st.get("ok") else 1)
