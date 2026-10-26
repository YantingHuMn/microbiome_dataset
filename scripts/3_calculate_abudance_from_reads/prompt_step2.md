You are resolving the sample metadata for ONE paper whose raw reads have already
been downloaded. You are running unattended — no human will answer questions, so
decide yourself, and where the evidence does not support a conclusion, say so
instead of guessing. A wrong label is far worse here than an admitted gap.

## Your paper
- Folder: {PAPER}
- DOI: {DOI}  (https://doi.org/{DOI})
- Title: {TITLE}
- Downloaded reads (READ ONLY, do not modify or delete): {READS_DIR}
- Original/supplementary files already on disk, if any (READ ONLY): {ORIG_DIR}
- Working dir for this task (write ONLY here): {WORK_DIR}
- Write your structured answer to exactly this path: {OUT_JSON}

## What already exists (inputs, not answers)
- `{WORK_DIR}/ena_runs.tsv` — every run this paper's FASTQs came from, with the
  ENA fields already retrieved: run, experiment, BioSample, study, library
  strategy/source/selection/layout, instrument, read count, scientific_name,
  host, isolation_source, sample_title, sample_alias.
- `{PAPER_FACTS}` — a human-curated facts file, if one exists. It is EVIDENCE
  you may cite and build on. It does NOT replace your own check: confirm its
  claims against the paper, and report anything it gets wrong or omits.

## The chain you must resolve
Knowing a BioSample accession is NOT knowing what the sample is. Resolve the
whole chain and say where each link came from:

    FASTQ file -> run (SRR/ERR/DRR) -> experiment -> BioSample
              -> the sample identifier used IN THE PAPER
              -> everything an analysis could condition on (see the field list)

The last two links are the ones that matter and the ones that are usually
missing from ENA. Look for them in, roughly this order:

1. the ENA/SRA sample attributes already in `ena_runs.tsv` (aliases and titles
   often encode group, subject, site or day)
2. the paper's Methods, figures, and figure legends
3. the supplementary material — tables of sample IDs, mapping files, metadata
   sheets. Download them if they are openly available.
4. the authors' code or data repository (GitHub/Zenodo/OSF/figshare) — mapping
   files and QIIME/phyloseq metadata files live there surprisingly often
5. the BioProject/BioSample records at NCBI, which sometimes carry attributes
   that the ENA mirror does not expose

Europe PMC usually has the open-access full text and the supplementary files:
  https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:%22{DOI}%22&format=json&resultType=core
  https://www.ebi.ac.uk/europepmc/webservices/rest/PMC<id>/fullTextXML
  https://www.ebi.ac.uk/europepmc/webservices/rest/PMC<id>/supplementaryFiles

## The per-sample fields to resolve

These two decide what the paper can be used for, so they must appear in every
sample's `fields`, even when their confidence is `unknown`:

{FIELD_TABLE_REQUIRED}

The rest are what downstream analyses condition on. Fill in every one the paper
supports; omit the ones it does not report:

{FIELD_TABLE_ANALYSIS}

`sample_type` and `body_site` are deliberately separate: a stool sample and a
colon biopsy can both be "gut", and a swab is not an anatomical site. Put the
material in `sample_type` and the anatomy in `body_site`. ENA's
`isolation_source` usually mixes the two, so split it rather than copying it
into both.

### Say paper-wide facts once

Most of these are constant across a whole paper. Put those in `sample_defaults`
with the same `{value, confidence, evidence}` shape and they apply to every
sample that does not override them. Put a field in a sample's own `fields` only
when it varies between samples or is genuinely sample-specific. Do not repeat a
constant 100 times — a truncated answer is worse than a compact one.

## Rules
- **Never invent a group label.** If you cannot tell which samples are cases and
  which are controls, that is the finding. Report it.
- **Separate what you confirmed from what you suspect.** Confidence is recorded
  **per field**, not per sample: the sample id can be certain while the group is
  only a guess, and that must be visible. Each field carries its own
  `confidence`:
  - `confirmed` — a specific document states it, and you cite where
  - `inferred` — a consistent, stated reasoning step from evidence (e.g. an
    alias pattern that a figure legend defines); say what the reasoning was
  - `unknown` — you could not establish it. Leave the value empty.
  Only `confirmed` will be treated downstream as fact.
- A sample count that matches the paper's N is supporting evidence, never on its
  own a mapping. "63 runs and the paper says 25+38" does not tell you which run
  is which.
- Mark extraction blanks, buffer controls, mock communities, dilution series and
  positive controls with `is_technical_control: true`. They are not biological
  samples.
- Do not modify anything under {READS_DIR} or {ORIG_DIR}.
- Do not download sequencing data. Supplementary tables and metadata files only.

## Output — write this JSON to {OUT_JSON} and nothing else to that path

```json
{
  "paper_folder": "{PAPER}",
  "doi": "{DOI}",
  "sample_id_scheme": {
    "description": "how the ENA alias/title maps to the paper's own sample ids, or why it cannot be mapped",
    "confidence": "confirmed|inferred|unknown",
    "evidence": "exact source: document, table/figure number, URL, file path"
  },
  "key_experiment": {
    "host": "", "body_site": "", "design": "", "disease": "",
    "group_names": [], "n_groups": 0,
    "confidence": "confirmed|inferred|unknown",
    "evidence": ""
  },
  "sample_defaults": {
    "_comment": "any of the analysis fields below that are constant across the paper; omit the rest",
{DEFAULTS_JSON}
  },
  "samples": [
    {
      "biosample_accession": "SAMN...",
      "run_accessions": ["SRR..."],
      "is_technical_control": false,
      "fields": {
{FIELDS_JSON}
      }
    }
  ],
  "unresolved": ["each thing you could not establish, and what would resolve it"],
  "sources_consulted": [
    {"type": "paper|supplementary|repository|ena|ncbi|other", "locator": "URL or path", "useful": true}
  ],
  "summary": "3-6 sentences: what is now known, what is not, and what the paper can and cannot be used for"
}
```

`paper_sample_id` and `experimental_group` must appear in every sample's
`fields` (use `"confidence": "unknown"` with an empty value if you could not
establish them). Every other field goes in `fields` when it varies between
samples, in `sample_defaults` when it is constant, and is simply omitted when
the paper does not report it.

## Scope: the downloaded samples only

`ena_runs.tsv` has an `in_manifest` column. Your answer must contain exactly one
entry in `samples` for each of the **{N_SAMPLES} BioSamples with
`in_manifest = yes`** — including the ones you could not resolve, because an
omitted sample is indistinguishable from an overlooked one.

Rows marked `in_manifest = no` belong to the same study but were not downloaded.
Read them freely — they often reveal the naming scheme, the full group
structure or the controls — but do NOT put them in `samples`. `sources_consulted` must list what you actually opened,
including sources that turned out to be useless.

Write the JSON file, then reply with a short plain-text summary. Do not print
the JSON in your reply.
