# Microbiome Dataset

This repository is used to collect, organize, and standardize microbiome datasets obtained from published studies.

## Data Types

The planned database may include:

- Viral abundance
- Bacterial abundance
- Disease status
- Human demographic and clinical metadata
- Study and publication information

## Repository Structure

- `scripts/find_papers/`: paper-and-dataset discovery pipeline (stages 0-6). See its
  [README](scripts/find_papers/README.md) for the full stage list. Stages 0-5 only
  search and list; stage 6 is the one stage that downloads and content-verifies files.
- `scripts/inst/`: cluster submission scripts.
- `metadata/`: Standardized study and sample metadata
- `docs/`: Documentation

## Requirements

stage0-5 use only the Python standard library. stage6 (the one stage that
opens actual spreadsheet/table files) needs `pandas`, `numpy`, `openpyxl`,
and `xlrd` -- see `scripts/find_papers/requirements.txt`. Install them into
whichever conda env `submit_all_steps.sh` activates (`virus` on the UNC
cluster) BEFORE running stage6:

```
pip install -r scripts/find_papers/requirements.txt
```

stage6 also checks for these itself at startup and exits immediately with
a clear message if any are missing, rather than continuing and silently
mis-scoring every `.xlsx`/`.xls` file it encounters as "no abundance
matrix found" (which is what happened before this check existed: a whole
cluster run logged nothing but per-file openpyxl import errors and never
actually opened any of that run's `.xlsx` supplementary tables).

## Data category

Every paper the pipeline keeps lands in one of two buckets:

- **`processed`** -- a matrix already exists somewhere public. `source_pipeline`
  (or `source_pipeline_hint` before content verification) names the tool that
  produced it, where recoverable: QIIME2, mothur, DADA2, Kraken2/Bracken,
  MetaPhlAn, HUMAnN, PICRUSt2, etc.
- **`raw_reads`** -- only sequencing reads are public; a 16S/metagenomics
  pipeline still has to be run to get an abundance matrix.

The content-confirmed, final split lives in two files under
`scripts/find_papers/results/stage6_verified/`:
`abundance_final.csv` (matrix in hand) and `needs_pipeline_or_review.csv`
(raw reads needing a pipeline, or ambiguous cases needing manual review --
distinguished inside that file by its `next_action` column).
