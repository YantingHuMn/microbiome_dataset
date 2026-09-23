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
