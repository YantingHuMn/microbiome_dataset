# Microbiome paper-to-data discovery pipeline

This directory contains a paper-and-dataset discovery pipeline for bacterial or
viral abundance matrices. Stages 0-5 are **listing-only**: they search paper
metadata, screen studies, extract data accessions/URLs, and list repository
file names -- never opening or downloading deposited data. Stage 6 is the one
explicit, separate stage that does download and open files, to turn a
filename-based guess into a content-verified fact.

Every record that survives the pipeline ends up in exactly one of two
buckets, tracked from stage5 onward:

- **`raw_reads`** -- only sequencing reads are public; someone still has to
  run a 16S/metagenomics pipeline to get an abundance matrix.
- **`processed`** -- a matrix already exists. `source_pipeline_hint`
  (stage5, text-mined) / `source_pipeline` (stage6, content-confirmed)
  records which tool produced it (QIIME2, mothur, DADA2, Kraken2/Bracken,
  MetaPhlAn, HUMAnN, PICRUSt2, ...) wherever that's recoverable, so a later
  harmonisation step knows the matrix's provenance instead of treating all
  "processed" data as interchangeable.

## Stages

1. `stage0_build_queries.py` creates reproducible bacterial/viral search queries.
2. `stage1_search_papers.py` searches Europe PMC and builds a deduplicated paper list.
3. `stage2_screen_papers.py` performs high-recall title/abstract screening.
4. `stage3_extract_data_links.py` reads open-access full text when available and extracts repository accessions (BioProject/SRA/GEO/MGnify/Zenodo/figshare/Dryad/OSF/GitHub), URLs, and data-availability text. Also records a `europepmc_supp` link for every PMCID, since the paper's own PMC supplementary bundle is itself a candidate data source that earlier versions of this pipeline never looked at.
5. `stage4_resolve_datasets.py` converts accessions into canonical landing/download locations and lists repository *file names* where public APIs allow it (Zenodo, figshare, Dryad, OSF, ENA, and now the PMC supplementary ZIP's member names, plus an MGnify reverse-lookup from BioProject accessions). It still never persists deposited file content to disk.
6. `stage5_classify_candidates.py` assigns paper- and dataset-level usability states from filename evidence alone, and splits `data_category` into `processed` / `raw_reads` / `unknown`. Datasets with a table-shaped file in a processed-data repository but no decisive filename match go to `needs_content_check` instead of being silently merged into `manual_review` -- that queue is exactly stage6's input.
7. `stage6_verify_abundance.py` downloads `abundance_ready` + `needs_content_check` datasets ONE AT A TIME, judges each table by structure (taxa-shaped row/column labels, numeric block, long-vs-wide), writes accepted rows to a long-format table, and deletes the raw file immediately after -- so peak disk/memory usage is bounded by a single dataset, not the whole corpus. This is the stage that turns a filename guess into `content_verified`.

## Principal outputs

Stage 0-5 (listing only):
- `papers_master.csv`: all unique discovered papers.
- `manual_whitelist.tsv`: manual additions that should always be included even when the strict keyword search misses them.
- `papers_screened.csv`: initial relevance and assay classification.
- `paper_data_links.csv`: one row per paper/accession or external data URL.
- `datasets_master.csv`: canonical repository records, listed file names, and potential download locations.
- `paper_candidates_final.csv` / `dataset_candidates_final.csv`: final paper-/dataset-level triage, now carrying `data_category` (`processed`/`raw_reads`/`unknown`) and `source_pipeline_hint`.
- `abundance_ready.csv`, `raw_reads_ready.csv`, `needs_content_check.csv`, `manual_review.csv`: actionable subsets.

Stage 6 (opens files, content-verified) -- the three deliverables:
- **`abundance_final.csv`**: papers with >=1 content-verified abundance matrix. `data_category=processed`, `source_pipeline` names the tool (QIIME2/mothur/DADA2/Kraken2-Bracken/MetaPhlAn/...) where recoverable. Nothing left to do with these.
- **`run_own_pipeline.csv`**: papers whose `next_action=run_own_pipeline`; only raw reads are public, so run your 16S/metagenomics pipeline.
- **`manual_content_review.csv`**: papers whose `next_action=manual_content_review`; a table-shaped file exists but neither the filename nor the automated content check settled whether it is an abundance matrix, so it must be opened manually.

Supporting/audit files (not the deliverable, but the trail behind it):
- `dataset_verification.csv`: one row per dataset actually opened -- `content_verified`, `n_tables_checked/accepted`, `data_category`, `pipeline_source`, and a JSON `note` with the per-file verdicts.
- `paper_verification.csv`: every column from `paper_candidates_final.csv` plus the stage6 verdict, before the three-way deliverable split. Cases with other actions such as `no_data_found` remain visible here rather than being mixed into either action queue.
- `abundance_long.tsv.gz`: full long-format table (`paper_id, dataset_id, sample_id, taxon, kingdom..species, value, value_type, pipeline_source, domain, source_file`) behind every row of `abundance_final.csv`.
- `matrices/<dataset_id>__matrix.tsv`: one taxa x sample matrix per verified dataset, built from the long table at the end of the run (never held in memory during downloading).

## Manual allowlist workflow

Use `Database/results/stage1_papers/manual_whitelist.tsv` for papers that are valid but missed by the automated query logic. Each run of `stage1_search_papers.py` loads this file and unions its records into `papers_master.csv` before writing the final CSV. To add a new paper later, append a new TSV row with the same columns as the master CSV and rerun the stage.
