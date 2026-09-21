# Microbiome paper-to-data discovery pipeline

This directory contains a **pre-download** pipeline for discovering papers and public datasets that may support bacterial or viral abundance matrices.  It searches paper metadata, screens studies, extracts data accessions and URLs, and records where data could be downloaded.

## Stages

1. `stage0_build_queries.py` creates reproducible bacterial/viral search queries.
2. `stage1_search_papers.py` searches Europe PMC and builds a deduplicated paper list.
3. `stage2_screen_papers.py` performs high-recall title/abstract screening.
4. `stage3_extract_data_links.py` reads open-access full text when available and extracts repository accessions, URLs, and data-availability text.
5. `stage4_resolve_datasets.py` converts accessions into canonical landing/download locations and retrieves lightweight repository metadata/file listings where public APIs allow it. It never downloads deposited data files.
6. `stage5_classify_candidates.py` assigns paper- and dataset-level usability states.


## Principal outputs

- `papers_master.csv`: all unique discovered papers.
- `papers_screened.csv`: initial relevance and assay classification.
- `paper_data_links.csv`: one row per paper/accession or external data URL.
- `datasets_master.csv`: canonical repository records and potential download locations.
- `paper_candidates_final.csv`: final paper-level triage.
- `dataset_candidates_final.csv`: final dataset-level triage.
- `abundance_ready.csv`, `raw_reads_ready.csv`, `manual_review.csv`: actionable subsets.
