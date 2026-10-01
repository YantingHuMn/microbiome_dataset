#!/bin/bash
#SBATCH --job-name=virusDB
#SBATCH --partition=small
#SBATCH --time=1-00:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=/hickory/proj/didonglab/dataset/virus/Database/results/logs/submit_all_steps.out
#SBATCH --error=/hickory/proj/didonglab/dataset/virus/Database/results/logs/submit_all_steps.err
#SBATCH --mail-user=yanting@unc.edu
#SBATCH --mail-type=END,FAIL
# #SBATCH --dependency=afterany:4605511
# #SBATCH --gres=gpu:1
# #SBATCH --array=0-3

WORKERS=${SLURM_CPUS_PER_TASK:-1}
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"

export PYTHONUNBUFFERED=1
# module purge
# module add r/4.4.0
# module load seurat/5.3.0-R4.4.0

if [ -d /nas/sycamore/home/yanting/.conda/envs/virus_syca ]; then
    source /nas/sycamore/apps/anaconda/2025.12-2/etc/profile.d/conda.sh
    conda activate virus_syca
else
    module load anaconda
    conda activate virus
fi
echo "conda env: ${CONDA_DEFAULT_ENV:-unknown} ($(command -v python3)) on $(hostname)"

cd /hickory/proj/didonglab/dataset/virus/yanting/microbiome_dataset
READ_DIR="/hickory/proj/didonglab/dataset/virus/yanting/microbiome_dataset/scripts/find_papers"
OUT_DIR="/hickory/proj/didonglab/dataset/virus/Database/results"

# python ${READ_DIR}/stage0_build_queries.py \
#     --out ${OUT_DIR}/stage0_queries/search_queries.tsv

# python "${READ_DIR}/stage1_search_papers.py" \
#     --out_dir "${OUT_DIR}" \
#     --workers "${WORKERS}"

# python "${READ_DIR}/stage2_screen_papers.py" \
#     --papers "${OUT_DIR}/stage1_papers/papers_master.csv" \
#     --outdir "${OUT_DIR}/stage2_screen"

# python "${READ_DIR}/stage3_extract_data_links.py" \
#     --papers "${OUT_DIR}/stage2_screen/papers_screened.csv" \
#     --outdir "${OUT_DIR}/stage3_links" \
#     --workers "${WORKERS}"

# python "${READ_DIR}/stage4_resolve_datasets.py" \
#     --links "${OUT_DIR}/stage3_links/paper_data_links.csv" \
#     --outdir "${OUT_DIR}/stage4_datasets" \
#     --workers "${WORKERS}" \
#     --skip-mgnify


# python "${READ_DIR}/stage5_classify_candidates.py" \
#     --papers "${OUT_DIR}/stage3_links/papers_with_data_text.csv" \
#     --datasets "${OUT_DIR}/stage4_datasets/datasets_master.csv" \
#     --outdir "${OUT_DIR}/stage5_final"

# python scripts/find_papers/augment_europepmc_supp.py fetch \
#     --master /hickory/proj/didonglab/dataset/virus/Database/results/stage4_datasets/datasets_master.csv \
#     --progress /hickory/proj/didonglab/dataset/virus/Database/results/stage4_datasets/augment_progress.jsonl \
#     --workers 8

# python scripts/find_papers/augment_europepmc_supp.py merge \
#     --master /hickory/proj/didonglab/dataset/virus/Database/results/stage4_datasets/datasets_master.csv \
#     --progress /hickory/proj/didonglab/dataset/virus/Database/results/stage4_datasets/augment_progress.jsonl

# python scripts/find_papers/stage5_classify_candidates.py \
#     --datasets "${OUT_DIR}/stage4_datasets/datasets_master.csv" \
#     --papers "${OUT_DIR}/stage3_links/papers_with_data_text.csv" \
#     --outdir "${OUT_DIR}/stage5_final"


# # stage6 downloads real files -- --workers here also multiplies peak
# # disk/memory (roughly WORKERS x --max-study-mb under --scratch), unlike
# # stages 1/3/4 where concurrency only overlaps wait time. Drop WORKERS (a
# # local override below) if --scratch is small or --mem needs headroom.
# python "${READ_DIR}/stage6_verify_abundance.py" \
#     --datasets "${OUT_DIR}/stage5_final/dataset_candidates_final.csv" \
#     --papers "${OUT_DIR}/stage5_final/paper_candidates_final.csv" \
#     --outdir "${OUT_DIR}/stage6_verified" \
#     --scratch "${OUT_DIR}/stage6_scratch" \
#     --max-study-mb 500 \
#     --workers "${WORKERS}"

# download + analysize (can pause)
python3 scripts/extract_abundance_matrix/stage7_extract_abundance_matrix.py extract \
    --abundance-ready /hickory/proj/didonglab/dataset/virus/Database/results/stage5_final/abundance_ready.csv \
    --datasets /hickory/proj/didonglab/dataset/virus/Database/results/stage5_final/dataset_candidates_final.csv \
    --data-dir /hickory/proj/didonglab/dataset/virus/Database/data/abundance_matrix_v2 \
    --progress /hickory/proj/didonglab/dataset/virus/Database/results/stage7_extract/progress_v2.jsonl \
    --scratch /tmp/${SLURM_JOB_ID}_mb_stage7 \
    --raw-store /hickory/proj/didonglab/dataset/virus/Database/data/raw_files \
    --workers "${WORKERS}" \
    --range "50"

# combine
python3 scripts/extract_abundance_matrix/stage7_extract_abundance_matrix.py build \
    --abundance-ready /hickory/proj/didonglab/dataset/virus/Database/results/stage5_final/abundance_ready.csv \
    --progress /hickory/proj/didonglab/dataset/virus/Database/results/stage7_extract/progress_v2.jsonl \
    --studies-tsv /hickory/proj/didonglab/dataset/virus/Database/metadata/studies_v2.tsv \
    --sample-tsv /hickory/proj/didonglab/dataset/virus/Database/metadata/sample_v2.tsv \
    --blocked-tsv /hickory/proj/didonglab/dataset/virus/Database/metadata/blocked_manual_download_v2.tsv \
    --table-manifest-tsv /hickory/proj/didonglab/dataset/virus/Database/metadata/table_manifest_v2.tsv \
    --needs-review-tsv /hickory/proj/didonglab/dataset/virus/Database/metadata/needs_review_tables_v2.tsv