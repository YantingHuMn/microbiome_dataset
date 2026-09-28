#!/bin/bash
#SBATCH --job-name=virusDB
#SBATCH --partition=interact
#SBATCH --time=8:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=/hickory/proj/didonglab/dataset/virus/Database/results/logs/submit_all_steps.out
#SBATCH --error=/hickory/proj/didonglab/dataset/virus/Database/results/logs/submit_all_steps.err
#SBATCH --mail-user=yanting@unc.edu
#SBATCH --mail-type=END,FAIL
# #SBATCH --dependency=afterok:2270381
# #SBATCH --gres=gpu:1
# #SBATCH --array=0-3

echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
# Force unbuffered stdout/stderr for every python invocation below. Without
# this, output redirected to a file (not a TTY) is block-buffered by
# default -- a script that's genuinely running can leave .out/.err looking
# completely empty for a long time even though it already printed its
# startup lines, because they're sitting in an unflushed buffer rather
# than on disk. This removes that ambiguity for good.
export PYTHONUNBUFFERED=1
# module purge
# module add r/4.4.0
# module load seurat/5.3.0-R4.4.0
#
# Environment activation is cluster-specific -- module systems and conda
# install paths differ between longleaf and sycamore. Detect by checking
# whether the sycamore conda install PATH exists, not by hostname --
# compute-node hostnames are usually generic (e.g. longleaf's own compute
# nodes are named "c151613", nothing like "longleaf"), so a job actually
# running on a sycamore compute node would NOT match "*sycamore*" even
# though the login node does. The shared /nas/sycamore filesystem, by
# contrast, is mounted identically on every node of that cluster
# (that's the whole point of a shared filesystem), so its existence is a
# reliable signal regardless of what the current node happens to be
# called. (A wrong/missing conda env here means every python3 call below
# silently runs in whatever env happened to be active instead -- e.g.
# check_dependencies() correctly caught this as a missing xlrd, but a
# less-guarded dependency could fail silently instead.)
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

WORKERS=8

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
    --data-dir /hickory/proj/didonglab/dataset/virus/Database/data \
    --progress /hickory/proj/didonglab/dataset/virus/Database/results/stage7_extract/progress.jsonl \
    --scratch /tmp/${SLURM_JOB_ID}_mb_stage7 \
    --workers 8

# combine
python3 scripts/extract_abundance_matrix/stage7_extract_abundance_matrix.py build \
    --abundance-ready /hickory/proj/didonglab/dataset/virus/Database/results/stage5_final/abundance_ready.csv \
    --progress /hickory/proj/didonglab/dataset/virus/Database/results/stage7_extract/progress.jsonl \
    --studies-tsv /hickory/proj/didonglab/dataset/virus/Database/metadata/studies.tsv \
    --sample-tsv /hickory/proj/didonglab/dataset/virus/Database/metadata/sample.tsv \
    --blocked-tsv /hickory/proj/didonglab/dataset/virus/Database/metadata/blocked_manual_download.tsv