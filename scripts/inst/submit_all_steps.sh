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
# module purge
# module add r/4.4.0
# module load seurat/5.3.0-R4.4.0
module load anaconda
conda activate virus

cd /hickory/proj/didonglab/dataset/virus/yanting/microbiome_dataset
READ_DIR="/hickory/proj/didonglab/dataset/virus/yanting/microbiome_dataset/scripts/find_papers"
OUT_DIR="/hickory/proj/didonglab/dataset/virus/Database/results"

# All stages below accept --workers N for concurrent, per-host rate-limited
# requests (see scripts/find_papers/_netutil.py). Turn this down (or to 1,
# the old sequential behavior) if any one API starts throttling/erroring --
# raising it further past ~8-16 mostly stops helping since most APIs' own
# per-host pacing becomes the bottleneck before your CPU/network does.
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

python "${READ_DIR}/stage4_resolve_datasets.py" \
    --links "${OUT_DIR}/stage3_links/paper_data_links.csv" \
    --outdir "${OUT_DIR}/stage4_datasets" \
    --workers "${WORKERS}"

python "${READ_DIR}/stage5_classify_candidates.py" \
    --papers "${OUT_DIR}/stage3_links/papers_with_data_text.csv" \
    --datasets "${OUT_DIR}/stage4_datasets/datasets_master.csv" \
    --outdir "${OUT_DIR}/stage5_final"

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
