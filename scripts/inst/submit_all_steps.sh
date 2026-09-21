#!/bin/bash
#SBATCH --job-name=virusDB
#SBATCH --partition=general
#SBATCH --time=1-00:00:00
#SBATCH --mem=32G
#SBATCH --output=/proj/didonglab/dataset/virus/Database/results/logs/submit_all_steps.out
#SBATCH --error=/proj/didonglab/dataset/virus/Database/results/logs/submit_all_steps.err
#SBATCH --mail-user=yanting@unc.edu
#SBATCH --mail-type=END,FAIL
# #SBATCH --dependency=afterok:65211362
# #SBATCH --gres=gpu:1
# #SBATCH --array=0-3

set -euo pipefail

# module purge
# module add r/4.4.0
# module load seurat/5.3.0-R4.4.0
module load anaconda
conda activate virus

cd /proj/didonglab/dataset/virus/yanting/microbiome_dataset

READ_DIR="/proj/didonglab/dataset/virus/yanting/microbiome_dataset/scripts/find_papers"
OUT_DIR="/proj/didonglab/dataset/virus/Database/results"

# python ${READ_DIR}/stage0_build_queries.py \
#     --out ${OUT_DIR}/stage0_queries/search_queries.tsv

python "${READ_DIR}/stage1_search_papers.py" \
    --out_dir "${OUT_DIR}"

python "${READ_DIR}/stage2_screen_papers.py" \
    --papers "${OUT_DIR}/stage1_papers/papers_master.csv" \
    --outdir "${OUT_DIR}/stage2_screen"

python "${READ_DIR}/stage3_extract_data_links.py" \
    --papers "${OUT_DIR}/stage2_screen/papers_screened.csv" \
    --outdir "${OUT_DIR}/stage3_links"

python "${READ_DIR}/stage4_resolve_datasets.py" \
    --links "${OUT_DIR}/stage3_links/paper_data_links.csv" \
    --outdir "${OUT_DIR}/stage4_datasets"

python "${READ_DIR}/stage5_classify_candidates.py" \
    --papers "${OUT_DIR}/stage3_links/papers_with_data_text.csv" \
    --datasets "${OUT_DIR}/stage4_datasets/datasets_master.csv" \
    --outdir "${OUT_DIR}/stage5_final"
