#!/bin/bash
#SBATCH --job-name=virusDB
#SBATCH --partition=small
#SBATCH --time=1-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --output=/hickory/proj/didonglab/dataset/virus/Database/results/logs/submit_specific_steps.out
#SBATCH --error=/hickory/proj/didonglab/dataset/virus/Database/results/logs/submit_specific_steps.err
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


# rewrite search abundance paper

