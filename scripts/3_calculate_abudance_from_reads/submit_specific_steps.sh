#!/bin/bash
#SBATCH --job-name=readsAbund
#SBATCH --partition=small
#SBATCH --time=6:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=/hickory/proj/didonglab/dataset/virus/Database/results/logs/reads_abundance1.out
#SBATCH --error=/hickory/proj/didonglab/dataset/virus/Database/results/logs/reads_abundance1.err
#SBATCH --mail-user=yanting@unc.edu
#SBATCH --mail-type=END,FAIL

# One paper at a time: download -> metadata -> abundance.
#
# Edit the block below, then just:    sbatch submit_specific_steps.sh
# Nothing needs to be set on the command line.
#
# Resuming is the default: steps already recorded done are skipped, and step 2
# reuses a previously validated metadata extraction rather than calling the
# model again. Set REFRESH_METADATA=1 to force a fresh extraction.
#
# Nothing here downloads a reference database. Step 3 refuses to run until its
# inputs are confirmed, so run it with PREFLIGHT=1 first. A preflight is
# recorded as `preflight_ok`, never as a finished abundance analysis.
#
# No credential is read, written or logged by this script: the model CLI uses
# its own existing login.

# ============================ EDIT THIS BLOCK ==============================
# --- which papers -----------------------------------------------------------
PAPER="Acuña2023,Bacci2021"   # one or more paper folders, comma or space separated
                              # (ignored if PAPERS_FILE is set)
PAPERS_FILE=""             # file with one paper folder per line, # for comments

# --- which steps ------------------------------------------------------------
STEPS="2,3"                  # 1=download  2=metadata  3=abundance, e.g. "1,2,3"
FORCE=""                   # 1 = redo steps already marked done
PREFLIGHT=""               # 1 = step 3 only checks inputs/tools, runs nothing

# --- step 2: resolving sample metadata from the paper -----------------------
REFRESH_METADATA=""        # 1 = call the model again even if an answer is cached
NO_LLM=""                  # 1 = ENA only. Cannot reach a full pass; debugging only
AGENT="claude"             # claude | codex
MODEL=""                   # e.g. sonnet / opus; empty = the CLI's own default
LLM_TIMEOUT=2400           # seconds per session before it is killed
LLM_RETRIES=1              # extra attempts after the first one fails
MAX_TURNS=120              # max turns per session (claude only)

# --- step 3: what to quantify ----------------------------------------------
# The ASSAY decides which branch runs; this only says what to get OUT of a
# shotgun library. It has no effect on a 16S paper.
#   bacteria | viruses | both
QUANTIFY="bacteria"

# --- step 3: 16S reference (installed) --------------------------------------
SILVA=/hickory/proj/didonglab/dataset/virus/Database/data/reference_db/dada2/silva_nr99_v138.2_toSpecies_trainset.fa.gz
SILVA_SPECIES=/hickory/proj/didonglab/dataset/virus/Database/data/reference_db/dada2/silva_v138.2_assignSpecies.fa.gz

# --- step 3: shotgun bacteria/archaea (NOT installed yet) -------------------
# See Database/data/reference_db/README.md for the verified install commands.
# Leave empty until sylph and its GTDB sketch exist; preflight will block and
# say so rather than running something wrong.
SYLPH_DB=""                # e.g. .../reference_db/sylph/v0.3-c200-gtdb-r214.syldb
SYLPH_TAX=""               # e.g. GTDB_r214; must match SYLPH_DB's release

# --- step 3: viral catalogue (NOT installed yet) ----------------------------
# The catalogue must match the sample environment. UHGV is HUMAN GUT ONLY.
VIRAL_INDEX=""             # bowtie2 index prefix, e.g. .../viral/uhgv/uhgv_votus
VIRAL_CATALOGUE_NAME="unnamed"

# --- step 3: host removal (optional, shotgun only) --------------------------
# Empty = not performed, and the output README says so rather than implying it.
HOST_INDEX=""              # bowtie2 index of the host genome
# ===========================================================================

set -uo pipefail
HERE=/hickory/proj/didonglab/dataset/virus/yanting/microbiome_dataset/scripts/3_calculate_abudance_from_reads

export PYTHONUNBUFFERED=1
if [ -d /nas/sycamore/home/yanting/.conda/envs/virus_syca ]; then
    source /nas/sycamore/apps/anaconda/2025.12-2/etc/profile.d/conda.sh
    conda activate virus_syca
else
    module load anaconda && conda activate virus
fi
# step 3 (DADA2) loads r/4.6.0 itself; dada2 and Biostrings are already in it.
echo "conda env: ${CONDA_DEFAULT_ENV:-unknown} on $(hostname), cpus=${SLURM_CPUS_PER_TASK:-?}"

ARGS=(--steps "$STEPS" --agent "$AGENT" --llm-timeout "$LLM_TIMEOUT"
      --llm-retries "$LLM_RETRIES" --max-turns "$MAX_TURNS")
[ -n "$MODEL" ]            && ARGS+=(--model "$MODEL")
[ -n "$FORCE" ]            && ARGS+=(--force)
[ -n "$PREFLIGHT" ]        && ARGS+=(--preflight)
[ -n "$REFRESH_METADATA" ] && ARGS+=(--refresh-metadata)
[ -n "$NO_LLM" ]           && ARGS+=(--no-llm)
[ -n "$SILVA" ]            && ARGS+=(--silva "$SILVA")
[ -n "$SILVA_SPECIES" ]    && ARGS+=(--silva-species "$SILVA_SPECIES")
[ -n "$QUANTIFY" ]         && ARGS+=(--quantify "$QUANTIFY")
[ -n "$SYLPH_DB" ]         && ARGS+=(--sylph-db "$SYLPH_DB")
[ -n "$SYLPH_TAX" ]        && ARGS+=(--sylph-tax "$SYLPH_TAX")
[ -n "$VIRAL_INDEX" ]      && ARGS+=(--viral-index "$VIRAL_INDEX")
[ -n "$VIRAL_CATALOGUE_NAME" ] && ARGS+=(--viral-catalogue-name "$VIRAL_CATALOGUE_NAME")
[ -n "$HOST_INDEX" ]       && ARGS+=(--host-index "$HOST_INDEX")

if [ -n "$PAPERS_FILE" ]; then
    ARGS+=(--papers-file "$PAPERS_FILE")
elif [ -n "$PAPER" ]; then
    # Split on commas and whitespace: "A,B" and "A B" both mean two papers.
    # Passing the whole string as one --paper would look for a folder literally
    # named "A,B" and fail only later, when step 2 cannot find its manifest.
    OLD_IFS=$IFS; IFS=', '
    for one in $PAPER; do
        [ -n "$one" ] && ARGS+=(--paper "$one")
    done
    IFS=$OLD_IFS
else
    echo "set PAPER or PAPERS_FILE in the block above"; exit 2
fi

python "$HERE/pipeline.py" "${ARGS[@]}"
