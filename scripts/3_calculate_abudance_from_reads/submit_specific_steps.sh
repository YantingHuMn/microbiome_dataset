#!/bin/bash
#SBATCH --job-name=readsAbund
#SBATCH --partition=small
#SBATCH --time=1-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --dependency=afterok:4950793
#SBATCH --output=/hickory/proj/didonglab/dataset/virus/Database/results/logs/reads_abundance.out
#SBATCH --error=/hickory/proj/didonglab/dataset/virus/Database/results/logs/reads_abundance.err
#SBATCH --mail-user=yanting@unc.edu
#SBATCH --mail-type=END,FAIL

# One paper at a time: download -> metadata -> abundance.

# ============================ EDIT THIS BLOCK ==============================
CLAIM_N="5"
CLAIM_SOURCE=/hickory/proj/didonglab/dataset/virus/Database/data/abudance_integration/1_2_status_raw_reads_only.csv
CLAIM_EXCLUDE="" # e.g. "risk_no_grouping.csv" to also skip those

PAPER=""
PAPERS_FILE=""  # default file path

STEPS="1,2,3" # 1=download  2=metadata  3=abundance, e.g. "1,2,3"
FORCE=""  # 1 = redo steps already marked done; blank means skip.
          # Keep this BLANK for a claim run: Bacci2021 already has step 1
          # and 2 done, and forcing would re-download and re-pay for them.
PREFLIGHT="" # 1 = step 3 only checks inputs/tools, runs nothing
ARCHIVE_OLD="1"  # 1 = move a previous run's dada2/shotgun/virome into
                 #     _archive/<utc>/ instead of overwriting it

# --- step 2: resolving sample metadata from the paper -----------------------
REFRESH_METADATA="" # 1 = call the model again even if an answer is cached
NO_LLM="" # 1 = ENA only. Cannot reach a full pass; debugging only
AGENT="codex"  # claude | codex
MODEL=""  # e.g. sonnet / opus; empty = the CLI's own default
LLM_TIMEOUT=2400 # seconds per session before it is killed
LLM_RETRIES=1 # extra attempts after the first one fails
MAX_TURNS=120 # max turns per session (claude only)

# --- step 3: what to quantify ----------------------------------------------
# The ASSAY decides which branch runs; this only says what to get OUT of a
# shotgun library. It has no effect on a 16S paper.
# bacteria | viruses | both
# "both" is the right default now that both databases are installed: it has no
# effect on a 16S paper, and for a bulk shotgun it gets you the viral table too
# instead of silently leaving it on the table. A VLP-enriched library still gets
# only the viral branch - counting bacteria in a library built to exclude them
# is blocked, not quietly done.
QUANTIFY="both"

# --- step 3: 16S reference (installed) --------------------------------------
SILVA=/hickory/proj/didonglab/dataset/virus/Database/data/reference_db/dada2/silva_nr99_v138.2_toSpecies_trainset.fa.gz
SILVA_SPECIES=/hickory/proj/didonglab/dataset/virus/Database/data/reference_db/dada2/silva_v138.2_assignSpecies.fa.gz

# --- step 3: shotgun bacteria/archaea (installed) --------------------------
# sylph 1.0.0 + GTDB r214 sketch (85,205 genomes). Leave empty to skip.
SYLPH_DB=/hickory/proj/didonglab/dataset/virus/Database/data/reference_db/sylph/v0.3-c1000-gtdb-r214.syldb
SYLPH_TAX=GTDB_r214        # must match SYLPH_DB's release

# --- step 3: viral catalogue (installed) ------------------------------------
# UHGV is a HUMAN GUT catalogue - 168,536 vOTUs. Leave VIRAL_INDEX empty to skip.
VIRAL_INDEX=/hickory/proj/didonglab/dataset/virus/Database/data/reference_db/viral/uhgv/uhgv_votus
VIRAL_CATALOGUE_NAME="UHGV vOTU representatives (human gut)"
# Which environment this catalogue is for. Step 3 compares it with the paper's
# host/body site and BLOCKS on a mismatch, so a mouse or plant library is never
# mapped against a human gut catalogue. Clear it to disable the check (you then
# own the consequences), or point VIRAL_INDEX at a catalogue for that
# environment and change this to match.
# Plain words matched as substrings, never a regex. Groups are separated by ';'
# and must ALL match; alternatives inside a group by '|'. Host and site are
# separate groups on purpose: with one combined list a MOUSE GUT sample would
# match "gut" and be mapped against a HUMAN GUT catalogue.
VIRAL_CATALOGUE_SCOPE="human|homo sapiens ; gut|stool|faec|fec|intestin|colon|rectum|caec|cecum"

# --- step 3: host removal (optional, shotgun only) --------------------------
# Empty = not performed, and the output README says so rather than implying it.
HOST_INDEX=""              # bowtie2 index of the host genome
# ===========================================================================

set -uo pipefail
HERE=/hickory/proj/didonglab/dataset/virus/yanting/microbiome_dataset/scripts/3_calculate_abudance_from_reads
OUT_ROOT=/hickory/proj/didonglab/dataset/virus/Database/data/reads_abundance

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
[ -n "$ARCHIVE_OLD" ]      && ARGS+=(--archive-old)
[ -n "$REFRESH_METADATA" ] && ARGS+=(--refresh-metadata)
[ -n "$NO_LLM" ]           && ARGS+=(--no-llm)
[ -n "$SILVA" ]            && ARGS+=(--silva "$SILVA")
[ -n "$SILVA_SPECIES" ]    && ARGS+=(--silva-species "$SILVA_SPECIES")
[ -n "$QUANTIFY" ]         && ARGS+=(--quantify "$QUANTIFY")
[ -n "$SYLPH_DB" ]         && ARGS+=(--sylph-db "$SYLPH_DB")
[ -n "$SYLPH_TAX" ]        && ARGS+=(--sylph-tax "$SYLPH_TAX")
[ -n "$VIRAL_INDEX" ]      && ARGS+=(--viral-index "$VIRAL_INDEX")
[ -n "$VIRAL_CATALOGUE_NAME" ] && ARGS+=(--viral-catalogue-name "$VIRAL_CATALOGUE_NAME")
[ -n "$VIRAL_CATALOGUE_SCOPE" ] && ARGS+=(--viral-catalogue-scope "$VIRAL_CATALOGUE_SCOPE")
[ -n "$HOST_INDEX" ]       && ARGS+=(--host-index "$HOST_INDEX")

ARGS+=(--release-claims)

if [ -n "$CLAIM_N" ]; then
    # Claim first, then run only what this job actually holds. Each paper is
    # released from the claim as soon as its own run ends, so the next
    # submission sees it as finished rather than as still being worked on.
    CLAIM_LIST="$OUT_ROOT/claims/${SLURM_JOB_ID:-local$$}.papers.txt"
    mkdir -p "$OUT_ROOT/claims"
    EXC=()
    for e in $CLAIM_EXCLUDE; do EXC+=(--exclude "$e"); done
    python "$HERE/claim_papers.py" --n "$CLAIM_N" --source "$CLAIM_SOURCE" \
        --out-root "$OUT_ROOT" --out "$CLAIM_LIST" "${EXC[@]}" || exit 1
    if [ ! -s "$CLAIM_LIST" ]; then
        echo "nothing left to claim - every paper is done, failed, unsupported, or held by another job"
        exit 0
    fi
    ARGS+=(--papers-file "$CLAIM_LIST")
elif [ -n "$PAPERS_FILE" ]; then
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
