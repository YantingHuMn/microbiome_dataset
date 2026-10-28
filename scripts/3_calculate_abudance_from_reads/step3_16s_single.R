#!/usr/bin/env Rscript
#
# DADA2 denoising for a SINGLE-END 16S paper. Called by step3_abundance.py.
#
# Why this is a separate script and not a flag on step3_16s.R: the paired script
# is built around the merge budget. It measures the amplicon length from R1/R2
# overlap, picks truncLen so the two mates still span it, and treats a collapsed
# merge as a hard failure. None of that exists here. There is one read, it is
# truncated on quality alone, and the ASV IS that fragment.
#
# What that costs, stated once so it travels with the data:
#   * the ASV covers only the 5' part of the amplicon, so two taxa that differ
#     only downstream of the cut are indistinguishable here
#   * taxonomic resolution is lower than the same region sequenced paired
#   * ASVs from this paper CANNOT be compared sequence-to-sequence with ASVs
#     from a paired run of the same region: they are different-length fragments
#
# truncLen is chosen from the data: expected-error profiles give the share of
# reads surviving maxEE at each possible cut, and the choice balances that
# against keeping the fragment long enough to classify. --trunc-len overrides.

suppressPackageStartupMessages({
  library(dada2); library(Biostrings); library(ShortRead)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1]]
}
int2 <- function(x) as.integer(strsplit(x, ",")[[1]])

paper        <- getarg("--paper")
reads_root   <- getarg("--reads-root")
out_dir      <- getarg("--out-dir")
region       <- getarg("--region", "unknown")
trim_left    <- int2(getarg("--trim-left", "0,0"))[1]
max_trunc    <- int2(getarg("--max-trunc", "0,0"))[1]
threads      <- as.integer(getarg("--threads", "4"))
silva        <- getarg("--silva")
silva_sp     <- getarg("--silva-species", NA)
max_ee       <- as.numeric(getarg("--max-ee", "2"))
trunc_manual <- if (is.null(getarg("--trunc-len"))) NULL else int2(getarg("--trunc-len"))[1]
probe_n      <- as.integer(getarg("--probe-samples", "6"))
probe_reads  <- as.integer(getarg("--probe-reads", "40000"))
min_reads    <- as.integer(getarg("--min-reads", "1000"))
min_len      <- as.integer(getarg("--min-asv-len", "200"))
min_err_bases<- as.numeric(getarg("--min-error-bases", "1e8"))
grouped      <- getarg("--grouped-analysis", "1") == "1"
amp_status   <- getarg("--amplicon-status", "unknown")
only_probe   <- "--probe-only" %in% args

stopifnot(!is.null(paper), !is.null(out_dir))
dir.create(file.path(out_dir, getarg("--out-sub", "dada2")), recursive = TRUE, showWarnings = FALSE)
out_sub <- getarg("--out-sub", "dada2")
res <- function(...) file.path(out_dir, out_sub, ...)
say <- function(...) cat(sprintf(...), "\n", sep = "")
fail <- function(...) { cat("FAILED: ", sprintf(...), "\n", sep = ""); quit(status = 1) }

# ---- samples and files -----------------------------------------------------
# The sample list is passed in, not assumed: a paper can hold several assay
# types, and feeding another branch's runs to this one would denoise shotgun or
# RNA reads as if they were this amplicon.
sfile <- getarg("--samples", file.path(out_dir, "samples_used.csv"))
if (!file.exists(sfile)) fail("%s missing; run step 3 preflight first", sfile)
meta <- read.csv(sfile, check.names = FALSE, fileEncoding = "UTF-8-BOM")
say("%s: %d biological samples (single-end)", paper, nrow(meta))

mapfile <- file.path(out_dir, "run_sample_mapping.csv")
if (!file.exists(mapfile)) fail("%s missing; run step 2 first", mapfile)
rmap <- read.csv(mapfile, check.names = FALSE, fileEncoding = "UTF-8-BOM")
rmap <- rmap[rmap$biosample_accession %in% meta$sample_id, , drop = FALSE]

layout <- unique(toupper(rmap$library_layout))
if (!identical(layout, "SINGLE"))
  fail(paste("this script is for single-end reads; run_sample_mapping.csv reports",
             "library_layout = %s. Use step3_16s.R for paired data."),
       paste(layout, collapse = "/"))

bad_merge <- meta$sample_id[meta$multi_run_handling ==
                            "MULTIPLE_EXPERIMENTS_do_not_merge_without_review"]
if (length(bad_merge))
  fail(paste("%d sample(s) have runs from different experiments and step 2 flagged",
             "them for review; pooling them would mix libraries: %s"),
       length(bad_merge), paste(head(bad_merge, 5), collapse = ", "))

cat_dir <- res("pooled_runs")
fq <- character(nrow(meta)); n_runs_used <- integer(nrow(meta))
for (i in seq_len(nrow(meta))) {
  sid <- meta$sample_id[i]
  f <- rmap$fastq_path[rmap$biosample_accession == sid]
  f <- sort(f[file.exists(f)])
  if (!length(f)) fail("%s: no FASTQ on disk from run_sample_mapping.csv", sid)
  n_runs_used[i] <- length(f)
  if (length(f) == 1) {
    fq[i] <- f
  } else {                      # every run is used, not just the first
    dir.create(cat_dir, showWarnings = FALSE, recursive = TRUE)
    dest <- file.path(cat_dir, sprintf("%s.fastq.gz", sid))
    if (!file.exists(dest)) {
      if (!file.create(dest)) fail("could not create %s", dest)
      if (!all(file.append(dest, f))) { unlink(dest); fail("could not pool runs for %s", sid) }
    }
    fq[i] <- dest
  }
}
names(fq) <- meta$sample_id
meta$n_runs_used <- n_runs_used
if (any(n_runs_used > 1))
  say("%d sample(s) sequenced across several runs; all runs pooled per sample",
      sum(n_runs_used > 1))

# ---- sequencing batches ----------------------------------------------------
flowcell_of <- function(path) {
  con <- gzfile(path, "r"); on.exit(close(con))
  h <- readLines(con, n = 1)
  m <- regmatches(h, regexpr("[A-Za-z0-9_-]+:[0-9]+:[A-Za-z0-9-]+:[0-9]+:", h))
  if (!length(m)) return(NA_character_)
  paste(strsplit(m, ":")[[1]][1:3], collapse = ":")
}
batch <- vapply(fq, function(f) tryCatch(flowcell_of(f), error = function(e) NA_character_), "")
if (anyNA(batch)) {
  say("batch: FASTQ headers carry no flow-cell id (stripped on submission); using")
  say("       ONE error model for all samples rather than inventing batches.")
  batch[] <- "unknown_single_model"
} else {
  say("batch: %d flow cell(s): %s", length(unique(batch)),
      paste(sprintf("%s (n=%d)", names(table(batch)), table(batch)), collapse = ", "))
}
meta$sequencing_batch <- batch

# ---- truncLen from expected-error profiles ---------------------------------
ee_profile <- function(files, n) {
  cum <- NULL; tot <- 0
  for (f in files) {
    fqs <- ShortRead::yield(ShortRead::FastqStreamer(f, n))
    q <- as(Biostrings::quality(fqs), "matrix")
    # A read shorter than the cut is discarded by filterAndTrim, so its missing
    # tail must score as a failure rather than as zero error.
    p <- 10^(-q / 10); p[is.na(p)] <- Inf
    ee <- t(apply(p, 1, cumsum))
    if (!is.null(cum)) {
      w <- min(ncol(cum), ncol(ee))
      cum <- rbind(cum[, seq_len(w), drop = FALSE], ee[, seq_len(w), drop = FALSE])
    } else cum <- ee
    tot <- tot + nrow(ee)
    if (tot >= n) break
  }
  cum
}
set.seed(1)
probe_ids <- unlist(lapply(split(meta$sample_id, meta$sequencing_batch), function(ids)
  ids[round(seq(1, length(ids), length.out = max(1, round(probe_n * length(ids) / nrow(meta)))))]))
probe_ids <- unique(probe_ids)
say("probe subset (%d samples): %s", length(probe_ids), paste(probe_ids, collapse = ", "))

ee <- ee_profile(fq[probe_ids], probe_reads)
read_len <- ncol(ee)
cap <- if (max_trunc > 0) min(max_trunc, read_len) else read_len
say("read length %d; truncLen capped at %d", read_len, cap)

if (!is.null(trunc_manual)) {
  trunc_len <- trunc_manual
  basis <- sprintf("MANUAL OVERRIDE via --trunc-len %d", trunc_len)
  grid <- NULL
} else {
  # There is no merge to score, so the trade-off is retention against fragment
  # length. Keep the longest cut that still holds a high share of reads: a
  # short fragment classifies poorly however many reads survive.
  cuts <- seq(cap, max(trim_left + min_len, 100), by = -5)
  keep <- vapply(cuts, function(L) if (L > ncol(ee)) NA_real_ else mean(ee[, L] <= max_ee), 0)
  grid <- data.frame(trunc_len = cuts, retained = keep,
                     asv_length = cuts - trim_left)
  grid <- grid[!is.na(grid$retained) & grid$asv_length >= min_len, ]
  if (!nrow(grid)) fail("no truncLen keeps an ASV of at least %d bp", min_len)
  best <- max(grid$retained)
  ok <- grid[grid$retained >= 0.95 * best, ]          # within 5% of the best
  trunc_len <- max(ok$trunc_len)                      # then the longest of those
  basis <- sprintf(paste("chosen from expected-error profiles on %d probe samples:",
                         "%.1f%% of reads pass maxEE=%g at this cut (best achievable",
                         "%.1f%%); the longest cut within 5%% of the best was taken so",
                         "the ASV stays as long as possible"),
                   length(probe_ids), 100 * grid$retained[grid$trunc_len == trunc_len],
                   max_ee, 100 * best)
  write.csv(grid, res("trunclen_scan.csv"), row.names = FALSE)
}
say("truncLen = %d  (ASV length %d bp)  [%s]", trunc_len, trunc_len - trim_left, basis)
if (only_probe) {
  writeLines(c(sprintf("truncLen\t%d", trunc_len), basis), res("trunclen_choice.txt"))
  say("--probe-only: stopping before the full run.")
  quit(status = 0)
}

# ---- filter ----------------------------------------------------------------
fdir <- res("filtered")
unlink(fdir, recursive = TRUE)        # never reuse a previous run's output
dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
filt <- file.path(fdir, paste0(meta$sample_id, ".fq.gz"))
ft <- filterAndTrim(fq, filt, truncLen = trunc_len, trimLeft = trim_left,
                    maxEE = max_ee, truncQ = 2, maxN = 0, rm.phix = TRUE,
                    compress = TRUE, multithread = threads, verbose = FALSE)
rownames(ft) <- meta$sample_id
kept <- file.exists(filt) & ft[, "reads.out"] >= min_reads
if (!any(kept)) fail("every sample fell below %d reads after filtering", min_reads)
if (any(!kept)) say(paste("%d sample(s) below %d reads after filtering, reported with",
                          "zero counts rather than dropped: %s"),
                    sum(!kept), min_reads,
                    paste(meta$sample_id[!kept], collapse = ", "))

# ---- error models, per real batch ------------------------------------------
bat <- meta$sequencing_batch[kept]
bases <- tapply(ft[kept, "reads.out"], bat, sum) * (trunc_len - trim_left)
big <- names(bases)[bases >= min_err_bases]
small <- setdiff(names(bases), big)
model_of <- if (length(small))
  ifelse(bat %in% big, bat, paste0("own_model_below_threshold:", bat)) else bat
if (length(small))
  say(paste("WARNING batch(es) %s hold under %.0e bases; they keep their OWN error",
            "model rather than being mixed into another flow cell's."),
      paste(small, collapse = ", "), min_err_bases)
for (b in names(bases))
  say("  batch %s: %d sample(s), %.2e bases", b, sum(bat == b), bases[[b]])

errs <- list()
for (b in unique(model_of)) {
  i <- which(model_of == b)
  errs[[b]] <- learnErrors(filt[kept][i], multithread = threads, verbose = FALSE)
}
pdf(res("error_models.pdf"))
for (b in names(errs)) print(plotErrors(errs[[b]], nominalQ = TRUE) + ggplot2::ggtitle(b))
invisible(dev.off())

# ---- denoise ---------------------------------------------------------------
dadas <- list()
for (b in unique(model_of)) {
  i <- which(model_of == b)
  d <- dada(filt[kept][i], err = errs[[b]], pool = "pseudo",
            multithread = threads, verbose = FALSE)
  if (length(i) == 1) d <- list(d)
  ids <- meta$sample_id[kept][i]
  names(d) <- ids
  dadas[ids] <- d
}
seqtab <- makeSequenceTable(dadas)
say("%d ASVs before chimera removal", ncol(seqtab))

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = threads, verbose = FALSE)
chim_reads <- sum(seqtab.nochim) / sum(seqtab)
if (chim_reads < 0.80)
  say(paste("WARNING chimera removal discarded %.1f%% of reads (healthy runs keep",
            ">90%%). Treat low-prevalence ASVs with care."), 100 * (1 - chim_reads))
say("chimera removal: %d of %d ASVs kept, %.1f%% of reads retained",
    ncol(seqtab.nochim), ncol(seqtab), 100 * chim_reads)

# ---- per-sample, per-step counts -------------------------------------------
getn <- function(x) vapply(x, function(y) sum(getUniques(y)), 0)
track <- data.frame(
  sample_id = meta$sample_id,
  raw = ft[, "reads.in"], filtered = ft[, "reads.out"],
  denoised = getn(dadas)[meta$sample_id],
  nonchimeric = rowSums(seqtab.nochim)[match(meta$sample_id, rownames(seqtab.nochim))],
  sequencing_batch = meta$sequencing_batch, error_model = "not_analysed",
  row.names = NULL, stringsAsFactors = FALSE)
track$error_model[match(meta$sample_id[kept], track$sample_id)] <- model_of
for (col in c("raw", "filtered", "denoised", "nonchimeric"))
  track[[col]][is.na(track[[col]])] <- 0
track$retention_filtered <- track$filtered / pmax(track$raw, 1)
track$retention_overall  <- track$nonchimeric / pmax(track$raw, 1)
write.csv(track, res("read_counts_per_step.csv"), row.names = FALSE)

analysed <- track$filtered >= min_reads
empty <- track$sample_id[analysed & track$nonchimeric == 0]
frac_empty <- length(empty) / max(sum(analysed), 1)
say("%d of %d analysed sample(s) ended empty", length(empty), sum(analysed))
if (frac_empty > 0.05)
  fail(paste("%.1f%% of analysed samples produced nothing (%d of %d). That is a",
             "systematic failure, not a few bad samples; see",
             "read_counts_per_step.csv."), 100 * frac_empty, length(empty), sum(analysed))

# ---- taxonomy --------------------------------------------------------------
taxa <- assignTaxonomy(seqtab.nochim, silva, multithread = threads,
                       minBoot = 50, tryRC = TRUE)
colnames(taxa)[colnames(taxa) == "Species"] <- "Species_naive_bayes"
if (!is.na(silva_sp) && file.exists(silva_sp)) {
  taxa <- addSpecies(taxa, silva_sp)
  colnames(taxa)[colnames(taxa) == "Species"] <- "Species_exact_match"
}

# ---- outputs ---------------------------------------------------------------
asv_id <- sprintf("ASV%04d", seq_len(ncol(seqtab.nochim)))
seqs <- colnames(seqtab.nochim)
writeXStringSet(setNames(DNAStringSet(seqs), asv_id), res("asv_sequences.fasta"))

counts <- seqtab.nochim; colnames(counts) <- asv_id
write.csv(data.frame(sample_id = rownames(counts), counts, check.names = FALSE),
          res("asv_counts.csv"), row.names = FALSE)
denom <- rowSums(counts)
if (any(denom == 0))
  fail("%d sample(s) have zero counts and cannot be put on a relative scale: %s",
       sum(denom == 0), paste(rownames(counts)[denom == 0], collapse = ", "))
rel <- sweep(counts, 1, denom, "/")
stopifnot(all(abs(rowSums(rel) - 1) < 1e-9))
write.csv(data.frame(sample_id = rownames(rel), rel, check.names = FALSE),
          res("asv_relative_abundance.csv"), row.names = FALSE)
saveRDS(list(seqtab = seqtab.nochim, asv_id = asv_id, taxa = taxa),
        res("seqtab_nochim.rds"))
write.csv(data.frame(asv_id, sequence = seqs, length = nchar(seqs),
                     as.data.frame(taxa, stringsAsFactors = FALSE), row.names = NULL),
          res("asv_taxonomy.csv"), row.names = FALSE)
meta$reads_in_matrix <- denom[match(meta$sample_id, rownames(counts))]
meta$reads_in_matrix[is.na(meta$reads_in_matrix)] <- 0
write.csv(meta, res("sample_metadata.csv"), row.names = FALSE)

writeLines(c(
  sprintf("paper: %s", paper),
  sprintf("region: %s", region),
  "library: SINGLE-END. Analysed with step3_16s_single.R, not the paired script.",
  if (amp_status == "verified")
    "primers: confirmed against the paper's Methods and re-checked in the reads."
  else paste("primers: MEASURED FROM THE READS ONLY - nobody has read this paper's",
             "Methods."),
  if (grouped) "grouped analysis: supported by the step-2 metadata gate."
  else paste("grouped analysis: NOT SUPPORTED - the step-2 gate passed this paper for",
             "quantification only; group labels were never confirmed."),
  sprintf("samples: %d, of which %d have reads in the matrix",
          nrow(meta), sum(meta$reads_in_matrix > 0)),
  sprintf("ASVs: %d", ncol(seqtab.nochim)),
  "",
  "parameters, and why:",
  sprintf("  trimLeft = %d   primer length (0 = reads arrived pre-trimmed)", trim_left),
  sprintf("  truncLen = %d   %s", trunc_len, basis),
  sprintf("  the ASV is therefore a %d bp fragment of the 5' end of the amplicon",
          trunc_len - trim_left),
  sprintf("  maxEE    = %g", max_ee),
  "  pool     = \"pseudo\"",
  sprintf("  error models: %d; batches from the flow-cell id in the FASTQ headers",
          length(unique(model_of))),
  sprintf("  chimeras: consensus; %d of %d ASVs kept, %.1f%% of reads retained%s",
          ncol(seqtab.nochim), ncol(seqtab), 100 * chim_reads,
          if (chim_reads < 0.80) "  <-- ABNORMAL" else ""),
  "",
  "WHAT SINGLE-END COSTS - read this before using the table:",
  "  The ASV covers only the 5' part of the amplicon. Two taxa that differ only",
  "  downstream of the cut collapse into one ASV here, so richness is",
  "  underestimated and taxonomic resolution is lower than the same region",
  "  sequenced paired.",
  "  These ASVs CANNOT be matched sequence-to-sequence against ASVs from a",
  "  paired run of the same region - they are fragments of different length.",
  "  Cross-paper merging must happen at a taxonomic rank (genus), never on the",
  "  ASV sequence.",
  "",
  "table orientation: one ROW per sample, one COLUMN per ASV.",
  "relative abundance denominator: each sample's own total over the ASVs in",
  "  asv_counts.csv (after chimera removal). Each ROW sums to 1. NOT raw reads;",
  "  per-step losses are in read_counts_per_step.csv."
), res("README.txt"))

say("done: %d ASVs x %d samples -> %s", ncol(seqtab.nochim), nrow(seqtab.nochim), res(""))
