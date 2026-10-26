#!/usr/bin/env Rscript
#
# DADA2 denoising for one paper. Called by step3_abundance.py, which has already
# verified the primers and the amplicon length against the paper and passes them
# in; nothing here is inferred from file names or hard-coded per paper.
#
# truncLen is CHOSEN FROM THE DATA, in three passes that all have to agree:
#   A  expected-error profiles of both mates (how much of each read survives
#      maxEE at each possible cut) - not a Q25 line, which says nothing about
#      how many reads actually pass a filter
#   B  the merge budget: a candidate is only considered if what is left after
#      trimLeft still spans the measured amplicon with minOverlap to spare
#   C  the real pipeline on a probe subset: filter -> denoise -> merge, scored
#      on reads that survive all the way to a merged pair
# --trunc-len overrides all of it (and is recorded as a manual override).
#
# Coordinates: truncLen is in RAW read bases. Verified against dada2's
# fastqPairedFilter: it sets startF = trimLeft+1 and endF = truncLen - trimLeft,
# applies narrow(start = startF) FIRST and narrow(end = endF) after, and keeps
# only reads whose ORIGINAL length is >= truncLen. So a read contributes bases
# trimLeft+1 .. truncLen and the retained length is truncLen - trimLeft. The
# merge budget below is written in those same raw units, which is why trimLeft
# appears in it.

suppressPackageStartupMessages({
  library(dada2); library(Biostrings); library(ShortRead)
})

# ---- arguments -------------------------------------------------------------
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
trim_left    <- int2(getarg("--trim-left", "0,0"))
max_trunc    <- int2(getarg("--max-trunc", "0,0"))
insert_p99   <- as.integer(getarg("--insert-p99"))
min_overlap  <- as.integer(getarg("--min-overlap", "12"))
threads      <- as.integer(getarg("--threads", "4"))
silva        <- getarg("--silva")
silva_sp     <- getarg("--silva-species", NA)
max_ee       <- as.numeric(int2(getarg("--max-ee", "2,5")))
trunc_manual <- if (is.null(getarg("--trunc-len"))) NULL else int2(getarg("--trunc-len"))
probe_n      <- as.integer(getarg("--probe-samples", "6"))   # pass B/C subset
probe_reads  <- as.integer(getarg("--probe-reads", "40000"))
only_probe   <- "--probe-only" %in% args        # validate before the full run
min_reads    <- as.integer(getarg("--min-reads", "1000"))
min_merge    <- as.numeric(getarg("--min-merge-frac", "0.30"))
min_err_bases<- as.numeric(getarg("--min-error-bases", "1e8"))
grouped     <- getarg("--grouped-analysis", "1") == "1"
amp_status  <- getarg("--amplicon-status", "unknown")

stopifnot(!is.null(paper), !is.null(out_dir), !is.na(insert_p99))
dir.create(file.path(out_dir, "dada2"), recursive = TRUE, showWarnings = FALSE)
res <- function(...) file.path(out_dir, "dada2", ...)
say <- function(...) cat(sprintf(...), "\n", sep = "")
# sprintf() the message first: say("FAILED: ", msg) would make the message an
# unused argument of the format string and print nothing but the prefix.
fail <- function(...) { cat("FAILED: ", sprintf(...), "\n", sep = ""); quit(status = 1) }

# ---- samples ---------------------------------------------------------------
# samples_used.csv is written by step3_abundance.py with the technical controls
# (blanks, dilution series) already removed, so they can never reach the matrix.
sfile <- file.path(out_dir, "samples_used.csv")
if (!file.exists(sfile)) fail("%s missing; run step 3 preflight first", sfile)
meta <- read.csv(sfile, check.names = FALSE, fileEncoding = "UTF-8-BOM")
say("%s: %d biological samples", paper, nrow(meta))

# The FASTQs come from run_sample_mapping.csv, which step 2 built from the
# download manifest plus ENA. Globbing the run directory instead would be a
# guess: a single-end run has no _2 file, and a fallback to "any .fastq.gz"
# silently hands the SAME file back as both mates.
mapfile <- file.path(out_dir, "run_sample_mapping.csv")
if (!file.exists(mapfile)) fail("%s missing; run step 2 first", mapfile)
rmap <- read.csv(mapfile, check.names = FALSE, fileEncoding = "UTF-8-BOM")
rmap <- rmap[rmap$biosample_accession %in% meta$sample_id, , drop = FALSE]

files_for <- function(sid, mate) {
  r <- rmap[rmap$biosample_accession == sid & rmap$mate == mate, , drop = FALSE]
  r <- r[order(r$run_accession), , drop = FALSE]
  r$fastq_path
}
layout <- unique(toupper(rmap$library_layout))
if (!identical(layout, "PAIRED"))
  fail(paste("this script is for paired reads; run_sample_mapping.csv reports",
             "library_layout = %s. A single-end paper needs its own branch, not",
             "this one."), paste(layout, collapse = "/"))

# Every run of a sample is used. Taking only the first would silently halve a
# sample sequenced across two lanes, and step 2 has already said whether the
# runs may be pooled at all.
bad_merge <- meta$sample_id[meta$multi_run_handling ==
                            "MULTIPLE_EXPERIMENTS_do_not_merge_without_review"]
if (length(bad_merge))
  fail(paste("%d sample(s) have runs from different experiments and step 2",
             "flagged them for review; pooling them here would mix libraries: %s"),
       length(bad_merge), paste(head(bad_merge, 5), collapse = ", "))

cat_dir <- file.path(out_dir, "dada2", "pooled_runs")
pool_runs <- function(sid, mate, paths) {
  if (length(paths) == 1) return(paths)
  dir.create(cat_dir, showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(cat_dir, sprintf("%s_R%d.fastq.gz", sid, mate))
  if (!file.exists(dest)) {               # gzip members concatenate as-is
    if (!file.create(dest)) fail("could not create %s", dest)
    if (!all(file.append(dest, paths))) {
      unlink(dest)                        # never leave a half-written pool
      fail("could not pool runs for %s", sid)
    }
  }
  dest
}

fwd <- rev <- character(nrow(meta))
n_runs_used <- integer(nrow(meta))
for (i in seq_len(nrow(meta))) {
  sid <- meta$sample_id[i]
  f <- files_for(sid, "R1"); r <- files_for(sid, "R2")
  if (!length(f) || !length(r))
    fail("%s: run_sample_mapping.csv lists %d R1 and %d R2 file(s)", sid, length(f), length(r))
  if (length(f) != length(r))
    fail("%s: %d R1 files but %d R2 files - the mates do not pair",
         sid, length(f), length(r))
  missing <- c(f, r)[!file.exists(c(f, r))]
  if (length(missing)) fail("%s: file(s) listed but not on disk: %s", sid,
                            paste(basename(missing), collapse = ", "))
  if (any(duplicated(c(f, r))))
    fail("%s: the same file is listed as both mates", sid)
  n_runs_used[i] <- length(f)
  fwd[i] <- pool_runs(sid, 1, f); rev[i] <- pool_runs(sid, 2, r)
}
names(fwd) <- names(rev) <- meta$sample_id
meta$n_runs_used <- n_runs_used
if (any(n_runs_used > 1))
  say("%d sample(s) sequenced across several runs; all runs pooled per sample (%d runs in total)",
      sum(n_runs_used > 1), sum(n_runs_used))  # counted as runs, both mates each

# Cheap sanity check before an expensive run: the mates must start with the
# same read. filterAndTrim below also runs with matchIDs = TRUE, which drops
# any pair whose ids disagree further in.
first_id <- function(path) {
  con <- gzfile(path, "r"); on.exit(close(con))
  sub("[ /].*$", "", sub("^@", "", readLines(con, n = 1)))
}
for (i in seq_len(nrow(meta))) {
  a <- tryCatch(first_id(fwd[i]), error = function(e) NA_character_)
  b <- tryCatch(first_id(rev[i]), error = function(e) NA_character_)
  if (is.na(a) || is.na(b) || a != b)
    fail("%s: R1 and R2 do not start with the same read id (%s vs %s)",
         meta$sample_id[i], a, b)
}
say("R1/R2 first read ids match for all %d samples", nrow(meta))

# ---- sequencing batches ----------------------------------------------------
# A batch is a real flow cell, read out of the FASTQ header. Arbitrary groupings
# of samples are NOT batches, so when the headers have been stripped (as SRA
# often does) every sample goes into one model and that is said out loud.
flowcell_of <- function(path) {
  # One line, not the whole file: readFastq() would load every read into memory
  # just to reach the first header.
  con <- gzfile(path, "r"); on.exit(close(con))
  h <- readLines(con, n = 1)
  m <- regmatches(h, regexpr("[A-Za-z0-9_-]+:[0-9]+:[A-Za-z0-9-]+:[0-9]+:", h))
  if (!length(m)) return(NA_character_)
  paste(strsplit(m, ":")[[1]][1:3], collapse = ":")
}
batch <- vapply(fwd, function(f) tryCatch(flowcell_of(f), error = function(e) NA_character_), "")
if (anyNA(batch)) {
  say("batch: FASTQ headers carry no instrument/run/flow-cell id (stripped on submission).")
  say("       Using ONE error model for all samples - inventing batches from sample")
  say("       names would be worse than admitting the batch structure is unknown.")
  batch[] <- "unknown_single_model"
} else {
  say("batch: %d flow cell(s) found in the read headers: %s", length(unique(batch)),
      paste(sprintf("%s (n=%d)", names(table(batch)), table(batch)), collapse = ", "))
}
meta$sequencing_batch <- batch

# ---- pass A: expected-error profiles ---------------------------------------
# For each possible cut, how many reads would survive maxEE? Computed exactly
# from the quality strings rather than read off a median-quality curve.
ee_profile <- function(files, n = 40000) {
  cum <- NULL; tot <- 0
  for (f in files) {
    fq <- ShortRead::yield(ShortRead::FastqStreamer(f, n))
    q <- as(Biostrings::quality(fq), "matrix")
    # A read shorter than the cut is DISCARDED by filterAndTrim, so its missing
    # tail must score as a failure. Leaving it as 0 error would make an
    # over-long truncLen look like the best option on a run of mixed lengths.
    p <- 10^(-q / 10); p[is.na(p)] <- Inf
    ee <- t(apply(p, 1, cumsum))
    if (!is.null(cum)) {                       # runs may differ in read length
      w <- min(ncol(cum), ncol(ee))
      cum <- rbind(cum[, seq_len(w), drop = FALSE], ee[, seq_len(w), drop = FALSE])
    } else cum <- ee
    tot <- tot + nrow(ee)
    if (tot >= n) break
  }
  cum
}
set.seed(1)
# Spread the probe over the flow cells rather than over the sample list, so the
# chosen truncLen is not tuned to whichever run happens to sort first.
probe_ids <- unlist(lapply(split(meta$sample_id, meta$sequencing_batch), function(ids)
  ids[round(seq(1, length(ids), length.out = max(1, round(probe_n * length(ids) / nrow(meta)))))]))
probe_ids <- unique(probe_ids)
say("probe subset (%d samples over %d batch(es)): %s", length(probe_ids),
    length(unique(meta$sequencing_batch[match(probe_ids, meta$sample_id)])),
    paste(probe_ids, collapse = ", "))

eeF <- ee_profile(fwd[probe_ids], probe_reads)
eeR <- ee_profile(rev[probe_ids], probe_reads)
pass_at <- function(ee, L, cap) if (L > ncol(ee)) NA else mean(ee[, L] <= cap)
readlen <- c(ncol(eeF), ncol(eeR))
cap <- pmin(max_trunc, readlen)
say("read length %d/%d; truncLen capped at %d/%d", readlen[1], readlen[2], cap[1], cap[2])

# ---- pass B: candidates that can still merge -------------------------------
need <- insert_p99 + min_overlap + sum(trim_left)
say("merge budget: truncLen[1]+truncLen[2] >= %d (insert p99 %d + overlap %d + trimLeft %d+%d)",
    need, insert_p99, min_overlap, trim_left[1], trim_left[2])
if (sum(cap) < need) fail("reads are too short: best possible sum %d < required %d",
                          sum(cap), need)

grid <- expand.grid(f = seq(cap[1], max(trim_left[1] + 50, cap[1] - 80), by = -10),
                    r = seq(cap[2], max(trim_left[2] + 50, cap[2] - 100), by = -10))
grid <- grid[grid$f + grid$r >= need, ]
grid$keepF <- vapply(grid$f, function(L) pass_at(eeF, L, max_ee[1]), 0)
grid$keepR <- vapply(grid$r, function(L) pass_at(eeR, L, max_ee[2]), 0)
grid$keep  <- grid$keepF * grid$keepR          # both mates must survive
grid$slack <- grid$f + grid$r - need
grid <- grid[order(-grid$keep, -grid$slack), ]
write.csv(grid, res("trunclen_pass_a_b.csv"), row.names = FALSE)
say("pass A/B: %d admissible candidates; best estimated filter retention %.1f%%",
    nrow(grid), 100 * grid$keep[1])

cands <- grid[seq_len(min(5, nrow(grid))), c("f", "r")]
# Keep one deliberately short-and-clean and one long-and-permissive candidate so
# pass C is choosing between real alternatives, not re-ranking near-identical ones.
cands <- unique(rbind(cands, grid[nrow(grid), c("f", "r")]))

# ---- pass C: the real pipeline on the probe subset -------------------------
tmp <- file.path(tempdir(), "probe"); dir.create(tmp, showWarnings = FALSE)
score_candidate <- function(f, r) {
  ff <- file.path(tmp, paste0(probe_ids, "_F.fq.gz"))
  rr <- file.path(tmp, paste0(probe_ids, "_R.fq.gz"))
  ft <- tryCatch(filterAndTrim(fwd[probe_ids], ff, rev[probe_ids], rr,
                               truncLen = c(f, r), trimLeft = trim_left,
                               maxEE = max_ee, truncQ = 2, maxN = 0, rm.phix = TRUE,
                               compress = TRUE, multithread = threads,
                               verbose = FALSE),
                 error = function(e) NULL)
  if (is.null(ft) || !any(file.exists(ff))) return(c(filt = 0, merged = 0))
  ok <- file.exists(ff) & file.exists(rr)
  eF <- learnErrors(ff[ok], multithread = threads, verbose = FALSE)
  eR <- learnErrors(rr[ok], multithread = threads, verbose = FALSE)
  dF <- dada(ff[ok], err = eF, multithread = threads, verbose = FALSE)
  dR <- dada(rr[ok], err = eR, multithread = threads, verbose = FALSE)
  mg <- mergePairs(dF, ff[ok], dR, rr[ok], minOverlap = min_overlap,
                   maxMismatch = 0, verbose = FALSE)
  nin <- sum(ft[, "reads.in"]); nmg <- sum(vapply(mg, function(m) sum(m$abundance), 0))
  unlink(c(ff, rr))
  c(filt = sum(ft[, "reads.out"]) / max(nin, 1), merged = nmg / max(nin, 1))
}

if (!is.null(trunc_manual)) {
  trunc_len <- trunc_manual
  basis <- sprintf("MANUAL OVERRIDE via --trunc-len %d,%d", trunc_len[1], trunc_len[2])
  if (sum(trunc_len) < need)
    say("WARNING: the manual truncLen sum %d is below the merge budget %d.",
        sum(trunc_len), need)
  scores <- NULL
} else {
  say("pass C: running filter -> denoise -> merge on the probe subset for %d candidates",
      nrow(cands))
  scores <- do.call(rbind, lapply(seq_len(nrow(cands)), function(i) {
    s <- score_candidate(cands$f[i], cands$r[i])
    say("  truncLen c(%d,%d): %.1f%% pass the filter, %.1f%% reach a merged pair",
        cands$f[i], cands$r[i], 100 * s["filt"], 100 * s["merged"])
    data.frame(f = cands$f[i], r = cands$r[i], filt = s["filt"], merged = s["merged"])
  }))
  scores <- scores[order(-scores$merged, -(scores$f + scores$r)), ]
  write.csv(scores, res("trunclen_pass_c.csv"), row.names = FALSE)
  if (scores$merged[1] <= 0) fail("no candidate truncLen merged any reads")
  trunc_len <- c(scores$f[1], scores$r[1])
  basis <- sprintf(paste("chosen by merged-read retention on %d probe samples:",
                         "%.1f%% of raw reads reach a merged pair (%.1f%% pass the filter)"),
                   length(probe_ids), 100 * scores$merged[1], 100 * scores$filt[1])
}
say("truncLen = c(%d, %d)  [%s]", trunc_len[1], trunc_len[2], basis)
if (only_probe) {
  writeLines(c(sprintf("truncLen\t%d,%d", trunc_len[1], trunc_len[2]), basis),
             res("trunclen_choice.txt"))
  say("--probe-only: stopping after validation, nothing was run on the full set.")
  quit(status = 0)
}

# ---- filter all samples ----------------------------------------------------
fdir <- file.path(out_dir, "dada2", "filtered")
# Clear it first. filterAndTrim writes no file for a sample whose reads are all
# removed, so a leftover file from an earlier run at the same path would be
# picked up by the file.exists() test below and denoised as if it were this
# run's output - with the previous run's truncLen baked into it.
unlink(fdir, recursive = TRUE)
dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
filtF <- file.path(fdir, paste0(meta$sample_id, "_F.fq.gz"))
filtR <- file.path(fdir, paste0(meta$sample_id, "_R.fq.gz"))
# matchIDs is deliberately OFF. Its auto-detection picks the whitespace field
# holding 6 colons, which for these SRA headers is
# "M06129:119:000000000-JL26M:1:1101:21802:1004/1" - the /1 and /2 suffixes mean
# no pair ever matches and the filter silently discards every read (verified:
# it removed 100% of reads on all six probe samples). Pairing is instead
# checked above from run_sample_mapping.csv and the mates' first read ids, and
# filterAndTrim itself errors on streams of unequal length.
ft <- filterAndTrim(fwd, filtF, rev, filtR, truncLen = trunc_len, trimLeft = trim_left,
                    maxEE = max_ee, truncQ = 2, maxN = 0, rm.phix = TRUE,
                    compress = TRUE, multithread = threads, verbose = FALSE)
rownames(ft) <- meta$sample_id
kept <- file.exists(filtF) & file.exists(filtR) & ft[, "reads.out"] >= min_reads
if (!any(kept)) fail("every sample fell below %d reads after filtering", min_reads)
if (any(!kept)) say(paste("%d sample(s) below %d reads after filtering, carried into",
                          "the report as empty rather than silently dropped: %s"),
                    sum(!kept), min_reads, paste(meta$sample_id[!kept], collapse = ", "))

# ---- error models, per real batch -----------------------------------------
# A batch with too little data gives a worse model than the pooled one, so small
# batches are pooled and the decision is recorded per sample.
bat <- meta$sequencing_batch[kept]
# Judged on sequenced bases, not on how many samples a batch happens to hold.
# learnErrors is called once PER MATE, so the budget it sees is reads x that
# mate's retained length - not the two mates added together. The forward read
# is the shorter of the two here, so it is the binding constraint. (Checked
# against learnErrors' own report: 284257 reads x 214 bp = 6.08e7, exactly the
# "total bases ... will be used" it printed.)
per_read <- min(trunc_len - trim_left)
bases <- tapply(ft[kept, "reads.out"], bat, sum) * per_read
big <- names(bases)[bases >= min_err_bases]
model_of <- ifelse(bat %in% big, bat, "pooled_small_batches")
small <- setdiff(names(bases), big)
if (length(small)) {
  pooled <- sum(bases[small])
  if (length(small) == 1 || pooled < min_err_bases) {
    # Pooling one small batch with nothing, or with other small batches that
    # together still fall short, does not buy a better model. Say so instead of
    # letting the label "pooled" imply it was fixed.
    model_of <- ifelse(bat %in% big, bat, paste0("own_model_below_threshold:", bat))
    say(paste("WARNING error models: batch(es) %s hold %.2e bases, under the %.0e",
              "needed for a well-estimated model. They keep their OWN model rather",
              "than being mixed into another flow cell's; treat ASVs seen only in",
              "those samples with care."),
        paste(small, collapse = ", "), pooled, min_err_bases)
  } else {
    say("error models: %s have enough data alone; %s are pooled to reach %.2e bases",
        paste(big, collapse = ", "), paste(small, collapse = ", "), pooled)
  }
}
for (b in names(bases))
  say("  batch %s: %d sample(s), %.2e bases -> model %s", b, sum(bat == b), bases[[b]],
      unique(model_of[bat == b]))
errF <- errR <- list()
for (b in unique(model_of)) {
  i <- which(model_of == b)
  errF[[b]] <- learnErrors(filtF[kept][i], multithread = threads, verbose = FALSE)
  errR[[b]] <- learnErrors(filtR[kept][i], multithread = threads, verbose = FALSE)
  say("error model %s: %d samples", b, length(i))
}
pdf(res("error_models.pdf")); for (b in names(errF)) {
  print(plotErrors(errF[[b]], nominalQ = TRUE) + ggplot2::ggtitle(paste("R1", b)))
  print(plotErrors(errR[[b]], nominalQ = TRUE) + ggplot2::ggtitle(paste("R2", b)))
}; invisible(dev.off())

# ---- denoise and merge -----------------------------------------------------
# pool = "pseudo": sensitivity close to full pooling for variants that are rare
# in any one sample but present across the study, without full pooling's cost.
dadaF <- dadaR <- mergers <- list()
for (b in unique(model_of)) {
  i <- which(model_of == b)
  dF <- dada(filtF[kept][i], err = errF[[b]], pool = "pseudo",
             multithread = threads, verbose = FALSE)
  dR <- dada(filtR[kept][i], err = errR[[b]], pool = "pseudo",
             multithread = threads, verbose = FALSE)
  if (length(i) == 1) { dF <- list(dF); dR <- list(dR) }
  mg <- mergePairs(dF, filtF[kept][i], dR, filtR[kept][i],
                   minOverlap = min_overlap, maxMismatch = 0, verbose = FALSE)
  if (length(i) == 1) mg <- list(mg)
  ids <- meta$sample_id[kept][i]
  names(dF) <- names(dR) <- names(mg) <- ids
  dadaF[ids] <- dF; dadaR[ids] <- dR; mergers[ids] <- mg
}

seqtab <- makeSequenceTable(mergers)
say("%d ASVs before length filtering and chimera removal", ncol(seqtab))

# ---- length filter ---------------------------------------------------------
# Keyed to the amplicon length MEASURED from these reads, with a tolerance wide
# enough for real taxon-to-taxon variation; everything dropped is written out.
lens <- nchar(colnames(seqtab))
lo <- insert_p99 - 60; hi <- insert_p99 + 20
in_range <- lens >= lo & lens <= hi
say("length filter: keeping ASVs of %d-%d bp (%d of %d ASVs, %.2f%% of reads)",
    lo, hi, sum(in_range), length(lens), 100 * sum(seqtab[, in_range]) / sum(seqtab))
write.csv(data.frame(length = as.integer(names(table(lens))),
                     n_asvs = as.integer(table(lens)),
                     reads = tapply(colSums(seqtab), lens, sum),
                     kept = as.integer(names(table(lens))) >= lo &
                            as.integer(names(table(lens))) <= hi),
          res("asv_length_distribution.csv"), row.names = FALSE)
seqtab <- seqtab[, in_range, drop = FALSE]

# ---- chimeras --------------------------------------------------------------
# "consensus" calls a sequence chimeric only when it is flagged in a majority of
# the samples it appears in, so a single odd sample cannot remove a real ASV.
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = threads, verbose = FALSE)
chim_reads <- sum(seqtab.nochim) / sum(seqtab)
# Losing a large share of READS to chimera removal is a warning about the data,
# not a routine step. It usually means the merge was marginal: pairs joined on
# a short, unverified overlap look chimeric downstream. The number travels into
# the README so nobody has to dig it out of a log.
if (chim_reads < 0.80)
  say(paste("WARNING chimera removal discarded %.1f%% of merged reads (healthy",
            "runs keep >90%%). The rare end of this ASV table is unreliable;",
            "check the merge overlap before using low-prevalence ASVs."),
      100 * (1 - chim_reads))
say("chimera removal: %d of %d ASVs kept, %.1f%% of reads retained",
    ncol(seqtab.nochim), ncol(seqtab), 100 * sum(seqtab.nochim) / sum(seqtab))

# ---- per-sample, per-step read counts --------------------------------------
getn <- function(x) vapply(x, function(y) sum(getUniques(y)), 0)
ids <- rownames(seqtab.nochim)
track <- data.frame(
  sample_id   = meta$sample_id,
  raw         = ft[, "reads.in"],
  filtered    = ft[, "reads.out"],
  denoised_f  = getn(dadaF)[meta$sample_id],
  denoised_r  = getn(dadaR)[meta$sample_id],
  merged      = getn(mergers)[meta$sample_id],
  nonchimeric = rowSums(seqtab.nochim)[match(meta$sample_id, ids)],
  sequencing_batch = meta$sequencing_batch,
  error_model = "not_analysed",
  row.names = NULL, stringsAsFactors = FALSE)
track$error_model[match(meta$sample_id[kept], track$sample_id)] <- model_of
# Only the read counts default to zero; blanking a text column to "0" would
# claim a sample had an error model called 0.
for (col in c("raw", "filtered", "denoised_f", "denoised_r", "merged", "nonchimeric"))
  track[[col]][is.na(track[[col]])] <- 0
track$retention_filtered   <- track$filtered    / pmax(track$raw, 1)
track$retention_merged     <- track$merged      / pmax(track$filtered, 1)
track$retention_overall    <- track$nonchimeric / pmax(track$raw, 1)
write.csv(track, res("read_counts_per_step.csv"), row.names = FALSE)

# ---- hard failures ---------------------------------------------------------
# An empty or collapsing merge is a result about the data, never something to
# paper over with a looser parameter.
# "Empty" means a sample that WAS denoised and produced nothing. A sample that
# fell below min_reads at the filter never entered denoising; it is already
# reported as dropped above, and counting it here confused "too few reads to
# start with" (ordinary) with "merging collapsed" (a real anomaly).
analysed <- track$filtered >= min_reads
empty <- track$sample_id[analysed & track$nonchimeric == 0]
weak  <- track$sample_id[analysed & track$retention_merged < min_merge]
med   <- median(track$retention_merged[analysed])
frac_empty <- length(empty) / max(sum(analysed), 1)
say("merge retention: median %.1f%%, %d sample(s) below %.0f%%, %d of %d analysed empty",
    100 * med, length(weak), 100 * min_merge, length(empty), sum(analysed))
# One bad sample in a hundred is a bad sample; a systematic collapse is a
# different thing and is the only reason to throw the whole paper away.
if (length(empty) && frac_empty <= 0.05 && med >= min_merge)
  say(paste("NOTE %d analysed sample(s) ended with zero reads (%.1f%% of those",
            "analysed): %s. They stay in read_counts_per_step.csv with zero",
            "counts and are absent from the matrix."),
      length(empty), 100 * frac_empty, paste(head(empty, 5), collapse = ", "))
if (med < min_merge || frac_empty > 0.05) {
  writeLines(c(sprintf("empty after merging: %s", paste(empty, collapse = ", ")),
               sprintf("below %.0f%% merge retention: %s", 100 * min_merge,
                       paste(weak, collapse = ", ")),
               sprintf("median merge retention: %.3f", med)),
             res("MERGE_FAILURE.txt"))
  fail(paste("merging is anomalous (median retention %.1f%%, %d empty samples of those analysed).",
             "Counts are in read_counts_per_step.csv; see MERGE_FAILURE.txt.",
             "Do NOT relax minOverlap to make this pass - check the amplicon",
             "length and the truncLen choice first."), 100 * med, length(empty))
}

# ---- taxonomy --------------------------------------------------------------
taxa <- assignTaxonomy(seqtab.nochim, silva, multithread = threads,
                       minBoot = 50, tryRC = TRUE)
# The trainset's own species column is a naive-Bayes call, which DADA2's own
# documentation warns is not reliable at species level for 16S. addSpecies is
# an exact full-length match. They are kept as two clearly named columns rather
# than one ambiguous "Species" and a silent "Species.1" duplicate.
colnames(taxa)[colnames(taxa) == "Species"] <- "Species_naive_bayes"
if (!is.na(silva_sp) && file.exists(silva_sp)) {
  taxa <- addSpecies(taxa, silva_sp)      # exact match only; never an alignment
  colnames(taxa)[colnames(taxa) == "Species"] <- "Species_exact_match"
}

# ---- outputs ---------------------------------------------------------------
asv_id <- sprintf("ASV%04d", seq_len(ncol(seqtab.nochim)))
seqs <- colnames(seqtab.nochim)
writeXStringSet(setNames(DNAStringSet(seqs), asv_id), res("asv_sequences.fasta"))

# One row per SAMPLE, one column per ASV - the orientation the rest of the
# database uses, and the one every R/Python stats package expects.
counts <- seqtab.nochim
colnames(counts) <- asv_id
write.csv(data.frame(sample_id = rownames(counts), counts, check.names = FALSE),
          res("asv_counts.csv"), row.names = FALSE)

# Relative abundance: each sample's counts divided by that SAMPLE's own total of
# the ASVs in this table (post length-filter, post-chimera). It is not divided
# by raw reads, so each ROW sums to exactly 1 and samples remain comparable;
# what each sample lost on the way is in read_counts_per_step.csv.
denom <- rowSums(counts)
# A zero-total row cannot be expressed as a relative abundance. Such a sample
# should already have been caught as "empty after merging"; if one survives to
# here the table would quietly break its own "every row sums to 1" contract.
if (any(denom == 0))
  fail("%d sample(s) have zero counts after chimera removal and cannot be put on a relative scale: %s",
       sum(denom == 0), paste(rownames(counts)[denom == 0], collapse = ", "))
rel <- sweep(counts, 1, denom, "/")
stopifnot(all(abs(rowSums(rel) - 1) < 1e-9))
write.csv(data.frame(sample_id = rownames(rel), rel, check.names = FALSE),
          res("asv_relative_abundance.csv"), row.names = FALSE)

# The count table as an R object, so the output stage can be redone without
# repeating the denoising run.
saveRDS(list(seqtab = seqtab.nochim, asv_id = asv_id, taxa = taxa),
        res("seqtab_nochim.rds"))

# An ASV is a sequence variant, not a species: two species can share an
# identical V3-V4 fragment and one genome can carry several distinct copies.
write.csv(data.frame(asv_id, sequence = seqs, length = nchar(seqs),
                     as.data.frame(taxa, stringsAsFactors = FALSE), row.names = NULL),
          res("asv_taxonomy.csv"), row.names = FALSE)

meta$reads_in_matrix <- denom[match(meta$sample_id, rownames(counts))]
meta$reads_in_matrix[is.na(meta$reads_in_matrix)] <- 0
write.csv(meta, res("sample_metadata.csv"), row.names = FALSE)

writeLines(c(
  sprintf("paper: %s", paper),
  sprintf("region: %s", region),
  if (amp_status == "verified")
    "primers: confirmed against the paper's Methods and re-checked in the reads."
  else paste("primers: MEASURED FROM THE READS ONLY - nobody has read this paper's",
             "Methods. The primer pair was identified by matching it in the reads and",
             "the amplicon length was measured from read-pair overlap, which is what",
             "the denoising parameters need; but it has not been cross-checked against",
             "what the authors say they did."),
  if (grouped) "grouped analysis: supported by the step-2 metadata gate."
  else paste("grouped analysis: NOT SUPPORTED. The step-2 gate passed this paper",
             "for quantification only - the per-sample matrices are sound, but the",
             "group labels were never confirmed, so no case/control or treatment",
             "comparison may be drawn from these tables."),
  sprintf("samples: %d biological samples, %d with reads in the matrix",
          nrow(meta), sum(meta$reads_in_matrix > 0)),
  sprintf("ASVs: %d", ncol(seqtab.nochim)),
  "",
  "parameters, and why:",
  sprintf("  trimLeft   = c(%d, %d)   primer lengths from the paper (0 = reads arrived pre-trimmed)",
          trim_left[1], trim_left[2]),
  sprintf("  truncLen   = c(%d, %d)   %s", trunc_len[1], trunc_len[2], basis),
  "  truncLen is in RAW read bases. filterAndTrim removes trimLeft first and then",
  "  truncates to truncLen - trimLeft, keeping only reads at least truncLen long,",
  sprintf("  so the retained biological length is %d / %d bp.",
          trunc_len[1] - trim_left[1], trunc_len[2] - trim_left[2]),
  sprintf("  merge budget: %d + %d = %d >= %d (insert p99 %d + overlap %d + trimLeft %d)",
          trunc_len[1], trunc_len[2], sum(trunc_len), need, insert_p99, min_overlap,
          sum(trim_left)),
  sprintf("  maxEE      = c(%g, %g)   R2 is allowed more error than R1 because the",
          max_ee[1], max_ee[2]),
  "               second mate is always the worse read; a single strict cutoff",
  "               would discard good pairs for R2's tail alone.",
  "  pool       = \"pseudo\"   recovers variants that are rare in any one sample",
  "               but present across the study, at a fraction of full pooling's cost.",
  sprintf("  error models: %d, one per sequencing batch with >= %.0e bases; smaller",
          length(unique(model_of)), min_err_bases),
  "               batches are pooled. Batches come from the flow-cell id in the",
  "               FASTQ headers, never from sample names.",
  "  chimeras   = removeBimeraDenovo(method = \"consensus\")   a sequence is removed",
  "               only if flagged in most of the samples it appears in.",
  "",
  "filters applied, in order:",
  sprintf("  1. maxN = 0, truncQ = 2, rm.phix = TRUE, maxEE = c(%g, %g)", max_ee[1], max_ee[2]),
  sprintf("  2. samples under %d reads after filtering are reported with zero counts", min_reads),
  sprintf("  3. ASVs outside %d-%d bp dropped (measured amplicon length +/- real variation)", lo, hi),
  sprintf("  4. consensus chimera removal: %d of %d ASVs kept, %.1f%% of merged reads retained%s",
          ncol(seqtab.nochim), ncol(seqtab), 100 * chim_reads,
          if (chim_reads < 0.80) "  <-- ABNORMAL, see warning above" else ""),
  "",
  "table orientation:",
  "  one ROW per sample, one COLUMN per ASV, first column sample_id.",
  "",
  "relative abundance denominator:",
  "  each sample's own total over the ASVs in asv_counts.csv (after the length",
  "  filter and chimera removal). Each ROW sums to 1. NOT raw reads: per-step",
  "  losses are in read_counts_per_step.csv.",
  "",
  "files: asv_counts.csv, asv_relative_abundance.csv, asv_sequences.fasta,",
  "       seqtab_nochim.rds,",
  "       asv_taxonomy.csv, sample_metadata.csv, read_counts_per_step.csv,",
  "       asv_length_distribution.csv, trunclen_pass_a_b.csv, trunclen_pass_c.csv,",
  "       error_models.pdf"
), res("README.txt"))

say("done: %d ASVs x %d samples -> %s", ncol(seqtab.nochim), nrow(seqtab.nochim), res(""))
