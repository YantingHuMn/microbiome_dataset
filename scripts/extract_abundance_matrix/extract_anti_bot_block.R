input_file <- "/hickory/proj/didonglab/dataset/virus/Database/metadata/blocked_manual_download.tsv"

df <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = "character"
)

stopifnot("reason" %in% names(df))
stopifnot("attempted_files" %in% names(df))

# 0：reason 和 dataset_id 都是空的 -- 前面 stage 压根没给这篇论文关联到数据源，
#     不是下载失败，需要你自己去确认这篇论文的数据到底存不存在、放在哪
is_no_dataset <- is.na(df$reason) | df$reason == ""

# 1：临时性问题，清掉重跑大概率解决
retry_pattern <- paste(
  "timeout", "BrokenPipeError", "Internal Server Error",
  "Service Temporarily Unavailable", "end-of-stream marker",
  "ConnectionResetError", "ConnectionError", "list_deposit_files error",
  sep = "|"
)
is_retry <- !is_no_dataset & !is.na(df$reason) &
  grepl(retry_pattern, df$reason, ignore.case = TRUE)

# 2：根本不是丰度数据（参考分类数据库 / 原始测序 reads），可以忽略
skip_pattern <- paste(
  "classifier\\.qza$",
  "_R[12][_.]", "_mpg_data_",
  "\\.fastq(\\.gz)?$", "\\.fq(\\.gz)?$",
  "_a[a-z]\\.gz$",
  sep = "|"
)
is_skip <- !is_no_dataset & !is_retry &
  !is.na(df$attempted_files) &
  grepl(skip_pattern, df$attempted_files, ignore.case = TRUE)

# 3：真正需要人工确认（反爬拦截、明确拒绝访问、文件格式解析不了）
# 原来的 anti_bot 规则并进这一类，再加上 Forbidden 和格式判断不了的情况
manual_pattern <- paste(
  "anti[ _-]?bot", "Forbidden", "format cannot be determined",
  sep = "|"
)
is_manual <- !is_no_dataset & !is_retry & !is_skip &
  !is.na(df$reason) &
  grepl(manual_pattern, df$reason, ignore.case = TRUE)

# 4：剩下没归类的
is_uncategorized <- !is_no_dataset & !is_retry & !is_skip & !is_manual

out_dir <- dirname(input_file)

write_group <- function(mask, filename) {
  write.table(
    df[mask, , drop = FALSE],
    file = file.path(out_dir, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

write_group(is_no_dataset, "blocked_0_no_dataset_linked.tsv")
write_group(is_retry, "blocked_1_retry_transient.tsv")
write_group(is_skip, "blocked_2_skip_not_real_data.tsv")
write_group(is_manual, "blocked_3_manual_check.tsv")
write_group(is_uncategorized, "blocked_4_uncategorized.tsv")

report <- function(mask, label) {
  cat(label, "rows:", sum(mask), "\n")
  cat(label, "unique paper:", length(unique(df[mask, , drop = FALSE]$paper_id)), "\n\n")
}

report(is_no_dataset,"0_no_dataset_linked")
report(is_retry, "1_retry_transient")
report(is_skip, "2_skip_not_real_data")
report(is_manual, "3_manual_check")
report(is_uncategorized, "4_uncategorized")

stopifnot(sum(is_no_dataset, is_retry, is_skip, is_manual, is_uncategorized) == nrow(df))
cat("sanity check passed -- every row classified exactly once, total:", nrow(df), "\n")