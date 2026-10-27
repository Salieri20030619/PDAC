#!/usr/bin/env Rscript
# ============================================================
# GEO series_matrix 解析 (手写, 不依赖 GEOquery)
# 
# 输入: 5 个 GEO 队列的 *_series_matrix.txt.gz
# 输出: 每个队列 2 个文件
#   {GSE}.expr.tsv.gz    - probe/gene × sample 表达矩阵 (原始 probe 级)
#   {GSE}.clinical.tsv   - sample × clinical fields 临床表
# 
# 注意: 表达矩阵这步保留原始 probe ID, gene symbol 映射在下一步 (3_probe2gene_geo.R)
#       因为 5 个队列平台不同, 需要分别下载对应 GPL 注释
# ============================================================

suppressPackageStartupMessages({
    library(data.table)
})

BULK <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk"

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# 通用 series_matrix 解析函数
# ============================================================
parse_series_matrix <- function(gse_id, gse_dir) {
    log_msg(sprintf("=== %s ===", gse_id))
    sm_file <- file.path(gse_dir, sprintf("%s_series_matrix.txt.gz", gse_id))
    
    if (!file.exists(sm_file)) {
        log_msg(sprintf("  ⚠️ 文件不存在: %s", sm_file))
        return(invisible(NULL))
    }
    
    # 读全部行
    log_msg("  读文件...")
    con <- gzfile(sm_file, "r")
    all_lines <- readLines(con, warn = FALSE)
    close(con)
    log_msg(sprintf("  总行数: %d", length(all_lines)))
    
    # 找表达数据起止位置
    begin_idx <- grep("^!series_matrix_table_begin", all_lines)
    end_idx   <- grep("^!series_matrix_table_end",   all_lines)
    
    if (length(begin_idx) == 0) {
        log_msg("  ⚠️ 找不到 !series_matrix_table_begin")
        return(invisible(NULL))
    }
    
    # ============================================================
    # 1. 解析元信息部分 (`!` 开头的行, 在表达数据之前)
    # ============================================================
    meta_lines <- all_lines[1:(begin_idx - 1)]
    
    # 只保留 !Sample_ 开头的行 (这些每个样本一列)
    sample_lines <- meta_lines[grepl("^!Sample_", meta_lines)]
    log_msg(sprintf("  !Sample_ 行数: %d", length(sample_lines)))
    
    # 解析每行: 第 1 列是字段名 (如 !Sample_geo_accession), 后续每列一个样本
    # 数据是用 tab 分隔, 带引号
    parse_line <- function(line) {
        # 去除引号, 按 tab 拆分
        parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
        # 去掉两端引号
        parts <- gsub('^"|"$', '', parts)
        parts
    }
    
    sample_data <- lapply(sample_lines, parse_line)
    
    # 第一个字段是字段名 (e.g. "!Sample_geo_accession")
    field_names <- sapply(sample_data, function(x) x[1])
    # 去掉 "!Sample_" 前缀作为列名
    field_names_clean <- sub("^!Sample_", "", field_names)
    
    # 检查每行的样本数是否一致
    n_samples_each <- sapply(sample_data, function(x) length(x) - 1)
    n_samples <- median(n_samples_each)
    log_msg(sprintf("  样本数: %d", n_samples))
    
    # 同名字段可能出现多次 (e.g. 多个 characteristics_ch1)
    # 把同名字段加索引: characteristics_ch1_1, characteristics_ch1_2, ...
    dedup_names <- character(length(field_names_clean))
    counts <- list()
    for (i in seq_along(field_names_clean)) {
        n <- field_names_clean[i]
        if (is.null(counts[[n]])) counts[[n]] <- 0
        counts[[n]] <- counts[[n]] + 1
        dedup_names[i] <- if (counts[[n]] == 1) n else sprintf("%s_%d", n, counts[[n]])
    }
    
    # 构建临床 data.table: 行=样本, 列=字段
    clin_mat <- matrix("", nrow = n_samples, ncol = length(dedup_names))
    colnames(clin_mat) <- dedup_names
    for (i in seq_along(sample_data)) {
        values <- sample_data[[i]][-1]  # 去掉字段名
        # 补齐到 n_samples 长度
        if (length(values) < n_samples) values <- c(values, rep("", n_samples - length(values)))
        if (length(values) > n_samples) values <- values[1:n_samples]
        clin_mat[, i] <- values
    }
    
    clin_dt <- as.data.table(clin_mat)
    
    # GSM ID 作为第 1 列
    if ("geo_accession" %in% colnames(clin_dt)) {
        setcolorder(clin_dt, c("geo_accession", setdiff(colnames(clin_dt), "geo_accession")))
    }
    
    log_msg(sprintf("  临床表: %d 样本 × %d 字段", nrow(clin_dt), ncol(clin_dt)))
    
    # ============================================================
    # 2. 解析表达数据部分
    # ============================================================
    if (length(end_idx) == 0) end_idx <- length(all_lines) + 1
    expr_lines <- all_lines[(begin_idx + 1):(end_idx - 1)]
    
    # 跳过空行
    expr_lines <- expr_lines[nchar(expr_lines) > 0]
    log_msg(sprintf("  表达行数 (含表头): %d", length(expr_lines)))
    
    if (length(expr_lines) < 2) {
        log_msg("  ⚠️ 表达数据为空")
        # 仍输出临床
        clin_file <- file.path(gse_dir, sprintf("%s.clinical.tsv", gse_id))
        fwrite(clin_dt, clin_file, sep = "\t")
        log_msg(sprintf("  ✓ %s", basename(clin_file)))
        return(invisible(list(clinical = clin_dt, expr = NULL)))
    }
    
    # 把表达数据写到临时文件, 用 fread 读
    tmp_file <- tempfile(fileext = ".tsv")
    writeLines(expr_lines, tmp_file)
    
    expr_dt <- fread(tmp_file, sep = "\t", header = TRUE, quote = '"', showProgress = FALSE)
    file.remove(tmp_file)
    
    # 第 1 列是 probe ID
    setnames(expr_dt, 1, "probe_id")
    log_msg(sprintf("  表达矩阵: %d probes × %d 样本", nrow(expr_dt), ncol(expr_dt) - 1))
    
    # ============================================================
    # 3. 写出
    # ============================================================
    expr_file <- file.path(gse_dir, sprintf("%s.expr.tsv.gz", gse_id))
    clin_file <- file.path(gse_dir, sprintf("%s.clinical.tsv", gse_id))
    
    fwrite(expr_dt, expr_file, sep = "\t", compress = "gzip")
    fwrite(clin_dt, clin_file, sep = "\t")
    
    log_msg(sprintf("  ✓ %s (%.1f MB)", basename(expr_file), file.size(expr_file) / 1024^2))
    log_msg(sprintf("  ✓ %s (%.1f KB)", basename(clin_file), file.size(clin_file) / 1024))
    
    # 预览临床字段
    log_msg("  临床字段名:")
    for (cn in colnames(clin_dt)) cat(sprintf("    %s\n", cn))
    
    # 预览有几个 characteristics 字段，每个的非空值前 5 个
    char_cols <- grep("characteristics", colnames(clin_dt), value = TRUE)
    if (length(char_cols) > 0) {
        log_msg("  characteristics 字段示例值:")
        for (cn in char_cols) {
            vals <- unique(clin_dt[[cn]])
            vals <- vals[vals != "" & !is.na(vals)]
            cat(sprintf("    %s: %s\n", cn, paste(head(vals, 3), collapse = " | ")))
        }
    }
    
    invisible(list(clinical = clin_dt, expr = expr_dt))
}

# ============================================================
# 主流程: 5 个队列
# ============================================================
gse_ids <- c("GSE71729", "GSE21501", "GSE85916", "GSE62452", "GSE57495")

results <- list()
for (gse in gse_ids) {
    results[[gse]] <- parse_series_matrix(gse, file.path(BULK, gse))
    cat("\n")
}

log_msg("=== 总结 ===")
for (gse in gse_ids) {
    r <- results[[gse]]
    if (is.null(r)) {
        cat(sprintf("  %s: 失败\n", gse))
    } else {
        cat(sprintf("  %s: %d 样本 × %d 临床字段, %s\n",
                   gse,
                   nrow(r$clinical),
                   ncol(r$clinical),
                   if (is.null(r$expr)) "无表达" else sprintf("%d probes", nrow(r$expr))))
    }
}
