#!/usr/bin/env Rscript
# ============================================================
# PACA-AU 数据整合 v2
# 用途: 合并 exp_array 为基因 × 样本矩阵, 合并 donor 等临床表
# 
# 输入: raw_per_patient/DOxxxxx/{donor,specimen,sample,exp_array,ssm_open}/part-*.gz
# 输出:
#   exp_array.PACA-AU.matrix.tsv.gz    基因 × 样本表达矩阵
#   exp_array.PACA-AU.long.tsv.gz      长格式备份（含完整元信息）
#   donor.PACA-AU.tsv.gz               临床
#   specimen.PACA-AU.tsv.gz            标本
#   sample.PACA-AU.tsv.gz              样本-平台
#   ssm_open.PACA-AU.tsv.gz            突变
# ============================================================

suppressPackageStartupMessages({
    library(data.table)
})

ROOT    <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/ICGC_PACA_AU"
RAW_DIR <- file.path(ROOT, "raw_per_patient")
OUT_DIR <- ROOT

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# 通用合并函数
# ============================================================
merge_data_type <- function(data_type, out_name) {
    log_msg(sprintf("=== %s ===", data_type))
    
    patient_dirs <- list.dirs(RAW_DIR, recursive = FALSE)
    all_files <- character(0)
    for (pdir in patient_dirs) {
        sub <- file.path(pdir, data_type)
        if (dir.exists(sub)) {
            all_files <- c(all_files, list.files(sub, "^part-.*\\.gz$", full.names = TRUE))
        }
    }
    log_msg(sprintf("  %d 个 part 文件", length(all_files)))
    if (length(all_files) == 0) return(NULL)
    
    # data.table fread 不带表头时返回 V1, V2, ...
    dt_list <- lapply(seq_along(all_files), function(i) {
        if (i %% 1000 == 0) log_msg(sprintf("    已读 %d / %d", i, length(all_files)))
        tryCatch(fread(all_files[i], header = FALSE, showProgress = FALSE),
                 error = function(e) NULL)
    })
    dt_list <- dt_list[!sapply(dt_list, is.null)]
    
    merged <- rbindlist(dt_list, fill = TRUE)
    log_msg(sprintf("  合并: %d 行 × %d 列", nrow(merged), ncol(merged)))
    
    out_file <- file.path(OUT_DIR, out_name)
    fwrite(merged, out_file, sep = "\t", compress = "gzip")
    log_msg(sprintf("  写出: %s (%.1f MB)", basename(out_file), file.size(out_file) / 1024^2))
    
    invisible(merged)
}

# ============================================================
# 1. donor / specimen / sample / ssm_open (无表头, 用 V1-Vn)
# ============================================================
donor    <- merge_data_type("donor",    "donor.PACA-AU.tsv.gz")
specimen <- merge_data_type("specimen", "specimen.PACA-AU.tsv.gz")
sample   <- merge_data_type("sample",   "sample.PACA-AU.tsv.gz")
ssm      <- merge_data_type("ssm_open", "ssm_open.PACA-AU.tsv.gz")

# ============================================================
# 2. exp_array (无表头) - 先存长格式
# ============================================================
exp_long <- merge_data_type("exp_array", "exp_array.PACA-AU.long.tsv.gz")

# ============================================================
# 3. exp_array 长→宽: 基因 × 样本矩阵
# 列名约定（根据原始数据看到的结构）:
#   V1=icgc_donor_id  V2=project  V3=icgc_specimen_id  V4=icgc_sample_id
#   V5=analyzed_sample_id  V6=analysis_id  V7=gene_model
#   V8=gene_id (ILMN_*)  V9=normalized_expression_value
# ============================================================
if (!is.null(exp_long)) {
    log_msg("=== exp_array 长 → 宽矩阵 ===")
    
    setnames(exp_long,
             old = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"),
             new = c("icgc_donor_id", "project_code", "icgc_specimen_id",
                     "icgc_sample_id", "analyzed_sample_id", "analysis_id",
                     "gene_model", "gene_id", "expr"))
    
    # 强制数值
    exp_long[, expr := as.numeric(expr)]
    
    log_msg(sprintf("  唯一样本数 (icgc_sample_id): %d", uniqueN(exp_long$icgc_sample_id)))
    log_msg(sprintf("  唯一基因数 (gene_id): %d", uniqueN(exp_long$gene_id)))
    log_msg(sprintf("  非空表达值数: %d / %d", sum(!is.na(exp_long$expr)), nrow(exp_long)))
    
    # 用 icgc_sample_id 做样本 ID（同一患者可能多个样本）
    # 如果同一基因×样本有多个值（多探针），取平均
    log_msg("  dcast 长→宽...")
    exp_mat <- dcast(
        exp_long,
        gene_id ~ icgc_sample_id,
        value.var = "expr",
        fun.aggregate = function(x) mean(x, na.rm = TRUE),
        fill = NA
    )
    
    log_msg(sprintf("  矩阵维度: %d 基因 × %d 样本", nrow(exp_mat), ncol(exp_mat) - 1))
    
    out_file <- file.path(OUT_DIR, "exp_array.PACA-AU.matrix.tsv.gz")
    fwrite(exp_mat, out_file, sep = "\t", compress = "gzip")
    log_msg(sprintf("  写出: %s (%.1f MB)", basename(out_file), file.size(out_file) / 1024^2))
    
    # 简要预览
    log_msg("  矩阵前 3 行 × 前 4 列预览:")
    print(exp_mat[1:3, 1:4])
}

log_msg("=== 完成 ===")
log_msg("输出文件:")
for (f in list.files(OUT_DIR, pattern = "PACA-AU\\.(matrix|tsv|long)\\.tsv\\.gz$", full.names = TRUE)) {
    cat(sprintf("  %s  (%.1f MB)\n", basename(f), file.size(f) / 1024^2))
}
