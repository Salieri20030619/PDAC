#!/usr/bin/env Rscript
# ============================================================
# TCGA-PAAD 临床表修复 + 精确样本量统计
#
# 问题: 原始 clinical 表 196 样本, 但表达只有 183
#       多出 13 个无表达数据的样本 (8 个 FZ 系列 patient: 7 个 01A + 6 个 11A)
#
# 修复: 只保留表达矩阵里的 183 个样本, 重做统计
# ============================================================

suppressPackageStartupMessages({ library(data.table) })

ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/TCGA_PAAD"
log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# 1. 读表达矩阵, 拿到 183 个样本 ID
# ============================================================
log_msg("=== 1. 读表达矩阵 ===")
expr_header <- fread(file.path(ROOT, "TCGA-PAAD.gene.tsv.gz"), nrows = 0)
expr_samples <- setdiff(colnames(expr_header), "gene")
log_msg(sprintf("  表达样本: %d", length(expr_samples)))

# ============================================================
# 2. 读原始 standardized clinical
# ============================================================
log_msg("=== 2. 读 standardized clinical ===")
std <- fread(file.path(ROOT, "TCGA-PAAD.clinical_std.tsv"))
log_msg(sprintf("  原始: %d 样本", nrow(std)))

# ============================================================
# 3. 过滤: 只保留有表达数据的样本
# ============================================================
log_msg("=== 3. 过滤到 183 个有表达数据的样本 ===")
std_with_expr <- std[sample_id %in% expr_samples]
log_msg(sprintf("  过滤后: %d 样本", nrow(std_with_expr)))

# 验证
missing_in_clin <- setdiff(expr_samples, std$sample_id)
log_msg(sprintf("  表达中但临床中无的样本: %d", length(missing_in_clin)))
if (length(missing_in_clin) > 0) {
    log_msg(sprintf("    %s", paste(missing_in_clin, collapse = ", ")))
}

# ============================================================
# 4. 精确统计
# ============================================================
log_msg("=== 4. 精确样本量统计 (基于 183 个表达样本) ===")
log_msg(sprintf("  总样本 (有表达): %d", nrow(std_with_expr)))

log_msg("  sample_type_simple 分布:")
print(table(std_with_expr$sample_type_simple, useNA = "ifany"))

log_msg(sprintf("  含 OS_time: %d", sum(!is.na(std_with_expr$os_time_month))))
log_msg(sprintf("  含 OS_status: %d", sum(!is.na(std_with_expr$os_status))))
log_msg(sprintf("  同时 time + status: %d",
                sum(!is.na(std_with_expr$os_time_month) & !is.na(std_with_expr$os_status))))

# Primary + 完整生存
primary_complete <- std_with_expr[sample_type_simple == "Primary" &
                                    !is.na(os_time_month) & !is.na(os_status)]
log_msg(sprintf("  Primary + 完整生存 (训练/验证可用): %d", nrow(primary_complete)))

# Primary + 完整生存 + 真正 PDAC 诊断 (curation 预览)
pdac_diagnoses <- c("Infiltrating duct carcinoma, NOS", "Adenocarcinoma, NOS")
primary_pdac <- primary_complete[primary_diagnosis %in% pdac_diagnoses]
log_msg(sprintf("  + 真正 PDAC 诊断 (排 NET/Mucinous/Undiff): %d", nrow(primary_pdac)))

# Primary + 完整生存 + PDAC 诊断 + 无 TCGA annotation flag
primary_pdac_clean <- primary_pdac[is.na(annotation_classification) | annotation_classification == ""]
log_msg(sprintf("  + 排 TCGA annotation flag (curated): %d", nrow(primary_pdac_clean)))

log_msg("  原始诊断分布 (curation 关键参考):")
print(sort(table(primary_complete$primary_diagnosis, useNA = "ifany"), decreasing = TRUE))

log_msg("  annotation_classification 分布:")
print(table(primary_complete$annotation_classification, useNA = "ifany"))

log_msg("  stage 分布 (filtered):")
print(table(primary_pdac_clean$stage, useNA = "ifany"))

# ============================================================
# 5. 输出更新后的 clinical_std (只保留 183 个有表达的样本)
# ============================================================
log_msg("=== 5. 输出更新文件 ===")

out_std <- file.path(ROOT, "TCGA-PAAD.clinical_std.tsv")
fwrite(std_with_expr, out_std, sep = "\t")
log_msg(sprintf("  ✓ %s (覆盖, 现 %d 样本)", basename(out_std), nrow(std_with_expr)))

# 额外: 记录被丢弃的 13 个样本
discarded <- std[!sample_id %in% expr_samples]
out_discarded <- file.path(ROOT, "TCGA-PAAD.clinical_discarded_no_expr.tsv")
fwrite(discarded, out_discarded, sep = "\t")
log_msg(sprintf("  ✓ %s (丢弃, 13 样本, 无 RNA-seq 数据)", basename(out_discarded)))

# 6. 关键样本量结论 (用于交接文档)
log_msg("=== 6. 关键样本量结论 (TCGA-PAAD) ===")
log_msg(sprintf("  原始 clinical 表: 196 样本"))
log_msg(sprintf("  原始 expression 矩阵: 183 样本"))
log_msg(sprintf("  交集 (实际可用): 183 样本"))
log_msg(sprintf("  其中 Primary + 完整生存: %d", nrow(primary_complete)))
log_msg(sprintf("  + 真 PDAC 诊断: %d", nrow(primary_pdac)))
log_msg(sprintf("  + 排 TCGA 异常 flag: %d (curated, 接近 Peran 2018 的 ~150)", nrow(primary_pdac_clean)))
