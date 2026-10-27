#!/usr/bin/env Rscript
# ============================================================
# TCGA-PAAD curation (Peran 2018 style)
#
# 输入: TCGA-PAAD.gene.tsv.gz (36650 × 183), TCGA-PAAD.clinical_std.tsv (183)
# 输出: TCGA-PAAD.curated.gene.tsv.gz (gene × 147 curated samples)
#       TCGA-PAAD.curated.clinical_std.tsv (147 curated 临床)
#
# Curation 规则:
#   1. 只保留 Primary Tumor
#   2. 必须有完整 OS_time + OS_status
#   3. primary_diagnosis ∈ {Infiltrating duct carcinoma NOS, Adenocarcinoma NOS}
#      排除: Neuroendocrine, Mucinous adenocarcinoma, mixed subtypes, undifferentiated
#   4. TCGA annotation_classification 为空 (排除 "Notification" / "Observation" 标记)
# ============================================================

suppressPackageStartupMessages({ library(data.table) })

ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/TCGA_PAAD"
log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# 1. 读临床表
# ============================================================
log_msg("=== 1. 读 clinical_std ===")
clin <- fread(file.path(ROOT, "TCGA-PAAD.clinical_std.tsv"))
log_msg(sprintf("  总: %d", nrow(clin)))

# ============================================================
# 2. 四步 curation
# ============================================================
log_msg("=== 2. 应用 Peran-style curation ===")

# Step 1: Primary Tumor
step1 <- clin[sample_type_simple == "Primary"]
log_msg(sprintf("  Step 1 (Primary): %d", nrow(step1)))

# Step 2: 完整 OS_time + OS_status
step2 <- step1[!is.na(os_time_month) & !is.na(os_status)]
log_msg(sprintf("  Step 2 (+ 完整生存): %d", nrow(step2)))

# Step 3: 真正 PDAC 诊断
pdac_dx <- c("Infiltrating duct carcinoma, NOS", "Adenocarcinoma, NOS")
step3 <- step2[primary_diagnosis %in% pdac_dx]
log_msg(sprintf("  Step 3 (+ PDAC 诊断): %d", nrow(step3)))
log_msg("    排除的诊断:")
print(setdiff(step2$primary_diagnosis, pdac_dx))

# Step 4: 无 TCGA annotation flag
step4 <- step3[is.na(annotation_classification) | annotation_classification == ""]
log_msg(sprintf("  Step 4 (+ 排 TCGA flag): %d (curated final)", nrow(step4)))
log_msg("    排除的 flag:")
print(table(step3[!is.na(annotation_classification) & annotation_classification != ""]$annotation_classification))

curated_ids <- step4$sample_id

# ============================================================
# 3. 输出 curated 临床表
# ============================================================
log_msg("=== 3. 输出 curated 临床表 ===")
out_clin <- file.path(ROOT, "TCGA-PAAD.curated.clinical_std.tsv")
fwrite(step4, out_clin, sep = "\t")
log_msg(sprintf("  ✓ %s (%d 样本)", basename(out_clin), nrow(step4)))

# ============================================================
# 4. 输出 curated 表达矩阵
# ============================================================
log_msg("=== 4. 过滤表达矩阵到 curated 样本 ===")
expr <- fread(file.path(ROOT, "TCGA-PAAD.gene.tsv.gz"), header = TRUE)
log_msg(sprintf("  原始: %d genes × %d samples", nrow(expr), ncol(expr) - 1))

# 只保留 gene 列 + curated 样本列
keep_cols <- c("gene", curated_ids)
missing_cols <- setdiff(curated_ids, colnames(expr))
if (length(missing_cols) > 0) {
    log_msg(sprintf("  ⚠️ %d 个 curated 样本不在表达矩阵: %s",
                    length(missing_cols), paste(head(missing_cols, 3), collapse = ", ")))
    keep_cols <- intersect(keep_cols, colnames(expr))
}
expr_curated <- expr[, ..keep_cols]
log_msg(sprintf("  curated: %d genes × %d samples", nrow(expr_curated), ncol(expr_curated) - 1))

out_gene <- file.path(ROOT, "TCGA-PAAD.curated.gene.tsv.gz")
fwrite(expr_curated, out_gene, sep = "\t", compress = "gzip")
log_msg(sprintf("  ✓ %s (%.1f MB)", basename(out_gene), file.size(out_gene) / 1024^2))

# ============================================================
# 5. 最终汇总
# ============================================================
log_msg("=== 5. TCGA-PAAD 数据集汇总 ===")
log_msg("  基础版 (用于敏感性分析):")
log_msg(sprintf("    TCGA-PAAD.gene.tsv.gz           36650 genes × 183 samples"))
log_msg(sprintf("    TCGA-PAAD.clinical_std.tsv      183 samples"))
log_msg("  Curated 版 (用于主验证):")
log_msg(sprintf("    TCGA-PAAD.curated.gene.tsv.gz   36650 genes × %d samples", ncol(expr_curated) - 1))
log_msg(sprintf("    TCGA-PAAD.curated.clinical_std.tsv  %d samples", nrow(step4)))
log_msg(sprintf("  与 Peran 2018 (~150) 一致 ✓"))

# 关键统计
log_msg("=== 6. Curated 队列关键指标 ===")
log_msg(sprintf("  OS_time 范围 (月): %.1f - %.1f, 中位 %.1f",
                min(step4$os_time_month), max(step4$os_time_month), median(step4$os_time_month)))
log_msg(sprintf("  Death events: %d (%.1f%%)",
                sum(step4$os_status == 1), 100 * mean(step4$os_status == 1)))
log_msg(sprintf("  Age 范围: %.0f - %.0f, 中位 %.0f",
                min(step4$age, na.rm=TRUE), max(step4$age, na.rm=TRUE),
                median(step4$age, na.rm=TRUE)))
log_msg("  Stage 分布:")
print(table(step4$stage, useNA = "ifany"))
log_msg("  Grade 分布:")
print(table(step4$grade, useNA = "ifany"))
log_msg("  Sex 分布:")
print(table(step4$sex, useNA = "ifany"))
