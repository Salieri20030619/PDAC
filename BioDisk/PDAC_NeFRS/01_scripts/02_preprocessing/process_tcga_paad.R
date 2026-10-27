#!/usr/bin/env Rscript
# ============================================================
# TCGA-PAAD 整理
#
# 输入:
#   TCGA-PAAD.star_tpm.tsv.gz       (TPM 表达, 60660 ENSG × 183 samples)
#   TCGA-PAAD.clinical.tsv.gz       (89 字段临床)
#   TCGA-PAAD.survival.tsv.gz       (Xena 已处理: sample, OS.time, OS, _PATIENT)
#
# 输出:
#   TCGA-PAAD.gene.tsv.gz           (gene-level TPM)
#   TCGA-PAAD.ensembl2symbol_map.tsv  映射记录
#   TCGA-PAAD.clinical_std.tsv      标准化临床
#   TCGA-PAAD.clinical_full.tsv     合并后完整临床备查
# ============================================================

suppressPackageStartupMessages({
    library(data.table)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
})

ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/TCGA_PAAD"
log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# 1. 读 TPM 表达矩阵
# ============================================================
log_msg("=== 1. 读 TPM 表达矩阵 ===")
expr <- fread(file.path(ROOT, "TCGA-PAAD.star_tpm.tsv.gz"),
              header = TRUE, showProgress = FALSE)
setnames(expr, 1, "ensembl_versioned")
log_msg(sprintf("  矩阵: %d ENSG × %d 样本", nrow(expr), ncol(expr) - 1))

# 去除版本号 (ENSG00000000003.15 → ENSG00000000003)
expr[, ensembl_id := sub("\\..*$", "", ensembl_versioned)]
log_msg(sprintf("  唯一 ENSG (去版本): %d", uniqueN(expr$ensembl_id)))

# ============================================================
# 2. ENSG → gene symbol 映射
# ============================================================
log_msg("=== 2. ENSG → gene symbol 映射 ===")

map <- AnnotationDbi::select(org.Hs.eg.db,
                              keys = expr$ensembl_id,
                              columns = c("SYMBOL", "ENTREZID"),
                              keytype = "ENSEMBL")
setDT(map)
log_msg(sprintf("  AnnotationDbi 返回: %d 行", nrow(map)))
log_msg(sprintf("  有 SYMBOL: %d (%.1f%%)",
                sum(!is.na(map$SYMBOL)),
                100 * mean(!is.na(map$SYMBOL))))

# 处理 1:many (一个 ENSG → 多个 SYMBOL)
# 策略: 取第一个 SYMBOL (主导 transcript)
map_unique <- map[!is.na(SYMBOL), .(SYMBOL = SYMBOL[1]), by = ENSEMBL]
setnames(map_unique, c("ensembl_id", "gene_symbol"))
log_msg(sprintf("  去重后 ENSG → SYMBOL 映射: %d", nrow(map_unique)))

# ============================================================
# 3. 应用映射 + 多 ENSG 聚合 (最大方差)
# ============================================================
log_msg("=== 3. 应用映射 + 聚合 (最大方差 ENSG) ===")

expr_mapped <- merge(expr, map_unique, by = "ensembl_id", all.x = FALSE)
log_msg(sprintf("  匹配上: %d / %d (%.1f%%)",
                nrow(expr_mapped), nrow(expr),
                100 * nrow(expr_mapped) / nrow(expr)))
log_msg(sprintf("  映射到 %d 个唯一 gene", uniqueN(expr_mapped$gene_symbol)))

# 样本列
sample_cols <- setdiff(colnames(expr_mapped), c("ensembl_versioned", "ensembl_id", "gene_symbol"))
log_msg(sprintf("  样本列数: %d", length(sample_cols)))

# 每基因 ENSG 数分布
ensg_per_gene <- expr_mapped[, .N, by = gene_symbol]
log_msg("  每基因 ENSG 数分布:")
log_msg(sprintf("    1 ENSG: %d 基因", sum(ensg_per_gene$N == 1)))
log_msg(sprintf("    2 ENSG: %d 基因", sum(ensg_per_gene$N == 2)))
log_msg(sprintf("    3+ ENSG: %d 基因", sum(ensg_per_gene$N >= 3)))

# 多 ENSG 聚合: 最大方差
log_msg("  计算方差, 选最大方差 ENSG...")
expr_mapped[, ensg_var := apply(.SD, 1, var, na.rm = TRUE), .SDcols = sample_cols]
selected <- expr_mapped[, .SD[which.max(ensg_var)], by = gene_symbol]
log_msg(sprintf("  聚合后唯一基因数: %d", nrow(selected)))

# 输出 gene-level matrix
gene_mat <- selected[, c("gene_symbol", sample_cols), with = FALSE]
setnames(gene_mat, "gene_symbol", "gene")

out_gene <- file.path(ROOT, "TCGA-PAAD.gene.tsv.gz")
fwrite(gene_mat, out_gene, sep = "\t", compress = "gzip")
log_msg(sprintf("  ✓ %s (%.1f MB)", basename(out_gene), file.size(out_gene) / 1024^2))

# 输出映射记录
map_record <- selected[, .(gene_symbol, selected_ensembl = ensembl_id, ensg_var)]
map_record <- merge(map_record, ensg_per_gene, by = "gene_symbol")
setnames(map_record, "N", "n_ensg_available")
out_map <- file.path(ROOT, "TCGA-PAAD.ensembl2symbol_map.tsv")
fwrite(map_record, out_map, sep = "\t")
log_msg(sprintf("  ✓ %s", basename(out_map)))

# ============================================================
# 4. 临床数据合并 + 标准化
# ============================================================
log_msg("=== 4. 临床数据合并 + 标准化 ===")

clin <- fread(file.path(ROOT, "TCGA-PAAD.clinical.tsv.gz"),
              header = TRUE, na.strings = c("", "NA", "not reported", "Not Reported", "Unknown"))
log_msg(sprintf("  临床表: %d 样本 × %d 字段", nrow(clin), ncol(clin)))

surv <- fread(file.path(ROOT, "TCGA-PAAD.survival.tsv.gz"),
              header = TRUE, na.strings = c("", "NA"))
log_msg(sprintf("  生存表: %d 样本 × %d 字段", nrow(surv), ncol(surv)))

# join clinical + survival, 以 sample 为 key
clin_full <- merge(clin, surv, by = "sample", all.x = TRUE)
log_msg(sprintf("  合并后: %d × %d", nrow(clin_full), ncol(clin_full)))

# 输出完整临床
out_full <- file.path(ROOT, "TCGA-PAAD.clinical_full.tsv")
fwrite(clin_full, out_full, sep = "\t")
log_msg(sprintf("  ✓ %s", basename(out_full)))

# 标准化字段 (与其他队列统一)
std <- data.table(
    sample_id     = clin_full$sample,
    patient_id    = clin_full$`_PATIENT`,
    
    # 生存
    os_time_day   = as.numeric(clin_full$OS.time),
    os_time_month = as.numeric(clin_full$OS.time) / 30.44,
    os_status     = as.integer(clin_full$OS),
    
    # sample type (从 TCGA barcode 第 14-15 位)
    # 01A/01B = Primary Tumor, 11A/11B = Solid Tissue Normal
    sample_type_code = substr(clin_full$sample, 14, 15),
    sample_type      = clin_full$sample_type.samples,
    
    # 人口学
    sex           = clin_full$gender.demographic,
    age           = as.numeric(clin_full$age_at_index.demographic),
    race          = clin_full$race.demographic,
    
    # 临床病理
    stage         = clin_full$ajcc_pathologic_stage.diagnoses,
    grade         = clin_full$tumor_grade.diagnoses,
    t_stage       = clin_full$ajcc_pathologic_t.diagnoses,
    n_stage       = clin_full$ajcc_pathologic_n.diagnoses,
    m_stage       = clin_full$ajcc_pathologic_m.diagnoses,
    
    # 诊断信息 (curation 用)
    primary_site  = clin_full$primary_site,
    primary_diagnosis = clin_full$primary_diagnosis.diagnoses,
    
    # 注释 (curation 用 — 这里包含 "should not be included" 这种关键信息)
    annotation_classification = clin_full$classification.annotations,
    annotation_category = clin_full$category.annotations,
    annotation_notes    = clin_full$notes.annotations,
    
    # 治疗
    prior_treatment = clin_full$prior_treatment.diagnoses
)

# sample_type_simple 简化
std[, sample_type_simple := fcase(
    sample_type_code %in% c("01", "01A", "01B", "01C"), "Primary",
    sample_type_code %in% c("11", "11A", "11B"), "Normal",
    sample_type_code %in% c("06", "06A"), "Metastasis",
    default = "Other"
)]

# 报告
log_msg("=== 5. 总样本统计 ===")
log_msg(sprintf("  总样本: %d", nrow(std)))
log_msg("  sample_type_simple 分布:")
print(table(std$sample_type_simple, useNA = "ifany"))

log_msg(sprintf("  含 OS_time: %d", sum(!is.na(std$os_time_month))))
log_msg(sprintf("  含 OS_status: %d", sum(!is.na(std$os_status))))
log_msg(sprintf("  同时有 time + status: %d",
                sum(!is.na(std$os_time_month) & !is.na(std$os_status))))

# Primary + 完整生存
primary_complete <- std[sample_type_simple == "Primary" &
                          !is.na(os_time_month) & !is.na(os_status)]
log_msg(sprintf("  Primary + 完整生存: %d", nrow(primary_complete)))

if (nrow(primary_complete) > 0) {
    log_msg(sprintf("    OS_time 范围 (月): %.1f - %.1f, 中位 %.1f",
                    min(primary_complete$os_time_month, na.rm = TRUE),
                    max(primary_complete$os_time_month, na.rm = TRUE),
                    median(primary_complete$os_time_month, na.rm = TRUE)))
    log_msg(sprintf("    Death events: %d (%.1f%%)",
                    sum(primary_complete$os_status == 1),
                    100 * mean(primary_complete$os_status == 1)))
    log_msg(sprintf("    age 范围: %.0f - %.0f, 中位 %.0f",
                    min(primary_complete$age, na.rm = TRUE),
                    max(primary_complete$age, na.rm = TRUE),
                    median(primary_complete$age, na.rm = TRUE)))
    log_msg("    stage 分布:")
    print(table(primary_complete$stage, useNA = "ifany"))
    log_msg("    primary_diagnosis 分布 (curation 关键信息):")
    print(sort(table(primary_complete$primary_diagnosis, useNA = "ifany"), decreasing = TRUE))
}

# 已知混杂样本预览 (annotation 字段非 NA 的)
log_msg("=== 6. 已知混杂样本预览 (来自 TCGA annotation 字段) ===")
flagged <- std[!is.na(annotation_classification) | !is.na(annotation_category)]
log_msg(sprintf("  含 annotation flag 的样本: %d", nrow(flagged)))
if (nrow(flagged) > 0) {
    log_msg("  flag 类别分布:")
    print(table(flagged$annotation_classification, useNA = "ifany"))
    log_msg("  category 分布:")
    print(table(flagged$annotation_category, useNA = "ifany"))
}

# 输出标准化临床
out_std <- file.path(ROOT, "TCGA-PAAD.clinical_std.tsv")
fwrite(std, out_std, sep = "\t")
log_msg(sprintf("=== 7. 输出 ==="))
log_msg(sprintf("  ✓ %s", basename(out_std)))

# 关键基因健康检查
log_msg("=== 8. 关键基因健康检查 ===")
check_genes <- c("KRAS", "TP53", "SMAD4", "CDKN2A", "GAPDH", "ACTB",
                 "GLIS3", "NEUROD1", "INSM1", "NFIB", "SOX4",
                 "GPX4", "SLC7A11", "ACSL4", "FTH1")
present <- check_genes %in% gene_mat$gene
for (g in check_genes) {
    cat(sprintf("    %-10s %s\n", g, if (g %in% gene_mat$gene) "✓" else "✗"))
}
log_msg(sprintf("  共 %d / %d 在矩阵中", sum(present), length(check_genes)))
