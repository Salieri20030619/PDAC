#!/usr/bin/env Rscript
# ============================================================
# PACA-AU 临床数据标准化
#
# 输入: donor / specimen / sample 三表 (无 header, ICGC ARGO 标准 schema)
# 表达矩阵: exp_array.PACA-AU.gene.tsv.gz (20744 genes × 269 SA*)
#
# 流程:
#   1. 给三表手动加 ICGC schema 列名
#   2. join: sample → specimen → donor
#   3. 标准化字段: sample_id / os_time_month / os_status / sample_type / age / sex / stage
#   4. 只保留表达矩阵里的 269 个 SA*
#
# 输出:
#   clinical_std.PACA-AU.tsv
#   clinical_full.PACA-AU.tsv  (含 join 后所有原始字段, 备查)
# ============================================================

suppressPackageStartupMessages({
    library(data.table)
})

ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/ICGC_PACA_AU"
log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# 1. ICGC ARGO Release 28 标准 schema
# ============================================================
# 注: 这是 ICGC 25K 项目 (含 PACA-AU) 的官方 schema
# 文档: https://docs.icgc-argo.org/docs/data-access/icgc-25k-data

donor_cols <- c(
    "icgc_donor_id",                                  # V1: DO*
    "project_code",                                   # V2: PACA-AU
    "study_donor_involved_in",                        # V3
    "submitted_donor_id",                             # V4: ICGC_*
    "donor_sex",                                      # V5: male/female
    "donor_vital_status",                             # V6: alive/deceased ← OS_status
    "disease_status_last_followup",                   # V7
    "donor_relapse_type",                             # V8
    "donor_age_at_diagnosis",                         # V9: 年龄
    "donor_age_at_enrollment",                        # V10
    "donor_age_at_last_followup",                     # V11
    "donor_relapse_interval",                         # V12
    "donor_diagnosis_icd10",                          # V13
    "donor_tumour_staging_system_at_diagnosis",       # V14
    "donor_tumour_stage_at_diagnosis",                # V15: stage
    "donor_tumour_stage_at_diagnosis_supplemental",   # V16
    "donor_survival_time",                            # V17: 天数 ← OS_time
    "donor_interval_of_last_followup",                # V18
    "prior_malignancy",                               # V19
    "cancer_type_prior_malignancy",                   # V20
    "cancer_history_first_degree_relative"            # V21
)

specimen_cols <- c(
    "icgc_specimen_id",                               # V1: SP*
    "project_code",                                   # V2
    "study_specimen_involved_in",                     # V3
    "submitted_specimen_id",                          # V4
    "icgc_donor_id",                                  # V5: DO*
    "submitted_donor_id",                             # V6
    "specimen_type",                                  # V7: Primary/Met/CellLine/Xenograft 等
    "specimen_type_other",                            # V8
    "specimen_interval",                              # V9
    "specimen_donor_treatment_type",                  # V10
    "specimen_donor_treatment_type_other",            # V11
    "specimen_processing",                            # V12
    "specimen_processing_other",                      # V13
    "specimen_storage",                               # V14
    "specimen_storage_other",                         # V15
    "tumour_confirmed",                               # V16
    "specimen_biobank",                               # V17
    "specimen_biobank_id",                            # V18
    "specimen_available",                             # V19
    "tumour_histological_type",                       # V20
    "tumour_grading_system",                          # V21
    "tumour_grade",                                   # V22
    "tumour_grade_supplemental",                      # V23
    "tumour_stage_system",                            # V24
    "tumour_stage",                                   # V25: 分期
    "tumour_stage_supplemental",                      # V26
    "digital_image_of_stained_section",               # V27
    "percentage_cellularity",                         # V28
    "level_of_cellularity"                            # V29
)

sample_cols <- c(
    "icgc_sample_id",                                 # V1: SA*
    "project_code",                                   # V2
    "submitted_sample_id",                            # V3
    "icgc_specimen_id",                               # V4: SP*
    "submitted_specimen_id",                          # V5
    "icgc_donor_id",                                  # V6: DO*
    "submitted_donor_id",                             # V7
    "analyzed_sample_interval",                       # V8
    "percentage_cellularity",                         # V9
    "level_of_cellularity",                           # V10
    "study"                                           # V11
)

# ============================================================
# 2. 读取 + 重命名三表
# ============================================================
log_msg("=== 1. 读取三表 + 加 schema 列名 ===")

donor <- fread(file.path(ROOT, "donor.PACA-AU.tsv.gz"),
               sep = "\t", header = FALSE, na.strings = c("", "NA"))
# 第 1 行是占位 V1-V21 (合并时遗留), 去掉
donor <- donor[!startsWith(V1, "V")]
setnames(donor, donor_cols)
log_msg(sprintf("  donor: %d × %d", nrow(donor), ncol(donor)))

specimen <- fread(file.path(ROOT, "specimen.PACA-AU.tsv.gz"),
                  sep = "\t", header = FALSE, na.strings = c("", "NA"), quote = "")
specimen <- specimen[!startsWith(V1, "V")]
setnames(specimen, specimen_cols)
log_msg(sprintf("  specimen: %d × %d", nrow(specimen), ncol(specimen)))

sample <- fread(file.path(ROOT, "sample.PACA-AU.tsv.gz"),
                sep = "\t", header = FALSE, na.strings = c("", "NA"))
sample <- sample[!startsWith(V1, "V")]
setnames(sample, sample_cols)
log_msg(sprintf("  sample: %d × %d", nrow(sample), ncol(sample)))

# ============================================================
# 3. 读表达矩阵, 提取 269 个 SA*
# ============================================================
log_msg("=== 2. 提取表达矩阵的 269 个样本 ===")
expr_header <- fread(file.path(ROOT, "exp_array.PACA-AU.gene.tsv.gz"),
                     nrows = 0)
expr_samples <- setdiff(colnames(expr_header), "gene")
log_msg(sprintf("  表达矩阵样本数 (SA*): %d", length(expr_samples)))
log_msg(sprintf("  示例: %s ...", paste(head(expr_samples, 3), collapse = ", ")))

# ============================================================
# 4. join: sample → specimen → donor
# ============================================================
log_msg("=== 3. 三表 join ===")

# sample 是细粒度 (1657 样本 / 461 donor), 以 sample 为中心 join
clin_full <- merge(sample, specimen, by = c("icgc_specimen_id", "icgc_donor_id", 
                                              "project_code"), 
                   all.x = TRUE, suffixes = c(".sample", ".specimen"))
log_msg(sprintf("  sample + specimen: %d × %d", nrow(clin_full), ncol(clin_full)))

clin_full <- merge(clin_full, donor, by = c("icgc_donor_id", "project_code"), 
                   all.x = TRUE)
log_msg(sprintf("  + donor: %d × %d", nrow(clin_full), ncol(clin_full)))

# ============================================================
# 5. 筛 269 个表达矩阵样本
# ============================================================
log_msg("=== 4. 筛 269 个表达矩阵样本 ===")
clin_269 <- clin_full[icgc_sample_id %in% expr_samples]
log_msg(sprintf("  匹配上: %d / %d", nrow(clin_269), length(expr_samples)))

missing <- setdiff(expr_samples, clin_269$icgc_sample_id)
if (length(missing) > 0) {
    log_msg(sprintf("  ⚠️ 未匹配上的 SA*: %d (%s)", 
                    length(missing), 
                    paste(head(missing, 3), collapse = ", ")))
}

# ============================================================
# 6. 标准化字段
# ============================================================
log_msg("=== 5. 字段标准化 ===")

std <- data.table(
    sample_id     = clin_269$icgc_sample_id,
    donor_id      = clin_269$icgc_donor_id,
    specimen_id   = clin_269$icgc_specimen_id,
    
    # 生存数据
    # ICGC donor_survival_time 单位是天, 转月 (× 1/30.44)
    os_time_month = as.numeric(clin_269$donor_survival_time) / 30.44,
    # alive=0, deceased=1
    os_status     = as.integer(ifelse(clin_269$donor_vital_status == "deceased", 1, 
                                ifelse(clin_269$donor_vital_status == "alive", 0, NA))),
    
    # specimen 类型 (curation 用)
    sample_type   = clin_269$specimen_type,
    
    # 人口学
    sex           = clin_269$donor_sex,
    age           = as.numeric(clin_269$donor_age_at_diagnosis),
    
    # 分期 (从 donor 表取, specimen 表也有但更不完整)
    stage         = clin_269$donor_tumour_stage_at_diagnosis,
    grade         = clin_269$tumour_grade,
    
    # 其他
    cellularity_pct = as.numeric(clin_269$percentage_cellularity.sample)
)

# 简化 sample_type 命名
std[, sample_type_simple := fcase(
    grepl("Primary tumour", sample_type, ignore.case = TRUE), "Primary",
    grepl("Metast", sample_type, ignore.case = TRUE), "Metastasis",
    grepl("Xenograft", sample_type, ignore.case = TRUE), "Xenograft",
    grepl("Cell line", sample_type, ignore.case = TRUE), "CellLine",
    grepl("Normal", sample_type, ignore.case = TRUE), "Normal",
    default = "Other"
)]

log_msg("  sample_type 原始分布:")
print(table(std$sample_type, useNA = "ifany"))

log_msg("  sample_type_simple 分布:")
print(table(std$sample_type_simple, useNA = "ifany"))

# ============================================================
# 7. 报告 + 输出
# ============================================================
log_msg("=== 6. 总样本统计 ===")
log_msg(sprintf("  总样本: %d", nrow(std)))
log_msg(sprintf("  含 OS_time: %d", sum(!is.na(std$os_time_month))))
log_msg(sprintf("  含 OS_status: %d", sum(!is.na(std$os_status))))
log_msg(sprintf("  同时有 time + status: %d", 
                sum(!is.na(std$os_time_month) & !is.na(std$os_status))))

# Primary + 完整生存 (训练集核心)
primary_complete <- std[sample_type_simple == "Primary" & 
                         !is.na(os_time_month) & !is.na(os_status)]
log_msg(sprintf("  Primary + 完整生存 (训练集核心): %d", nrow(primary_complete)))

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
}

# 输出
out_std <- file.path(ROOT, "clinical_std.PACA-AU.tsv")
out_full <- file.path(ROOT, "clinical_full.PACA-AU.tsv")
fwrite(std, out_std, sep = "\t")
fwrite(clin_269, out_full, sep = "\t")
log_msg(sprintf("=== 输出 ==="))
log_msg(sprintf("  ✓ %s", basename(out_std)))
log_msg(sprintf("  ✓ %s", basename(out_full)))
