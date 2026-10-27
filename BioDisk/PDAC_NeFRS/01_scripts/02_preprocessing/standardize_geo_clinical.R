#!/usr/bin/env Rscript
# ============================================================
# 5 个 GEO 队列临床数据标准化
#
# 输出每个队列一个统一格式的临床表:
#   {GSE}.clinical_std.tsv
# 标准字段:
#   sample_id     - GSM ID
#   sample_type   - Primary/Met/Normal/CellLine
#   os_time_month - 生存时间 (月, GSE85916 已 ×12)
#   os_status     - 1=死亡, 0=删失
#   stage         - AJCC 分期 (有就有，没就 NA)
#   grade         - 病理分级 (G1-G4)
#   t_stage / n_stage - T/N 分期 (有就有)
#   subtype       - 论文已分类 (如 Classical/Basal-like)
# ============================================================

suppressPackageStartupMessages({ library(data.table) })

BULK <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk"
log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# 工具函数: 在所有 characteristics_* 字段里找匹配前缀的值
# 例如 find_by_prefix(row, "survival_months") -> 返回数值
# ============================================================
find_by_prefix <- function(row, prefix, char_cols) {
    for (cn in char_cols) {
        v <- row[[cn]]
        if (is.na(v) || v == "") next
        if (grepl(paste0("^\\s*", prefix, "\\s*:"), v, ignore.case = TRUE)) {
            # 去除前缀和冒号, 返回值
            val <- sub(paste0("^\\s*", prefix, "\\s*:\\s*"), "", v, ignore.case = TRUE)
            return(trimws(val))
        }
    }
    return(NA_character_)
}

# 应用于整个 data.table 的某个前缀
extract_field <- function(clin, prefix) {
    char_cols <- grep("^characteristics", colnames(clin), value = TRUE)
    sapply(seq_len(nrow(clin)), function(i) {
        row <- as.list(clin[i, ])
        find_by_prefix(row, prefix, char_cols)
    })
}

# ============================================================
# 1. GSE71729 (Moffitt 2015, Agilent dual-color)
# ============================================================
process_GSE71729 <- function() {
    log_msg("=== GSE71729 ===")
    clin <- fread(file.path(BULK, "GSE71729/GSE71729.clinical.tsv"))
    
    std <- data.table(
        sample_id     = clin$geo_accession,
        os_time_month = as.numeric(extract_field(clin, "survival_months")),
        os_status     = as.integer(extract_field(clin, "death_event_1death_0censor")),
        sample_type   = extract_field(clin, "tissue type"),
        cell_line     = extract_field(clin, "cell line/tissue"),
        tumor_subtype = extract_field(clin, "tumor_subtype_0na_1classical_2basal"),
        stroma_subtype = extract_field(clin, "stroma_subtype_0na_1low_2normal_3activated")
    )
    
    # tumor_subtype 编码转可读: 0=NA, 1=Classical, 2=Basal-like
    std[, tumor_subtype := c("0" = NA_character_, "1" = "Classical", "2" = "Basal-like")[tumor_subtype]]
    std[, stroma_subtype := c("0" = NA_character_, "1" = "Low", "2" = "Normal", "3" = "Activated")[stroma_subtype]]
    
    # sample_type 缺失时从 cell_line 推断 (CellLine 类样本没 tissue type)
    std[is.na(sample_type), sample_type := ifelse(grepl("CellLine|BXPC3|Capan|HPAC|HPAFII|Hs766T|HUPT3|MiaPaca|Panc|HPNE|SW1990|T3M4", cell_line, ignore.case = TRUE),
                                                   "CellLine", NA_character_)]
    
    log_msg(sprintf("  总: %d 样本", nrow(std)))
    log_msg("  sample_type 分布:")
    print(table(std$sample_type, useNA = "ifany"))
    
    # 用于生存分析的子集
    valid <- std[!is.na(os_time_month) & !is.na(os_status)]
    log_msg(sprintf("  含完整生存: %d 样本", nrow(valid)))
    log_msg("    其中 sample_type:")
    print(table(valid$sample_type, useNA = "ifany"))
    
    return(std)
}

# ============================================================
# 2. GSE21501 (Stratford 2010, Agilent dual-color)
# ============================================================
process_GSE21501 <- function() {
    log_msg("=== GSE21501 ===")
    clin <- fread(file.path(BULK, "GSE21501/GSE21501.clinical.tsv"))
    
    std <- data.table(
        sample_id     = clin$geo_accession,
        os_time_month = as.numeric(extract_field(clin, "os time")),
        os_status     = as.integer(extract_field(clin, "os event")),
        sample_type   = "Primary",  # 全是 PDAC primary, 但有 30 个 mix/normal 需要确认
        tissue        = extract_field(clin, "tissue"),
        tumor_type    = extract_field(clin, "tumor type"),
        t_stage       = extract_field(clin, "t stage"),
        n_stage       = extract_field(clin, "n stage"),
        risk_group    = extract_field(clin, "risk group")
    )
    
    # 没有 OS 数据的样本可能是 reference/normal mix
    std[is.na(os_time_month), sample_type := "Non-tumor"]
    
    log_msg(sprintf("  总: %d 样本", nrow(std)))
    log_msg(sprintf("  含完整生存: %d", nrow(std[!is.na(os_time_month) & !is.na(os_status)])))
    log_msg("  T stage 分布:")
    print(table(std$t_stage, useNA = "ifany"))
    
    return(std)
}

# ============================================================
# 3. GSE85916 (Kirby 2016, Affymetrix HuGene 2.0)
# ============================================================
process_GSE85916 <- function() {
    log_msg("=== GSE85916 ===")
    clin <- fread(file.path(BULK, "GSE85916/GSE85916.clinical.tsv"))
    
    std <- data.table(
        sample_id     = clin$geo_accession,
        # os.year 单位是年, 转月
        os_time_month = as.numeric(extract_field(clin, "os.year")) * 12,
        os_status     = as.integer(extract_field(clin, "death")),
        sample_type   = "Primary"  # 全是 PDAC primary
    )
    
    log_msg(sprintf("  总: %d 样本", nrow(std)))
    log_msg(sprintf("  含完整生存: %d", nrow(std[!is.na(os_time_month) & !is.na(os_status)])))
    log_msg(sprintf("  OS_time 范围 (月): %.1f - %.1f", min(std$os_time_month, na.rm = TRUE),
                                                          max(std$os_time_month, na.rm = TRUE)))
    
    return(std)
}

# ============================================================
# 4. GSE62452 (Yang 2016, Affymetrix HuGene 1.0)
# ============================================================
process_GSE62452 <- function() {
    log_msg("=== GSE62452 ===")
    clin <- fread(file.path(BULK, "GSE62452/GSE62452.clinical.tsv"))
    
    std <- data.table(
        sample_id     = clin$geo_accession,
        os_time_month = as.numeric(extract_field(clin, "survival months")),
        # "?" 转 NA
        os_status_raw = extract_field(clin, "survival status"),
        tissue        = extract_field(clin, "tissue"),
        grade         = extract_field(clin, "grading"),
        stage         = extract_field(clin, "Stage")
    )
    
    std[, os_status := as.integer(ifelse(os_status_raw %in% c("0", "1"), os_status_raw, NA))]
    std[, os_status_raw := NULL]
    
    # tumor vs adjacent normal
    std[, sample_type := ifelse(grepl("tumor", tissue, ignore.case = TRUE) & 
                                !grepl("non-tumor|adjacent", tissue, ignore.case = TRUE),
                                "Primary", "Normal")]
    
    log_msg(sprintf("  总: %d 样本", nrow(std)))
    log_msg("  sample_type 分布:")
    print(table(std$sample_type, useNA = "ifany"))
    log_msg(sprintf("  含完整生存 (含 ? 转 NA): %d", 
                    nrow(std[!is.na(os_time_month) & !is.na(os_status)])))
    log_msg("  stage 分布:")
    print(table(std$stage, useNA = "ifany"))
    
    return(std)
}

# ============================================================
# 5. GSE57495 (Chen 2015, Affymetrix HuGene 1.0)
# ============================================================
process_GSE57495 <- function() {
    log_msg("=== GSE57495 ===")
    clin <- fread(file.path(BULK, "GSE57495/GSE57495.clinical.tsv"))
    
    std <- data.table(
        sample_id     = clin$geo_accession,
        os_time_month = as.numeric(extract_field(clin, "overall survival \\(month\\)")),
        # DEAD=1, ALIVE=0
        vital_raw     = extract_field(clin, "vital.status"),
        stage         = extract_field(clin, "Stage"),
        sample_type   = "Primary"
    )
    std[, os_status := as.integer(ifelse(vital_raw == "DEAD", 1, 
                                  ifelse(vital_raw == "ALIVE", 0, NA)))]
    std[, vital_raw := NULL]
    
    log_msg(sprintf("  总: %d 样本", nrow(std)))
    log_msg(sprintf("  含完整生存: %d", nrow(std[!is.na(os_time_month) & !is.na(os_status)])))
    log_msg("  stage 分布:")
    print(table(std$stage, useNA = "ifany"))
    
    return(std)
}

# ============================================================
# 主流程
# ============================================================
results <- list(
    GSE71729 = process_GSE71729(),
    GSE21501 = process_GSE21501(),
    GSE85916 = process_GSE85916(),
    GSE62452 = process_GSE62452(),
    GSE57495 = process_GSE57495()
)

# 写出
cat("\n")
log_msg("=== 写出标准化临床表 ===")
for (gse in names(results)) {
    out <- file.path(BULK, gse, sprintf("%s.clinical_std.tsv", gse))
    fwrite(results[[gse]], out, sep = "\t")
    log_msg(sprintf("  ✓ %s", basename(out)))
}

# 总结
cat("\n")
log_msg("=== 最终可用样本量 ===")
for (gse in names(results)) {
    r <- results[[gse]]
    valid <- r[!is.na(os_time_month) & !is.na(os_status) & sample_type == "Primary"]
    log_msg(sprintf("  %s: %d 总 → %d 含完整生存的 Primary", gse, nrow(r), nrow(valid)))
}
