#!/usr/bin/env Rscript
# ============================================================
# W5 Step 1: TCGA-LIHC (HCC) + TCGA-COAD/READ (CRC) 预处理
# 照搬 process_tcga_paad.R 约定: star_tpm → ENSG去版本 → org.Hs.eg.db
#   → SYMBOL (1:many取第一) → 多ENSG按最大方差聚合 → gene-level log2(TPM+1)
# 临床: clinical.tsv.gz + survival.tsv.gz 合并 → 标准化字段 (与8队列一致)
# 输出 (不z-score, 与8队列 .gene.tsv.gz 同格式):
#   HCC/HCC.gene.tsv.gz + HCC.clinical_std.tsv
#   CRC/CRC.gene.tsv.gz + CRC.clinical_std.tsv   (COAD+READ合并)
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(AnnotationDbi); library(org.Hs.eg.db)
})
PAN <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/pancancer"
log_msg <- function(m) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), m))
gcol <- function(dt, nm) if (nm %in% names(dt)) dt[[nm]] else NA   # 字段缺失保护

# ---- 缓存 ENSG→SYMBOL 映射 (避免重复查 org.Hs.eg.db) ----
ENSMAP <- NULL
ens2sym <- function(ens_ids) {
  m <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(ens_ids),
                             columns = "SYMBOL", keytype = "ENSEMBL")
  setDT(m); m <- m[!is.na(SYMBOL), .(SYMBOL = SYMBOL[1]), by = ENSEMBL]
  setnames(m, c("ensembl_id", "gene_symbol")); m
}

# ---- 1) star_tpm → gene-level (最大方差聚合) ----
read_genelevel <- function(tpm_path) {
  log_msg(sprintf("  读 %s", basename(tpm_path)))
  expr <- fread(tpm_path, header = TRUE, showProgress = FALSE)
  setnames(expr, 1, "ensembl_versioned")
  expr[, ensembl_id := sub("\\..*$", "", ensembl_versioned)]
  map <- ens2sym(expr$ensembl_id)
  expr <- merge(expr, map, by = "ensembl_id", all.x = FALSE)
  scols <- setdiff(names(expr), c("ensembl_versioned", "ensembl_id", "gene_symbol"))
  log_msg(sprintf("    映射上 %d ENSG → %d 唯一基因 × %d 样本",
                  nrow(expr), uniqueN(expr$gene_symbol), length(scols)))
  expr[, v := apply(.SD, 1, var, na.rm = TRUE), .SDcols = scols]
  sel <- expr[, .SD[which.max(v)], by = gene_symbol]
  gm <- sel[, c("gene_symbol", scols), with = FALSE]; setnames(gm, "gene_symbol", "gene")
  gm
}

# ---- 2) 临床标准化 (照搬 PAAD std 字段) ----
build_std <- function(clin_path, surv_path) {
  clin <- fread(clin_path, na.strings = c("", "NA", "not reported", "Not Reported", "Unknown"))
  surv <- fread(surv_path, na.strings = c("", "NA"))
  cf <- merge(clin, surv, by = "sample", all.x = TRUE)
  std <- data.table(
    sample_id      = cf$sample,
    patient_id     = gcol(cf, "_PATIENT"),
    os_time_day    = as.numeric(gcol(cf, "OS.time")),
    os_time_month  = as.numeric(gcol(cf, "OS.time")) / 30.44,
    os_status      = as.integer(gcol(cf, "OS")),
    sample_type_code = substr(cf$sample, 14, 15),
    sample_type    = gcol(cf, "sample_type.samples"),
    sex            = gcol(cf, "gender.demographic"),
    age            = as.numeric(gcol(cf, "age_at_index.demographic")),
    stage          = gcol(cf, "ajcc_pathologic_stage.diagnoses"),
    grade          = gcol(cf, "tumor_grade.diagnoses"),
    primary_site   = gcol(cf, "primary_site"),
    primary_diagnosis = gcol(cf, "primary_diagnosis.diagnoses"),
    annotation_classification = gcol(cf, "classification.annotations")
  )
  std[, sample_type_simple := fcase(
    sample_type_code %in% c("01","01A","01B","01C"), "Primary",
    sample_type_code %in% c("11","11A","11B"), "Normal",
    sample_type_code %in% c("06","06A","02","02A"), "Metastasis_or_Recurrent",
    default = "Other")]
  std
}

# ---- 3) 轻度 curation: Primary + 完整生存 (+报告诊断分布,诊断过滤可调) ----
curate <- function(std, label, keep_dx = NULL) {
  s1 <- std[sample_type_simple == "Primary"]
  s2 <- s1[!is.na(os_time_month) & !is.na(os_status) & os_time_month > 0]
  log_msg(sprintf("  [%s] Primary %d → +完整生存 %d", label, nrow(s1), nrow(s2)))
  log_msg(sprintf("  [%s] primary_diagnosis 分布:", label))
  print(sort(table(s2$primary_diagnosis, useNA = "ifany"), decreasing = TRUE))
  if (!is.null(keep_dx)) {
    s3 <- s2[primary_diagnosis %in% keep_dx]
    log_msg(sprintf("  [%s] +诊断过滤 %d (排除 %d)", label, nrow(s3), nrow(s2) - nrow(s3)))
    s3
  } else s2
}

# ============================================================
# HCC = TCGA-LIHC
# ============================================================
log_msg("===== HCC (TCGA-LIHC) =====")
hcc_gm  <- read_genelevel(file.path(PAN, "HCC/TCGA-LIHC.star_tpm.tsv.gz"))
hcc_std <- build_std(file.path(PAN, "HCC/TCGA-LIHC.clinical.tsv.gz"),
                     file.path(PAN, "HCC/TCGA-LIHC.survival.tsv.gz"))
# HCC 诊断: 先 NULL 看全分布; 跑完按输出把非HCC(如胆管癌)加进排除即可
hcc_keep <- c("Hepatocellular carcinoma, NOS",
              "Hepatocellular carcinoma, fibrolamellar",
              "Hepatocellular carcinoma, clear cell type",
              "Hepatocellular carcinoma, spindle cell variant",
              "Hepatocellular carcinoma, sarcomatoid")
hcc_cur <- curate(hcc_std, "HCC", keep_dx = hcc_keep)
hcc_ids <- intersect(hcc_cur$sample_id, names(hcc_gm))
fwrite(hcc_gm[, c("gene", hcc_ids), with = FALSE],
       file.path(PAN, "HCC/HCC.gene.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(hcc_cur[sample_id %in% hcc_ids],
       file.path(PAN, "HCC/HCC.clinical_std.tsv"), sep = "\t")
log_msg(sprintf("  ✓ HCC: %d 基因 × %d 样本 (events %d)",
                nrow(hcc_gm), length(hcc_ids), sum(hcc_cur[sample_id %in% hcc_ids]$os_status)))

# ============================================================
# CRC = TCGA-COAD + TCGA-READ (合并)
# ============================================================
log_msg("===== CRC (TCGA-COAD + TCGA-READ) =====")
coad_gm <- read_genelevel(file.path(PAN, "CRC/TCGA-COAD.star_tpm.tsv.gz"))
read_gm <- read_genelevel(file.path(PAN, "CRC/TCGA-READ.star_tpm.tsv.gz"))
common  <- intersect(coad_gm$gene, read_gm$gene)
log_msg(sprintf("  COAD %d 基因, READ %d 基因, 共有 %d → 合并",
                nrow(coad_gm), nrow(read_gm), length(common)))
crc_gm <- merge(coad_gm[gene %in% common], read_gm[gene %in% common], by = "gene")
coad_std <- build_std(file.path(PAN, "CRC/TCGA-COAD.clinical.tsv.gz"),
                      file.path(PAN, "CRC/TCGA-COAD.survival.tsv.gz"))
read_std <- build_std(file.path(PAN, "CRC/TCGA-READ.clinical.tsv.gz"),
                      file.path(PAN, "CRC/TCGA-READ.survival.tsv.gz"))
crc_std <- rbind(coad_std, read_std)
crc_keep <- c("Adenocarcinoma, NOS", "Mucinous adenocarcinoma",
              "Adenocarcinoma with mixed subtypes", "Papillary adenocarcinoma, NOS")
crc_cur <- curate(crc_std, "CRC", keep_dx = crc_keep)
crc_ids <- intersect(crc_cur$sample_id, names(crc_gm))
fwrite(crc_gm[, c("gene", crc_ids), with = FALSE],
       file.path(PAN, "CRC/CRC.gene.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(crc_cur[sample_id %in% crc_ids],
       file.path(PAN, "CRC/CRC.clinical_std.tsv"), sep = "\t")
log_msg(sprintf("  ✓ CRC: %d 基因 × %d 样本 (events %d)",
                nrow(crc_gm), length(crc_ids), sum(crc_cur[sample_id %in% crc_ids]$os_status)))

# ---- 关键基因健康检查 (NeFRS hub + 持家) ----
log_msg("===== 关键基因健康检查 =====")
chk <- c("NRP1","EGFR","BMPR1B","GRN","CDC25A","REST","WWTR1","GAPDH","ACTB")
for (nm in c("HCC","CRC")) {
  g <- if (nm == "HCC") hcc_gm$gene else crc_gm$gene
  log_msg(sprintf("  [%s] %s", nm, paste(sprintf("%s%s", chk, ifelse(chk %in% g, "✓","✗")), collapse=" ")))
}
log_msg("Done.")
