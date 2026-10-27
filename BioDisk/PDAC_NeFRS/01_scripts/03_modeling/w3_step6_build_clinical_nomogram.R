# W3 Day 4 Step A: 重新构造 3 队列临床表 (含完整 age/sex/stage/grade)
# 从原始 clinical_full 重新读, 不依赖之前不完整的 clinical_std

suppressPackageStartupMessages({
  library(data.table); library(survival); library(glmnet)
})

RAW   <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk"
W2_IN <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/02_w2_inputs"
MIME  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime"
OUT   <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/05_nomogram"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

cat("=== W3 Day 4 Step A: Rebuild clinical tables ===\n")

# === 加载基础数据 ===
cohort_data <- readRDS(file.path(W2_IN, "cohort_data_8.rds"))
mime_input  <- readRDS(file.path(MIME, "mime_input_list.rds"))
res         <- readRDS(file.path(MIME, "mime_results.rds"))
mime_genes  <- colnames(mime_input[[1]])[-(1:3)]
ml <- res$ml.res[["Enet[α=0.1]"]]

# 计算 NeFRS
compute_rs <- function(cn) {
  co <- cohort_data[[cn]]
  X <- matrix(0, nrow=ncol(co$expression), ncol=length(mime_genes),
              dimnames=list(colnames(co$expression), mime_genes))
  avail <- intersect(mime_genes, rownames(co$expression))
  for (g in avail) X[, g] <- co$expression[g, ]
  rs <- as.numeric(predict(ml, newx=X, s="lambda.min"))
  names(rs) <- rownames(X)
  rs
}

# === Helper: TNM → AJCC stage (PDAC 2010 AJCC 7th) ===
tnm_to_ajcc <- function(tnm_str) {
  if (is.na(tnm_str) || tnm_str == "" || tnm_str == "TXNXMX") return(NA)
  
  # 提取 T/N/M (去掉 a/b/c sub-stage 修饰)
  t_val <- sub(".*?T(\\d|X).*", "\\1", tnm_str)
  n_val <- sub(".*?N(\\d|X).*", "\\1", tnm_str)
  m_val <- sub(".*?M(\\d|X).*", "\\1", tnm_str)
  
  # M1 → IV
  if (m_val == "1") return("IV")
  # T4 → III
  if (t_val == "4") return("III")
  # TX/NX → NA
  if (t_val == "X") return(NA)
  
  # T1N0 → IA, T1N1 → IIB
  # T2N0 → IB, T2N1 → IIB
  # T3N0 → IIA, T3N1 → IIB
  if (t_val == "1" && n_val == "0") return("IA")
  if (t_val == "2" && n_val == "0") return("IB")
  if (t_val == "3" && n_val == "0") return("IIA")
  if (t_val %in% c("1","2","3") && n_val == "1") return("IIB")
  if (t_val %in% c("1","2","3") && n_val == "X") return("II")  # 不确定 N → II 大类
  
  return(NA)
}

# 简化 4 类
simplify_stage <- function(ajcc) {
  if (is.na(ajcc)) return(NA)
  if (ajcc %in% c("I", "IA", "IB")) return("I")
  if (ajcc == "IIA") return("IIA")
  if (ajcc == "IIB" || ajcc == "II") return("IIB")
  if (ajcc %in% c("III", "IV")) return("III-IV")
  return(NA)
}

# === 1. PACA-AU clinical_full ===
cat("\n--- PACA-AU ---\n")
paca_full <- fread(file.path(RAW, "ICGC_PACA_AU/clinical_full.PACA-AU.tsv"))
paca_std  <- fread(file.path(RAW, "ICGC_PACA_AU/clinical_std.PACA-AU.tsv"))
cat(sprintf("  clinical_full: %d rows, %d cols\n", nrow(paca_full), ncol(paca_full)))
cat(sprintf("  clinical_std: %d rows\n", nrow(paca_std)))

# clinical_full 用 icgc_donor_id, clinical_std 用 sample_id
# 通过 std 表的 donor_id 关联
paca_full_keys <- paca_full[, .(donor_id = icgc_donor_id,
                                  age = as.numeric(donor_age_at_diagnosis),
                                  sex_raw = donor_sex,
                                  stage_raw = tumour_stage,
                                  grade_raw = tumour_grade)]

paca_clin <- merge(paca_std[, .(sample_id, donor_id, os_time_month, os_status)],
                    unique(paca_full_keys, by="donor_id"),
                    by="donor_id", all.x=TRUE)

# 加 NeFRS
paca_rs <- compute_rs("PACA_AU")
paca_clin$NeFRS <- paca_rs[paca_clin$sample_id]

# 协变量统一化
paca_clin[, Sex := factor(ifelse(sex_raw == "male", "Male", "Female"),
                           levels=c("Female","Male"))]
paca_clin[, Stage_AJCC := factor(sapply(stage_raw, function(s) simplify_stage(tnm_to_ajcc(s))),
                                  levels=c("I","IIA","IIB","III-IV"))]
paca_clin[, Grade := factor(
  ifelse(grepl("^1", grade_raw), "G1",
  ifelse(grepl("^2", grade_raw), "G2",
  ifelse(grepl("^[34]", grade_raw), "G3-G4", NA))),
  levels=c("G1","G2","G3-G4"))]
paca_clin[, age := as.numeric(age)]

cat("\n  Coverage:\n")
cat(sprintf("    NeFRS: %d/%d (%.0f%%)\n", sum(!is.na(paca_clin$NeFRS)), nrow(paca_clin),
            100*sum(!is.na(paca_clin$NeFRS))/nrow(paca_clin)))
cat(sprintf("    age:   %d/%d (%.0f%%)\n", sum(!is.na(paca_clin$age)), nrow(paca_clin),
            100*sum(!is.na(paca_clin$age))/nrow(paca_clin)))
cat(sprintf("    Sex:   %d/%d (%.0f%%)\n", sum(!is.na(paca_clin$Sex)), nrow(paca_clin),
            100*sum(!is.na(paca_clin$Sex))/nrow(paca_clin)))
cat(sprintf("    Stage: %d/%d (%.0f%%)\n", sum(!is.na(paca_clin$Stage_AJCC)), nrow(paca_clin),
            100*sum(!is.na(paca_clin$Stage_AJCC))/nrow(paca_clin)))
cat(sprintf("    Grade: %d/%d (%.0f%%)\n", sum(!is.na(paca_clin$Grade)), nrow(paca_clin),
            100*sum(!is.na(paca_clin$Grade))/nrow(paca_clin)))

cat("\n  Distribution:\n")
cat("  Sex:");        print(table(paca_clin$Sex, useNA="ifany"))
cat("  Stage_AJCC:"); print(table(paca_clin$Stage_AJCC, useNA="ifany"))
cat("  Grade:");      print(table(paca_clin$Grade, useNA="ifany"))
cat(sprintf("  age range: %.0f - %.0f, median %.0f\n",
            min(paca_clin$age, na.rm=TRUE), max(paca_clin$age, na.rm=TRUE),
            median(paca_clin$age, na.rm=TRUE)))

# === 2. TCGA clinical_full ===
cat("\n--- TCGA ---\n")
tcga_full <- fread(file.path(RAW, "TCGA_PAAD/TCGA-PAAD.clinical_full.tsv"))
tcga_curated <- fread(file.path(RAW, "TCGA_PAAD/TCGA-PAAD.curated.clinical_std.tsv"))
cat(sprintf("  clinical_full: %d rows\n", nrow(tcga_full)))
cat(sprintf("  curated: %d rows\n", nrow(tcga_curated)))

# clinical_full 用 'sample' (TCGA-XX-XXXX-01A)
# curated 用 sample_id (TCGA-XX-XXXX-01A)
tcga_full_keys <- tcga_full[, .(sample_id = sample,
                                  age = as.numeric(age_at_index.demographic),
                                  sex_raw = gender.demographic,
                                  stage_raw = ajcc_pathologic_stage.diagnoses)]
# Note: tumor_grade.diagnoses is all empty, skip grade for TCGA

tcga_clin <- merge(tcga_curated[, .(sample_id, os_time_month, os_status)],
                    unique(tcga_full_keys, by="sample_id"),
                    by="sample_id", all.x=TRUE)

tcga_rs <- compute_rs("TCGA_PAAD")
tcga_clin$NeFRS <- tcga_rs[tcga_clin$sample_id]

tcga_clin[, Sex := factor(ifelse(sex_raw == "male", "Male", "Female"),
                           levels=c("Female","Male"))]
# TCGA stage: "Stage IA" → "I", "Stage IIA" → "IIA", etc
tcga_clin[, Stage_simple := sub("^Stage ", "", stage_raw)]
tcga_clin[, Stage_AJCC := factor(
  ifelse(Stage_simple %in% c("I","IA","IB"), "I",
  ifelse(Stage_simple == "IIA", "IIA",
  ifelse(Stage_simple == "IIB", "IIB",
  ifelse(Stage_simple %in% c("III","IV"), "III-IV", NA)))),
  levels=c("I","IIA","IIB","III-IV"))]
tcga_clin[, Grade := factor(NA, levels=c("G1","G2","G3-G4"))]  # TCGA grade 空, 全 NA

cat("\n  Coverage:\n")
cat(sprintf("    NeFRS: %d/%d\n", sum(!is.na(tcga_clin$NeFRS)), nrow(tcga_clin)))
cat(sprintf("    age:   %d/%d\n", sum(!is.na(tcga_clin$age)), nrow(tcga_clin)))
cat(sprintf("    Stage: %d/%d\n", sum(!is.na(tcga_clin$Stage_AJCC)), nrow(tcga_clin)))
cat(sprintf("    Grade: 0/%d (TCGA tumor_grade field empty in GDC)\n", nrow(tcga_clin)))

cat("\n  Distribution:\n")
cat("  Sex:");        print(table(tcga_clin$Sex, useNA="ifany"))
cat("  Stage_AJCC:"); print(table(tcga_clin$Stage_AJCC, useNA="ifany"))

# === 3. E-MTAB-6134 ===
cat("\n--- E-MTAB-6134 ---\n")
emtab <- fread(file.path(RAW, "E_MTAB_6134/clinical_std.E-MTAB-6134.tsv"))
cat(sprintf("  clinical_std: %d rows\n", nrow(emtab)))

emtab_clin <- emtab[, .(sample_id, os_time_month = OS_time, os_status = OS_status,
                         sex_raw = sex, stage_raw = stage_TNM, grade_raw = grade)]
emtab_rs <- compute_rs("E_MTAB_6134")
emtab_clin$NeFRS <- emtab_rs[emtab_clin$sample_id]

emtab_clin[, age := NA_real_]  # E-MTAB 无 age
emtab_clin[, Sex := factor(ifelse(sex_raw == "male", "Male", "Female"),
                            levels=c("Female","Male"))]
emtab_clin[, Stage_AJCC := factor(sapply(stage_raw, function(s) simplify_stage(tnm_to_ajcc(s))),
                                   levels=c("I","IIA","IIB","III-IV"))]
emtab_clin[, Grade := factor(
  ifelse(grade_raw == "tumour grading G1", "G1",
  ifelse(grade_raw == "tumour grading G2", "G2",
  ifelse(grade_raw == "tumour grading G3", "G3-G4", NA))),
  levels=c("G1","G2","G3-G4"))]

cat("\n  Coverage:\n")
cat(sprintf("    NeFRS: %d/%d\n", sum(!is.na(emtab_clin$NeFRS)), nrow(emtab_clin)))
cat(sprintf("    age:   0/%d (E-MTAB-6134 真无 patient age, sdrf 验证)\n", nrow(emtab_clin)))
cat(sprintf("    Stage: %d/%d\n", sum(!is.na(emtab_clin$Stage_AJCC)), nrow(emtab_clin)))
cat(sprintf("    Grade: %d/%d\n", sum(!is.na(emtab_clin$Grade)), nrow(emtab_clin)))

cat("\n  Distribution:\n")
cat("  Sex:");        print(table(emtab_clin$Sex, useNA="ifany"))
cat("  Stage_AJCC:"); print(table(emtab_clin$Stage_AJCC, useNA="ifany"))
cat("  Grade:");      print(table(emtab_clin$Grade, useNA="ifany"))

# === 保存 ===
out_cols <- c("sample_id","os_time_month","os_status","NeFRS",
              "age","Sex","Stage_AJCC","Grade")
clin_list <- list(
  PACA_AU = paca_clin[, ..out_cols],
  TCGA_PAAD = tcga_clin[, ..out_cols],
  E_MTAB_6134 = emtab_clin[, ..out_cols]
)

saveRDS(clin_list, file.path(OUT, "clinical_with_nefrs_3cohorts.rds"))
for (cn in names(clin_list)) {
  fwrite(clin_list[[cn]], file.path(OUT, sprintf("clinical_%s.tsv", cn)), sep="\t")
}

cat("\n=== Step A 完成 ===\n")
cat(sprintf("Saved: %s/clinical_with_nefrs_3cohorts.rds\n", OUT))
cat("3 队列样本量:\n")
for (cn in names(clin_list)) {
  d <- clin_list[[cn]]
  cat(sprintf("  %-15s: n=%d, NeFRS=%d, age=%d, Stage=%d, Grade=%d\n",
              cn, nrow(d),
              sum(!is.na(d$NeFRS)),
              sum(!is.na(d$age)),
              sum(!is.na(d$Stage_AJCC)),
              sum(!is.na(d$Grade))))
}
