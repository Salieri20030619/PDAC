# 接续: 从已存的 UCell scores 跑 NeFRS vs Hwang 17 相关性
suppressPackageStartupMessages({ library(data.table) })

OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"

cat("=== Path B Correlations ===\n")
ucell_scores <- readRDS(file.path(OUT, "w4_ucell_nefrs_plus_hwang17.rds"))
annot <- readRDS(file.path(OUT, "w4_annot_with_ucell.rds"))

cat(sprintf("UCell: %d x %d\n", nrow(ucell_scores), ncol(ucell_scores)))
cat(sprintf("Annot: %d x %d\n", nrow(annot), ncol(annot)))

# 合并 — 注意 annot 的列名带空格, 用 backtick
score_df <- data.table(barcode = rownames(ucell_scores))
for (col in colnames(ucell_scores)) {
  score_df[[col]] <- as.numeric(ucell_scores[, col])
}

# 把 annot 关键列重命名再 merge
annot_clean <- copy(annot)
setnames(annot_clean, "Level 1 Annotation", "Level1")
setnames(annot_clean, "Level 2 Annotation", "Level2")
setnames(annot_clean, "Level 3 Annotation", "Level3")

# 用 annot_clean 里的旧 NeFRS_all 重命名避免冲突
annot_clean[, c("NeFRS_all","NeFRS_risk","NeFRS_protective") := NULL]

# 直接 merge
m <- merge(annot_clean[, .(barcode, Level1, Level3, treatment_status)],
            score_df, by="barcode", sort=FALSE)
cat(sprintf("Merged: %d cells\n", nrow(m)))

nefrs_col <- "NeFRS_all_UCell"
hwang_cols <- setdiff(colnames(ucell_scores),
                       c("NeFRS_all_UCell","NeFRS_risk_UCell","NeFRS_protective_UCell"))

# === A. 全 224k cells ===
cat("\n=== Spearman: NeFRS_all vs Hwang 17 programs (ALL 224k cells) ===\n")
result_all <- data.table()
for (col in hwang_cols) {
  rho <- cor(m[[nefrs_col]], m[[col]], method="spearman", use="pairwise.complete.obs")
  result_all <- rbind(result_all, data.table(program=col, rho=round(rho, 3)))
}
result_all <- result_all[order(-rho)]
print(result_all)

# === B. Malignant only ===
cat("\n=== Spearman: NeFRS_all vs Hwang (Malignant only) ===\n")
malig <- m[Level1 == "Epithelial (malignant)"]
cat(sprintf("n = %d\n", nrow(malig)))
result_malig <- data.table()
for (col in hwang_cols) {
  rho <- cor(malig[[nefrs_col]], malig[[col]], method="spearman", use="pairwise.complete.obs")
  result_malig <- rbind(result_malig, data.table(program=col, rho=round(rho, 3)))
}
result_malig <- result_malig[order(-rho)]
print(result_malig)

# === C. Untreated Malignant ===
cat("\n=== Spearman: NeFRS_all vs Hwang (Untreated Malignant only) ===\n")
ut_malig <- m[Level1 == "Epithelial (malignant)" & treatment_status == "Untreated"]
cat(sprintf("n = %d\n", nrow(ut_malig)))
result_ut <- data.table()
for (col in hwang_cols) {
  rho <- cor(ut_malig[[nefrs_col]], ut_malig[[col]], method="spearman", use="pairwise.complete.obs")
  result_ut <- rbind(result_ut, data.table(program=col, rho=round(rho, 3)))
}
result_ut <- result_ut[order(-rho)]
print(result_ut)

# === D. 综合表 ===
cat("\n=== 综合相关性 (3 子集) ===\n")
combined <- merge(result_all[, .(program, rho_all=rho)],
                   result_malig[, .(program, rho_malig=rho)], by="program")
combined <- merge(combined, result_ut[, .(program, rho_unt_malig=rho)], by="program")
combined <- combined[order(-rho_malig)]
print(combined)

fwrite(combined, file.path(OUT, "w4_nefrs_vs_hwang17_correlation.tsv"), sep="\t")
cat("\n=== Done ===\n")

# === E. NeFRS_risk vs NRP 单独详细 (Hwang NRP 是论文核心) ===
cat("\n=== 重点对比: NeFRS_risk vs Hwang NRP (3 子集) ===\n")
nrp_col <- "mal_lineage_Neural-like_progenitor_UCell"
for (sub_name in c("ALL", "Malignant", "Untreated_Malignant")) {
  if (sub_name == "ALL") sub <- m
  else if (sub_name == "Malignant") sub <- malig
  else sub <- ut_malig
  
  rho_all  <- cor(sub$NeFRS_all_UCell,        sub[[nrp_col]], method="spearman")
  rho_risk <- cor(sub$NeFRS_risk_UCell,       sub[[nrp_col]], method="spearman")
  rho_prot <- cor(sub$NeFRS_protective_UCell, sub[[nrp_col]], method="spearman")
  
  cat(sprintf("\n--- %s (n=%d) ---\n", sub_name, nrow(sub)))
  cat(sprintf("  NeFRS_all        vs NRP: rho = %+.3f\n", rho_all))
  cat(sprintf("  NeFRS_risk (β>0) vs NRP: rho = %+.3f\n", rho_risk))
  cat(sprintf("  NeFRS_protective vs NRP: rho = %+.3f\n", rho_prot))
}

cat("\n=== 终极 done ===\n")
