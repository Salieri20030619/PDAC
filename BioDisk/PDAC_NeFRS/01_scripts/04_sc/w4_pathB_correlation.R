# 接续: UCell scores 已存, 计算相关性
suppressPackageStartupMessages({ library(data.table); library(Matrix) })

OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"

cat("=== NeFRS vs Hwang 17 programs Spearman correlation ===\n")

# 加载 UCell scores
ucell_scores <- readRDS(file.path(OUT, "w4_ucell_nefrs_plus_hwang17.rds"))
cat(sprintf("UCell scores: %d x %d\n", nrow(ucell_scores), ncol(ucell_scores)))

# 加载 annot (注意它有空格列名)
annot <- readRDS(file.path(OUT, "w4_annot_with_ucell.rds"))
annot <- as.data.table(annot)
# 重命名空格列
setnames(annot, "Level 1 Annotation", "Level1", skip_absent=TRUE)
setnames(annot, "Level 2 Annotation", "Level2", skip_absent=TRUE)
setnames(annot, "Level 3 Annotation", "Level3", skip_absent=TRUE)
cat(sprintf("annot cols: %s\n", paste(colnames(annot)[1:10], collapse=", ")))

# 合并 UCell scores
score_df <- data.table(barcode = rownames(ucell_scores))
for (col in colnames(ucell_scores)) {
  score_df[[col]] <- ucell_scores[, col]
}
annot2 <- merge(annot[, .(barcode, Level1, Level3, treatment_status)],
                 score_df, by="barcode", sort=FALSE)
cat(sprintf("Merged: %d cells\n", nrow(annot2)))

# 关键列
nefrs_col <- "NeFRS_all_UCell"
hwang_cols <- setdiff(colnames(ucell_scores),
                       c("NeFRS_all_UCell","NeFRS_risk_UCell","NeFRS_protective_UCell"))

# Subset
malig <- annot2[Level1 == "Epithelial (malignant)"]
ut_malig <- annot2[Level1 == "Epithelial (malignant)" & treatment_status == "Untreated"]
cat(sprintf("\nAll cells: %d\n", nrow(annot2)))
cat(sprintf("Malignant cells: %d\n", nrow(malig)))
cat(sprintf("Untreated malignant: %d\n", nrow(ut_malig)))

# === 主表: 3 subset × 3 NeFRS × 17 Hwang ===
cor_results <- data.table()
for (nefrs in c("NeFRS_all_UCell", "NeFRS_risk_UCell", "NeFRS_protective_UCell")) {
  for (h in hwang_cols) {
    rho_all <- cor(annot2[[nefrs]], annot2[[h]],
                    method="spearman", use="pairwise.complete.obs")
    rho_malig <- cor(malig[[nefrs]], malig[[h]],
                      method="spearman", use="pairwise.complete.obs")
    rho_ut_malig <- cor(ut_malig[[nefrs]], ut_malig[[h]],
                         method="spearman", use="pairwise.complete.obs")
    cor_results <- rbind(cor_results, data.table(
      NeFRS_type = nefrs,
      Hwang_program = h,
      rho_all = round(rho_all, 3),
      rho_malig = round(rho_malig, 3),
      rho_ut_malig = round(rho_ut_malig, 3)
    ))
  }
}

# 主关注: NeFRS_all 跟 17 程序 (按 rho_malig 排序)
cat("\n=== Main: NeFRS_all vs Hwang 17 (sorted by rho in Malignant) ===\n")
main_table <- cor_results[NeFRS_type == "NeFRS_all_UCell"][order(-rho_malig)]
print(main_table[, .(Hwang_program, rho_all, rho_malig, rho_ut_malig)])

cat("\n=== Risk-only NeFRS (β>0 27 基因) vs Hwang 17 ===\n")
risk_table <- cor_results[NeFRS_type == "NeFRS_risk_UCell"][order(-rho_malig)]
print(risk_table[, .(Hwang_program, rho_all, rho_malig, rho_ut_malig)])

cat("\n=== Protective-only NeFRS (β<0 26 基因) vs Hwang 17 ===\n")
prot_table <- cor_results[NeFRS_type == "NeFRS_protective_UCell"][order(-rho_malig)]
print(prot_table[, .(Hwang_program, rho_all, rho_malig, rho_ut_malig)])

# === 保存 ===
fwrite(cor_results, file.path(OUT, "w4_nefrs_vs_hwang17_correlation.tsv"), sep="\t")
cat(sprintf("\nSaved: %s/w4_nefrs_vs_hwang17_correlation.tsv\n", OUT))

# === 关键问题: NeFRS vs 真实 NRP 在 Malignant 内的 rho ===
nrp_col <- "mal_lineage_Neural-like_progenitor_UCell"
ut_malig_v <- ut_malig[!is.na(get(nefrs_col)) & !is.na(get(nrp_col))]
cat(sprintf("\n=== ⭐ 关键: NeFRS_all vs Hwang NRP 200 ===\n"))
cat(sprintf("In all cells:       rho = %+.3f (n=%d)\n",
            cor(annot2[[nefrs_col]], annot2[[nrp_col]], method="spearman"),
            nrow(annot2)))
cat(sprintf("In Malignant only:  rho = %+.3f (n=%d)\n",
            cor(malig[[nefrs_col]], malig[[nrp_col]], method="spearman"),
            nrow(malig)))
cat(sprintf("In Untreated Mal:   rho = %+.3f (n=%d)\n",
            cor(ut_malig[[nefrs_col]], ut_malig[[nrp_col]], method="spearman"),
            nrow(ut_malig)))

# === 跨 cell type 内, NeFRS 跟 NRP 在每个 Level 3 cell type 的 rho ===
cat("\n=== NeFRS vs NRP rho by Level 3 cell type (top 15 cell types) ===\n")
rho_by_l3 <- annot2[, .(rho = cor(get(nefrs_col), get(nrp_col), method="spearman"),
                         n = .N), by=Level3][n >= 100][order(-rho)]
print(head(rho_by_l3, 15))

cat("\n=== Done ===\n")
