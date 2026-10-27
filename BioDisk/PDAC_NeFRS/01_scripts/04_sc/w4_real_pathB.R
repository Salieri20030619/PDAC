# 真实路径 B: NeFRS vs Hwang 17 programs UCell 相关性分析
# 不再用 W1 错误的"Hwang NRP 49", 用真实 Hwang Sup Table 2 的 200 基因列表
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(UCell)
})

OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
HWANG_FILE <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/hwang2022_supp/hwang2022_programs_clean.tsv"

cat("=== Path B Real: NeFRS vs Hwang programs ===\n")

# === 1. 读 Hwang 17 程序 ===
hwang <- fread(HWANG_FILE)
cat(sprintf("Hwang programs file: %d rows\n", nrow(hwang)))
cat(sprintf("Categories: %s\n", paste(unique(hwang$category), collapse=" | ")))
cat(sprintf("Programs (17): %s\n", paste(unique(hwang$program), collapse=" | ")))

# 构建 17 个程序的基因 list
hwang_progs <- split(hwang$gene, hwang$program)
# 为防 Adhesive 名字冲突, 用 category+program 命名
hwang$key <- paste0(
  ifelse(hwang$category == "Malignant state programs", "mal_state",
         ifelse(hwang$category == "Malignant lineage programs", "mal_lineage", "CAF")),
  "_", gsub(" ", "_", hwang$program)
)
prog_list <- split(hwang$gene, hwang$key)
cat(sprintf("\n=== Programs to score (n=%d) ===\n", length(prog_list)))
for (k in names(prog_list)) {
  cat(sprintf("  %s: %d genes\n", k, length(prog_list[[k]])))
}

# === 2. 加载 X (sparse) ===
cat("\nLoading X_sparse_full.rds...\n")
X <- readRDS(file.path(OUT, "w4_X_sparse_full.rds"))
gene_names <- colnames(X)
cell_barcodes <- rownames(X)
cat(sprintf("X: %d cells x %d genes\n", nrow(X), ncol(X)))

# 检查每个程序在 GSE202051 里有多少基因
cat("\n=== Hwang program gene availability in GSE202051 ===\n")
for (k in names(prog_list)) {
  avail <- intersect(prog_list[[k]], gene_names)
  cat(sprintf("  %-40s : %d / %d (%.0f%%)\n",
              k, length(avail), length(prog_list[[k]]),
              100*length(avail)/length(prog_list[[k]])))
}

# === 3. 转置 X ===
cat("\nTransposing...\n")
X_t <- t(X)

# === 4. 加上 NeFRS 3 signatures ===
coef_df <- fread("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime/nefrs_v3_enet_coefficients.tsv")
nefrs_signatures <- list(
  NeFRS_all = coef_df$gene,
  NeFRS_risk = coef_df[coef > 0]$gene,
  NeFRS_protective = coef_df[coef < 0]$gene
)

# 合并所有 signatures (3 NeFRS + 17 Hwang)
all_signatures <- c(nefrs_signatures, prog_list)
cat(sprintf("\nTotal signatures to score: %d (3 NeFRS + 17 Hwang)\n", length(all_signatures)))

# === 5. UCell ===
cat("\nRunning UCell (estimated 5-15 min for 20 signatures)...\n")
t0 <- Sys.time()
ucell_scores <- ScoreSignatures_UCell(matrix=X_t, features=all_signatures,
                                       ncores=8, chunk.size=5000)
cat(sprintf("Done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units="mins"))))
cat(sprintf("Output: %d x %d\n", nrow(ucell_scores), ncol(ucell_scores)))
cat(sprintf("Columns: %s\n", paste(colnames(ucell_scores), collapse=" | ")))

# === 6. 保存 ===
saveRDS(ucell_scores, file.path(OUT, "w4_ucell_nefrs_plus_hwang17.rds"))
cat("Saved: w4_ucell_nefrs_plus_hwang17.rds\n")

# === 7. NeFRS vs Hwang programs 相关性 (Spearman) ===
cat("\n=== Spearman correlation: NeFRS_all vs Hwang 17 programs ===\n")
# 在 Malignant 细胞内算 (Hwang programs 主要针对 Malignant)
annot <- readRDS(file.path(OUT, "w4_annot_with_ucell.rds"))
cat(sprintf("Annot rows: %d\n", nrow(annot)))

# 合并 ucell 到 annot
score_df <- data.table(barcode = rownames(ucell_scores))
for (col in colnames(ucell_scores)) {
  score_df[[col]] <- ucell_scores[, col]
}
annot2 <- merge(annot[, .(barcode, Level1=`Level 1 Annotation`, Level3=`Level 3 Annotation`, treatment_status)],
                 score_df, by="barcode", sort=FALSE)
cat(sprintf("Merged: %d cells\n", nrow(annot2)))

# 总体 Spearman: NeFRS_all_UCell vs 17 Hwang programs
nefrs_col <- "NeFRS_all_UCell"
hwang_cols <- setdiff(colnames(ucell_scores), c("NeFRS_all_UCell","NeFRS_risk_UCell","NeFRS_protective_UCell"))
cat(sprintf("\nNeFRS_all_UCell vs 17 Hwang programs (all 224k cells):\n"))
for (col in hwang_cols) {
  cor_all <- cor(annot2[[nefrs_col]], annot2[[col]], method="spearman", use="pairwise.complete.obs")
  cat(sprintf("  %-40s : rho = %+.3f\n", col, cor_all))
}

# Malignant 内 Spearman
cat("\n=== Within Malignant cells only ===\n")
malig <- annot2[Level1 == "Epithelial (malignant)"]
cat(sprintf("Malignant cells: %d\n", nrow(malig)))
for (col in hwang_cols) {
  cor_m <- cor(malig[[nefrs_col]], malig[[col]], method="spearman", use="pairwise.complete.obs")
  cat(sprintf("  %-40s : rho = %+.3f\n", col, cor_m))
}

# Untreated Malignant 内 (避开治疗影响)
cat("\n=== Within Untreated Malignant cells only ===\n")
ut_malig <- annot2[Level1 == "Epithelial (malignant)" & treatment_status == "Untreated"]
cat(sprintf("Untreated malignant: %d\n", nrow(ut_malig)))
for (col in hwang_cols) {
  cor_u <- cor(ut_malig[[nefrs_col]], ut_malig[[col]], method="spearman", use="pairwise.complete.obs")
  cat(sprintf("  %-40s : rho = %+.3f\n", col, cor_u))
}

# === 8. 保存相关性结果 ===
cor_results <- data.table()
for (col in hwang_cols) {
  cor_all <- cor(annot2[[nefrs_col]], annot2[[col]], method="spearman", use="pairwise.complete.obs")
  cor_m <- cor(malig[[nefrs_col]], malig[[col]], method="spearman", use="pairwise.complete.obs")
  cor_u <- cor(ut_malig[[nefrs_col]], ut_malig[[col]], method="spearman", use="pairwise.complete.obs")
  cor_results <- rbind(cor_results, data.table(
    program = col,
    rho_all = round(cor_all, 3),
    rho_malig = round(cor_m, 3),
    rho_unt_malig = round(cor_u, 3)
  ))
}
cor_results <- cor_results[order(-rho_malig)]
print(cor_results)
fwrite(cor_results, file.path(OUT, "w4_nefrs_vs_hwang17_correlation.tsv"), sep="\t")
cat("\nSaved: w4_nefrs_vs_hwang17_correlation.tsv\n")

cat("\n=== Done ===\n")
