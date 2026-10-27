# W4 Task 3: UCell NeFRS module score + 53 个基因签名
suppressPackageStartupMessages({
  library(UCell); library(Matrix); library(data.table)
})

OUT  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
FIG  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"

cat("=== W4 Task 3: UCell NeFRS module score ===\n")
t_total <- Sys.time()

# === 1. 加载数据 ===
cat("\nLoading X (sparse full matrix)...\n")
X_full <- readRDS(file.path(OUT, "w4_X_sparse_full.rds"))
cat(sprintf("Loaded: %d cells x %d genes\n", nrow(X_full), ncol(X_full)))

meta <- readRDS(file.path(OUT, "w4_metadata_only.rds"))
annot_df <- meta$annot

# === 2. UCell 要求基因 x 细胞 ( reversed) ===
# X_full 现在是 cells x genes, 需要转置
cat("\nTransposing matrix (cells x genes -> genes x cells)...\n")
t0 <- Sys.time()
X_t <- t(X_full)
cat(sprintf("Transposed in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units="mins"))))
cat(sprintf("Class: %s, dim: %d x %d\n", class(X_t), nrow(X_t), ncol(X_t)))

# === 3. NeFRS 53 基因签名 ===
NEFRS_GENES <- c("FGF13","MAP1S","REST","NRP1","CDC25A","KIF21A","TOR2A","LGMN",
                  "WWTR1","DTNBP1","DCX","EREG","ISCA2","FADS1","KRT18","BMPR1B",
                  "SOCS2","EGFR","CP","SLC25A46","MECP2","FBLN1","NQO1","ETHE1",
                  "EPHB6","VPS13A","NOTCH1","HES1","CX3CR1","SPATA2","SIRT2",
                  "NR2F2","FOXO4","VIPR1","MCF2L","STXBP1","BOC","MYCN","GABRB3",
                  "ATOH1","APLP1","FLRT3","WNT7A","ITGB1","IFRD1","EFNB2",
                  "SEMA3A","CAPRIN2","GALNT14","GJA1","HELLS","WNT7B","BRD7")

# 同时拆分 NeFRS 为正/负方向(根据 coef 符号)
# 读 v2 公式
coef_path <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime/nefrs_v3_enet_coefficients.tsv"
coef_df <- fread(coef_path)
nefrs_risk <- coef_df[coef > 0]$gene  # risk 方向 (HR>1)
nefrs_prot <- coef_df[coef < 0]$gene  # protective (HR<1)
cat(sprintf("\nNeFRS risk genes (β>0): %d\n", length(nefrs_risk)))
cat(sprintf("NeFRS protective genes (β<0): %d\n", length(nefrs_prot)))

# === 4. UCell module score (双签名: risk + protective) ===
signatures <- list(
  NeFRS_all = NEFRS_GENES,
  NeFRS_risk = nefrs_risk,
  NeFRS_protective = nefrs_prot
)

cat("\nRunning UCell ScoreSignatures_UCell...\n")
cat("(224k cells x 53 genes, expected 5-15 min)\n")
t0 <- Sys.time()
ucell_scores <- ScoreSignatures_UCell(matrix=X_t, features=signatures, ncores=8, chunk.size=5000)
cat(sprintf("UCell completed in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units="mins"))))
cat(sprintf("Output dim: %d x %d\n", nrow(ucell_scores), ncol(ucell_scores)))
cat(sprintf("Output columns: %s\n", paste(colnames(ucell_scores), collapse=", ")))

# === 5. 合并到 annot ===
cat("\nMerging UCell scores into annot...\n")
score_df <- as.data.frame(ucell_scores)
score_df$barcode <- rownames(score_df)
annot_df <- merge(annot_df, score_df, by="barcode", sort=FALSE, all.x=TRUE)

# === 6. 总览 ===
cat("\n=== UCell score distribution (all cells) ===\n")
for (s in colnames(ucell_scores)) {
  v <- ucell_scores[, s]
  cat(sprintf("  %s: median=%.4f, mean=%.4f, range=[%.4f, %.4f]\n",
              s, median(v, na.rm=TRUE), mean(v, na.rm=TRUE),
              min(v, na.rm=TRUE), max(v, na.rm=TRUE)))
}

cat("\n=== NeFRS_all score by Level 1 Annotation ===\n")
score_by_l1 <- annot_df[!is.na(annot_df$NeFRS_all_UCell), 
                         .(median_score = median(NeFRS_all_UCell, na.rm=TRUE),
                           n = .N),
                         by = `Level 1 Annotation`]
score_by_l1 <- score_by_l1[order(-median_score)]
print(score_by_l1)

cat("\n=== NeFRS_all score by Level 3 Annotation (top 20) ===\n")
score_by_l3 <- annot_df[!is.na(annot_df$NeFRS_all_UCell), 
                         .(median_score = median(NeFRS_all_UCell, na.rm=TRUE),
                           n = .N),
                         by = `Level 3 Annotation`]
score_by_l3 <- score_by_l3[order(-median_score)][1:20]
print(score_by_l3)

cat("\n=== NeFRS_all score by treatment_status ===\n")
score_by_treat <- annot_df[!is.na(annot_df$NeFRS_all_UCell), 
                            .(median_score = median(NeFRS_all_UCell, na.rm=TRUE),
                              n = .N),
                            by = treatment_status]
score_by_treat <- score_by_treat[order(-median_score)]
print(score_by_treat)

# === 7. 保存 ===
saveRDS(annot_df, file.path(OUT, "w4_annot_with_ucell.rds"))
fwrite(annot_df, file.path(OUT, "w4_annot_with_ucell.tsv"), sep="\t")
saveRDS(ucell_scores, file.path(OUT, "w4_ucell_scores_matrix.rds"))

cat(sprintf("\n=== Task 3 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t_total, units="mins"))))
cat("Files:\n")
cat(sprintf("  %s/w4_annot_with_ucell.rds (含 NeFRS_all/risk/protective UCell)\n", OUT))
cat(sprintf("  %s/w4_annot_with_ucell.tsv\n", OUT))
