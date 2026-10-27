# W4 Task 3 (拆分 - 只算 UCell + 保存)
suppressPackageStartupMessages({
  library(UCell); library(Matrix); library(data.table)
})

OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"

cat("=== UCell only — 计算 + 保存 ===\n")
t_total <- Sys.time()

cat("Loading X_sparse_full...\n")
X_full <- readRDS(file.path(OUT, "w4_X_sparse_full.rds"))
X_t <- t(X_full)
cat(sprintf("Matrix: %d x %d (genes x cells)\n", nrow(X_t), ncol(X_t)))

# NeFRS 53 + risk/protective 拆分
NEFRS_GENES <- c("FGF13","MAP1S","REST","NRP1","CDC25A","KIF21A","TOR2A","LGMN",
                  "WWTR1","DTNBP1","DCX","EREG","ISCA2","FADS1","KRT18","BMPR1B",
                  "SOCS2","EGFR","CP","SLC25A46","MECP2","FBLN1","NQO1","ETHE1",
                  "EPHB6","VPS13A","NOTCH1","HES1","CX3CR1","SPATA2","SIRT2",
                  "NR2F2","FOXO4","VIPR1","MCF2L","STXBP1","BOC","MYCN","GABRB3",
                  "ATOH1","APLP1","FLRT3","WNT7A","ITGB1","IFRD1","EFNB2",
                  "SEMA3A","CAPRIN2","GALNT14","GJA1","HELLS","WNT7B","BRD7")
coef_df <- fread("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime/nefrs_v3_enet_coefficients.tsv")
signatures <- list(
  NeFRS_all = NEFRS_GENES,
  NeFRS_risk = coef_df[coef > 0]$gene,
  NeFRS_protective = coef_df[coef < 0]$gene
)

cat("Running UCell...\n")
ucell_scores <- ScoreSignatures_UCell(matrix=X_t, features=signatures,
                                       ncores=8, chunk.size=5000)
cat(sprintf("UCell done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t_total, units="mins"))))

# 立刻保存!
saveRDS(ucell_scores, file.path(OUT, "w4_ucell_scores_matrix.rds"))
cat(sprintf("Saved: %s/w4_ucell_scores_matrix.rds\n", OUT))
cat(sprintf("Dim: %d x %d\n", nrow(ucell_scores), ncol(ucell_scores)))
cat(sprintf("Columns: %s\n", paste(colnames(ucell_scores), collapse=", ")))

cat(sprintf("\nDone in %.1f min\n",
            as.numeric(difftime(Sys.time(), t_total, units="mins"))))
