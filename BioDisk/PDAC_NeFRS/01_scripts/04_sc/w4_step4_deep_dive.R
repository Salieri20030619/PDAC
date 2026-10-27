# W4 Task 4: 关键假设验证 - Ductal atypical + treatment 细致分析
suppressPackageStartupMessages({ library(data.table) })

OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
annot <- readRDS(file.path(OUT, "w4_annot_with_ucell.rds"))

cat("=== Task 4: Deep dive ===\n")

# === 1. Ductal atypical 在不同 treatment status 中的分布 ===
cat("\n=== Ductal (atypical) by treatment status ===\n")
dap <- annot[Level3 == "Ductal (atypical)"]
print(dap[, .(median_NeFRS = median(NeFRS_all, na.rm=TRUE),
              median_risk  = median(NeFRS_risk, na.rm=TRUE),
              n = .N),
          by = treatment_status][order(-median_NeFRS)])

# === 2. Malignant 在不同 treatment 中, 看不同细胞数量分布 ===
cat("\n=== Malignant cells distribution by treatment ===\n")
malig <- annot[Level1 == "Epithelial (malignant)"]
print(malig[, .N, by=treatment_status][order(-N)])

# === 3. Ductal atypical 比例: CRT 后是否富集? (Hwang 核心论点) ===
cat("\n=== Cell composition ratios by treatment ===\n")
comp <- annot[Level3 %in% c("Ductal (atypical)", "Malignant", "Ductal", "Schwann", "myCAF")]
comp_pct <- comp[, .N, by=.(treatment_status, Level3)]
comp_total <- annot[, .(total = .N), by=treatment_status]
comp_pct <- merge(comp_pct, comp_total, by="treatment_status")
comp_pct[, pct := round(100*N/total, 2)]
comp_pct <- comp_pct[order(treatment_status, -pct)]
print(comp_pct)

# === 4. 重要: 看 Untreated 患者内, NeFRS 高 vs 低的细胞类型分布 ===
cat("\n=== Untreated cells, NeFRS_all top 10% vs bottom 10% ===\n")
unt <- annot[treatment_status == "Untreated" & !is.na(NeFRS_all)]
cutoff_top <- quantile(unt$NeFRS_all, 0.90, na.rm=TRUE)
cutoff_bot <- quantile(unt$NeFRS_all, 0.10, na.rm=TRUE)
cat(sprintf("Untreated 90%% cutoff: %.4f\n", cutoff_top))
cat(sprintf("Untreated 10%% cutoff: %.4f\n", cutoff_bot))

top_cells <- unt[NeFRS_all > cutoff_top]
bot_cells <- unt[NeFRS_all < cutoff_bot]

cat("\n=== Untreated Top 10% NeFRS - cell composition (Level3) ===\n")
top_comp <- top_cells[, .(n=.N), by=Level3][order(-n)][1:15]
top_comp[, pct := round(100*n/nrow(top_cells), 2)]
print(top_comp)

cat("\n=== Untreated Bot 10% NeFRS - cell composition (Level3) ===\n")
bot_comp <- bot_cells[, .(n=.N), by=Level3][order(-n)][1:15]
bot_comp[, pct := round(100*n/nrow(bot_cells), 2)]
print(bot_comp)

# === 5. CRT 治疗后的 Ductal atypical 是否仍高于 Untreated Ductal atypical ===
cat("\n=== Ductal (atypical) Untreated vs CRT NeFRS comparison ===\n")
dap_uvt <- annot[Level3 == "Ductal (atypical)" & treatment_status %in% c("Untreated","CRT","CRTl","CRTx")]
print(dap_uvt[, .(median = median(NeFRS_all, na.rm=TRUE),
                   mean = mean(NeFRS_all, na.rm=TRUE),
                   q75 = quantile(NeFRS_all, 0.75, na.rm=TRUE),
                   n = .N), by=treatment_status][order(-median)])

# === 6. 关键 hub 基因(NRP1, SEMA3A 等) 单基因表达 by cell type ===
# 这里需要表达矩阵
cat("\n=== Loading X for hub gene check ===\n")
X <- readRDS(file.path(OUT, "w4_X_sparse_full.rds"))
gene_names <- colnames(X)
hub_genes <- c("NRP1", "SEMA3A", "REST", "EGFR", "EPHB6", "MYCN")
hub_in <- intersect(hub_genes, gene_names)
cat(sprintf("Hub genes available: %s\n", paste(hub_in, collapse=", ")))

# 取 hub 基因 sparse 列
hub_idx <- match(hub_in, gene_names)
hub_mat <- as.matrix(X[, hub_idx, drop=FALSE])
cat(sprintf("Hub matrix: %d x %d\n", nrow(hub_mat), ncol(hub_mat)))

# 加到 annot
hub_df <- data.table(barcode = rownames(hub_mat))
for (g in hub_in) {
  hub_df[[g]] <- hub_mat[, g]
}
annot_hub <- merge(annot, hub_df, by="barcode", sort=FALSE)

# 看 hub 基因在 Top cell types 的表达 (% non-zero, mean expression)
cat("\n=== Hub gene expression % nonzero by Level 3 ===\n")
for (g in hub_in) {
  cat(sprintf("\n--- %s ---\n", g))
  hub_summary <- annot_hub[, .(pct_nonzero = round(100*sum(get(g) > 0, na.rm=TRUE)/.N, 2),
                                mean_expr = round(mean(get(g), na.rm=TRUE), 4),
                                n = .N),
                            by = Level3][order(-pct_nonzero)][1:10]
  print(hub_summary)
}

cat("\n=== Done ===\n")
