# W4 Task 3b: 接续 — 把 annot_df 转 data.table 后做分组分析
suppressPackageStartupMessages({
  library(data.table); library(Matrix)
})

OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"

cat("=== W4 Task 3b: Group analyses ===\n")

# 读 UCell scores + metadata
ucell_scores <- readRDS(file.path(OUT, "w4_ucell_scores_matrix.rds"))
meta <- readRDS(file.path(OUT, "w4_metadata_only.rds"))

cat(sprintf("UCell scores: %d x %d\n", nrow(ucell_scores), ncol(ucell_scores)))

# 合并 (UCell scores 的 rownames = barcode)
score_df <- data.table(barcode = rownames(ucell_scores),
                        NeFRS_all = ucell_scores[, "NeFRS_all_UCell"],
                        NeFRS_risk = ucell_scores[, "NeFRS_risk_UCell"],
                        NeFRS_protective = ucell_scores[, "NeFRS_protective_UCell"])

annot <- as.data.table(meta$annot)
cat(sprintf("annot dim: %d x %d\n", nrow(annot), ncol(annot)))
annot <- merge(annot, score_df, by="barcode", sort=FALSE, all.x=TRUE)
cat(sprintf("After merge: %d x %d\n", nrow(annot), ncol(annot)))

# 重命名 Level 列 (空格转下划线方便后续)
setnames(annot, "Level 1 Annotation", "Level1")
setnames(annot, "Level 2 Annotation", "Level2")
setnames(annot, "Level 3 Annotation", "Level3")

# === 总览 ===
cat("\n=== UCell score distribution (all cells) ===\n")
for (s in c("NeFRS_all","NeFRS_risk","NeFRS_protective")) {
  v <- annot[[s]]
  cat(sprintf("  %s: median=%.4f, mean=%.4f, IQR=[%.4f, %.4f]\n",
              s, median(v, na.rm=TRUE), mean(v, na.rm=TRUE),
              quantile(v, 0.25, na.rm=TRUE), quantile(v, 0.75, na.rm=TRUE)))
}

# === 1. NeFRS by Level 1 Annotation ===
cat("\n=== NeFRS_all by Level 1 (broad celltype) ===\n")
score_l1 <- annot[!is.na(NeFRS_all),
                   .(median_score = median(NeFRS_all, na.rm=TRUE),
                     mean_score = mean(NeFRS_all, na.rm=TRUE),
                     n = .N),
                   by = Level1][order(-median_score)]
print(score_l1)

# === 2. NeFRS by Level 3 (细粒度) ===
cat("\n=== NeFRS_all by Level 3 (细粒度) — Top 20 ===\n")
score_l3 <- annot[!is.na(NeFRS_all),
                   .(median_score = median(NeFRS_all, na.rm=TRUE),
                     mean_score = mean(NeFRS_all, na.rm=TRUE),
                     n = .N),
                   by = Level3][order(-median_score)]
print(head(score_l3, 20))

# === 3. ⭐ NeFRS_risk (NeFRS 的 risk 方向) by Level 3 ===
cat("\n=== NeFRS_risk (β>0 27 基因) by Level 3 — Top 20 ===\n")
score_l3_risk <- annot[!is.na(NeFRS_risk),
                        .(median_risk = median(NeFRS_risk, na.rm=TRUE),
                          median_prot = median(NeFRS_protective, na.rm=TRUE),
                          n = .N),
                        by = Level3][order(-median_risk)]
print(head(score_l3_risk, 20))

# === 4. ⭐⭐⭐ Treatment status comparison (Hwang 论文核心) ===
cat("\n=== NeFRS by treatment_status (Hwang 论文核心) ===\n")
score_treat <- annot[!is.na(NeFRS_all),
                      .(median_all = median(NeFRS_all, na.rm=TRUE),
                        median_risk = median(NeFRS_risk, na.rm=TRUE),
                        median_prot = median(NeFRS_protective, na.rm=TRUE),
                        mean_all = mean(NeFRS_all, na.rm=TRUE),
                        n = .N),
                      by = treatment_status][order(-median_all)]
print(score_treat)

# === 5. ⭐⭐⭐ Untreated vs CRT — 直接对比 ===
cat("\n=== Direct comparison: Untreated vs CRT (in malignant cells only) ===\n")
malig_sub <- annot[Level1 == "Epithelial (malignant)" & 
                    treatment_status %in% c("Untreated","CRT") & 
                    !is.na(NeFRS_all)]
cat(sprintf("Malignant cells with Untreated or CRT: %d\n", nrow(malig_sub)))
cat("Untreated malignant n:", sum(malig_sub$treatment_status=="Untreated"), "\n")
cat("CRT malignant n:", sum(malig_sub$treatment_status=="CRT"), "\n")

if (nrow(malig_sub) > 100) {
  # Wilcoxon test
  wt <- wilcox.test(malig_sub[treatment_status=="Untreated"]$NeFRS_all,
                     malig_sub[treatment_status=="CRT"]$NeFRS_all)
  cat(sprintf("\nWilcoxon NeFRS_all: W = %.0f, p = %.2e\n", wt$statistic, wt$p.value))
  
  cat("\nMalignant NeFRS_all distribution:\n")
  print(malig_sub[, .(median = median(NeFRS_all, na.rm=TRUE),
                       mean = mean(NeFRS_all, na.rm=TRUE),
                       q75 = quantile(NeFRS_all, 0.75, na.rm=TRUE),
                       n = .N), by=treatment_status])
  
  # 同样测 NeFRS_risk
  wt_risk <- wilcox.test(malig_sub[treatment_status=="Untreated"]$NeFRS_risk,
                          malig_sub[treatment_status=="CRT"]$NeFRS_risk)
  cat(sprintf("\nWilcoxon NeFRS_risk: W = %.0f, p = %.2e\n", wt_risk$statistic, wt_risk$p.value))
  
  cat("\nMalignant NeFRS_risk distribution:\n")
  print(malig_sub[, .(median = median(NeFRS_risk, na.rm=TRUE),
                       mean = mean(NeFRS_risk, na.rm=TRUE),
                       q75 = quantile(NeFRS_risk, 0.75, na.rm=TRUE),
                       n = .N), by=treatment_status])
}

# === 6. NeFRS 在 Schwann 细胞中 (神经直接相关) ===
cat("\n=== NeFRS in Schwann + neurons ===\n")
nerve_cells <- annot[Level3 %in% c("Schwann", "Intra-pancreatic neurons") &
                      !is.na(NeFRS_all)]
cat(sprintf("Nerve-related cells: %d\n", nrow(nerve_cells)))
print(nerve_cells[, .(median_all = median(NeFRS_all, na.rm=TRUE),
                       median_risk = median(NeFRS_risk, na.rm=TRUE),
                       n = .N), by=Level3])

# === 保存 ===
fwrite(annot, file.path(OUT, "w4_annot_with_ucell.tsv"), sep="\t")
saveRDS(annot, file.path(OUT, "w4_annot_with_ucell.rds"))

# 保存关键统计结果
saveRDS(list(
  by_level1 = score_l1,
  by_level3 = score_l3,
  by_level3_risk = score_l3_risk,
  by_treatment = score_treat,
  malig_comparison = if(exists("wt")) list(wilcox_all=wt, wilcox_risk=wt_risk) else NULL
), file.path(OUT, "w4_score_summaries.rds"))

cat("\n=== Done ===\n")
cat("Files saved:\n")
cat(sprintf("  %s/w4_annot_with_ucell.rds (annot + 3 NeFRS scores per cell)\n", OUT))
cat(sprintf("  %s/w4_score_summaries.rds (5 分组统计)\n", OUT))
