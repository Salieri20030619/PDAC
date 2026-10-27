# 严肃诊断: 验证我之前的几个断言对不对
suppressPackageStartupMessages({
  library(data.table); library(Matrix)
})

OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
annot <- readRDS(file.path(OUT, "w4_annot_with_ucell.rds"))
X <- readRDS(file.path(OUT, "w4_X_sparse_full.rds"))

cat("=== 严肃诊断 1: NeFRS 53 基因, 哪几个在 scRNA 里有信号 ===\n")
NEFRS <- c("FGF13","MAP1S","REST","NRP1","CDC25A","KIF21A","TOR2A","LGMN",
            "WWTR1","DTNBP1","DCX","EREG","ISCA2","FADS1","KRT18","BMPR1B",
            "SOCS2","EGFR","CP","SLC25A46","MECP2","FBLN1","NQO1","ETHE1",
            "EPHB6","VPS13A","NOTCH1","HES1","CX3CR1","SPATA2","SIRT2",
            "NR2F2","FOXO4","VIPR1","MCF2L","STXBP1","BOC","MYCN","GABRB3",
            "ATOH1","APLP1","FLRT3","WNT7A","ITGB1","IFRD1","EFNB2",
            "SEMA3A","CAPRIN2","GALNT14","GJA1","HELLS","WNT7B","BRD7")

# 每个基因在所有 224k 细胞中的 % nonzero
nefrs_idx <- match(NEFRS, colnames(X))
pct_nonzero <- numeric(53)
mean_expr <- numeric(53)
for (i in seq_along(NEFRS)) {
  if (is.na(nefrs_idx[i])) {
    pct_nonzero[i] <- NA; mean_expr[i] <- NA
  } else {
    v <- X[, nefrs_idx[i]]
    pct_nonzero[i] <- 100 * sum(v > 0) / length(v)
    mean_expr[i] <- mean(v)
  }
}

# 读 coef 看 risk/protective
coef_df <- fread("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime/nefrs_v3_enet_coefficients.tsv")
diag_df <- data.table(gene = NEFRS, pct_nonzero = round(pct_nonzero, 2),
                       mean_expr = round(mean_expr, 4))
diag_df <- merge(diag_df, coef_df[, .(gene, coef, priority)], by="gene", all.x=TRUE)
diag_df <- diag_df[order(-pct_nonzero)]

cat("\n=== NeFRS 53 基因 in scRNA: 表达率排序 ===\n")
print(diag_df)

cat(sprintf("\n>50%% expressed: %d / 53\n", sum(diag_df$pct_nonzero > 50, na.rm=TRUE)))
cat(sprintf(">10%% expressed: %d / 53\n", sum(diag_df$pct_nonzero > 10, na.rm=TRUE)))
cat(sprintf("<1%%  expressed: %d / 53 (基本 dropout)\n",
            sum(diag_df$pct_nonzero < 1, na.rm=TRUE)))

# 哪些是 risk hub 但 dropout 严重
hub_dropout <- diag_df[priority %in% c("P1_TopPriority","P2_HighPriority") & pct_nonzero < 5]
cat("\n=== Hub 基因 dropout 严重的 ===\n")
print(hub_dropout)

cat("\n=== 严肃诊断 2: NeFRS 高分到底是哪些细胞? ===\n")
# 全数据集 Top 10% NeFRS, 看 Level3 分布跟 baseline 比 (FET)
overall_n <- annot[, .N, by=Level3]
overall_n[, baseline_pct := round(100*N/sum(N), 2)]

cutoff_top10 <- quantile(annot$NeFRS_all, 0.90, na.rm=TRUE)
top_n <- annot[NeFRS_all > cutoff_top10, .N, by=Level3]
top_n[, top10_pct := round(100*N/sum(N), 2)]

enrich <- merge(overall_n[, .(Level3, baseline_pct)],
                 top_n[, .(Level3, top10_pct)], by="Level3", all=TRUE)
enrich[is.na(top10_pct), top10_pct := 0]
enrich[, fold_enrichment := round(top10_pct / baseline_pct, 2)]
cat("\n=== Top 10% NeFRS vs 全 dataset (fold enrichment) ===\n")
print(enrich[order(-fold_enrichment)][1:20])

cat("\n=== 严肃诊断 3: 仅 Untreated 内, Ductal atypical 是否优势? ===\n")
# 之前看的是全 dataset top 10%, 但全 dataset 已经混合 treated/untreated
# 严格做法: 只在 Untreated 内分位数
unt <- annot[treatment_status == "Untreated" & !is.na(NeFRS_all)]
cat(sprintf("Untreated cells: %d\n", nrow(unt)))

q90 <- quantile(unt$NeFRS_all, 0.90, na.rm=TRUE)
unt_top <- unt[NeFRS_all > q90]
unt_top_pct <- unt_top[, .N, by=Level3]
unt_top_pct[, pct := round(100*N/sum(N), 2)]
unt_all_pct <- unt[, .N, by=Level3]
unt_all_pct[, baseline_pct := round(100*N/sum(N), 2)]
unt_enrich <- merge(unt_all_pct[, .(Level3, baseline_pct)],
                     unt_top_pct[, .(Level3, top_pct=pct)], by="Level3", all=TRUE)
unt_enrich[is.na(top_pct), top_pct := 0]
unt_enrich[, fold := round(top_pct / baseline_pct, 2)]
cat("\n=== Untreated 内 Top 10% NeFRS Level3 enrichment ===\n")
print(unt_enrich[order(-fold)][1:15])

cat("\n=== 严肃诊断 4: CRT 治疗类型分布 ===\n")
# 看 treatment_status 跟 pid (patient) / sampleid 的关系
# 同一个 patient 可能有多个 treatment 阶段?
t_sample <- annot[, .(n_cells = .N, n_samples = uniqueN(sampleid)), by=treatment_status]
print(t_sample)

cat("\n=== treatment by sampleid count ===\n")
ts_table <- annot[, .(treatment_status = unique(treatment_status)), by=sampleid]
print(table(ts_table$treatment_status))

cat("\n=== 严肃诊断 5: Malignant 内 NeFRS 分布 (避开 cell type bias) ===\n")
# 只看 Malignant 细胞内不同 treatment 的 NeFRS
malig <- annot[Level1 == "Epithelial (malignant)" & !is.na(NeFRS_all)]
cat(sprintf("Malignant cells: %d\n", nrow(malig)))
malig_treat <- malig[, .(median = median(NeFRS_all, na.rm=TRUE),
                          mean = mean(NeFRS_all, na.rm=TRUE),
                          q25 = quantile(NeFRS_all, 0.25, na.rm=TRUE),
                          q75 = quantile(NeFRS_all, 0.75, na.rm=TRUE),
                          n = .N), by=treatment_status][order(-median)]
print(malig_treat)

cat("\n=== Done ===\n")
