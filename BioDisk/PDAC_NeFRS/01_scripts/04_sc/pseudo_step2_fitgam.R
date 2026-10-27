# 拟时序后半段: 全上皮轨迹(原作UMAP) + fitGAM + associationTest
# 降采样~20k保留稀有类型, 控制内存
suppressPackageStartupMessages({
  library(data.table); library(slingshot); library(SingleCellExperiment)
  library(tradeSeq); library(Matrix); library(BiocParallel)
})
SC  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/09_pseudotime"
MIME<- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime"

cat("=== fitGAM: 全上皮轨迹 (原作UMAP) ===\n"); t0 <- Sys.time()

X <- readRDS(file.path(SC,"w4_X_sparse_full.rds"))
annot <- as.data.table(readRDS(file.path(SC,"w4_annot_with_ucell.rds")))
meta <- readRDS(file.path(SC,"w4_metadata_only.rds"))
coef <- fread(file.path(MIME,"nefrs_v3_enet_coefficients.tsv"))

# 上皮 + 原作 UMAP
epi_types <- c("Acinar","ADM","Ductal","Ductal (atypical)","Malignant")
um <- data.table(barcode=meta$cell_barcodes, UMAP_1=meta$umap[,1], UMAP_2=meta$umap[,2])
um <- merge(um, annot[,.(barcode,Level2)], by="barcode")
epi <- um[Level2 %in% epi_types]
cat(sprintf("上皮: %d\n", nrow(epi)))

# === C 折中分层采样 ~3万 (每段充足, 恶性不过度主导) ===
# 目的: fitGAM 拟合轨迹各阶段曲线, 需每段足够细胞; 非丰度估计, 故不保持原比例
# 文档需注明此采样策略
set.seed(42)
caps <- c("Acinar"=99999, "ADM"=99999, "Ductal"=8000,
          "Ductal (atypical)"=99999, "Malignant"=8000)
samp_idx <- c()
for (ct in epi_types) {
  ci <- which(epi$Level2==ct)
  cap <- caps[[ct]]
  if (length(ci) > cap) ci <- sample(ci, cap)
  samp_idx <- c(samp_idx, ci)
}
epi_s <- epi[samp_idx]
cat(sprintf("分层采样后: %d (各类型 cap: Ductal/Malignant 8k, 其余全留)
", nrow(epi_s)))
print(table(epi_s$Level2))

# 表达矩阵 (原始 counts 给 fitGAM, 它要 counts)
# X 是 log-norm, fitGAM 理想要 raw counts; 但只有 log-norm, 用 expm1 近似还原
idx <- match(epi_s$barcode, rownames(X))
Xs <- X[idx,]                       # cells x genes, log-norm
cat(sprintf("Xs: %d x %d\n", nrow(Xs), ncol(Xs)))
rm(X); gc()

# slingshot on 原作 UMAP
sce <- SingleCellExperiment(assays=list(logcounts=t(as.matrix(Xs))))
reducedDim(sce,"UMAP") <- as.matrix(epi_s[,.(UMAP_1,UMAP_2)])
colData(sce)$ct <- factor(epi_s$Level2)
sce <- slingshot(sce, clusterLabels="ct", reducedDim="UMAP",
                 start.clus="Acinar", end.clus="Malignant")
cat("Lineages:\n"); print(slingLineages(SlingshotDataSet(sce)))
gc()

# === fitGAM (counts = expm1 还原, 控制 workers=2 防撑爆) ===
cat("\n[fitGAM] 拟合 (genes x cells, workers=2)...\n")
counts_mat <- as.matrix(expm1(t(as.matrix(Xs))))  # genes x cells, 近似 counts
# 只对表达基因拟合 (减负担): 检出>5% 的基因
det <- rowMeans(counts_mat > 0)
keep_g <- which(det > 0.05)
cat(sprintf("拟合基因: %d (检出>5%%)\n", length(keep_g)))
counts_fit <- counts_mat[keep_g,]
rm(counts_mat); gc()

pt_mat <- slingPseudotime(sce)
w_mat <- slingCurveWeights(sce)

# 过滤: 至少在一条 lineage 有非 NA pseudotime 的细胞
# slingPseudotime 用 na=FALSE 把 NA 填成该 lineage 内值, 但更稳的是过滤全 NA 细胞
valid_cell <- rowSums(!is.na(pt_mat)) > 0
cat(sprintf("有效细胞 (非全NA pseudotime): %d / %d\n", sum(valid_cell), length(valid_cell)))
# slingshot 的 weights 已处理分配; fitGAM 要求 pt 无 NA -> 用 slingPseudotime(na=FALSE)
pt_mat <- slingPseudotime(sce, na=FALSE)
w_mat <- slingCurveWeights(sce)
counts_fit <- counts_fit[, valid_cell]
pt_mat <- pt_mat[valid_cell,, drop=FALSE]
w_mat <- w_mat[valid_cell,, drop=FALSE]
cat(sprintf("fitGAM 输入: %d genes x %d cells\n", nrow(counts_fit), ncol(counts_fit)))

set.seed(42)
t_gam <- Sys.time()
sce_gam <- fitGAM(counts=counts_fit, pseudotime=pt_mat, cellWeights=w_mat,
                  nknots=6, parallel=TRUE, BPPARAM=MulticoreParam(workers=4))
cat(sprintf("  fitGAM 完成 (%.1f min)\n", as.numeric(difftime(Sys.time(),t_gam,units="mins"))))
gc()

# === associationTest: 找轨迹驱动基因 ===
cat("\n[associationTest] 全基因扫...\n")
assoc <- associationTest(sce_gam)
assoc$gene <- rownames(assoc)
assoc <- as.data.table(assoc)[order(-waldStat)]
assoc[, padj := p.adjust(pvalue, "BH")]
cat(sprintf("显著随轨迹变化基因 (padj<0.05): %d\n", sum(assoc$padj<0.05, na.rm=TRUE)))
cat("\nTop 20 轨迹驱动基因:\n")
print(head(assoc[, .(gene, waldStat, padj)], 20))

# NeFRS 成员在其中的排名
assoc[, is_nefrs := gene %in% coef$gene]
cat(sprintf("\nNeFRS 成员中显著的: %d / %d\n",
            sum(assoc$is_nefrs & assoc$padj<0.05, na.rm=TRUE), sum(assoc$is_nefrs)))
cat("NeFRS 成员的轨迹关联 (top):\n")
print(head(assoc[is_nefrs==TRUE, .(gene, waldStat, padj)], 15))

# === 保存 ===
saveRDS(sce_gam, file.path(OUT,"pseudo_fitgam_sce.rds"))
fwrite(assoc, file.path(OUT,"pseudo_associationTest.tsv"), sep="\t")
saveRDS(list(sce=sce, epi_s=epi_s), file.path(OUT,"pseudo_allepi_origUMAP.rds"))
cat(sprintf("\n=== fitGAM 完成, 总 %.1f min ===\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
gc()
