# scTenifoldKnk 小规模测试: NRP1, 恶性细胞, 确认跑通
suppressPackageStartupMessages({
  library(Seurat); library(scTenifoldKnk); library(data.table); library(Matrix)
})
options(bitmapType="cairo")
SC  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"

cat("=== scTenifoldKnk 测试 (NRP1, 恶性细胞) ===\n")
t0 <- Sys.time()

X <- readRDS(file.path(SC,"w4_X_sparse_full.rds"))   # cells x genes, log-norm
annot <- as.data.table(readRDS(file.path(SC,"w4_annot_with_ucell.rds")))

# 恶性细胞
mal <- annot[Level2=="Malignant"]$barcode
idx <- match(mal, rownames(X)); idx <- idx[!is.na(idx)]
Xm <- X[idx,]
cat(sprintf("恶性细胞: %d, 基因: %d\n", nrow(Xm), ncol(Xm)))
rm(X); gc()

# scTenifoldKnk 要 counts (genes x cells); X 是 log-norm, expm1 近似还原
# 先测试: 用 HVG 5000 + NRP1 降维 (全 2 万基因太大)
gene <- "NRP1"
cat(sprintf("KO 基因: %s, 在数据中: %s\n", gene, gene %in% colnames(Xm)))

# HVG: 用方差选 top 3000 (测试用小一点)
gene_var <- apply(Xm, 2, var)
hvg <- names(sort(gene_var, decreasing=TRUE))[1:3000]
genes_use <- unique(c(gene, hvg))
genes_use <- genes_use[genes_use %in% colnames(Xm)]
cat(sprintf("用基因数 (HVG3000+target): %d\n", length(genes_use)))

# counts matrix: genes x cells, NRP1 放第一行
countMatrix <- as.matrix(t(expm1(Xm[, genes_use])))  # genes x cells
countMatrix <- round(countMatrix)  # scTenifoldKnk 要整数 counts
if (which(rownames(countMatrix)==gene) != 1) {
  countMatrix <- countMatrix[c(gene, setdiff(rownames(countMatrix),gene)), ]
}
cat(sprintf("矩阵: %d genes x %d cells, 目标基因位置: %d\n",
            nrow(countMatrix), ncol(countMatrix), which(rownames(countMatrix)==gene)))
rm(Xm); gc()

# 跑 scTenifoldKnk
cat("\n开始 scTenifoldKnk (测试, 可能 20-40 min)...\n")
res <- tryCatch({
  scTenifoldKnk(countMatrix=countMatrix, gKO=gene, qc=FALSE)
}, error=function(e){cat("FAIL:", conditionMessage(e),"\n"); NULL})

if (!is.null(res)) {
  cat(sprintf("\n✅ 成功! 耗时 %.1f min\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
  cat("diffRegulation 前 10:\n")
  print(head(res$diffRegulation, 10))
  saveRDS(res, "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/10_knockout/Results/TEST_NRP1.rds")
  cat("测试结果存了\n")
} else {
  cat("\n❌ 测试失败, 看错误\n")
}
