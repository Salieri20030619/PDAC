# 用真实 counts 跑 NRP1 KO, 对比之前伪 counts 的扰动数 (81)
suppressPackageStartupMessages({library(scTenifoldKnk); library(Matrix); library(data.table)})
options(bitmapType="cairo")
KO <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/10_knockout"

cat("=== 真实 counts NRP1 KO 测试 ===\n")
t0 <- Sys.time()
counts <- readRDS(file.path(KO,"malignant_real_counts.rds"))  # cells x genes, 真 counts
cat(sprintf("恶性 counts: %d x %d\n", nrow(counts), ncol(counts)))

gene <- "NRP1"
# HVG: 真 counts 上用方差选 (跟之前同口径但用真 counts)
gene_var <- apply(counts, 2, var)
hvg <- names(sort(gene_var, decreasing=TRUE))[1:3000]
genes_use <- unique(c(gene, hvg)); genes_use <- genes_use[genes_use %in% colnames(counts)]
cat(sprintf("基因数: %d\n", length(genes_use)))

# genes x cells, NRP1 第一行, 真 counts (不用 expm1!)
cm <- as.matrix(t(counts[, genes_use]))
if (which(rownames(cm)==gene)!=1) cm <- cm[c(gene,setdiff(rownames(cm),gene)),]
cat(sprintf("矩阵: %d x %d, 目标基因位置 %d, 最大值 %.0f\n",
            nrow(cm), ncol(cm), which(rownames(cm)==gene), max(cm)))
rm(counts); gc()

cat("\n跑 scTenifoldKnk (真 counts)...\n")
res <- tryCatch(scTenifoldKnk(countMatrix=cm, gKO=gene, qc=FALSE),
                error=function(e){cat("FAIL:",conditionMessage(e),"\n"); NULL})

if (!is.null(res)) {
  dr <- res$diffRegulation
  cat(sprintf("\n=== 真实 counts NRP1 KO 结果 ===\n"))
  cat(sprintf("耗时 %.1f min\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
  cat(sprintf("P<0.05: %d, adjP<0.05: %d\n", sum(dr$p.value<0.05), sum(dr$p.adj<0.05)))
  cat(sprintf("对比: 伪counts 时 P<0.05=114, adjP<0.05=81\n\n"))
  cat("Top 10 扰动 (by FC):\n")
  print(head(dr[order(-dr$FC),][dr$gene!=gene,][,c("gene","Z","FC","p.adj")],10))
  saveRDS(res, file.path(KO,"Results/TEST_NRP1_realcounts.rds"))
}
