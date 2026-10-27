# fitGAM 出图: NeFRS 基因沿轨迹平滑曲线 + 轨迹驱动基因
suppressPackageStartupMessages({
  library(data.table); library(tradeSeq); library(ggplot2); library(SingleCellExperiment)
})
OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/09_pseudotime"
FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures/Fig8_pseudotime"
MIME<- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime"

cat("=== fitGAM 出图 ===\n")
sce_gam <- readRDS(file.path(OUT,"pseudo_fitgam_sce.rds"))
assoc <- fread(file.path(OUT,"pseudo_associationTest.tsv"))
coef <- fread(file.path(MIME,"nefrs_v3_enet_coefficients.tsv"))

# NeFRS 成员在轨迹驱动中的情况
assoc[, is_nefrs := gene %in% coef$gene]
nefrs_assoc <- assoc[is_nefrs==TRUE][order(-waldStat)]
cat(sprintf("NeFRS 成员 %d 个, 沿轨迹显著 (padj<0.05): %d\n",
            nrow(nefrs_assoc), sum(nefrs_assoc$padj<0.05, na.rm=TRUE)))
cat("\nNeFRS 成员轨迹关联 top 15:\n")
print(head(nefrs_assoc[, .(gene, waldStat, meanLogFC, padj)], 15))

# === Fig8E: NeFRS top hub 基因沿 pseudotime 平滑曲线 ===
# 选 fitGAM 拟合过的 NeFRS 基因里 waldStat 最高的几个
hub_candidates <- c("NRP1","EGFR","BMPR1B","KRT18","CDC25A","WNT7B")
fitted_genes <- rownames(sce_gam)
hubs <- intersect(hub_candidates, fitted_genes)
hubs <- intersect(hubs, nefrs_assoc$gene)  # 确保在 assoc 里
cat(sprintf("\n出图 hub 基因: %s\n", paste(hubs, collapse=", ")))

plots <- list()
for (g in hubs) {
  p <- tryCatch(
    plotSmoothers(sce_gam, counts=NULL, gene=g) +
      ggtitle(g) + theme_bw(base_size=10) +
      theme(legend.position="none", plot.title=element_text(size=11, face="bold")),
    error=function(e){cat(sprintf("  %s 出图失败: %s\n", g, conditionMessage(e))); NULL})
  if (!is.null(p)) plots[[g]] <- p
}

if (length(plots)>0) {
  library(patchwork)
  combined <- wrap_plots(plots, ncol=3) +
    plot_annotation(title="Fig 8E. NeFRS hub genes along Acinar->Malignant trajectory (fitGAM)",
                    subtitle="tradeSeq smoothed expression; 29581 epithelial cells (stratified sampling)")
  ggsave(file.path(FIG,"Fig8E_nefrs_genes_smoothers.pdf"), combined, width=12, height=7)
  cat("Fig8E 保存\n")
}

# === Fig8F: 轨迹驱动基因 top 30 (waldStat) 表 + NeFRS 标记 ===
top30 <- head(assoc[order(-waldStat)], 30)
top30[, label := ifelse(is_nefrs, paste0(gene,"*"), gene)]
cat("\nTop 30 轨迹驱动基因 (* = NeFRS 成员):\n")
print(top30[, .(gene, waldStat, meanLogFC, is_nefrs)])

fwrite(nefrs_assoc, file.path(OUT,"nefrs_members_trajectory_assoc.tsv"), sep="\t")
cat("\n存 nefrs_members_trajectory_assoc.tsv\n")
