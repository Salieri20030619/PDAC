# B: 恶性细胞内 NeFRS 高/低对比 (回避轨迹, 直接对比)
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(Matrix); library(ggplot2)
})
SC  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/09_pseudotime"
FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures/Fig8_pseudotime"

cat("=== B: 恶性细胞 NeFRS 高低对比 ===\n")
t0 <- Sys.time()

X <- readRDS(file.path(SC, "w4_X_sparse_full.rds"))
annot <- as.data.table(readRDS(file.path(SC, "w4_annot_with_ucell.rds")))

mal_cells <- annot[Level2=="Malignant"]$barcode
idx <- match(mal_cells, rownames(X)); idx <- idx[!is.na(idx)]
X_m <- X[idx, ]
annot_m <- annot[match(rownames(X_m), barcode)]
cat(sprintf("恶性细胞: %d\n", nrow(X_m)))
rm(X); gc()

# NeFRS 高/低三分 (上1/3 vs 下1/3)
nefrs <- annot_m$NeFRS_all
q33 <- quantile(nefrs, 1/3); q67 <- quantile(nefrs, 2/3)
grp <- ifelse(nefrs >= q67, "High", ifelse(nefrs <= q33, "Low", "Mid"))
cat(sprintf("NeFRS 分组: High %d / Mid %d / Low %d\n",
            sum(grp=="High"), sum(grp=="Mid"), sum(grp=="Low")))
cat(sprintf("High NeFRS median: %.3f, Low: %.3f\n",
            median(nefrs[grp=="High"]), median(nefrs[grp=="Low"])))

# 建 Seurat (修正版)
expr_t <- t(X_m)
assay_obj <- CreateAssayObject(counts=expr_t)
seu <- CreateSeuratObject(counts=assay_obj)
seu[["RNA"]]$data <- expr_t
seu$nefrs_grp <- grp
seu$pid <- annot_m$pid
rm(X_m, expr_t, assay_obj); gc()

# 差异表达 High vs Low
cat("\n[DE] High vs Low NeFRS 恶性细胞...\n")
Idents(seu) <- "nefrs_grp"
seu_hl <- subset(seu, idents=c("High","Low"))
de <- FindMarkers(seu_hl, ident.1="High", ident.2="Low",
                  logfc.threshold=0.25, min.pct=0.1, verbose=FALSE)
de$gene <- rownames(de)
de <- as.data.table(de)[order(-avg_log2FC)]
cat(sprintf("DE 基因 (p_adj<0.05): %d\n", sum(de$p_val_adj<0.05)))
cat("\nTop 15 上调 (High NeFRS):\n")
print(head(de[avg_log2FC>0, .(gene, avg_log2FC, p_val_adj, pct.1, pct.2)], 15))
cat("\nTop 15 下调:\n")
print(head(de[avg_log2FC<0][order(avg_log2FC), .(gene, avg_log2FC, p_val_adj)], 15))

fwrite(de, file.path(OUT, "malignant_NeFRS_hilo_DE.tsv"), sep="\t")

# Hwang 17 programs 在 High vs Low 的差异 (复用之前的 program score?)
# 简化: 看几个关键通路 marker 基因
cat("\n[通路] High NeFRS 恶性细胞特征基因检查...\n")
prog_markers <- list(
  Stress_TNF = c("NFKB1","TNF","JUN","FOS","ATF3","EGR1"),
  EMT = c("VIM","ZEB1","SNAI2","CDH2","FN1"),
  Cycling = c("MKI67","TOP2A","CCNB1","CDK1"),
  Neural = c("NRP1","NCAM1","STMN2","TUBB3")
)
data_mat <- GetAssayData(seu_hl, layer="data")
for (prog in names(prog_markers)) {
  gs <- intersect(prog_markers[[prog]], rownames(data_mat))
  if (length(gs)==0) next
  hi <- rowMeans(data_mat[gs, seu_hl$nefrs_grp=="High", drop=FALSE])
  lo <- rowMeans(data_mat[gs, seu_hl$nefrs_grp=="Low", drop=FALSE])
  cat(sprintf("  %s: High %.3f vs Low %.3f (diff %.3f)\n",
              prog, mean(hi), mean(lo), mean(hi)-mean(lo)))
}

# 火山图
de$sig <- "ns"
de$sig[de$p_val_adj<0.05 & de$avg_log2FC>0.5] <- "Up (High)"
de$sig[de$p_val_adj<0.05 & de$avg_log2FC< -0.5] <- "Down"
library(ggrepel)
top_lab <- rbind(head(de[sig=="Up (High)"],15), head(de[sig=="Down"][order(avg_log2FC)],10))
p <- ggplot(de, aes(avg_log2FC, -log10(p_val_adj+1e-300), color=sig)) +
  geom_point(size=0.6, alpha=0.5) +
  scale_color_manual(values=c("ns"="gray70","Up (High)"="#C0392B","Down"="#2980B9")) +
  geom_text_repel(data=top_lab, aes(label=gene), size=2.5, max.overlaps=30, color="black") +
  labs(title="Fig 8D. DE genes: High vs Low NeFRS malignant cells",
       x="log2FC (High/Low)", y="-log10(adj p)") +
  theme_bw(base_size=12)
ggsave(file.path(FIG,"Fig8D_malignant_hilo_volcano.pdf"), p, width=9, height=7)

cat(sprintf("\n=== B 完成, %.1f min ===\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
gc()
