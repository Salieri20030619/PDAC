# 拟时序段1: slingshot 轨迹 (全上皮谱系, 起点 Acinar) + NeFRS 沿 pseudotime
# 不含 fitGAM (段2单独跑), 内存安全
suppressPackageStartupMessages({
  library(Seurat); library(slingshot); library(SingleCellExperiment)
  library(data.table); library(Matrix); library(ggplot2); library(viridis)
})

SC  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/09_pseudotime"
FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures/Fig8_pseudotime"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

cat("=== 拟时序段1: slingshot ===\n")
t0 <- Sys.time()

# === 1. 读数据, 抽上皮谱系 ===
cat("[1] 读 X + annot...\n")
X <- readRDS(file.path(SC, "w4_X_sparse_full.rds"))   # cells x genes
annot <- as.data.table(readRDS(file.path(SC, "w4_annot_with_ucell.rds")))

epi_types <- c("Acinar","ADM","Ductal","Ductal (atypical)","Malignant")
epi_cells <- annot[Level2 %in% epi_types]$barcode
cat(sprintf("上皮谱系: %d cells\n", length(epi_cells)))

# 对齐 X 的行 (cells)
idx <- match(epi_cells, rownames(X))
idx <- idx[!is.na(idx)]
X_epi <- X[idx, ]
annot_epi <- annot[match(rownames(X_epi), barcode)]
cat(sprintf("X_epi: %d cells x %d genes\n", nrow(X_epi), ncol(X_epi)))
rm(X); gc()

cat("\n各类型细胞数:\n")
print(table(annot_epi$Level2))

# === 2. 建 Seurat (X 已是 log-norm) ===
# 注: Seurat 5.4.0 CreateSeuratObject(counts=矩阵) 有 bug, 先建 assay 再组装
cat("\n[2] 建 Seurat + PCA + UMAP...\n")
expr_t <- t(X_epi)  # genes x cells
assay_obj <- CreateAssayObject(counts=expr_t)
seu <- CreateSeuratObject(counts=assay_obj)
# X 已 log-norm: 把 data layer 设为 X (而非重新 NormalizeData)
seu[["RNA"]]$data <- expr_t
seu$Level2 <- annot_epi$Level2
seu$NeFRS_all <- annot_epi$NeFRS_all
seu$NeFRS_risk <- annot_epi$NeFRS_risk
rm(X_epi, expr_t, assay_obj); gc()

seu <- FindVariableFeatures(seu, nfeatures=2000, verbose=FALSE)
seu <- ScaleData(seu, verbose=FALSE)
seu <- RunPCA(seu, npcs=30, verbose=FALSE)
seu <- RunUMAP(seu, dims=1:20, verbose=FALSE)
cat(sprintf("  降维完成 (%.1f min)\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
gc()

# === 3. slingshot (起点 Acinar) ===
cat("\n[3] slingshot (start.clus = Acinar)...\n")
sce <- as.SingleCellExperiment(seu)
colData(sce)$Level2 <- factor(seu$Level2)
reducedDim(sce, "UMAP") <- Embeddings(seu, "umap")

sce <- slingshot(sce, clusterLabels="Level2", reducedDim="UMAP",
                 start.clus="Acinar")
pt <- slingPseudotime(sce)
cat(sprintf("  推断 %d 条 lineage\n", ncol(pt)))
print(apply(pt, 2, function(x) range(x, na.rm=TRUE)))

# lineage 结构
cat("\nLineage 细胞类型顺序:\n")
sds <- SlingshotDataSet(sce)
print(slingLineages(sds))

# === 4. NeFRS 沿 pseudotime (核心结果, 不需 fitGAM) ===
cat("\n[4] NeFRS 沿 pseudotime...\n")
pt_main <- pt[, 1]  # 主 lineage
df <- data.frame(
  pseudotime = pt_main,
  NeFRS = seu$NeFRS_all,
  NeFRS_risk = seu$NeFRS_risk,
  celltype = seu$Level2,
  UMAP_1 = Embeddings(seu,"umap")[,1],
  UMAP_2 = Embeddings(seu,"umap")[,2]
)
df <- df[!is.na(df$pseudotime), ]

# 按 pseudotime 分 bin, 看 NeFRS 趋势
df$pt_bin <- cut(df$pseudotime, breaks=20, labels=FALSE)
bin_sum <- aggregate(cbind(NeFRS, NeFRS_risk) ~ pt_bin, df, median)
cat("NeFRS 沿 pseudotime (20 bins):\n")
print(head(bin_sum, 20))

# 相关性
cor_nefrs <- cor(df$pseudotime, df$NeFRS, method="spearman", use="complete.obs")
cat(sprintf("\nNeFRS vs pseudotime Spearman rho = %.3f\n", cor_nefrs))

# === 5. 保存 + 图 ===
saveRDS(sce, file.path(OUT, "pseudo_sce_slingshot.rds"))
saveRDS(seu, file.path(OUT, "pseudo_seu_epi.rds"))
fwrite(df, file.path(OUT, "pseudo_nefrs_along_pt.tsv"), sep="\t")

# Fig8A: UMAP by celltype + 轨迹
ct_cols <- c("Acinar"="#27AE60","ADM"="#16A085","Ductal"="#2980B9",
             "Ductal (atypical)"="#E67E22","Malignant"="#C0392B")
crv <- slingCurves(sce)
crv_df <- do.call(rbind, lapply(seq_along(crv), function(i){
  p <- crv[[i]]$s[crv[[i]]$ord,,drop=FALSE]
  data.frame(UMAP_1=p[,1], UMAP_2=p[,2], Lineage=paste0("L",i))
}))

p1 <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color=celltype), size=0.4, alpha=0.6) +
  scale_color_manual(values=ct_cols) +
  geom_path(data=crv_df, aes(UMAP_1,UMAP_2,group=Lineage), color="black", linewidth=1) +
  labs(title="Fig 8A. Epithelial trajectory (Acinar to Malignant)", color="Cell type") +
  theme_bw(base_size=12)
ggsave(file.path(FIG,"Fig8A_trajectory_celltype.pdf"), p1, width=8, height=6)

# Fig8B: UMAP by NeFRS
p2 <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color=NeFRS), size=0.4, alpha=0.7) +
  scale_color_viridis_c(option="magma") +
  geom_path(data=crv_df, aes(UMAP_1,UMAP_2,group=Lineage), color="black", linewidth=1) +
  labs(title="Fig 8B. NeFRS along trajectory", color="NeFRS") +
  theme_bw(base_size=12)
ggsave(file.path(FIG,"Fig8B_trajectory_NeFRS.pdf"), p2, width=8, height=6)

# Fig8C: NeFRS vs pseudotime (核心图)
p3 <- ggplot(df, aes(pseudotime, NeFRS)) +
  geom_point(aes(color=celltype), size=0.3, alpha=0.3) +
  scale_color_manual(values=ct_cols) +
  geom_smooth(method="loess", color="black", linewidth=1.2, se=TRUE) +
  labs(title=sprintf("Fig 8C. NeFRS along pseudotime (rho=%.3f)", cor_nefrs),
       x="Pseudotime (Acinar to Malignant)", y="NeFRS score", color="Cell type") +
  theme_bw(base_size=12)
ggsave(file.path(FIG,"Fig8C_NeFRS_vs_pseudotime.pdf"), p3, width=8, height=6)

cat("\n=== 段1 完成 ===\n")
cat(sprintf("耗时 %.1f min\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
cat("输出: pseudo_sce_slingshot.rds, pseudo_seu_epi.rds, Fig8A/B/C\n")
gc()
