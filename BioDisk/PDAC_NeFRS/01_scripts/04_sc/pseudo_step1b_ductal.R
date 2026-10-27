# 拟时序 A: 导管谱系局部轨迹 (Ductal → Ductal atypical), 带 Harmony 整合
# 生物学依据: 正常导管→异型导管是连续过渡, NeFRS atypical>Ductal 有梯度
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(slingshot); library(SingleCellExperiment)
  library(data.table); library(Matrix); library(ggplot2); library(viridis)
})

SC  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/09_pseudotime"
FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures/Fig8_pseudotime"

cat("=== 拟时序 A: 导管谱系轨迹 (Harmony 整合) ===\n")
t0 <- Sys.time()

# === 1. 抽导管谱系 ===
X <- readRDS(file.path(SC, "w4_X_sparse_full.rds"))
annot <- as.data.table(readRDS(file.path(SC, "w4_annot_with_ucell.rds")))

duct_types <- c("Ductal","Ductal (atypical)")
duct_cells <- annot[Level2 %in% duct_types]$barcode
idx <- match(duct_cells, rownames(X)); idx <- idx[!is.na(idx)]
X_d <- X[idx, ]
annot_d <- annot[match(rownames(X_d), barcode)]
cat(sprintf("导管谱系: %d cells\n", nrow(X_d)))
print(table(annot_d$Level2))
cat(sprintf("patient 数: %d\n", length(unique(annot_d$pid))))
rm(X); gc()

# === 2. 建 Seurat (修正版建对象) + Harmony ===
cat("\n[2] Seurat + Harmony 整合 (by pid)...\n")
expr_t <- t(X_d)
assay_obj <- CreateAssayObject(counts=expr_t)
seu <- CreateSeuratObject(counts=assay_obj)
seu[["RNA"]]$data <- expr_t   # X 已 log-norm
seu$Level2 <- annot_d$Level2
seu$NeFRS_all <- annot_d$NeFRS_all
seu$NeFRS_risk <- annot_d$NeFRS_risk
seu$pid <- annot_d$pid
rm(X_d, expr_t, assay_obj); gc()

seu <- FindVariableFeatures(seu, nfeatures=2000, verbose=FALSE)
seu <- ScaleData(seu, verbose=FALSE)
seu <- RunPCA(seu, npcs=30, verbose=FALSE)
seu <- RunHarmony(seu, group.by.vars="pid", dims.use=1:30, verbose=FALSE)
seu <- RunUMAP(seu, reduction="harmony", dims=1:20, verbose=FALSE)
cat(sprintf("  整合+降维完成 (%.1f min)\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
gc()

# === 3. slingshot (起点 Ductal, 终点 atypical) ===
cat("\n[3] slingshot (start=Ductal, end=Ductal atypical)...\n")
sce <- as.SingleCellExperiment(seu)
colData(sce)$Level2 <- factor(seu$Level2)
reducedDim(sce, "UMAP") <- Embeddings(seu, "umap")
sce <- slingshot(sce, clusterLabels="Level2", reducedDim="UMAP",
                 start.clus="Ductal", end.clus="Ductal (atypical)")
pt <- slingPseudotime(sce)
cat(sprintf("  %d lineage\n", ncol(pt)))
print(slingLineages(SlingshotDataSet(sce)))

# === 4. NeFRS 沿 pseudotime ===
pt_main <- pt[, 1]
df <- data.frame(pseudotime=pt_main, NeFRS=seu$NeFRS_all, NeFRS_risk=seu$NeFRS_risk,
                 celltype=seu$Level2,
                 UMAP_1=Embeddings(seu,"umap")[,1], UMAP_2=Embeddings(seu,"umap")[,2])
df <- df[!is.na(df$pseudotime), ]
cor_n <- cor(df$pseudotime, df$NeFRS, method="spearman", use="complete.obs")
cat(sprintf("\nNeFRS vs pseudotime rho = %.3f\n", cor_n))

df$pt_bin <- cut(df$pseudotime, breaks=20, labels=FALSE)
bin_sum <- aggregate(cbind(NeFRS,NeFRS_risk) ~ pt_bin, df, median)
cat("NeFRS 沿 pseudotime (bins):\n"); print(bin_sum)

# === 5. 保存 + 图 (se=FALSE 避免内存爆) ===
saveRDS(sce, file.path(OUT, "pseudo_ductal_sce.rds"))
saveRDS(seu, file.path(OUT, "pseudo_ductal_seu.rds"))
fwrite(df, file.path(OUT, "pseudo_ductal_nefrs.tsv"), sep="\t")

ct_cols <- c("Ductal"="#2980B9","Ductal (atypical)"="#E67E22")
crv <- slingCurves(sce)
crv_df <- do.call(rbind, lapply(seq_along(crv), function(i){
  p<-crv[[i]]$s[crv[[i]]$ord,,drop=FALSE]; data.frame(UMAP_1=p[,1],UMAP_2=p[,2],L=paste0("L",i))}))

p1 <- ggplot(df, aes(UMAP_1,UMAP_2)) +
  geom_point(aes(color=celltype), size=0.5, alpha=0.6) +
  scale_color_manual(values=ct_cols) +
  geom_path(data=crv_df, aes(UMAP_1,UMAP_2,group=L), color="black", linewidth=1) +
  labs(title="Fig 8A2. Ductal lineage trajectory (Harmony integrated)", color="Cell type") +
  theme_bw(base_size=12)
ggsave(file.path(FIG,"Fig8A2_ductal_celltype.pdf"), p1, width=8, height=6)

p2 <- ggplot(df, aes(UMAP_1,UMAP_2)) +
  geom_point(aes(color=NeFRS), size=0.5, alpha=0.7) +
  scale_color_viridis_c(option="magma") +
  geom_path(data=crv_df, aes(UMAP_1,UMAP_2,group=L), color="black", linewidth=1) +
  labs(title="Fig 8B2. NeFRS along ductal trajectory", color="NeFRS") +
  theme_bw(base_size=12)
ggsave(file.path(FIG,"Fig8B2_ductal_NeFRS.pdf"), p2, width=8, height=6)

# Fig8C2: NeFRS vs pt, 用 bin 均值线 (不用 loess se, 避免内存爆)
bin_df <- df
bin_df$pt_bin20 <- cut(bin_df$pseudotime, breaks=20)
binmean <- aggregate(NeFRS ~ pt_bin20, bin_df, median)
binmean$pt_mid <- aggregate(pseudotime ~ pt_bin20, bin_df, median)$pseudotime
p3 <- ggplot(df, aes(pseudotime, NeFRS)) +
  geom_point(aes(color=celltype), size=0.3, alpha=0.3) +
  scale_color_manual(values=ct_cols) +
  geom_line(data=binmean, aes(pt_mid, NeFRS), color="black", linewidth=1.2, inherit.aes=FALSE) +
  geom_point(data=binmean, aes(pt_mid, NeFRS), color="black", size=2, inherit.aes=FALSE) +
  labs(title=sprintf("Fig 8C2. NeFRS along ductal pseudotime (rho=%.3f)", cor_n),
       x="Pseudotime (Ductal to atypical)", y="NeFRS", color="Cell type") +
  theme_bw(base_size=12)
ggsave(file.path(FIG,"Fig8C2_ductal_NeFRS_pt.pdf"), p3, width=8, height=6)

cat(sprintf("\n=== A 完成, %.1f min ===\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
gc()
