# Werba Step 2: 标准化 + Harmony 整合 + 聚类 + major cell type 注释 + NeFRS 评分
# 目标: 复现 Hwang Fig6A 的 NeFRS by cell type 分布 (独立验证)
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(Matrix); library(data.table); library(dplyr)
})

OUT  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/07_werba"
MIME <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime"

cat("=== Werba Step 2: Integrate + Annotate + NeFRS ===\n")
t_total <- Sys.time()

seu <- readRDS(file.path(OUT, "werba_seu_qc.rds"))
cat(sprintf("Loaded: %d genes x %d cells\n", nrow(seu), ncol(seu)))

# === 1. 标准化 + HVG + PCA ===
cat("\n[1] Normalize + HVG + Scale + PCA...\n")
t0 <- Sys.time()
seu <- NormalizeData(seu, verbose=FALSE)
seu <- FindVariableFeatures(seu, nfeatures=2000, verbose=FALSE)
seu <- ScaleData(seu, verbose=FALSE)
seu <- RunPCA(seu, npcs=30, verbose=FALSE)
cat(sprintf("  PCA done (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))

# === 2. Harmony 整合 (去样本批次) ===
cat("\n[2] Harmony integration (by sample)...\n")
t0 <- Sys.time()
seu <- RunHarmony(seu, group.by.vars="sample", dims.use=1:30, verbose=FALSE)
cat(sprintf("  Harmony done (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))

# === 3. 聚类 ===
cat("\n[3] Clustering (Louvain, Harmony embedding)...\n")
t0 <- Sys.time()
seu <- FindNeighbors(seu, reduction="harmony", dims=1:30, verbose=FALSE)
seu <- FindClusters(seu, resolution=0.5, verbose=FALSE)
cat(sprintf("  %d clusters (%.1f min)\n", length(unique(seu$seurat_clusters)),
            as.numeric(difftime(Sys.time(), t0, units="mins"))))

# === 4. Major cell type 注释 (marker module score) ===
cat("\n[4] Cell type annotation (marker-based)...\n")
markers <- list(
  Malignant_Ductal = c("KRT19","KRT18","KRT8","EPCAM","CEACAM5","TSPAN8","MUC1","SOX9"),
  CAF              = c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRB","ACTA2","FAP"),
  Endothelial      = c("PECAM1","VWF","CDH5","CLDN5","PLVAP"),
  Tcell            = c("CD3D","CD3E","CD8A","CD4","IL7R","CCL5"),
  Bcell_Plasma     = c("CD79A","CD79B","MS4A1","MZB1","IGHG1","JCHAIN"),
  Myeloid          = c("CD68","CD14","LYZ","ITGAM","AIF1","C1QA","C1QB"),
  Acinar           = c("PRSS1","CPA1","CTRB1","CTRB2","CELA3A","REG1A"),
  Endocrine        = c("INS","GCG","SST","PPY","CHGA","CHGB"),
  Mast             = c("CPA3","TPSAB1","TPSB2","MS4A2"),
  Schwann          = c("SOX10","S100B","PLP1","MPZ","CDH19")
)
# 只保留存在的基因
markers <- lapply(markers, function(g) intersect(g, rownames(seu)))
for (ct in names(markers)) {
  seu <- AddModuleScore(seu, features=list(markers[[ct]]), name=paste0("ms_", ct), seed=42)
}
# module score 列名 (Seurat 加了 "1" 后缀)
ms_cols <- paste0("ms_", names(markers), "1")

# 每个 cluster 的平均 module score → argmax 注释
clust_scores <- seu@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(ms_cols), mean))
clust_scores <- as.data.frame(clust_scores)
rownames(clust_scores) <- clust_scores$seurat_clusters

# argmax 分配 cell type
clust_anno <- apply(clust_scores[, ms_cols], 1, function(x) {
  names(markers)[which.max(x)]
})
names(clust_anno) <- clust_scores$seurat_clusters

cat("\n=== Cluster → cell type assignment ===\n")
for (cl in names(clust_anno)) {
  n_cl <- sum(seu$seurat_clusters == cl)
  cat(sprintf("  Cluster %s (n=%d) → %s\n", cl, n_cl, clust_anno[cl]))
}

seu$cell_type <- unname(clust_anno[as.character(seu$seurat_clusters)])

cat("\n=== Cell type distribution ===\n")
print(table(seu$cell_type))

# === 5. NeFRS 评分 (sum coef × expr) ===
cat("\n[5] NeFRS score (sum coef × expr)...\n")
coef_df <- fread(file.path(MIME, "nefrs_v3_enet_coefficients.tsv"))
avail <- intersect(coef_df$gene, rownames(seu))
cat(sprintf("NeFRS genes in Werba: %d / %d\n", length(avail), nrow(coef_df)))

# 用 normalized data (logNorm) 算 score
expr <- GetAssayData(seu, layer="data")[avail, ]
coef_v <- coef_df$coef[match(avail, coef_df$gene)]
nefrs_score <- as.numeric(t(expr) %*% coef_v)
seu$NeFRS <- nefrs_score
cat(sprintf("NeFRS score: range [%.3f, %.3f], median %.3f\n",
            min(nefrs_score), max(nefrs_score), median(nefrs_score)))

# === 6. NeFRS by cell type (对照 Hwang Fig6A) ===
cat("\n=== NeFRS by cell type (median, descending) ===\n")
nefrs_by_ct <- seu@meta.data %>%
  group_by(cell_type) %>%
  summarise(n=n(), median=median(NeFRS), mean=mean(NeFRS),
            p90=quantile(NeFRS, 0.9)) %>%
  arrange(desc(median)) %>%
  as.data.frame()
print(nefrs_by_ct)

# === 7. 保存 ===
fwrite(nefrs_by_ct, file.path(OUT, "werba_nefrs_by_celltype.tsv"), sep="\t")
saveRDS(seu, file.path(OUT, "werba_seu_annotated.rds"))
# 轻量 metadata
meta_out <- seu@meta.data[, c("sample","seurat_clusters","cell_type","NeFRS",
                               "nFeature_RNA","nCount_RNA","percent.mt")]
fwrite(data.table(barcode=rownames(meta_out), meta_out),
       file.path(OUT, "werba_metadata.tsv"), sep="\t")

cat("\n=== Saved ===\n")
cat("  werba_seu_annotated.rds\n")
cat("  werba_nefrs_by_celltype.tsv\n")
cat("  werba_metadata.tsv\n")
cat(sprintf("=== Step 2 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t_total, units="mins"))))
