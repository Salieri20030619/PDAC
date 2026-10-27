# W4 Task 2: 提取核心数据 + 完整 label 解码
suppressPackageStartupMessages({
  library(rhdf5); library(data.table); library(Matrix)
})

H5AD <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/scRNA/GSE202051/GSE202051_totaldata-final-toshare.h5ad"
OUT  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

cat("=== W4 Task 2: Extract essential data ===\n")
t_total <- Sys.time()

# === 1. 基因名 + cell barcode ===
gene_names <- h5read(H5AD, "/var/_index")
cell_barcodes <- h5read(H5AD, "/obs/_index")
cat(sprintf("Cells: %d, Genes: %d\n", length(cell_barcodes), length(gene_names)))

# === 2. UMAP / Harmony ===
cat("Reading UMAP + Harmony embeddings...\n")
umap <- t(h5read(H5AD, "/obsm/X_umap"))
colnames(umap) <- c("UMAP_1", "UMAP_2")
harmony <- t(h5read(H5AD, "/obsm/X_harmony"))

# === 3. 读取所有 categorical 字段 + 解码 ===
cat("\n--- Reading categorical annotations ---\n")
categorical_cols <- c("Level 1 Annotation", "Level 2 Annotation", "Level 3 Annotation",
                      "broad_celltypes", "celltype_infercnv", "celltypes",
                      "cnv_leiden", "leiden", "new_celltypes",
                      "treatment_status", "new_treatment", "response",
                      "sampleid", "pid", "batch")

annot_df <- data.frame(barcode = cell_barcodes, stringsAsFactors=FALSE)

for (col in categorical_cols) {
  # 读 integer codes
  codes <- tryCatch(h5read(H5AD, paste0("/obs/", col)),
                     error=function(e) NULL)
  if (is.null(codes)) {
    cat(sprintf("  ✗ %s: codes missing\n", col)); next
  }
  
  # 读 string labels
  labels <- tryCatch(h5read(H5AD, paste0("/obs/__categories/", col)),
                      error=function(e) NULL)
  if (is.null(labels)) {
    annot_df[[col]] <- codes  # 没 label, 用数字
    cat(sprintf("  - %s: %d codes (no labels)\n", col, length(unique(codes))))
    next
  }
  
  # 解码: codes 是 0-indexed integer, labels 是 string vector
  decoded <- labels[codes + 1]
  annot_df[[col]] <- decoded
  
  # 报告 top 5 类别
  top5 <- sort(table(decoded), decreasing=TRUE)
  cat(sprintf("  ✓ %s (%d cats): %s\n", col, length(labels),
              paste(head(names(top5), 5), collapse=" / ")))
}

# 看 treatment_status 详细分布
cat("\n=== treatment_status full distribution ===\n")
print(table(annot_df$treatment_status, useNA="ifany"))

# 看 Level 3 详细分布 (Hwang NRP 应该在里面)
cat("\n=== Level 3 Annotation full distribution ===\n")
print(sort(table(annot_df$`Level 3 Annotation`), decreasing=TRUE))

# 看 new_celltypes 详细分布
cat("\n=== new_celltypes full distribution ===\n")
print(sort(table(annot_df$new_celltypes), decreasing=TRUE))

# 看 broad_celltypes
cat("\n=== broad_celltypes ===\n")
print(sort(table(annot_df$broad_celltypes), decreasing=TRUE))

# === 4. 保存 metadata-only (轻量, 不含表达矩阵) ===
saveRDS(list(
  cell_barcodes = cell_barcodes,
  gene_names = gene_names,
  umap = umap,
  annot = annot_df
), file.path(OUT, "w4_metadata_only.rds"))
fwrite(annot_df, file.path(OUT, "w4_cell_annotations.tsv"), sep="\t")
cat(sprintf("\nMetadata saved (without expression): %s/w4_metadata_only.rds\n", OUT))

# === 5. 读 NeFRS 53 基因表达矩阵 ===
NEFRS_GENES <- c("FGF13","MAP1S","REST","NRP1","CDC25A","KIF21A","TOR2A","LGMN",
                  "WWTR1","DTNBP1","DCX","EREG","ISCA2","FADS1","KRT18","BMPR1B",
                  "SOCS2","EGFR","CP","SLC25A46","MECP2","FBLN1","NQO1","ETHE1",
                  "EPHB6","VPS13A","NOTCH1","HES1","CX3CR1","SPATA2","SIRT2",
                  "NR2F2","FOXO4","VIPR1","MCF2L","STXBP1","BOC","MYCN","GABRB3",
                  "ATOH1","APLP1","FLRT3","WNT7A","ITGB1","IFRD1","EFNB2",
                  "SEMA3A","CAPRIN2","GALNT14","GJA1","HELLS","WNT7B","BRD7")

nefrs_idx <- which(gene_names %in% NEFRS_GENES)
nefrs_avail <- gene_names[nefrs_idx]
cat(sprintf("\n=== NeFRS 53 → GSE202051 ===\n"))
cat(sprintf("Found: %d / 53\n", length(nefrs_idx)))
missing <- setdiff(NEFRS_GENES, nefrs_avail)
if (length(missing) > 0) cat(sprintf("Missing: %s\n", paste(missing, collapse=", ")))

# === 6. 读 sparse 矩阵 ===
cat("\nReading X (sparse matrix, may take 1-3 min)...\n")
t0 <- Sys.time()
data_X    <- h5read(H5AD, "/X/data")
indices_X <- h5read(H5AD, "/X/indices")
indptr_X  <- h5read(H5AD, "/X/indptr")
cat(sprintf("X arrays loaded: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units="mins"))))
cat(sprintf("  data length: %d, indices length: %d, indptr length: %d\n",
            length(data_X), length(indices_X), length(indptr_X)))

# anndata CSR: shape = (n_obs, n_var)
# R Matrix sparseMatrix(j, p, x, dims, index1=FALSE) 构造 CSC
# CSR 的 indices 是 column indices (gene), indptr 是 row offsets (cell)
# 转 CSC 后, R 的 dgCMatrix 行=cell, 列=gene
n_obs <- length(indptr_X) - 1
n_var <- length(gene_names)

cat(sprintf("Building sparse matrix %d cells x %d genes...\n", n_obs, n_var))
t0 <- Sys.time()
X_sparse <- sparseMatrix(
  j = as.integer(indices_X) + 1L,
  p = as.integer(indptr_X),
  x = as.numeric(data_X),
  dims = c(n_obs, n_var),
  index1 = FALSE,
  repr = "C"
)
cat(sprintf("Sparse built in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units="mins"))))

rownames(X_sparse) <- cell_barcodes
colnames(X_sparse) <- gene_names

# 提取 NeFRS 子矩阵
nefrs_mat <- X_sparse[, nefrs_idx, drop=FALSE]
cat(sprintf("\nNeFRS sub-matrix: %d cells x %d genes\n", nrow(nefrs_mat), ncol(nefrs_mat)))

# === 7. 保存 ===
# 把 NeFRS 矩阵转回 cells x genes 但用 sparse 保留 (节省空间)
saveRDS(nefrs_mat, file.path(OUT, "w4_nefrs_expression.rds"))

# 同时保存全 X 也存一份 sparse (5-10 GB, 用于后续 module score)
saveRDS(X_sparse, file.path(OUT, "w4_X_sparse_full.rds"))

cat(sprintf("\n=== Task 2 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t_total, units="mins"))))
cat("Files:\n")
cat(sprintf("  %s/w4_metadata_only.rds\n", OUT))
cat(sprintf("  %s/w4_cell_annotations.tsv\n", OUT))
cat(sprintf("  %s/w4_nefrs_expression.rds\n", OUT))
cat(sprintf("  %s/w4_X_sparse_full.rds\n", OUT))

cat("\nMemory:\n")
gc()
