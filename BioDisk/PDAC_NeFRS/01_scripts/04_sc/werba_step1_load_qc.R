# Werba 2023 (GSE205013) Step 1: 读入 27 样本 10x + 合并 + QC
# Seurat v5 API
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table)
})

RAW <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/scRNA/GSE205013_Werba"
OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/07_werba"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

cat("=== Werba Step 1: Load + Merge + QC ===\n")
t_total <- Sys.time()

# === 1. 找所有样本 (GSM_Pxx 前缀) ===
mtx_files <- list.files(RAW, pattern="_matrix\\.mtx\\.gz$")
samples <- gsub("_matrix\\.mtx\\.gz$", "", mtx_files)
cat(sprintf("Found %d samples\n", length(samples)))

# === 2. 逐样本读入 (10x 三件套, 手动读避免命名问题) ===
read_one <- function(prefix) {
  bc   <- file.path(RAW, paste0(prefix, "_barcodes.tsv.gz"))
  feat <- file.path(RAW, paste0(prefix, "_features.tsv.gz"))
  mtx  <- file.path(RAW, paste0(prefix, "_matrix.mtx.gz"))

  m <- readMM(mtx)
  barcodes <- readLines(bc)
  features <- fread(feat, header=FALSE)
  # 第 2 列 = gene symbol
  rownames(m) <- make.unique(features$V2)
  colnames(m) <- paste0(prefix, "|", barcodes)  # 防 barcode 跨样本冲突
  m
}

cat("\nReading samples...\n")
mat_list <- list()
for (s in samples) {
  t0 <- Sys.time()
  mat_list[[s]] <- read_one(s)
  cat(sprintf("  %s: %d genes x %d cells (%.1fs)\n",
              s, nrow(mat_list[[s]]), ncol(mat_list[[s]]),
              as.numeric(difftime(Sys.time(), t0, units="secs"))))
}

# === 3. 合并成一个大矩阵 (所有样本基因一致 36601) ===
cat("\nMerging all samples...\n")
# 确认基因一致
gene_sets <- lapply(mat_list, rownames)
all_same <- all(sapply(gene_sets[-1], function(g) identical(g, gene_sets[[1]])))
cat(sprintf("All samples same gene order: %s\n", all_same))

if (all_same) {
  big_mat <- do.call(cbind, mat_list)
} else {
  # 取交集对齐
  common <- Reduce(intersect, gene_sets)
  cat(sprintf("Using %d common genes\n", length(common)))
  big_mat <- do.call(cbind, lapply(mat_list, function(m) m[common, ]))
}
cat(sprintf("Merged: %d genes x %d cells\n", nrow(big_mat), ncol(big_mat)))
rm(mat_list); gc()

# === 4. 建 Seurat 对象 ===
cat("\nCreating Seurat object...\n")
sample_id <- sub("\\|.*$", "", colnames(big_mat))
seu <- CreateSeuratObject(counts=big_mat, min.cells=3, min.features=200,
                          meta.data=data.frame(row.names=colnames(big_mat),
                                               sample=sample_id))
cat(sprintf("Seurat: %d genes x %d cells (after min.cells=3, min.features=200)\n",
            nrow(seu), ncol(seu)))
rm(big_mat); gc()

# === 5. QC 指标 ===
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")
cat("\n=== QC metrics (before filter) ===\n")
cat(sprintf("nFeature: median=%.0f, range=[%d, %d]\n",
            median(seu$nFeature_RNA), min(seu$nFeature_RNA), max(seu$nFeature_RNA)))
cat(sprintf("nCount: median=%.0f\n", median(seu$nCount_RNA)))
cat(sprintf("percent.mt: median=%.1f%%, >20%%: %d cells\n",
            median(seu$percent.mt), sum(seu$percent.mt > 20)))

# === 6. QC 过滤 ===
n_before <- ncol(seu)
seu <- subset(seu, subset = nFeature_RNA >= 200 & nFeature_RNA <= 7500 & percent.mt < 20)
cat(sprintf("\nAfter QC filter: %d cells (removed %d, %.1f%%)\n",
            ncol(seu), n_before - ncol(seu), 100*(n_before-ncol(seu))/n_before))

# 每样本细胞数
cat("\n=== Cells per sample (after QC) ===\n")
print(table(seu$sample))

# === 7. 保存 ===
saveRDS(seu, file.path(OUT, "werba_seu_qc.rds"))
cat(sprintf("\nSaved: werba_seu_qc.rds\n"))
cat(sprintf("=== Step 1 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t_total, units="mins"))))
