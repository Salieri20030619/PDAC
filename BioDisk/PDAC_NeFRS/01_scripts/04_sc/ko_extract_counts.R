# 提取 Hwang h5ad 真实 counts (/layers/counts), 恶性细胞
suppressPackageStartupMessages({library(rhdf5); library(Matrix); library(data.table)})
H5AD <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/scRNA/GSE202051/GSE202051_totaldata-final-toshare.h5ad"
SC   <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
OUT  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/10_knockout"

cat("=== 提取真实 counts (/layers/counts) ===\n")
t0 <- Sys.time()

# 基因名 + barcode
gene_names <- h5read(H5AD, "/var/_index")
cell_barcodes <- h5read(H5AD, "/obs/_index")
cat(sprintf("h5ad: %d 基因 x %d 细胞\n", length(gene_names), length(cell_barcodes)))

# /layers/counts 是 CSR (cells x genes, anndata 默认)
data_c <- h5read(H5AD, "/layers/counts/data")
indices_c <- h5read(H5AD, "/layers/counts/indices")
indptr_c <- h5read(H5AD, "/layers/counts/indptr")
cat(sprintf("counts nnz: %d\n", length(data_c)))

# 检查是否整数 (确认是 counts)
cat(sprintf("counts 前几值: %s\n", paste(head(data_c,8), collapse=", ")))
cat(sprintf("是否整数: %s\n", all(data_c[1:1000]==round(data_c[1:1000]))))

# 构建 sparse (CSR: indptr 长度 = n_cells+1)
n_cells <- length(cell_barcodes); n_genes <- length(gene_names)
counts <- sparseMatrix(j=as.integer(indices_c)+1L, p=as.numeric(indptr_c), x=as.numeric(data_c),
                       dims=c(n_cells, n_genes),
                       dimnames=list(cell_barcodes, gene_names))
cat(sprintf("counts 矩阵: %d x %d\n", nrow(counts), ncol(counts)))
cat(sprintf("最大值: %.0f (真 counts 应该是大整数)\n", max(counts)))

# 恶性细胞子集
annot <- as.data.table(readRDS(file.path(SC,"w4_annot_with_ucell.rds")))
mal <- annot[Level2=="Malignant"]$barcode
idx <- match(mal, rownames(counts)); idx <- idx[!is.na(idx)]
counts_mal <- counts[idx,]
cat(sprintf("\n恶性细胞 counts: %d x %d\n", nrow(counts_mal), ncol(counts_mal)))

saveRDS(counts_mal, file.path(OUT, "malignant_real_counts.rds"))
cat(sprintf("存 malignant_real_counts.rds, %.1f min\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
H5closeAll()
