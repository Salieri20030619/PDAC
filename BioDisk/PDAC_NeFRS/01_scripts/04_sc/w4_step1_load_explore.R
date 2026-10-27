# W4 Task 1: GSE202051 加载 + 探索 (用 anndataR 替代 zellkonverter)
suppressPackageStartupMessages({
  library(anndataR); library(SingleCellExperiment)
  library(HDF5Array); library(DelayedArray); library(data.table)
})

H5AD <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/scRNA/GSE202051/GSE202051_totaldata-final-toshare.h5ad"
OUT  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

cat("=== W4 Task 1: Load + Explore GSE202051 (anndataR) ===\n")
cat(sprintf("File: %s\n", H5AD))
cat(sprintf("Size: %.2f GB\n", file.size(H5AD)/1e9))

# === 用 anndataR 加载 ===
cat("\n--- Loading h5ad via anndataR ---\n")
t0 <- Sys.time()

# anndataR 的 read_h5ad 返回 anndataR 对象, 然后转 SCE
ad <- tryCatch({
  read_h5ad(H5AD, mode="r")  # mode="r" 是 backed (不全读入内存)
}, error = function(e) {
  cat("backed mode failed:", conditionMessage(e), "\n")
  cat("trying in-memory mode...\n")
  read_h5ad(H5AD)
})

elapsed <- difftime(Sys.time(), t0, units="mins")
cat(sprintf("Loaded in %.1f minutes\n", as.numeric(elapsed)))
cat(sprintf("Class: %s\n", paste(class(ad), collapse=", ")))

# 看 ad 有哪些 slot
cat("\n--- anndataR object structure ---\n")
print(ad)

# === 基本统计 ===
cat("\n--- Basic stats ---\n")
n_obs <- ad$n_obs
n_var <- ad$n_var
cat(sprintf("Cells (obs): %d\n", n_obs))
cat(sprintf("Genes (var): %d\n", n_var))

# === Cell metadata (obs) ===
cat("\n--- obs (cell metadata) columns ---\n")
obs_df <- as.data.frame(ad$obs)
cat(sprintf("Total columns: %d\n", ncol(obs_df)))
print(colnames(obs_df))

cat("\n--- obs column unique values ---\n")
for (col in colnames(obs_df)) {
  v <- obs_df[[col]]
  if (is.factor(v) || is.character(v)) {
    uv <- unique(as.character(v))
    cat(sprintf("  %-30s: %d unique (%s)\n",
                col, length(uv),
                paste(head(uv, 5), collapse=" / ")))
  } else if (is.numeric(v)) {
    cat(sprintf("  %-30s: numeric, range [%.2f, %.2f]\n",
                col, min(v, na.rm=TRUE), max(v, na.rm=TRUE)))
  } else {
    cat(sprintf("  %-30s: %s\n", col, class(v)[1]))
  }
}

# === 看是否有 cell type / NRP 注释 ===
cat("\n--- Potential cell type / NRP annotations ---\n")
candidates <- grep("cell|type|cluster|leiden|louvain|nrp|annot|state|class",
                    colnames(obs_df), value=TRUE, ignore.case=TRUE)
cat(sprintf("Candidates: %s\n", paste(candidates, collapse=", ")))
for (col in candidates) {
  cat(sprintf("\n  === %s ===\n", col))
  tab <- sort(table(obs_df[[col]]), decreasing=TRUE)
  print(head(tab, 15))
}

# === Reduced dims (obsm 在 anndata 里) ===
cat("\n--- obsm (reduced dims) keys ---\n")
obsm_keys <- tryCatch(names(ad$obsm), error=function(e) character(0))
print(obsm_keys)
for (k in obsm_keys) {
  m <- ad$obsm[[k]]
  cat(sprintf("  %s: dim %d x %d\n", k, nrow(m), ncol(m)))
}

# === uns (unstructured) ===
cat("\n--- uns keys ---\n")
uns_keys <- tryCatch(names(ad$uns), error=function(e) character(0))
print(head(uns_keys, 30))

# === 基因名 ===
cat("\n--- Gene names sample ---\n")
gene_names <- rownames(ad$var)
if (is.null(gene_names)) gene_names <- ad$var[[1]]
cat(sprintf("First 10: %s\n", paste(head(gene_names, 10), collapse=", ")))

# === NeFRS 53 基因匹配 ===
NEFRS_GENES <- c("FGF13","MAP1S","REST","NRP1","CDC25A","KIF21A","TOR2A","LGMN",
                  "WWTR1","DTNBP1","DCX","EREG","ISCA2","FADS1","KRT18","BMPR1B",
                  "SOCS2","EGFR","CP","SLC25A46","MECP2","FBLN1","NQO1","ETHE1",
                  "EPHB6","VPS13A","NOTCH1","HES1","CX3CR1","SPATA2","SIRT2",
                  "NR2F2","FOXO4","VIPR1","MCF2L","STXBP1","BOC","MYCN","GABRB3",
                  "ATOH1","APLP1","FLRT3","WNT7A","ITGB1","IFRD1","EFNB2",
                  "SEMA3A","CAPRIN2","GALNT14","GJA1","HELLS","WNT7B","BRD7")
nefrs_in <- intersect(NEFRS_GENES, gene_names)
nefrs_miss <- setdiff(NEFRS_GENES, gene_names)
cat(sprintf("\n=== NeFRS 53 genes in GSE202051 ===\n"))
cat(sprintf("Found: %d / 53 (%.0f%%)\n", length(nefrs_in), 100*length(nefrs_in)/53))
if (length(nefrs_miss) > 0) {
  cat(sprintf("Missing: %s\n", paste(nefrs_miss, collapse=", ")))
}

# === 保存 metadata ===
saveRDS(list(
  n_cells = n_obs, n_genes = n_var,
  obs_columns = colnames(obs_df),
  obsm_keys = obsm_keys,
  candidates = candidates,
  nefrs_in = nefrs_in, nefrs_miss = nefrs_miss
), file.path(OUT, "w4_step1_metadata.rds"))

# 同时保存 obs (cell-level metadata, 不大, ~10-50 MB)
fwrite(obs_df, file.path(OUT, "w4_obs_metadata.tsv"), sep="\t")
cat(sprintf("\nSaved obs to: %s/w4_obs_metadata.tsv\n", OUT))

cat("\n=== Task 1 done ===\n")
