# GSE199102 空间 GeoMx: NeFRS 在 Epithelial/CAF/Immune 区的空间分布
# 直接用作者处理好的 RDS (Q3 归一化表达 + ROI metadata), 跳过 GeomxTools
suppressPackageStartupMessages({ library(data.table) })

SPAT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/spatial/GSE199102"
OUT  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/08_spatial"
MIME <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

cat("=== GSE199102 Spatial: NeFRS by segment ===\n")

# === 1. 读作者 RDS (已 gunzip 到 /tmp/test.RDS) ===
rds_path <- "/tmp/test.RDS"
if (!file.exists(rds_path)) {
  system(sprintf("cp %s/GSE199102_dfs_Q321-6.RDS.gz /tmp/t.RDS.gz && gunzip -f /tmp/t.RDS.gz", SPAT))
  rds_path <- "/tmp/t.RDS"
}
dfs <- readRDS(rds_path)
cat(sprintf("List with %d elements\n", length(dfs)))

# [[2]] = Q3 归一化表达 (gene × ROI), [[3]] = ROI metadata
norm_expr <- dfs[[2]]
meta <- dfs[[3]]
cat(sprintf("Norm expr: %d genes x %d ROI\n", nrow(norm_expr), ncol(norm_expr)-1))
cat(sprintf("Metadata: %d ROI x %d cols\n", nrow(meta), ncol(meta)))

# 第一列是 Gene 名
gene_col <- colnames(norm_expr)[1]
cat(sprintf("Gene 列名: %s\n", gene_col))
genes <- norm_expr[[gene_col]]

# 表达矩阵 (去掉 gene 列, 转 matrix)
expr_mat <- as.matrix(norm_expr[, -1])
rownames(expr_mat) <- genes
cat(sprintf("expr_mat: %d x %d\n", nrow(expr_mat), ncol(expr_mat)))

# === 2. NeFRS 系数 ===
coef_df <- fread(file.path(MIME, "nefrs_v3_enet_coefficients.tsv"))
avail <- intersect(coef_df$gene, rownames(expr_mat))
cat(sprintf("\nNeFRS genes in spatial: %d / %d\n", length(avail), nrow(coef_df)))
missing <- setdiff(coef_df$gene, rownames(expr_mat))
if (length(missing)>0) cat(sprintf("Missing: %s\n", paste(missing, collapse=", ")))

# === 3. 算每个 ROI 的 NeFRS score ===
# GeoMx Q3 归一化值通常 log 之前, 先 log2 转换 (跟 bulk 一致尺度)
expr_log <- log2(expr_mat[avail, ] + 1)
coef_v <- coef_df$coef[match(avail, coef_df$gene)]
nefrs_roi <- as.numeric(t(expr_log) %*% coef_v)
names(nefrs_roi) <- colnames(expr_mat)
cat(sprintf("NeFRS ROI score: range [%.3f, %.3f], median %.3f\n",
            min(nefrs_roi), max(nefrs_roi), median(nefrs_roi)))

# === 4. 对齐 metadata (Sample_ID) ===
meta_dt <- as.data.table(meta)
# expr 列名 vs meta Sample_ID 对齐
meta_dt[, nefrs := nefrs_roi[match(Sample_ID, names(nefrs_roi))]]
cat(sprintf("\nMatched: %d / %d ROI\n", sum(!is.na(meta_dt$nefrs)), nrow(meta_dt)))

# === 5. NeFRS by Segment (核心结果) ===
cat("\n=== NeFRS by Segment (Epithelial / CAF / Immune) ===\n")
seg_sum <- meta_dt[!is.na(nefrs), .(
  n = .N,
  median = round(median(nefrs), 4),
  mean = round(mean(nefrs), 4),
  sd = round(sd(nefrs), 4)
), by = segment][order(-median)]
print(seg_sum)

# 统计检验: Epithelial vs CAF vs Immune (Kruskal-Wallis)
kw <- kruskal.test(nefrs ~ as.factor(segment), data = meta_dt[!is.na(nefrs)])
cat(sprintf("\nKruskal-Wallis p = %.3e\n", kw$p.value))
# 两两 Wilcoxon
segs <- unique(meta_dt$segment)
cat("\n两两 Wilcoxon:\n")
for (i in 1:(length(segs)-1)) for (j in (i+1):length(segs)) {
  a <- meta_dt[segment==segs[i] & !is.na(nefrs)]$nefrs
  b <- meta_dt[segment==segs[j] & !is.na(nefrs)]$nefrs
  w <- wilcox.test(a, b)
  cat(sprintf("  %s vs %s: p=%.3e (median %.3f vs %.3f)\n",
              segs[i], segs[j], w$p.value, median(a), median(b)))
}

# === 6. NeFRS by treatment (Untreated vs CRT) ===
cat("\n=== NeFRS by treatment ===\n")
tx_sum <- meta_dt[!is.na(nefrs), .(n=.N, median=round(median(nefrs),4)), by=treatment][order(-median)]
print(tx_sum)

# === 7. NeFRS by segment × treatment ===
cat("\n=== NeFRS by segment × treatment ===\n")
st_sum <- meta_dt[!is.na(nefrs), .(n=.N, median=round(median(nefrs),4)),
                  by=.(segment, treatment)][order(segment, -median)]
print(st_sum)

# === 8. 保存 ===
out_meta <- meta_dt[, .(Sample_ID, Patient, segment, treatment,
                        X_coordinate, Y_coordinate,
                        Epithelial_proportion, CAF_proportion, Immune_proportion,
                        nefrs)]
fwrite(out_meta, file.path(OUT, "spatial_nefrs_roi.tsv"), sep="\t")
fwrite(seg_sum, file.path(OUT, "spatial_nefrs_by_segment.tsv"), sep="\t")
saveRDS(list(meta=meta_dt, nefrs_roi=nefrs_roi, seg_sum=seg_sum),
        file.path(OUT, "spatial_nefrs_result.rds"))
cat("\n=== Saved ===\n")
cat("  spatial_nefrs_roi.tsv (639 ROI + NeFRS + 坐标)\n")
cat("  spatial_nefrs_by_segment.tsv\n")
cat("  spatial_nefrs_result.rds\n")
