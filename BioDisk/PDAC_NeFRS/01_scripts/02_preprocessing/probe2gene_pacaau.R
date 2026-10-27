#!/usr/bin/env Rscript
# ============================================================
# PACA-AU probe → gene symbol 映射 + 多 probe 聚合
#
# 输入:
#   exp_array.PACA-AU.matrix.tsv.gz  (47265 ILMN probes × 269 samples)
#   GPL10558.annot.gz                (NCBI GEO 官方注释 48137 probes)
#
# 输出:
#   exp_array.PACA-AU.gene.tsv.gz    (gene symbol × 269 samples)
#   probe2gene_map.tsv               (probe→gene 映射记录, 含选择信息)
#   probe2gene_summary.tsv           (统计摘要)
#
# 聚合策略:
#   同 gene 多 probe 取最大方差 probe (业界标准)
# ============================================================

suppressPackageStartupMessages({
    library(data.table)
})

ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/ICGC_PACA_AU"
EXP_FILE     <- file.path(ROOT, "exp_array.PACA-AU.matrix.tsv.gz")
ANNOT_FILE   <- file.path(ROOT, "GPL10558.annot.gz")
OUT_GENE     <- file.path(ROOT, "exp_array.PACA-AU.gene.tsv.gz")
OUT_MAP      <- file.path(ROOT, "probe2gene_map.tsv")
OUT_SUMMARY  <- file.path(ROOT, "probe2gene_summary.tsv")

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# 1. 读取 GPL10558 注释 (跳过 header，从 !platform_table_begin 后开始)
# ============================================================
log_msg("=== 读取 GPL10558 注释 ===")

# annot 文件结构：!platform_table_begin 是第 29 行，第 30 行是列名
# 我们用 skip = 29 让 fread 从列名行开始
annot <- fread(ANNOT_FILE,
               sep = "\t",
               skip = "ID\tGene title",   # 自动定位到列名行
               header = TRUE,
               quote = "",
               fill = TRUE,
               showProgress = FALSE)

log_msg(sprintf("  原始注释: %d 行 × %d 列", nrow(annot), ncol(annot)))
log_msg(sprintf("  列名: %s", paste(colnames(annot)[1:5], collapse = " | ")))

# 只保留 ID 和 Gene symbol 两列
annot <- annot[, .(probe_id = ID, gene_symbol = `Gene symbol`)]

# 过滤：去掉空 symbol 和 control probe
annot <- annot[!is.na(probe_id) & probe_id != ""]
annot[, gene_symbol := trimws(gene_symbol)]
annot_valid <- annot[!is.na(gene_symbol) & gene_symbol != ""]

log_msg(sprintf("  总 probe 数: %d", nrow(annot)))
log_msg(sprintf("  含 gene symbol 的 probe: %d (%.1f%%)",
                nrow(annot_valid),
                100 * nrow(annot_valid) / nrow(annot)))
log_msg(sprintf("  唯一 gene symbol 数: %d",
                length(unique(annot_valid$gene_symbol))))

# ============================================================
# 2. 读取 PACA-AU 表达矩阵
# ============================================================
log_msg("=== 读取 PACA-AU 表达矩阵 ===")
exp_mat <- fread(EXP_FILE, header = TRUE, showProgress = FALSE)
log_msg(sprintf("  矩阵: %d 行 × %d 列", nrow(exp_mat), ncol(exp_mat)))
log_msg(sprintf("  样本数: %d", ncol(exp_mat) - 1))

# 处理可能的 V7/V8 残留行 (之前发现的列名混入数据)
exp_mat <- exp_mat[gene_id != "V8" & !is.na(gene_id)]
setnames(exp_mat, "gene_id", "probe_id")
log_msg(sprintf("  过滤后: %d probes × %d 样本", nrow(exp_mat), ncol(exp_mat) - 1))

# ============================================================
# 3. probe → gene 映射
# ============================================================
log_msg("=== probe → gene 映射 ===")

# 用 inner join
mapped <- merge(exp_mat, annot_valid, by = "probe_id", all.x = FALSE)

n_total      <- nrow(exp_mat)
n_mapped     <- nrow(mapped)
n_unmapped   <- n_total - n_mapped
log_msg(sprintf("  matrix 中总 probes: %d", n_total))
log_msg(sprintf("  成功映射到 gene: %d (%.1f%%)", n_mapped, 100*n_mapped/n_total))
log_msg(sprintf("  未映射 (control/missing): %d", n_unmapped))
log_msg(sprintf("  映射到的唯一 gene 数: %d", uniqueN(mapped$gene_symbol)))

# 看每个 gene 对应多少 probes 的分布
probe_per_gene <- mapped[, .N, by = gene_symbol]
log_msg("  每基因对应 probe 数分布:")
log_msg(sprintf("    1 probe:     %d 基因", sum(probe_per_gene$N == 1)))
log_msg(sprintf("    2 probes:    %d 基因", sum(probe_per_gene$N == 2)))
log_msg(sprintf("    3 probes:    %d 基因", sum(probe_per_gene$N == 3)))
log_msg(sprintf("    4+ probes:   %d 基因", sum(probe_per_gene$N >= 4)))

# ============================================================
# 4. 多 probe 聚合：取最大方差 probe (业界标准)
# ============================================================
log_msg("=== 多 probe 聚合 (最大方差 probe) ===")

# 提取样本列
sample_cols <- setdiff(colnames(exp_mat), "probe_id")
log_msg(sprintf("  269 样本列已识别"))

# 计算每个 probe 的方差
log_msg("  计算每 probe 方差...")
mapped[, probe_var := apply(.SD, 1, var, na.rm = TRUE), .SDcols = sample_cols]

# 每个基因选方差最大的 probe
log_msg("  选择最大方差 probe...")
selected <- mapped[, .SD[which.max(probe_var)], by = gene_symbol]
log_msg(sprintf("  聚合后唯一基因数: %d", nrow(selected)))

# ============================================================
# 5. 输出
# ============================================================
log_msg("=== 输出文件 ===")

# 5.1 gene-level expression matrix
gene_mat <- selected[, c("gene_symbol", sample_cols), with = FALSE]
setnames(gene_mat, "gene_symbol", "gene")
fwrite(gene_mat, OUT_GENE, sep = "\t", compress = "gzip")
log_msg(sprintf("  ✓ %s (%.1f MB)", basename(OUT_GENE), file.size(OUT_GENE) / 1024^2))

# 5.2 probe→gene 映射记录 (含选择的 probe 信息, 方便回溯)
map_record <- selected[, .(gene_symbol, selected_probe = probe_id, probe_var)]
# 加上每基因的总 probe 数
map_record <- merge(map_record, probe_per_gene, by = "gene_symbol")
setnames(map_record, "N", "n_probes_available")
fwrite(map_record, OUT_MAP, sep = "\t")
log_msg(sprintf("  ✓ %s", basename(OUT_MAP)))

# 5.3 统计摘要
summary_lines <- c(
    sprintf("PACA-AU probe → gene mapping summary"),
    sprintf("Generated: %s", Sys.time()),
    sprintf(""),
    sprintf("Input matrix:    %d probes × %d samples", n_total, length(sample_cols)),
    sprintf("Annotation:      GPL10558 (Illumina HumanHT-12 V4.0, %d probes total)", nrow(annot)),
    sprintf("                 %d probes with gene symbol (%.1f%%)", 
            nrow(annot_valid), 100*nrow(annot_valid)/nrow(annot)),
    sprintf(""),
    sprintf("Mapping result:"),
    sprintf("  probes mapped:       %d / %d (%.1f%%)", n_mapped, n_total, 100*n_mapped/n_total),
    sprintf("  probes unmapped:     %d", n_unmapped),
    sprintf("  unique genes:        %d", uniqueN(mapped$gene_symbol)),
    sprintf(""),
    sprintf("Probes per gene distribution:"),
    sprintf("  1 probe:     %d genes", sum(probe_per_gene$N == 1)),
    sprintf("  2 probes:    %d genes", sum(probe_per_gene$N == 2)),
    sprintf("  3 probes:    %d genes", sum(probe_per_gene$N == 3)),
    sprintf("  4+ probes:   %d genes", sum(probe_per_gene$N >= 4)),
    sprintf(""),
    sprintf("Aggregation strategy: maximum-variance probe per gene"),
    sprintf("Final gene-level matrix: %d genes × %d samples", nrow(selected), length(sample_cols)),
    sprintf(""),
    sprintf("Output files:"),
    sprintf("  %s", basename(OUT_GENE)),
    sprintf("  %s", basename(OUT_MAP)),
    sprintf("  %s", basename(OUT_SUMMARY))
)
writeLines(summary_lines, OUT_SUMMARY)
log_msg(sprintf("  ✓ %s", basename(OUT_SUMMARY)))

# 预览
log_msg("=== 矩阵预览 ===")
log_msg("前 5 行 × 前 4 列:")
print(gene_mat[1:5, 1:4])

# 关键基因健康检查
log_msg("=== 关键基因健康检查 ===")
check_genes <- c("KRAS", "TP53", "SMAD4", "CDKN2A", "GAPDH", "ACTB",
                 "GLIS3", "NEUROD1", "INSM1", "NFIB", "SOX4",  # NRP
                 "GPX4", "SLC7A11", "ACSL4", "FTH1")             # 铁死亡
present <- check_genes %in% gene_mat$gene
log_msg(sprintf("  PDAC 关键基因 + NRP + 铁死亡 (%d 个) 检查:", length(check_genes)))
for (g in check_genes) {
    status <- if (g %in% gene_mat$gene) "✓" else "✗ MISSING"
    log_msg(sprintf("    %-10s %s", g, status))
}
log_msg(sprintf("  共 %d / %d 在矩阵中", sum(present), length(check_genes)))

log_msg("=== 完成 ===")
