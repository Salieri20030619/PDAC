#!/usr/bin/env Rscript
# 探索 FerrDb v2 (ferrdbhs from ferroviz)

suppressPackageStartupMessages({
    library(ferroviz)
    library(data.table)
})

data(ferrdbhs)
setDT(ferrdbhs)

cat("=== 1. 基本统计 ===\n")
cat(sprintf("总行数: %d\n", nrow(ferrdbhs)))
cat(sprintf("唯一基因数: %d\n", uniqueN(ferrdbhs$hs.gene)))

cat("\n=== 2. class 分布 ===\n")
print(table(ferrdbhs$class, useNA = "ifany"))

cat("\n=== 3. 是否有重复基因 (一个基因多 class)? ===\n")
dup_count <- ferrdbhs[, .N, by = hs.gene][N > 1]
cat(sprintf("出现多次的基因数: %d\n", nrow(dup_count)))
if (nrow(dup_count) > 0) {
    cat("前 10 个多 class 基因:\n")
    examples <- ferrdbhs[hs.gene %in% head(dup_count$hs.gene, 10)]
    print(examples[order(hs.gene)])
}

cat("\n=== 4. 关键基因健康检查 (是否在数据库里?) ===\n")
key_ferroptosis <- c("GPX4", "SLC7A11", "ACSL4", "FTH1", "FTL", "TFRC", 
                      "NCOA4", "VDAC2", "VDAC3", "AIFM2", "FSP1", "NRF2",
                      "KEAP1", "NFE2L2", "HMOX1", "PTGS2", "ALOX15")
for (g in key_ferroptosis) {
    in_db <- g %in% ferrdbhs$hs.gene
    if (in_db) {
        cls <- ferrdbhs[hs.gene == g, paste(class, collapse=", ")]
        cat(sprintf("  %-10s ✓ (%s)\n", g, cls))
    } else {
        cat(sprintf("  %-10s ✗\n", g))
    }
}

cat("\n=== 5. 是否覆盖 PACA-AU 训练表达矩阵 ===\n")
expr_file <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/ICGC_PACA_AU/exp_array.PACA-AU.gene.tsv.gz"
paca_genes <- fread(expr_file, select = 1, header = TRUE)
paca_gene_list <- paca_genes[[1]]
cat(sprintf("PACA-AU 训练基因数: %d\n", length(paca_gene_list)))

unique_ferrdb_genes <- unique(ferrdbhs$hs.gene)
in_paca <- unique_ferrdb_genes %in% paca_gene_list
cat(sprintf("FerrDb 基因在 PACA-AU 中: %d / %d (%.1f%%)\n",
            sum(in_paca), length(unique_ferrdb_genes), 100 * mean(in_paca)))

# 不在 PACA-AU 中的基因是哪些类型?
missing_genes <- unique_ferrdb_genes[!in_paca]
cat(sprintf("\n不在 PACA-AU 中的 %d 个 FerrDb 基因前 30 个:\n", length(missing_genes)))
print(head(missing_genes, 30))

# 这些缺失基因有多少是 lncRNA / 假基因 / circRNA?
# 简单判断: 名字以 LINC / RP11 / AC0 / AL0 / MIR 开头
non_coding_pattern <- grepl("^(LINC|RP11|AC0|AL0|MIR|SNHG|HOTAIR|XIST|LOC)", missing_genes)
cat(sprintf("\n推断为 lncRNA/伪基因/未定基因: %d / %d (%.1f%%)\n",
            sum(non_coding_pattern), length(missing_genes),
            100 * mean(non_coding_pattern)))

cat("\n=== 6. 跨 6 个队列的覆盖情况 ===\n")
cohorts <- list(
    "PACA-AU" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/ICGC_PACA_AU/exp_array.PACA-AU.gene.tsv.gz",
    "GSE71729" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/GSE71729/GSE71729.gene.tsv.gz",
    "GSE21501" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/GSE21501/GSE21501.gene.tsv.gz",
    "GSE85916" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/GSE85916/GSE85916.gene.tsv.gz",
    "GSE62452" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/GSE62452/GSE62452.gene.tsv.gz",
    "GSE57495" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/GSE57495/GSE57495.gene.tsv.gz",
    "TCGA_curated" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/TCGA_PAAD/TCGA-PAAD.curated.gene.tsv.gz"
)

for (name in names(cohorts)) {
    g <- fread(cohorts[[name]], select = 1, header = TRUE)[[1]]
    in_cohort <- unique_ferrdb_genes %in% g
    cat(sprintf("  %-15s: %d 基因  →  FerrDb 覆盖 %d / %d (%.1f%%)\n",
                name, length(g), sum(in_cohort), length(unique_ferrdb_genes),
                100 * mean(in_cohort)))
}

# 7. 6 个 bulk 队列共有 FerrDb 基因数
cat("\n=== 7. 6 个 bulk 队列共有 FerrDb 基因 (LASSO 用得到的) ===\n")
all_cohort_genes <- list()
for (name in names(cohorts)) {
    all_cohort_genes[[name]] <- fread(cohorts[[name]], select = 1, header = TRUE)[[1]]
}
common_genes <- Reduce(intersect, all_cohort_genes)
cat(sprintf("6 队列共有基因数: %d\n", length(common_genes)))

ferrdb_in_common <- intersect(unique_ferrdb_genes, common_genes)
cat(sprintf("FerrDb 在 6 队列共有基因中: %d\n", length(ferrdb_in_common)))

cat("\n=== 8. 保存结果 ===\n")
# 保存 FerrDb v2 完整列表 (含 class)
out_dir <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/01_gene_pool/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 标准化为我们的基因集 data.frame 格式 (复用 AML 模板结构)
ferrdb_std <- ferrdbhs[, .(
    GeneSymbol = hs.gene,
    GeneCategory = "Ferroptosis",
    SubCategory = class,
    Source = "FerrDb_v2"
)]

# 一基因多 class 的处理: 合并到一行
ferrdb_dedup <- ferrdb_std[, .(
    GeneCategory = "Ferroptosis",
    SubCategory = paste(unique(SubCategory), collapse = "|"),
    Source = "FerrDb_v2"
), by = GeneSymbol]

cat(sprintf("去重后唯一基因: %d\n", nrow(ferrdb_dedup)))
cat("class 合并后分布 (前 15 种):\n")
print(head(sort(table(ferrdb_dedup$SubCategory), decreasing=TRUE), 15))

# 输出
fwrite(ferrdb_dedup, file.path(out_dir, "01_ferrdb_v2_genes.tsv"), sep = "\t")
fwrite(data.table(gene = ferrdb_in_common), 
       file.path(out_dir, "01_ferrdb_v2_in_all_cohorts.tsv"), sep = "\t")
cat(sprintf("\n✓ 已保存: %s\n", file.path(out_dir, "01_ferrdb_v2_genes.tsv")))
cat(sprintf("✓ 已保存: %s\n", file.path(out_dir, "01_ferrdb_v2_in_all_cohorts.tsv")))
