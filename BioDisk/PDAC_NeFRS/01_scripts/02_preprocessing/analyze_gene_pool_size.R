#!/usr/bin/env Rscript
# 分析候选基因池大小, 帮助决定聚焦阈值

suppressPackageStartupMessages({ library(data.table) })

POOL <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/01_gene_pool"

ferrdb <- fread(file.path(POOL, "01_ferrdb_v2_genes.tsv"))
nrp <- fread(file.path(POOL, "02_nrp_signature_genes.tsv"))
neural <- fread(file.path(POOL, "03_neural_dev_genes.tsv"))

# 6 队列共有基因 (用作覆盖率分母)
cohorts <- list(
    "PACA-AU" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/ICGC_PACA_AU/exp_array.PACA-AU.gene.tsv.gz",
    "GSE71729" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/GSE71729/GSE71729.gene.tsv.gz",
    "GSE21501" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/GSE21501/GSE21501.gene.tsv.gz",
    "GSE85916" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/GSE85916/GSE85916.gene.tsv.gz",
    "GSE62452" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/GSE62452/GSE62452.gene.tsv.gz",
    "GSE57495" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/GSE57495/GSE57495.gene.tsv.gz",
    "TCGA_curated" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/TCGA_PAAD/TCGA-PAAD.curated.gene.tsv.gz"
)

all_cohort_genes <- list()
for (name in names(cohorts)) {
    all_cohort_genes[[name]] <- fread(cohorts[[name]], select = 1, header = TRUE)[[1]]
}
common_genes <- Reduce(intersect, all_cohort_genes)

cat("=== A. 三个源基因数 (raw) ===\n")
cat(sprintf("  FerrDb v2: %d\n", nrow(ferrdb)))
cat(sprintf("  NRP signature: %d\n", nrow(nrp)))
cat(sprintf("  神经发育 GO: %d\n", nrow(neural)))

cat("\n=== B. 三源交叉重叠 ===\n")
ferrdb_set <- ferrdb$GeneSymbol
nrp_set <- nrp$GeneSymbol
neural_set <- neural$GeneSymbol

cat(sprintf("  FerrDb ∩ NRP: %d\n", length(intersect(ferrdb_set, nrp_set))))
cat(sprintf("  FerrDb ∩ Neural: %d\n", length(intersect(ferrdb_set, neural_set))))
cat(sprintf("  NRP ∩ Neural: %d\n", length(intersect(nrp_set, neural_set))))
cat(sprintf("  三源都在: %d\n", length(Reduce(intersect, list(ferrdb_set, nrp_set, neural_set)))))

# 最有意思: 同时在 FerrDb 和神经发育里的基因 (= 神经×铁死亡交叉点)
ferrdb_neural <- intersect(ferrdb_set, neural_set)
cat(sprintf("\n=== C. 神经×铁死亡交叉基因 (FerrDb ∩ Neural): %d ===\n", length(ferrdb_neural)))
cat("这些是 NeFRS 最核心的候选基因 (理论上 signature 应该重点选这些):\n")
# 加上 FerrDb class 信息
intersect_dt <- ferrdb[GeneSymbol %in% ferrdb_neural]
cat("class 分布:\n")
print(table(intersect_dt$SubCategory))
cat("\n前 20 个 (含 class):\n")
print(intersect_dt[1:20, .(GeneSymbol, SubCategory)])

cat("\n=== D. 神经发育 GO 不同阈值下的基因数 ===\n")
for (thresh in c(1, 2, 3, 4, 5, 6, 7, 8, 10)) {
    n_total <- sum(neural$n_GO_terms >= thresh)
    n_in_cohorts <- sum(neural[n_GO_terms >= thresh]$GeneSymbol %in% common_genes)
    cat(sprintf("  GO term ≥ %2d : %d 基因 (6 队列共有: %d)\n", thresh, n_total, n_in_cohorts))
}

cat("\n=== E. 不同候选池策略评估 ===\n")
strategies <- list(
    "A: 最保守 (FerrDb + NRP, 无 GO)" = unique(c(ferrdb_set, nrp_set)),
    "B: 加 Neural GO ≥ 7" = unique(c(ferrdb_set, nrp_set, neural[n_GO_terms >= 7]$GeneSymbol)),
    "C: 加 Neural GO ≥ 5" = unique(c(ferrdb_set, nrp_set, neural[n_GO_terms >= 5]$GeneSymbol)),
    "D: 加 Neural GO ≥ 3" = unique(c(ferrdb_set, nrp_set, neural[n_GO_terms >= 3]$GeneSymbol)),
    "E: 加 Neural 全部 GO" = unique(c(ferrdb_set, nrp_set, neural_set))
)

for (name in names(strategies)) {
    genes <- strategies[[name]]
    n_total <- length(genes)
    n_in_cohorts <- sum(genes %in% common_genes)
    # 与训练集样本数比 (267)
    ratio <- n_in_cohorts / 267
    cat(sprintf("  %-40s : %4d 基因 (6 队列共有: %4d, 基因/样本=%.1f)\n",
                name, n_total, n_in_cohorts, ratio))
}

cat("\n=== F. 推荐 ===\n")
cat("LASSO Cox 经验法则: 基因数 ≤ 5-10 × 样本数 = 1300-2600\n")
cat("AIDPS 用 32 个 CPG (consensus prognostic genes), 极保守\n")
cat("AML 模板用 80 基因, 中保守\n")
cat("我们目标: 800-2000 候选基因, LASSO + 单变量 Cox 预筛后落到 ~200-300, 最终 ~30 signature\n")
