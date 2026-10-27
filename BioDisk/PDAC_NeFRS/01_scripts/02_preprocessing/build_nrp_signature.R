#!/usr/bin/env Rscript
# Source 2: Hwang 2022 Nature Genetics NRP (Neural-like progenitor) signature
# 真实版本: 从 Hwang 2022 Sup Table 2 提取的 top 200 cNMF-weighted NRP genes
#
# 来源:
#   Hwang WL et al. Nat Genet 2022;54:1178-1191
#   "Single-nucleus and spatial transcriptome profiling of pancreatic cancer
#    identifies multicellular dynamics associated with neoadjuvant treatment"
#   PMID: 35902743
#   Supplementary Table 2 列 14 "Neural-like progenitor"
#
# 历史说明:
#   原 W1 版本 (61 基因 手工 curated) 已存档至:
#   _archive_w1_nrp_pre_fix/02_nrp_signature_W1_curated.tsv
#   原版基因为通用神经发育 markers, 跟 Hwang Sup Table 2 重叠仅 5 个
#   本次替换 = 用论文官方 top 200 cNMF NRP weighted genes 全量替换

suppressPackageStartupMessages({
  library(data.table)
})

NRP_FILE <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/hwang2022_supp/hwang_NRP_200.txt"
PROGS_FILE <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/hwang2022_supp/hwang2022_programs_clean.tsv"

cat("=== Building NRP signature from Hwang 2022 Sup Table 2 ===\n")

progs <- fread(PROGS_FILE)
nrp_with_rank <- progs[program == "Neural-like progenitor",
                        .(GeneSymbol = gene, rank = rank)]
stopifnot(nrow(nrp_with_rank) == 200)

nrp_signature <- data.table(
  GeneSymbol = nrp_with_rank$GeneSymbol,
  GeneCategory = "NRP",
  SubCategory = sprintf("Hwang2022_NRP_top200_rank%d", nrp_with_rank$rank),
  Source = "Hwang2022_NatGenet_SupTable2"
)

# 保存
out_dir <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/01_gene_pool/"
out_file <- file.path(out_dir, "02_nrp_signature_genes.tsv")
fwrite(nrp_signature, out_file, sep = "\t")

cat(sprintf("✓ Saved: %s (%d genes)\n", out_file, nrow(nrp_signature)))
cat("\nFirst 10 genes by rank (highest cNMF weight first):\n")
print(head(nrp_signature, 10))

cat("\nLast 10 genes by rank:\n")
print(tail(nrp_signature, 10))
