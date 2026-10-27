# W3 Day 3 Step 5: Figure 2A - 53 基因的 4 源来源 UpSet 图
suppressPackageStartupMessages({
  library(data.table); library(UpSetR)
})

MIME  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime"
GP_IN <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/01_gene_pool"
FIG   <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"

# 读 NeFRS 53 基因 + 4 源来源信息
nz <- fread(file.path(MIME, "nefrs_v3_enet_coefficients.tsv"))
cat(sprintf("NeFRS 53 genes loaded\n"))

# 看 Sources 列格式
cat("\nFirst 5 Sources entries:\n")
print(head(nz[, .(gene, Sources, n_sources)], 5))

# 解析 Sources 列 (用 | 分隔的多源标记)
parse_sources <- function(s) {
  if (is.na(s)) return(character(0))
  unlist(strsplit(s, "\\|"))
}

# 4 源标准化命名
source_map <- c(
  "FerrDb_v2"                        = "Ferroptosis",
  "Hwang2022_NatGenet"               = "Hwang_NRP",
  "GO_BP_high_specificity"           = "Neural_GO",
  "Literature_review_PDAC_PNI"       = "PDAC_PNI"
)

# 构造 binary membership matrix (53 基因 × 4 源)
upset_df <- data.frame(gene = nz$gene)
for (src_raw in names(source_map)) {
  src_name <- source_map[[src_raw]]
  upset_df[[src_name]] <- as.integer(sapply(nz$Sources, function(s) {
    src_raw %in% parse_sources(s)
  }))
}

cat("\nUpSet membership matrix:\n")
print(head(upset_df, 10))

cat("\nSource counts in 53 NeFRS genes:\n")
print(colSums(upset_df[, -1]))

# 保存
fwrite(upset_df, file.path(FIG, "Fig2A_membership_matrix.tsv"), sep="\t")

# === Figure 2A: UpSet plot ===
outpath <- file.path(FIG, "Fig2A_NeFRS_4sources_UpSet.pdf")
pdf(outpath, width=11, height=7, pointsize=11)

upset(
  upset_df,
  sets = c("Ferroptosis", "Hwang_NRP", "Neural_GO", "PDAC_PNI"),
  order.by = "freq",
  keep.order = TRUE,
  number.angles = 0,
  point.size = 3.5,
  line.size = 1.5,
  mainbar.y.label = "Number of NeFRS genes",
  sets.x.label = "Genes per source",
  text.scale = c(1.4, 1.2, 1.3, 1.0, 1.3, 1.2),  # title, ylab, xlab, ticks, set_size_ylab, set_count
  matrix.color = "#1F4E79",
  main.bar.color = "#0072B2",
  sets.bar.color = "#C00000",
  shade.color = "gray85",
  shade.alpha = 0.4,
  matrix.dot.alpha = 0.5
)

dev.off()
cat(sprintf("\nSaved: %s (%d bytes)\n", outpath, file.info(outpath)$size))

# === 同时输出每个交集的基因详情 (Figure 2A 配套表格) ===
cat("\n=== Intersection details ===\n")

# 用 upset 风格的二进制 vector 标识 4 源 membership
nz$mem_string <- apply(upset_df[, -1], 1, paste, collapse="")
intersection_summary <- nz[, .(genes = paste(gene, collapse=", "),
                                 count = .N,
                                 priority_list = paste(unique(priority), collapse=";")),
                            by = mem_string]
intersection_summary <- intersection_summary[order(-count)]

cat("\nMembership pattern (Ferro|Hwang|GO|PNI) -> count:\n")
print(intersection_summary)
fwrite(intersection_summary,
       file.path(FIG, "Fig2A_intersection_details.tsv"), sep="\t")

cat("\nDone.\n")
