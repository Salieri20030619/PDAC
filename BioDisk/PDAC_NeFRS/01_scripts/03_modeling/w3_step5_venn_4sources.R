# W3 Day 3 Step 5 (方案B): 4 圈 Venn 图 - NeFRS 53 基因 4 源来源
suppressPackageStartupMessages({
  library(data.table); library(VennDiagram); library(grid)
})

MIME <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime"
FIG  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"

cat("=== W3 Step 5: 4-source Venn diagram ===\n")

# 读 NeFRS 53 基因 + 4 源
nz <- fread(file.path(MIME, "nefrs_v3_enet_coefficients.tsv"))
cat(sprintf("Total Enet 53 genes loaded\n"))

# 4 源标准化命名
source_map <- c(
  "FerrDb_v2"                  = "Ferroptosis",
  "Hwang2022_NatGenet_SupTable2"         = "Hwang_NRP",
  "GO_BP_high_specificity"     = "Neural_GO",
  "Literature_review_PDAC_PNI" = "PDAC_PNI"
)

# 解析 Sources, 构造 4 个集合 (gene list per source)
parse_sources <- function(s) {
  if (is.na(s)) return(character(0))
  unlist(strsplit(s, "\\|"))
}

sets <- list()
for (src_raw in names(source_map)) {
  src_name <- source_map[[src_raw]]
  sets[[src_name]] <- nz$gene[sapply(nz$Sources, function(s) {
    src_raw %in% parse_sources(s)
  })]
}

cat("\n=== Set sizes ===\n")
for (n in names(sets)) {
  cat(sprintf("  %-15s: %d genes\n", n, length(sets[[n]])))
}

# 计算所有 4-set 交集大小 (sanity check 用)
cat("\n=== Pairwise intersections ===\n")
for (i in 1:3) {
  for (j in (i+1):4) {
    n1 <- names(sets)[i]; n2 <- names(sets)[j]
    n_inter <- length(intersect(sets[[i]], sets[[j]]))
    cat(sprintf("  %s ∩ %s: %d\n", n1, n2, n_inter))
  }
}

cat(sprintf("\n=== Triple intersections ===\n"))
for (i in 1:2) {
  for (j in (i+1):3) {
    for (k in (j+1):4) {
      n_inter <- length(Reduce(intersect, sets[c(i,j,k)]))
      cat(sprintf("  %s ∩ %s ∩ %s: %d\n",
                  names(sets)[i], names(sets)[j], names(sets)[k], n_inter))
    }
  }
}
cat(sprintf("\n  All 4: %d\n", length(Reduce(intersect, sets))))

# === 绘制 4-set Venn 图 ===
outpath <- file.path(FIG, "Fig2A_NeFRS_4sources_Venn.pdf")

# VennDiagram 4-set 默认是椭圆型 ("ellipse" 布局),
# 颜色: 4 个柔和但有区分度的色调
venn_colors <- c("#E63946",   # Ferroptosis 红
                 "#F4A261",   # Hwang_NRP 橙
                 "#2A9D8F",   # Neural_GO 青绿
                 "#264653")   # PDAC_PNI 深蓝

# 关键: VennDiagram::venn.diagram 写文件,不要先 pdf()
# 否则会有空白页问题
venn.plot <- venn.diagram(
  x = sets,
  category.names = c("Ferroptosis\n(FerrDb v2)",
                     "Hwang NRP\n(Nature Genetics 2022)",
                     "Neural Development\n(GO BP)",
                     "PDAC PNI\n(Literature)"),
  filename = NULL,
  
  # 整体外观
  lwd = 2.5,
  lty = 1,
  fill = venn_colors,
  alpha = 0.45,
  
  # 集合标签
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.col = venn_colors,
  cat.default.pos = "outer",
  cat.dist = c(0.22, 0.22, 0.12, 0.12),
  cat.pos = c(-15, 15, -15, 15),
  
  # 交集内数字
  cex = 1.4,
  fontface = "bold",
  
  # 标题
  main = sprintf("NeFRS v3 (%d genes): 4-Source Biological Prior", nrow(nz)),
  main.cex = 1.6,
  main.fontface = "bold",
  main.pos = c(0.5, 1.03),
  
  # 副标题: hub 提示
  sub = sprintf("Total %d genes; P1 hub: %s; P2 hub: %s", nrow(nz), paste(nz[priority=="P1_TopPriority"]$gene, collapse=", "), paste(nz[priority=="P2_HighPriority"]$gene, collapse=", ")),
  sub.cex = 0.95,
  sub.fontface = "italic",
  sub.col = "gray30",
  sub.pos = c(0.5, -0.02),
  
  margin = 0.12,
  
  # 关闭日志输出
  disable.logging = TRUE
)

# 用 cairo_pdf 保证 Unicode + 渲染质量
cairo_pdf(outpath, width=11, height=10, pointsize=11)
grid.draw(venn.plot)
dev.off()

cat(sprintf("\nSaved: %s (%d bytes)\n", outpath, file.info(outpath)$size))

# === 同时输出 Hub 标注表 (论文 Table 2 用) ===
hub_df <- nz[priority %in% c("P1_TopPriority", "P2_HighPriority"),
             .(gene, coef=round(coef, 4), priority,
               n_sources, Sources, rank)]
setorder(hub_df, rank)
fwrite(hub_df, file.path(FIG, "Fig2A_Hub_genes_annotation.tsv"), sep="\t")

cat("\n=== Hub gene annotation (saved to Fig2A_Hub_genes_annotation.tsv) ===\n")
print(hub_df)

cat("\nDone.\n")
