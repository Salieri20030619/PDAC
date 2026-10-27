# W4 Fig6: NeFRS 广义 signature 单细胞特征
#   Fig6A: NeFRS_v3 score by cell type (展示 stromal-malignant 分布)
#   Fig6B: NeFRS_v3 vs 17 Hwang programs heatmap (展示多通路关联)
# 新叙事: NeFRS 是"神经-应激-基质"综合 signature, 不是 NRP-focused
suppressPackageStartupMessages({
  library(data.table)
})

SC  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"

cat("=== W4 Fig6: NeFRS broad signature single-cell ===\n")

# ============================================================
# Fig6A: NeFRS_v3 score by cell type
# ============================================================
ct <- fread(file.path(SC, "w4_v3_nefrs69_by_celltype.tsv"))
# 按 median 降序 (已经是), 取前 15
ct <- ct[order(-median)]
ct <- ct[n >= 100]  # 至少 100 cells

# 颜色: 按细胞大类着色
# malignant/ductal = 红系, CAF/stroma = 橙系, vascular = 蓝系, neural = 绿系, immune = 灰
ct[, category := fcase(
  Level3 %in% c("Malignant", "Ductal (atypical)", "Ductal", "ADM"), "Epithelial/Malignant",
  Level3 %in% c("myCAF", "CAF"), "CAF",
  Level3 %in% c("Vascular", "Pericyte", "Vascular smooth muscle", "Lymphatic"), "Vascular/Mural",
  Level3 %in% c("Schwann"), "Neural",
  Level3 %in% c("Alpha", "Hormone-negative neuroendocrine"), "Endocrine",
  default = "Other"
)]
cat_colors <- c(
  "Epithelial/Malignant" = "#C0392B",
  "CAF"                   = "#E67E22",
  "Vascular/Mural"        = "#2980B9",
  "Neural"                = "#27AE60",
  "Endocrine"             = "#8E44AD",
  "Other"                 = "#95A5A6"
)
ct[, col := cat_colors[category]]

outA <- file.path(FIG, "Fig6A_NeFRS_by_celltype.pdf")
pdf(outA, width=9, height=7, pointsize=11)
par(mar=c(11, 4.5, 4, 1.5))

bp <- barplot(ct$median, names.arg=NA, col=ct$col, border="gray30",
              ylim=c(0, max(ct$p90)*1.05),
              ylab="NeFRS score (median, sum coef×expr)",
              main="Fig 6A. NeFRS enriched at stromal-malignant interface\n(GSE202051, 224k cells)",
              cex.main=1.15, cex.lab=1.05)

# p90 误差棒 (median → p90)
arrows(bp, ct$median, bp, ct$p90, length=0.03, angle=90, col="gray50", lwd=1)

# x 轴标签 (倾斜)
text(bp, par("usr")[3] - 0.012, labels=ct$Level3, srt=45, adj=1, xpd=TRUE, cex=0.72)

# n 标在柱顶
text(bp, ct$p90 + 0.012, labels=sprintf("n=%d", ct$n), cex=0.6, col="gray40", srt=90, adj=0)

# legend
legend("topright", legend=names(cat_colors), fill=cat_colors,
       border="gray30", bty="n", cex=0.85, title="Cell category")

# 副标题
mtext("Top 15 cell types by NeFRS; immune (Lymphoid/Myeloid) rank lowest (not shown)",
      side=1, line=9.5, cex=0.75, col="gray40")
dev.off()
cat(sprintf("Fig6A saved: %d bytes\n", file.info(outA)$size))

# ============================================================
# Fig6B: NeFRS_v3 vs 17 Hwang programs heatmap
# ============================================================
cor_df <- fread(file.path(SC, "w4_v3_nefrs69_vs_hwang17_correlation.tsv"))
# clean program names
cor_df[, prog_clean := gsub("_UCell$", "", program)]
cor_df[, prog_clean := gsub("mal_state_", "Malig: ", prog_clean)]
cor_df[, prog_clean := gsub("mal_lineage_", "Lineage: ", prog_clean)]
cor_df[, prog_clean := gsub("CAF_", "CAF: ", prog_clean)]
cor_df[, prog_clean := gsub("_", " ", prog_clean)]

# 按 rho_malig 降序
cor_df <- cor_df[order(rho_malig)]  # 升序, barplot 从下往上

# 颜色: 程序类别
cor_df[, ptype := fcase(
  grepl("^Malig:", prog_clean), "Malignant state",
  grepl("^Lineage:", prog_clean), "Malignant lineage",
  grepl("^CAF:", prog_clean), "CAF program",
  default = "Other"
)]
ptype_col <- c("Malignant state"="#C0392B", "Malignant lineage"="#E67E22", "CAF program"="#2980B9")
cor_df[, col := ptype_col[ptype]]

# 高亮 NRP + TNF-NFkB
cor_df[, highlight := grepl("Neural-like progenitor|TNF-NFkB", prog_clean)]

outB <- file.path(FIG, "Fig6B_NeFRS_vs_Hwang17_programs.pdf")
pdf(outB, width=9, height=8, pointsize=11)
par(mar=c(4.5, 14, 4, 2))

bp <- barplot(cor_df$rho_malig, horiz=TRUE, names.arg=NA,
              col=cor_df$col, border=ifelse(cor_df$highlight, "black", "gray50"),
              xlim=c(0, max(cor_df$rho_malig)*1.45),
              xlab="Spearman rho (NeFRS vs program, malignant cells)",
              main="Fig 6B. NeFRS associates with multiple programs\nStress + CAF + Cycling dominate, not NRP alone",
              cex.main=1.1, cex.lab=1.0,
              lwd=ifelse(cor_df$highlight, 1.6, 0.5))

# y 轴程序名
text(par("usr")[1] - 0.005, bp, labels=cor_df$prog_clean, adj=1, xpd=TRUE,
     cex=0.78, font=ifelse(cor_df$highlight, 2, 1))

# rho 值标在柱尾
text(cor_df$rho_malig + 0.005, bp, labels=sprintf("%.3f", cor_df$rho_malig),
     adj=0, cex=0.7, col="gray20",
     font=ifelse(cor_df$highlight, 2, 1))

# NRP 标星
nrp_idx <- which(grepl("Neural-like progenitor", cor_df$prog_clean))
if (length(nrp_idx) > 0) {
  text(cor_df$rho_malig[nrp_idx] + 0.04, bp[nrp_idx], "<- NRP (Hwang focus)",
       adj=0, cex=0.72, col="#C0392B", font=2)
}
tnf_idx <- which(grepl("TNF-NFkB", cor_df$prog_clean))
if (length(tnf_idx) > 0) {
  text(cor_df$rho_malig[tnf_idx] + 0.04, bp[tnf_idx], "<- strongest",
       adj=0, cex=0.72, col="#C0392B", font=2)
}

legend("bottomright", legend=names(ptype_col), fill=ptype_col,
       border="gray30", bty="n", cex=0.85, title="Program type")
dev.off()
cat(sprintf("Fig6B saved: %d bytes\n", file.info(outB)$size))

cat("\n=== Done ===\n")
cat("Fig6A: NeFRS by cell type (stromal-malignant 富集)\n")
cat("Fig6B: NeFRS vs 17 Hwang programs (多通路关联, 应激主导)\n")
