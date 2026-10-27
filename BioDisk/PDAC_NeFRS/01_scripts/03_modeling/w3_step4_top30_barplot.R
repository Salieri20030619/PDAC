# W3 Day 3 Step 4: NeFRS Top 30 基因 bar plot (论文 Figure 2C)
suppressPackageStartupMessages({
  library(data.table)
})

MIME <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime"
FIG  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"

# 读 53 基因系数
nz <- fread(file.path(MIME, "nefrs_v3_enet_coefficients.tsv"))
cat(sprintf("Total Enet 53 genes loaded: %d\n", nrow(nz)))

# 按 |coef| 排序, 取 Top 30
nz_top <- nz[order(-abs(coef))][1:30]
nz_top[, abs_coef := abs(coef)]

# 颜色按 priority 分
nz_top$col <- ifelse(nz_top$priority == "P1_TopPriority", "#C00000",
              ifelse(nz_top$priority == "P2_HighPriority", "#FFA500", "#888888"))

# 方向标记
nz_top$direction <- ifelse(nz_top$coef > 0, "Risk", "Protective")

# 反转使最大的在顶部
nz_top <- nz_top[order(abs_coef)]
nz_top$y <- seq_len(nrow(nz_top))

cat("\nTop 30 NeFRS genes:\n")
print(nz_top[order(-abs_coef), .(rank, gene, coef, priority, direction)])

# === Bar plot ===
outpath <- file.path(FIG, "Fig2C_NeFRS_Top30_barplot.pdf")
pdf(outpath, width=9, height=10, pointsize=10)
par(mar=c(4.5, 7, 3.5, 1.5))

# x range
x_range <- range(nz_top$coef) * 1.15
plot(NA, xlim=x_range, ylim=c(0.5, nrow(nz_top)+0.5),
     xlab="NeFRS Coefficient (beta)", ylab="",
     main="NeFRS v2 Top 30 Contributing Genes",
     cex.main=1.3, cex.lab=1.1, las=1, yaxt="n")

# 网格
abline(v=0, lty=1, lwd=1.5, col="gray40")
abline(v=seq(-0.15, 0.15, 0.05), lty=3, lwd=0.4, col="gray85")

# 画 bar
for (i in seq_len(nrow(nz_top))) {
  r <- nz_top[i, ]
  rect(min(0, r$coef), r$y - 0.35,
       max(0, r$coef), r$y + 0.35,
       col=r$col, border="white", lwd=0.5)
}

# y 轴基因名
axis(2, at=nz_top$y, labels=nz_top$gene, las=2, cex.axis=0.85,
     col.axis="black", lwd=0, lwd.ticks=0)

# 标 * 在 P1 hub 基因前
for (i in seq_len(nrow(nz_top))) {
  r <- nz_top[i, ]
  if (r$priority == "P1_TopPriority") {
    mtext("*", side=2, at=r$y, line=-0.3, las=2, cex=1.3, col="#C00000", font=2)
  }
}

# 在每个 bar 上标 coef 数值
for (i in seq_len(nrow(nz_top))) {
  r <- nz_top[i, ]
  label_x <- ifelse(r$coef > 0, r$coef + 0.005, r$coef - 0.005)
  label_adj <- ifelse(r$coef > 0, 0, 1)
  text(label_x, r$y, sprintf("%.3f", r$coef),
       adj=label_adj, cex=0.7, col="gray30")
}

# Legend
legend("bottomright",
       legend=c("P1 TopPriority (3 sources, *)",
                "P2 HighPriority (2 sources)",
                "P3 Standard (1 source)"),
       fill=c("#C00000", "#FFA500", "#888888"),
       border=NA, bty="n", cex=0.95)

# 方向注释 (放图内顶部, 避免和坐标重叠)
text(min(x_range)*0.95, nrow(nz_top)+0.3, "<- Protective", cex=0.9, col="#0072B2", adj=0)
text(max(x_range)*0.95, nrow(nz_top)+0.3, "Risk ->",       cex=0.9, col="#C00000", adj=1)

dev.off()

cat(sprintf("\nSaved: %s (%d bytes)\n", outpath, file.info(outpath)$size))

# 同时保存数据表
fwrite(nz_top[order(-abs_coef), .(rank, gene, coef, priority, n_sources, direction)],
       file.path(FIG, "Fig2C_NeFRS_Top30_data.tsv"), sep="\t")

cat("\nDone.\n")
