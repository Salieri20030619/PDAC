# W3 Day 3 Step 4: NeFRS 全 53 基因 bar plot, Hub 高亮 (论文 Figure 2C)
suppressPackageStartupMessages({
  library(data.table)
})

MIME <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime"
FIG  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"

nz <- fread(file.path(MIME, "nefrs_v3_enet_coefficients.tsv"))
cat(sprintf("Total Enet %d genes (v3 rebuild): EXTRA\n", nrow(nz)))

# 按 coef 值排序 (waterfall plot 风格)
nz <- nz[order(coef)]   # protective 在上, risk 在下
nz$y <- seq_len(nrow(nz))

# 颜色: P1 红, P2 橙, P3 浅灰
nz$col <- ifelse(nz$priority == "P1_TopPriority", "#C00000",
          ifelse(nz$priority == "P2_HighPriority", "#FF8C00", "#B0B0B0"))

# Hub 基因 (P1 + P2) 加粗显示
nz$is_hub <- nz$priority %in% c("P1_TopPriority", "P2_HighPriority")

n_p1 <- sum(nz$priority == "P1_TopPriority")
n_p2 <- sum(nz$priority == "P2_HighPriority")
n_p3 <- sum(nz$priority == "P3_Standard")
cat(sprintf("P1: %d, P2: %d, P3: %d, total: %d\n", n_p1, n_p2, n_p3, nrow(nz)))

# === Bar plot ===
outpath <- file.path(FIG, "Fig2C_NeFRS_all53_barplot.pdf")
pdf(outpath, width=10, height=14, pointsize=10)
par(mar=c(4.5, 7, 3.5, 1.5))

x_range <- range(nz$coef) * 1.18
plot(NA, xlim=x_range, ylim=c(0.5, nrow(nz)+1.5),
     xlab="NeFRS Coefficient (beta)", ylab="",
     main=sprintf("NeFRS v3 Signature: All %d Genes", nrow(nz)),
     cex.main=1.35, cex.lab=1.1, las=1, yaxt="n")

# 网格
abline(v=0, lty=1, lwd=1.5, col="gray40")
abline(v=seq(-0.15, 0.15, 0.05), lty=3, lwd=0.4, col="gray88")

# 画 bar
for (i in seq_len(nrow(nz))) {
  r <- nz[i, ]
  bar_height <- ifelse(r$is_hub, 0.42, 0.32)  # hub bar 稍粗
  rect(min(0, r$coef), r$y - bar_height,
       max(0, r$coef), r$y + bar_height,
       col=r$col,
       border=ifelse(r$is_hub, "black", "white"),
       lwd=ifelse(r$is_hub, 0.8, 0.3))
}

# y 轴基因名 (hub 加粗)
for (i in seq_len(nrow(nz))) {
  r <- nz[i, ]
  axis(2, at=r$y, labels=r$gene, las=2,
       cex.axis=ifelse(r$is_hub, 0.85, 0.72),
       font=ifelse(r$is_hub, 2, 1),
       col.axis=ifelse(r$is_hub, "black", "gray35"),
       lwd=0, lwd.ticks=0)
}

# 标 ★ 在 P1 hub 旁
for (i in seq_len(nrow(nz))) {
  r <- nz[i, ]
  if (r$priority == "P1_TopPriority") {
    mtext("*", side=2, at=r$y, line=-0.4, las=2, cex=1.4, col="#C00000", font=2)
  }
}

# 在 hub bar 上标 coef
for (i in seq_len(nrow(nz))) {
  r <- nz[i, ]
  if (r$is_hub) {
    label_x <- ifelse(r$coef > 0, r$coef + 0.004, r$coef - 0.004)
    label_adj <- ifelse(r$coef > 0, 0, 1)
    text(label_x, r$y, sprintf("%.3f", r$coef),
         adj=label_adj, cex=0.75, col="gray20", font=2)
  }
}

# Legend
legend("bottomright",
       legend=c(sprintf("P1 TopPriority (* / %d genes, 3 sources)", n_p1),
                sprintf("P2 HighPriority (%d genes, 2 sources)", n_p2),
                sprintf("P3 Standard (%d genes, 1 source)", n_p3)),
       fill=c("#C00000", "#FF8C00", "#B0B0B0"),
       border=c("black", "black", NA),
       bty="n", cex=0.95)

# 方向注释
text(min(x_range)*0.97, nrow(nz)+1.1, "<- Protective",
     cex=1.0, col="#0072B2", adj=0, font=2)
text(max(x_range)*0.97, nrow(nz)+1.1, "Risk ->",
     cex=1.0, col="#C00000", adj=1, font=2)

# 副标题: hub 数量统计
mtext(sprintf("8 Hub genes (* / orange) from 4-source biological prior; %d single-source P3 genes",
              n_p3),
      side=3, line=0.3, cex=0.9, col="gray30")

dev.off()
cat(sprintf("\nSaved: %s (%d bytes)\n", outpath, file.info(outpath)$size))

# 保存完整 53 基因数据 (按 coef 排序)
fwrite(nz[order(-abs(coef)), .(rank, gene, coef, priority, n_sources)],
       file.path(FIG, "Fig2C_NeFRS_all53_data.tsv"), sep="\t")

cat("Done.\n")
