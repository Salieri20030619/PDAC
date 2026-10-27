# W3 Day 3 Step 2: Forest Plot (publication quality, v2)
suppressPackageStartupMessages({ library(data.table) })

FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"
df <- fread(file.path(FIG, "Fig3_KM_summary.tsv"))

# 排序
role_order <- c("TRAIN", "SEL", "HO")
df$role_factor <- factor(df$role, levels=role_order)
df <- df[order(df$role_factor, -HR), ]
df$y <- nrow(df):1

df$col <- ifelse(df$role == "TRAIN", "#C00000",
          ifelse(df$role == "HO", "#0072B2", "#666666"))
df$cex <- sqrt(df$n) / sqrt(max(df$n)) * 1.8 + 0.6

fmt_p <- function(p) {
  if (p < 1e-15) return("< 1e-15")
  if (p < 1e-4) return(sprintf("%.1e", p))
  if (p < 0.001) return(sprintf("%.4f", p))
  return(sprintf("%.3f", p))
}
df$p_label <- sapply(df$KM_p, fmt_p)
df$hr_label <- sprintf("%.2f (%.2f-%.2f)", df$HR, df$HR_lo, df$HR_hi)

outpath <- file.path(FIG, "Fig4_Forest_8cohort.pdf")
pdf(outpath, width=11, height=6.5, pointsize=10)

# 左 60%, 右 40%
layout(matrix(c(1, 2), nrow=1), widths=c(6, 4))

# === Panel 1: 左侧文字表 ===
par(mar=c(4.5, 0.5, 3.5, 0.2))
plot.new()
plot.window(xlim=c(0, 10), ylim=c(-0.6, nrow(df) + 1.2))

# 主标题
mtext("NeFRS v2 Prognostic Performance Across 8 PDAC Cohorts",
      side=3, line=1.5, at=10.5, font=2, cex=1.2, adj=0.5)

col_x <- c(cohort=0.1, role=2.6, n=3.8, c=5.3, hr=6.8, p=9.7)

text(col_x["cohort"], nrow(df)+0.7, "Cohort",      adj=0, font=2, cex=0.95)
text(col_x["role"],   nrow(df)+0.7, "Role",        adj=0, font=2, cex=0.95)
text(col_x["n"],      nrow(df)+0.7, "N (events)",  adj=0, font=2, cex=0.95)
text(col_x["c"],      nrow(df)+0.7, "C-index",     adj=0, font=2, cex=0.95)
text(col_x["hr"],     nrow(df)+0.7, "HR (95% CI)", adj=0, font=2, cex=0.95)
text(col_x["p"],      nrow(df)+0.7, "P value",     adj=1, font=2, cex=0.95)
abline(h=nrow(df)+0.4, lwd=0.8, col="gray50")

for (i in seq_len(nrow(df))) {
  r <- df[i, ]; y <- r$y
  text(col_x["cohort"], y, r$cohort,
       adj=0, cex=0.95, col=r$col, font=ifelse(r$role=="HO", 2, 1))
  text(col_x["role"],   y, r$role,
       adj=0, cex=0.88, col=r$col, font=ifelse(r$role=="HO", 2, 1))
  text(col_x["n"],      y, sprintf("%d", r$n),       adj=0, cex=0.92)
  text(col_x["c"],      y, sprintf("%.3f", r$c_index), adj=0, cex=0.92, font=2)
  text(col_x["hr"],     y, r$hr_label,               adj=0, cex=0.92)
  text(col_x["p"],      y, r$p_label,                adj=1, cex=0.92,
       font=ifelse(r$KM_p < 0.05, 2, 1),
       col=ifelse(r$KM_p < 0.05, "#C00000", "#666666"))
}

# === 修复 2: 左下角图例 (3 种角色颜色) ===
legend_y <- -0.3
points(col_x["cohort"]+0.05, legend_y, pch=15, cex=1.5, col="#C00000"); text(col_x["cohort"]+0.3, legend_y, "TRAIN", adj=0, cex=0.9, col="#C00000")
points(2.4, legend_y, pch=15, cex=1.5, col="#666666"); text(2.6, legend_y, "SEL-VAL", adj=0, cex=0.9, col="#666666")
points(4.4, legend_y, pch=15, cex=1.5, col="#0072B2"); text(4.6, legend_y, "HOLDOUT", adj=0, cex=0.9, col="#0072B2", font=2)

# === Panel 2: 右侧 forest ===
par(mar=c(4.5, 0.5, 3.5, 1.2))

x_lim <- c(0.4, 7)
x_breaks <- c(0.5, 1, 2, 4, 6)

plot.new()
plot.window(xlim=log(x_lim), ylim=c(-0.6, nrow(df) + 1.2))

text(log(2), nrow(df)+0.7, "Forest Plot (HR ± 95% CI)", adj=0.5, font=2, cex=0.95)
abline(h=nrow(df)+0.4, lwd=0.8, col="gray50")

# 参考线
abline(v=log(1), lwd=1.5, lty=2, col="gray50")
abline(v=log(c(0.5, 2, 4, 6)), lwd=0.4, lty=3, col="gray85")

# 8 个 HR 点 + CI
for (i in seq_len(nrow(df))) {
  r <- df[i, ]; y <- r$y
  segments(log(r$HR_lo), y, log(r$HR_hi), y, lwd=2.2, col=r$col)
  segments(log(r$HR_lo), y-0.15, log(r$HR_lo), y+0.15, lwd=2.2, col=r$col)
  segments(log(r$HR_hi), y-0.15, log(r$HR_hi), y+0.15, lwd=2.2, col=r$col)
  points(log(r$HR), y, pch=15, cex=r$cex, col=r$col)
}

# X 轴
axis(1, at=log(x_breaks), labels=x_breaks, cex.axis=0.9, padj=-0.5)
mtext("Hazard Ratio (log scale)", side=1, line=2.0, cex=0.95)

# === 修复 1: Protective/Risky 标注放在图内空白区 (y=-0.3, 在数据下方但在 X 轴上方) ===
text(log(0.65), -0.3, "<- Protective", cex=0.85, col="#0072B2", adj=0.5)
text(log(3.5),  -0.3, "Risky ->",      cex=0.85, col="#C00000", adj=0.5)

dev.off()

cat(sprintf("Saved: %s (%d bytes)\n", outpath, file.info(outpath)$size))
