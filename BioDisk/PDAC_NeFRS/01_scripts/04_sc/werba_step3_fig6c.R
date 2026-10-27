suppressPackageStartupMessages({ library(data.table) })
OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/07_werba"
FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"

ct <- fread(file.path(OUT, "werba_nefrs_by_celltype.tsv"))
ct <- ct[order(-median)]

ct[, category := fcase(
  cell_type %in% c("Malignant_Ductal","Acinar"), "Epithelial/Malignant",
  cell_type == "CAF", "CAF",
  cell_type == "Endothelial", "Vascular",
  cell_type == "Endocrine", "Endocrine",
  cell_type %in% c("Tcell","Bcell_Plasma","Myeloid","Mast"), "Immune",
  default="Other"
)]
cat_colors <- c("Epithelial/Malignant"="#C0392B","CAF"="#E67E22",
                "Vascular"="#2980B9","Endocrine"="#8E44AD",
                "Immune"="#7F8C8D","Other"="#BDC3C7")
ct[, col := cat_colors[category]]

outC <- file.path(FIG, "Fig6C_Werba_NeFRS_by_celltype.pdf")
pdf(outC, width=8, height=6, pointsize=11)
par(mar=c(8, 4.5, 4, 1.5))
ylim_lo <- min(ct$median)*1.3; ylim_hi <- max(ct$p90)*1.05
bp <- barplot(ct$median, names.arg=NA, col=ct$col, border="gray30",
              ylim=c(ylim_lo, ylim_hi),
              ylab="NeFRS score (median, sum coef x expr)",
              main="Fig 6C. Werba 2023 independent validation (scRNA, 239k cells)\nMalignant + Vascular high, Immune lowest - consistent with Hwang",
              cex.main=0.92, cex.lab=1.0)
abline(h=0, lty=2, col="gray50")
arrows(bp, ct$median, bp, ct$p90, length=0.03, angle=90, col="gray50", lwd=1)
text(bp, ylim_lo - 0.004, labels=ct$cell_type, srt=45, adj=1, xpd=TRUE, cex=0.78)
text(bp, ct$p90+0.006, labels=sprintf("n=%d",ct$n), cex=0.6, col="gray40", srt=90, adj=0)
legend("topright", legend=names(cat_colors), fill=cat_colors, border="gray30",
       bty="n", cex=0.72, title="Cell category")
dev.off()
cat(sprintf("Fig6C saved: %d bytes\n", file.info(outC)$size))
