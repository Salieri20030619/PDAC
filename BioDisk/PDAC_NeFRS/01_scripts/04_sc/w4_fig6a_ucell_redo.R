# Fig6A 重做: UCell 评分 (去深度), 替换受深度污染的 sum-coef 版本
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
SC  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"

dt <- fread(file.path(SC,"w4_nefrs_ucell_allcells.tsv"))
cat(sprintf("细胞数: %d\n", nrow(dt)))

# 按 Level1 汇总 (UCell)
by1 <- dt[, .(n=.N, median=median(nefrs_ucell), mean=mean(nefrs_ucell),
              p25=quantile(nefrs_ucell,.25), p75=quantile(nefrs_ucell,.75)),
          by=Level1][order(-median)]
cat("\n=== UCell NeFRS by Level1 ===\n"); print(by1)

# 细胞类别上色
by1[, category := fcase(
  grepl("malignant", Level1, ignore.case=TRUE), "Epithelial (malignant)",
  grepl("non-malignant", Level1, ignore.case=TRUE), "Epithelial (normal)",
  grepl("fibroblast", Level1, ignore.case=TRUE), "CAF",
  grepl("Endothelial|Vascular|Pericyte", Level1), "Vascular/Mural",
  grepl("Schwann|neuron", Level1, ignore.case=TRUE), "Neural",
  grepl("Endocrine", Level1), "Endocrine",
  grepl("Lymphoid|Myeloid", Level1), "Immune",
  default="Other")]
cat_cols <- c("Epithelial (malignant)"="#C0392B","Epithelial (normal)"="#E74C3C",
              "CAF"="#E67E22","Vascular/Mural"="#2980B9","Neural"="#27AE60",
              "Endocrine"="#8E44AD","Immune"="#7F8C8D","Other"="#BDC3C7")
by1[, Level1 := factor(Level1, levels=rev(by1$Level1))]

p <- ggplot(by1, aes(x=Level1, y=median, fill=category)) +
  geom_col(width=0.7, color="gray30") +
  geom_errorbar(aes(ymin=p25, ymax=p75), width=0.3, color="gray50") +
  geom_text(aes(label=paste0("n=",n)), hjust=-0.1, size=2.6, color="gray40") +
  coord_flip() +
  scale_fill_manual(values=cat_cols) +
  geom_hline(yintercept=0, lty=2, color="gray50") +
  labs(title="Fig 6A. NeFRS (UCell) by cell type — GSE202051 snRNA, 224k cells",
       subtitle="UCell rank-based scoring (depth-robust, rho=0.17 vs sum-coef 0.42); malignant high, immune lowest",
       y="NeFRS UCell score (median, IQR)", x="", fill="Cell category") +
  theme_bw(base_size=11) +
  theme(plot.title=element_text(size=10.5), plot.subtitle=element_text(size=8, color="gray40"))
ggsave(file.path(FIG,"Fig6A_NeFRS_by_celltype.pdf"), p, width=9.5, height=6)
cat("\nFig6A (UCell 版) 保存\n")

fwrite(by1, file.path(SC,"w4_fig6a_ucell_by_celltype.tsv"), sep="\t")
