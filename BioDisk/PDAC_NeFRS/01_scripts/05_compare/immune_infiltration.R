# 免疫浸润: ssGSEA 28 免疫细胞 (mmc3) + NeFRS 高低对比 + 相关性
suppressPackageStartupMessages({
  library(data.table); library(GSVA); library(ggplot2); library(pheatmap)
  library(RColorBrewer); library(reshape2)
})
options(bitmapType="cairo")
ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
OUT  <- file.path(ROOT,"02_results/14_immune")
FIG  <- file.path(ROOT,"03_figures/Fig13_immune")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE); dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

cat("=== 免疫浸润 ssGSEA ===\n")

# 1. mmc3 28 免疫细胞 gene set
gs <- fread(file.path(ROOT,"00_raw_data/mmc3.txt"))
setnames(gs, c("gene","cell_type","immunity"))
gene_list <- split(gs$gene, gs$cell_type)
cat(sprintf("mmc3: %d 免疫细胞类型, %d 基因\n", length(gene_list), nrow(gs)))

# 2. NeFRS 系数
nefrs <- fread(file.path(ROOT,"02_results/04_mime/nefrs_v3_enet_coefficients.tsv"))
nefrs_coef <- setNames(nefrs$coef, nefrs$gene)

# 3. 队列 (主用 PACA-AU 训练 + TCGA + E-MTAB 大队列, 多队列验证免疫模式)
cohorts <- list(
  PACA_AU   = list(expr="00_raw_data/bulk/ICGC_PACA_AU/exp_array.PACA-AU.gene.tsv.gz",
                   clin="00_raw_data/bulk/ICGC_PACA_AU/clinical_std.PACA-AU.tsv"),
  TCGA_PAAD = list(expr="00_raw_data/bulk/TCGA_PAAD/TCGA-PAAD.curated.gene.tsv.gz",
                   clin="00_raw_data/bulk/TCGA_PAAD/TCGA-PAAD.curated.clinical_std.tsv"),
  GSE71729  = list(expr="00_raw_data/bulk/GSE71729/GSE71729.gene.tsv.gz",
                   clin="00_raw_data/bulk/GSE71729/GSE71729.clinical_std.tsv")
)

all_diff <- data.table(); all_cor <- data.table()
for (cn in names(cohorts)) {
  cat(sprintf("\n--- %s ---\n", cn))
  ef <- file.path(ROOT, cohorts[[cn]]$expr); cf <- file.path(ROOT, cohorts[[cn]]$clin)
  expr <- as.data.frame(fread(ef)); rownames(expr) <- expr[[1]]; expr[[1]] <- NULL
  expr <- as.matrix(expr); rownames(expr) <- toupper(rownames(expr))
  cat(sprintf("  %d 基因 x %d 样本\n", nrow(expr), ncol(expr)))

  # ssGSEA
  gl <- lapply(gene_list, function(g) intersect(toupper(g), rownames(expr)))
  gl <- gl[sapply(gl, length)>=3]
  param <- ssgseaParam(expr, gl, normalize=TRUE)
  scores <- gsva(param, verbose=FALSE)   # 免疫细胞 x 样本
  fwrite(as.data.frame(scores), file.path(OUT, sprintf("ssgsea_%s.csv",cn)), row.names=TRUE)
  cat(sprintf("  ssGSEA: %d 免疫细胞\n", nrow(scores)))

  # NeFRS risk (z-score × coef)
  ng <- intersect(names(nefrs_coef), rownames(expr))
  em_z <- t(scale(t(expr[ng,,drop=FALSE]))); em_z[is.na(em_z)] <- 0
  risk <- colSums(em_z * nefrs_coef[ng])
  risk <- risk[colnames(scores)]
  grp <- ifelse(risk>=median(risk,na.rm=T),"High","Low")

  # 高低差异 (Wilcoxon)
  for (ic in rownames(scores)) {
    h <- scores[ic, grp=="High"]; l <- scores[ic, grp=="Low"]
    w <- wilcox.test(h,l)
    all_diff <- rbind(all_diff, data.table(cohort=cn, cell=ic,
                      high_mean=mean(h), low_mean=mean(l), diff=mean(h)-mean(l), p=w$p.value))
  }
  # NeFRS-免疫细胞相关
  for (ic in rownames(scores)) {
    ct <- cor.test(risk, scores[ic,], method="spearman")
    all_cor <- rbind(all_cor, data.table(cohort=cn, cell=ic, rho=ct$estimate, p=ct$p.value))
  }
}
all_diff[, padj := p.adjust(p,"BH"), by=cohort]
all_cor[, padj := p.adjust(p,"BH"), by=cohort]
fwrite(all_diff, file.path(OUT,"immune_highlow_diff.tsv"), sep="\t")
fwrite(all_cor, file.path(OUT,"immune_nefrs_cor.tsv"), sep="\t")

# === 汇总: PACA-AU 主队列 ===
cat("\n=== PACA-AU NeFRS 高低免疫差异 (top) ===\n")
pa <- all_diff[cohort=="PACA_AU"][order(diff)]
cat("高 NeFRS 降低的免疫细胞 (diff<0, top):\n")
print(head(pa[,.(cell,diff,padj)],8))
cat("\n高 NeFRS 升高的免疫细胞 (diff>0, top):\n")
print(head(pa[order(-diff),.(cell,diff,padj)],8))

# === 图: NeFRS-免疫相关性 (跨队列一致性) ===
cor_wide <- dcast(all_cor, cell~cohort, value.var="rho")
cor_mat <- as.matrix(cor_wide[,-1]); rownames(cor_mat) <- cor_wide$cell
cor_mat[is.na(cor_mat)] <- 0
ord <- order(rowMeans(cor_mat))
pheatmap(cor_mat[ord,], cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
         breaks=seq(-0.6,0.6,length.out=101),
         main="Fig13A. NeFRS vs immune cell infiltration (Spearman rho)",
         fontsize_row=8, filename=file.path(FIG,"Fig13A_immune_cor_heatmap.pdf"),
         width=5, height=8)
cat("\nFig13A 保存\n")

# 图: PACA-AU 高低差异 bar
pa[, cell := factor(cell, levels=cell[order(diff)])]
pa[, sig := ifelse(padj<0.05,"*","")]
p2 <- ggplot(pa, aes(cell, diff, fill=diff>0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values=c("TRUE"="#C0392B","FALSE"="#2980B9"), guide="none") +
  geom_text(aes(label=sig), hjust=ifelse(pa$diff>0,-0.2,1.2), size=4) +
  geom_hline(yintercept=0) +
  labs(title="Fig13B. Immune infiltration: High vs Low NeFRS (PACA-AU)",
       subtitle="* padj<0.05; red=higher in high-NeFRS, blue=lower",
       x="", y="ssGSEA diff (High - Low NeFRS)") +
  theme_bw(base_size=10) + theme(plot.subtitle=element_text(size=8,color="grey40"))
ggsave(file.path(FIG,"Fig13B_PACA_highlow.pdf"), p2, width=8, height=9)
cat("Fig13B 保存\n")
cat("\n=== 完成 ===\n")
