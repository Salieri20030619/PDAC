# 药敏: NeFRS 高低组 IC50 差异 (3队列) → 跨队列一致的敏感/耐药药
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})
options(bitmapType="cairo")
ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
DRUG <- file.path(ROOT,"02_results/16_drug"); FIG <- file.path(ROOT,"03_figures/Fig15_drug")
dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

nefrs <- fread(file.path(ROOT,"02_results/04_mime/nefrs_v3_enet_coefficients.tsv"))
nefrs_coef <- setNames(nefrs$coef, nefrs$gene)

cohorts <- list(
  PACA_AU   = list(drug="PACA_AU_DrugPred.csv",  expr="00_raw_data/bulk/ICGC_PACA_AU/exp_array.PACA-AU.gene.tsv.gz"),
  TCGA_PAAD = list(drug="TCGA_PAAD_DrugPred.csv",expr="00_raw_data/bulk/TCGA_PAAD/TCGA-PAAD.curated.gene.tsv.gz"),
  GSE71729  = list(drug="GSE71729_DrugPred.csv", expr="00_raw_data/bulk/GSE71729/GSE71729.gene.tsv.gz")
)

all_diff <- data.table()
for (cn in names(cohorts)) {
  cat(sprintf("\n--- %s ---\n", cn))
  # 预测 IC50
  ic50 <- as.data.frame(fread(file.path(DRUG, cohorts[[cn]]$drug)))
  rownames(ic50) <- ic50[[1]]; ic50[[1]] <- NULL; ic50 <- as.matrix(ic50)
  # NeFRS risk
  expr <- as.data.frame(fread(file.path(ROOT, cohorts[[cn]]$expr)))
  rownames(expr) <- expr[[1]]; expr[[1]] <- NULL; expr <- as.matrix(expr); rownames(expr) <- toupper(rownames(expr))
  ng <- intersect(names(nefrs_coef), rownames(expr))
  em_z <- t(scale(t(expr[ng,,drop=FALSE]))); em_z[is.na(em_z)] <- 0
  risk <- colSums(em_z * nefrs_coef[ng])
  # 对齐
  common <- intersect(rownames(ic50), names(risk))
  ic50 <- ic50[common,]; risk_c <- risk[common]
  grp <- ifelse(risk_c>=median(risk_c),"High","Low")
  cat(sprintf("  %d 样本 (High %d/Low %d), %d 药\n", length(common), sum(grp=="High"), sum(grp=="Low"), ncol(ic50)))

  # 每药 高低组 IC50 差异 (Wilcoxon); IC50 低=敏感
  for (d in colnames(ic50)) {
    h <- ic50[grp=="High",d]; l <- ic50[grp=="Low",d]
    w <- wilcox.test(h,l)
    # log2 fold change (High/Low), 负=高NeFRS更敏感(IC50低)
    log2fc <- log2(median(h)/median(l))
    all_diff <- rbind(all_diff, data.table(cohort=cn, drug=d,
                      high_med=median(h), low_med=median(l), log2fc=log2fc, p=w$p.value))
  }
}
all_diff[, padj := p.adjust(p,"BH"), by=cohort]
fwrite(all_diff, file.path(DRUG,"drug_highlow_diff.tsv"), sep="\t")

# 跨队列一致: 3队列都显著(padj<0.05)且 log2fc 同向
cat("\n=== 跨队列一致的差异药物 ===\n")
fc_wide <- dcast(all_diff, drug~cohort, value.var="log2fc"); setDT(fc_wide)
p_wide  <- dcast(all_diff, drug~cohort, value.var="padj"); setDT(p_wide)
cc <- intersect(c("PACA_AU","TCGA_PAAD","GSE71729"), colnames(fc_wide))
fc_wide <- fc_wide[complete.cases(fc_wide[,..cc])]
fc_wide[, mean_log2fc := rowMeans(.SD), .SDcols=cc]
fc_wide[, consistent := apply(.SD,1,function(x) all(sign(x)==sign(x[1]))), .SDcols=cc]
# 合并 padj (3队列都<0.05)
psig <- p_wide[, .(drug, allsig=(get(cc[1])<0.05 & get(cc[2])<0.05 & get(cc[3])<0.05))]
fc_wide <- merge(fc_wide, psig, by="drug")
robust <- fc_wide[consistent==TRUE & allsig==TRUE]
fwrite(robust, file.path(DRUG,"drug_consistent.tsv"), sep="\t")

cat(sprintf("3队列都显著且同向: %d 药\n", nrow(robust)))
cat("\n高 NeFRS 更敏感 (mean log2fc<0, IC50低, top10):\n")
print(head(robust[mean_log2fc<0][order(mean_log2fc),.(drug,mean_log2fc)],10))
cat("\n高 NeFRS 更耐药 (mean log2fc>0, top10):\n")
print(head(robust[mean_log2fc>0][order(-mean_log2fc),.(drug,mean_log2fc)],10))

# 图: 跨队列一致 top 敏感/耐药药
top_sens <- head(robust[mean_log2fc<0][order(mean_log2fc)],15)
top_res  <- head(robust[mean_log2fc>0][order(-mean_log2fc)],15)
pdt <- rbind(top_sens, top_res)
pdt[, drug_clean := gsub("_[0-9]+$","",drug)]
pdt[, drug_clean := factor(drug_clean, levels=drug_clean[order(mean_log2fc)])]
pdt[, dir := ifelse(mean_log2fc<0,"High-NeFRS sensitive","High-NeFRS resistant")]
p <- ggplot(pdt, aes(mean_log2fc, drug_clean, fill=dir)) +
  geom_col() +
  scale_fill_manual(values=c("High-NeFRS sensitive"="#C0392B","High-NeFRS resistant"="#2471A3")) +
  geom_vline(xintercept=0) +
  labs(title="Fig15. NeFRS-associated drug sensitivity (GDSC2, consistent in 3 cohorts)",
       subtitle="Negative = lower IC50 in high-NeFRS (more sensitive)",
       x="Mean log2(IC50 High/Low)", y="", fill="") +
  theme_bw(base_size=9) + theme(legend.position="bottom",
        plot.subtitle=element_text(size=8,color="grey40"))
ggsave(file.path(FIG,"Fig15A_drug_sensitivity.pdf"), p, width=9, height=9)
cat("\nFig15A 保存\n=== 完成 ===\n")
