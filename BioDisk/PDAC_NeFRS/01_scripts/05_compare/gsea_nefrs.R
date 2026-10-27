# GSEA: NeFRS 相关基因排序 → clusterProfiler GSEA, 3 队列串行 (跑完释放再下一个)
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db)
  library(enrichplot); library(ggplot2)
})
options(bitmapType="cairo")
ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
OUT  <- file.path(ROOT,"02_results/15_gsea"); FIG <- file.path(ROOT,"03_figures/Fig14_gsea")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE); dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

nefrs <- fread(file.path(ROOT,"02_results/04_mime/nefrs_v3_enet_coefficients.tsv"))
nefrs_coef <- setNames(nefrs$coef, nefrs$gene)

cohorts <- list(
  PACA_AU   = "00_raw_data/bulk/ICGC_PACA_AU/exp_array.PACA-AU.gene.tsv.gz",
  TCGA_PAAD = "00_raw_data/bulk/TCGA_PAAD/TCGA-PAAD.curated.gene.tsv.gz",
  GSE71729  = "00_raw_data/bulk/GSE71729/GSE71729.gene.tsv.gz"
)

# 单队列 GSEA 函数 (跑完返回结果, 局部变量自动释放)
run_gsea_cohort <- function(cn, ef) {
  cat(sprintf("\n========== %s ==========\n", cn))
  expr <- as.data.frame(fread(file.path(ROOT,ef)))
  rownames(expr) <- expr[[1]]; expr[[1]] <- NULL; expr <- as.matrix(expr)
  rownames(expr) <- toupper(rownames(expr))
  cat(sprintf("  %d 基因 x %d 样本\n", nrow(expr), ncol(expr)))

  ng <- intersect(names(nefrs_coef), rownames(expr))
  em_z <- t(scale(t(expr[ng,,drop=FALSE]))); em_z[is.na(em_z)] <- 0
  risk <- colSums(em_z * nefrs_coef[ng])

  cor_vec <- apply(expr, 1, function(g) if(sd(g)==0) NA else cor(g, risk, method="spearman"))
  cor_vec <- cor_vec[!is.na(cor_vec)]

  sym2ent <- bitr(names(cor_vec), "SYMBOL", "ENTREZID", org.Hs.eg.db)
  sym2ent <- sym2ent[!duplicated(sym2ent$SYMBOL),]
  gl <- cor_vec[sym2ent$SYMBOL]; names(gl) <- sym2ent$ENTREZID
  gl <- sort(gl[!is.na(gl)], decreasing=TRUE)
  cat(sprintf("  排序基因: %d\n", length(gl)))

  set.seed(1)
  gse_go <- gseGO(gl, OrgDb=org.Hs.eg.db, ont="BP", minGSSize=10, maxGSSize=500,
                  pvalueCutoff=0.05, verbose=FALSE, eps=0)
  gse_go_s <- setReadable(gse_go, org.Hs.eg.db, "ENTREZID")
  fwrite(as.data.table(gse_go_s@result), file.path(OUT, sprintf("gsea_GOBP_%s.tsv",cn)), sep="\t")
  cat(sprintf("  GO-BP 显著: %d\n", nrow(gse_go@result)))

  gse_kegg <- tryCatch(gseKEGG(gl, organism="hsa", minGSSize=10, maxGSSize=500,
                               pvalueCutoff=0.05, verbose=FALSE, eps=0),
                       error=function(e){cat("  KEGG失败\n");NULL})
  if(!is.null(gse_kegg)){
    gk <- setReadable(gse_kegg, org.Hs.eg.db, "ENTREZID")
    fwrite(as.data.table(gk@result), file.path(OUT, sprintf("gsea_KEGG_%s.tsv",cn)), sep="\t")
    cat(sprintf("  KEGG 显著: %d\n", nrow(gse_kegg@result)))
  }

  go_dt <- as.data.table(gse_go@result)
  cat("  高 NeFRS 上调 top5:\n")
  print(head(go_dt[NES>0][order(-NES),.(Description, NES, p.adjust)],5))
  cat("  高 NeFRS 下调 top5:\n")
  print(head(go_dt[NES<0][order(NES),.(Description, NES, p.adjust)],5))

  go_dt[, cohort:=cn]
  return(go_dt[,.(Description, NES, p.adjust, setSize, cohort)])
}

# 串行跑 3 队列
all_go <- list()
for (cn in names(cohorts)) {
  all_go[[cn]] <- run_gsea_cohort(cn, cohorts[[cn]])
  gc()  # 每队列跑完释放
}
all_go_dt <- rbindlist(all_go)
fwrite(all_go_dt, file.path(OUT,"gsea_GOBP_3cohorts.tsv"), sep="\t")

# 跨队列一致的通路 (3队列都显著且 NES 同向)
cat("\n\n========== 跨队列一致通路 ==========\n")
nes_wide <- dcast(all_go_dt, Description~cohort, value.var="NES")
setDT(nes_wide)
cohort_cols <- intersect(c("PACA_AU","TCGA_PAAD","GSE71729"), colnames(nes_wide))
nes_wide <- nes_wide[rowSums(!is.na(nes_wide[,..cohort_cols]))==length(cohort_cols)]  # 3队列都有
nes_wide[, mean_NES := rowMeans(.SD,na.rm=T), .SDcols=cohort_cols]
nes_wide[, consistent := apply(.SD, 1, function(x) all(sign(x)==sign(x[1]))), .SDcols=cohort_cols]
nes_cons <- nes_wide[consistent==TRUE]
fwrite(nes_cons, file.path(OUT,"gsea_consistent_3cohorts.tsv"), sep="\t")
cat(sprintf("3队列都显著且同向的通路: %d\n", nrow(nes_cons)))
cat("\n上调 (mean NES>0) top10:\n")
print(head(nes_cons[mean_NES>0][order(-mean_NES),.(Description, mean_NES)],10))
cat("\n下调 (mean NES<0) top10:\n")
print(head(nes_cons[mean_NES<0][order(mean_NES),.(Description, mean_NES)],10))

# 图: 跨队列一致 top 通路
top_up <- head(nes_cons[mean_NES>0][order(-mean_NES)],10)
top_dn <- head(nes_cons[mean_NES<0][order(mean_NES)],10)
pdt <- rbind(top_up, top_dn)
pdt[, Description := factor(Description, levels=Description[order(mean_NES)])]
pdt[, dir := ifelse(mean_NES>0,"Up in high-NeFRS","Down in high-NeFRS")]
p <- ggplot(pdt, aes(mean_NES, Description, fill=dir)) +
  geom_col() +
  scale_fill_manual(values=c("Up in high-NeFRS"="#C0392B","Down in high-NeFRS"="#2471A3")) +
  geom_vline(xintercept=0) +
  labs(title="Fig14A. GSEA GO-BP: NeFRS-correlated pathways (consistent in 3 cohorts)",
       subtitle="Mean NES across PACA-AU / TCGA / GSE71729", x="Mean NES", y="", fill="") +
  theme_bw(base_size=9) + theme(legend.position="bottom",
        plot.subtitle=element_text(size=8,color="grey40"))
ggsave(file.path(FIG,"Fig14A_GSEA_GOBP_consistent.pdf"), p, width=10, height=8)
cat("\nFig14A 保存\n=== 完成 ===\n")
