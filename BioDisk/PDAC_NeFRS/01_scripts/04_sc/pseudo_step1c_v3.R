# B v3: 恶性 NeFRS 高低 DE (presto, 快)
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(presto); library(ggplot2); library(ggrepel)
})
SC  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/09_pseudotime"
FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures/Fig8_pseudotime"
cat("=== B v3: presto DE ===\n"); t0 <- Sys.time()

X <- readRDS(file.path(SC,"w4_X_sparse_full.rds"))
annot <- as.data.table(readRDS(file.path(SC,"w4_annot_with_ucell.rds")))
mal <- annot[Level2=="Malignant"]$barcode
idx <- match(mal, rownames(X)); idx <- idx[!is.na(idx)]
Xm <- X[idx,]; am <- annot[match(rownames(Xm),barcode)]; rm(X); gc()

nefrs <- am$NeFRS_all
q33<-quantile(nefrs,1/3); q67<-quantile(nefrs,2/3)
grp <- ifelse(nefrs>=q67,"High",ifelse(nefrs<=q33,"Low","Mid"))
keep <- grp %in% c("High","Low")
Xk <- Xm[keep,]; gk <- grp[keep]
cat(sprintf("High %d vs Low %d\n", sum(gk=="High"), sum(gk=="Low")))

# presto: genes x cells
cat("[presto wilcoxauc]...\n")
res <- wilcoxauc(t(Xk), gk)   # 需要 genes x cells
res <- as.data.table(res)[group=="High"][order(-logFC)]
res[, padj := p.adjust(pval,"BH")]
cat(sprintf("DE 显著 padj<0.05 & |logFC|>0.25: %d\n",
            sum(res$padj<0.05 & abs(res$logFC)>0.25, na.rm=TRUE)))
cat("\nTop 15 上调 (High):\n")
print(head(res[logFC>0 & padj<0.05, .(feature,logFC,auc,padj,pct_in,pct_out)],15))
cat("\nTop 15 下调:\n")
print(head(res[logFC<0 & padj<0.05][order(logFC), .(feature,logFC,padj)],15))
fwrite(res, file.path(OUT,"malignant_NeFRS_hilo_DE.tsv"), sep="\t")

# 通路 marker
cat("\n[通路] High vs Low:\n")
progs <- list(Stress_TNF=c("NFKB1","JUN","FOS","ATF3","EGR1","JUNB"),
              EMT=c("VIM","ZEB1","SNAI2","FN1","TGFBI"),
              Cycling=c("MKI67","TOP2A","CCNB1","CDK1","PCNA"),
              Neural=c("NRP1","NCAM1","STMN2","TUBB3","GAP43"))
Xhi<-Xm[grp=="High",]; Xlo<-Xm[grp=="Low",]
for(p in names(progs)){
  gs<-intersect(progs[[p]],colnames(Xm)); if(length(gs)==0)next
  cat(sprintf("  %-12s High %.3f vs Low %.3f (%+.3f)\n", p,
      mean(Matrix::colMeans(Xhi[,gs,drop=FALSE])),
      mean(Matrix::colMeans(Xlo[,gs,drop=FALSE])),
      mean(Matrix::colMeans(Xhi[,gs,drop=FALSE]))-mean(Matrix::colMeans(Xlo[,gs,drop=FALSE]))))
}

# 火山图
res$sig<-"ns"
res$sig[res$padj<0.05 & res$logFC>0.5]<-"Up(High)"
res$sig[res$padj<0.05 & res$logFC< -0.5]<-"Down"
lab<-rbind(head(res[sig=="Up(High)"],15),head(res[sig=="Down"][order(logFC)],10))
res$nlp<- -log10(res$padj+1e-300)
pp<-ggplot(res,aes(logFC,nlp,color=sig))+geom_point(size=0.6,alpha=0.5)+
  scale_color_manual(values=c("ns"="gray70","Up(High)"="#C0392B","Down"="#2980B9"))+
  geom_text_repel(data=lab,aes(label=feature),size=2.5,max.overlaps=30,color="black")+
  labs(title="Fig 8D. DE: High vs Low NeFRS malignant cells",x="logFC(High/Low)",y="-log10(adj p)")+
  theme_bw(base_size=12)
ggsave(file.path(FIG,"Fig8D_malignant_hilo_volcano.pdf"),pp,width=9,height=7)
cat(sprintf("\n=== B 完成 %.1f min ===\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
