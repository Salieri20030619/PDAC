#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; GSE<-file.path(ROOT,"02_results/15_gsea"); PFIG<-file.path(ROOT,"03_figures/panels")
KG<-fread(file.path(GSE,"gsea_KEGG_PACA_AU.tsv"))[order(-abs(NES))][1:14]
abbr<-c(
 "Endocrine and other factor-regulated calcium reabsorption"="Endocrine Ca reabsorption",
 "Maturity onset diabetes of the young"="MODY (diabetes)",
 "Ribosome biogenesis in eukaryotes"="Ribosome biogenesis",
 "Th1 and Th2 cell differentiation"="Th1/Th2 differentiation",
 "Nucleocytoplasmic transport"="Nucleocytoplasmic transport",
 "Nucleotide excision repair"="Nucleotide excision repair")
KG[,Desc:=Description]
for(k in names(abbr)) KG[Desc==k, Desc:=abbr[[k]]]
KG[nchar(Desc)>26, Desc:=paste0(substr(Desc,1,24),"\u2026")]
KG[,Desc:=factor(Desc,levels=rev(Desc))]
pG<-ggplot(KG,aes(NES,Desc,fill=NES>0))+geom_vline(xintercept=0,linewidth=0.3,colour="grey70")+geom_col(width=0.72)+
  scale_fill_manual(values=c("TRUE"=COL_HIGH,"FALSE"=COL_LOW),guide="none")+
  labs(x="NES (PACA-AU)",y=NULL,title="GSEA KEGG (PACA-AU)")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=5.2),
    plot.margin=margin(2,3,2,2))
save_panel(pG,file.path(PFIG,"MF7_G_kegg"),w_mm=74,h_mm=54)
cat("G 重做完成\n")
