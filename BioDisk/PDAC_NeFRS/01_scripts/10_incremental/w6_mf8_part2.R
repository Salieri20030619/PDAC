#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; DRG<-file.path(ROOT,"02_results/16_drug"); SIG<-file.path(ROOT,"02_results/13_sig_compare"); PFIG<-file.path(ROOT,"03_figures/panels")
clean<-function(x) gsub("_[0-9]+$","",x)

## E re-render: label left of line, no clip
SC<-fread(file.path(SIG,"sig_compare_all_cohorts.tsv"))[cohort=="PACA_AU"]
NR<-fread(file.path(SIG,"nefrs_rank_summary.tsv")); nefrs_c<-NR[cohort=="PACA_AU"]$NeFRS_cindex
pE<-ggplot(SC,aes(cindex))+geom_histogram(bins=24,fill="grey70",colour="white",linewidth=0.2)+
  geom_vline(xintercept=nefrs_c,colour=COL_HIGH,linewidth=0.8)+
  annotate("text",x=nefrs_c-0.006,y=Inf,label="NeFRS 0.754 (1/84)",colour=COL_HIGH,hjust=1,vjust=1.5,size=5/.pt,fontface="bold")+
  labs(x="C-index (PACA-AU)",y="Published signatures (n)",title="NeFRS vs 84 published signatures")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"))
save_panel(pE,file.path(PFIG,"MF8_E_dist"),w_mm=72,h_mm=50)

## F drug response consistency across 3 cohorts (top sensitive + resistant)
DH<-fread(file.path(DRG,"drug_highlow_diff.tsv"))
agg<-DH[,.(mlfc=mean(log2fc)),by=drug]
sel<-c(agg[order(mlfc)][1:7]$drug, agg[order(-mlfc)][1:7]$drug)
DHt<-DH[drug %in% sel]; DHt[,drugc:=clean(drug)]
ord<-DHt[,.(m=mean(log2fc)),by=drugc][order(m)]$drugc
DHt[,drugc:=factor(drugc,levels=ord)]; DHt[,cohort:=factor(cohort,levels=c("PACA_AU","TCGA_PAAD","GSE71729"))]
pF<-ggplot(DHt,aes(cohort,drugc,fill=log2fc))+geom_tile(colour="white",linewidth=0.3)+
  scale_fill_gradient2(low=COL_HIGH,mid="white",high=COL_LOW,midpoint=0,name="log2FC")+
  labs(x=NULL,y=NULL,title="Drug response (3-cohort consistency)")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=5),
    axis.text.x=element_text(size=5,angle=30,hjust=1),legend.key.width=unit(2,"mm"),legend.key.height=unit(3,"mm"),
    legend.title=element_text(size=5.5),legend.text=element_text(size=4.8),panel.grid=element_blank())
save_panel(pF,file.path(PFIG,"MF8_F_drugheat"),w_mm=66,h_mm=50)
cat("E修 + F 完成\n")
