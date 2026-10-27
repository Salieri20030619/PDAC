#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; DRG<-file.path(ROOT,"02_results/16_drug"); PFIG<-file.path(ROOT,"03_figures/panels")
DP<-fread(file.path(DRG,"PACA_AU_DrugPred.csv")); setnames(DP,1,"sample")
cl<-readRDS(file.path(ROOT,"02_results/05_nomogram/clinical_with_nefrs_3cohorts.rds"))
paca<-as.data.table(if(is.list(cl)&&!is.data.frame(cl)) cl[[grep("PACA",names(cl))[1]]] else cl)
if("cohort" %in% colnames(paca)) paca<-paca[grepl("PACA",cohort)]
paca<-paca[!is.na(NeFRS)]
paca[,grp:=ifelse(NeFRS>median(NeFRS,na.rm=TRUE),"High","Low")]
DPl<-melt(DP,id.vars="sample",variable.name="drug",value.name="pred"); DPl[,drugc:=gsub("_[0-9]+$","",drug)]
cand<-c("BI-2536","Dasatinib","MK-2206","AZD5363")
DPc<-merge(DPl[drugc %in% cand],paca[,.(sample=sample_id,grp)],by="sample")
DPc<-DPc[grp %in% c("Low","High")]                          # only two groups
cat("匹配样本:",length(unique(DPc$sample)),"| 组:",paste(names(table(DPc$grp)),table(DPc[!duplicated(sample)]$grp),collapse=" "),"\n")
DPc[,grp:=factor(grp,levels=c("Low","High"))]; DPc[,drugc:=factor(drugc,levels=cand)]
pv<-DPc[,.(p=wilcox.test(pred~grp)$p.value),by=drugc]
lab<-pv[,.(drugc,lbl=ifelse(p<1e-4,"****",ifelse(p<1e-3,"***",ifelse(p<1e-2,"**","*"))))]
pG<-ggplot(DPc,aes(grp,pred,fill=grp))+geom_boxplot(outlier.size=0.25,linewidth=0.3,width=0.6)+
  facet_wrap(~drugc,scales="free_y",nrow=1)+
  geom_text(data=lab,aes(x=1.5,y=Inf,label=lbl),inherit.aes=FALSE,vjust=1.3,size=6/.pt)+
  scale_fill_manual(values=c(High=COL_HIGH,Low=COL_LOW),guide="none")+
  labs(x=NULL,y="Predicted response (lower = sensitive)",title="Candidate drugs (high vs low NeFRS, PACA-AU)")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),strip.text=element_text(size=5.8,face="bold"),
    panel.spacing=unit(2.4,"mm"),axis.text.x=element_text(size=5.6))
save_panel(pG,file.path(PFIG,"MF8_G_box"),w_mm=120,h_mm=52)
cat("G 去NA 完成\n")
