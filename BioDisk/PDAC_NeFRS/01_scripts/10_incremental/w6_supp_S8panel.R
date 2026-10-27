#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; SP<-file.path(ROOT,"02_results/08_spatial"); PFIG<-file.path(ROOT,"03_figures/panels")
MO<-fread(file.path(SP,"spatial_moran.tsv")); MO[,Patient:=factor(Patient,levels=Patient[order(MoranI)])]
mm<-mean(MO$MoranI)
pM<-ggplot(MO,aes(MoranI,Patient))+geom_col(fill="grey62",width=0.72)+
  geom_vline(xintercept=0,colour="grey30",linewidth=0.3)+
  geom_vline(xintercept=mm,colour=COL_HIGH,linetype=2,linewidth=0.5)+
  annotate("text",x=mm,y=Inf,label=sprintf("mean = %.2f",mm),colour=COL_HIGH,hjust=-0.05,vjust=1.6,size=5/.pt)+
  labs(x="Moran's I (epithelial NeFRS)",y=NULL,title="Spatial autocorrelation per patient")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=4.6))
save_panel(pM,file.path(PFIG,"S8_moran"),w_mm=82,h_mm=72)
RO<-fread(file.path(SP,"spatial_nefrs_roi.tsv")); RO<-RO[segment %in% c("Epithelial","CAF","Immune")]
RO[,trt2:=factor(ifelse(treatment=="Untreated","Untreated","Treated/CRT"),levels=c("Untreated","Treated/CRT"))]
RO[,segment:=factor(segment,levels=c("Epithelial","CAF","Immune"))]
pv<-RO[,.(p=wilcox.test(nefrs~trt2)$p.value),by=segment]; print(pv)
pT<-ggplot(RO,aes(segment,nefrs,fill=trt2))+geom_boxplot(outlier.size=0.2,linewidth=0.3,width=0.66)+
  scale_fill_manual(values=c("#9E9E9E",COL_HIGH),name=NULL)+
  labs(x=NULL,y="NeFRS",title="NeFRS by segment: pre vs post CRT (n.s.)")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),legend.position="top",legend.key.size=unit(2.6,"mm"),legend.text=element_text(size=5.4))
save_panel(pT,file.path(PFIG,"S8_treat"),w_mm=92,h_mm=72)
cat("S8 新panel完成\n")
