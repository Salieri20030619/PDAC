#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
SC<-file.path(ROOT,"02_results/06_sc"); PS<-file.path(ROOT,"02_results/09_pseudotime")
SPA<-file.path(ROOT,"02_results/08_spatial"); PFIG<-file.path(ROOT,"03_figures/panels")

PT<-fread(file.path(PS,"pseudo_nefrs_along_pt.tsv"))
set.seed(42); PTs<-PT[sample(.N, min(25000,.N))]
ctlev<-c("Malignant","Ductal (atypical)","Ductal","ADM","Acinar")
ctcol<-c("Malignant"=COL_HIGH,"Ductal (atypical)"="#E69F00","Ductal"="#56B4E9","ADM"="#009E73","Acinar"="#999999")
PTs[,celltype:=factor(celltype,levels=ctlev)]
uth<-function() theme_void(base_size=6.5)+theme(
  legend.key.size=unit(2.4,"mm"),legend.text=element_text(size=5.4),legend.title=element_text(size=6),
  plot.title=element_text(size=7,face="bold",hjust=0.5),plot.margin=margin(1,1,1,1))

## A UMAP by celltype
pA<-ggplot(PTs,aes(UMAP_1,UMAP_2,colour=celltype))+geom_point(size=0.12,alpha=0.5,stroke=0)+
  scale_colour_manual(values=ctcol,name=NULL)+
  guides(colour=guide_legend(override.aes=list(size=1.4,alpha=1),ncol=2))+
  labs(title="Epithelial cell type")+uth()+theme(legend.position="bottom")
save_panel(pA,file.path(PFIG,"MF5_A_umap_celltype"),w_mm=58,h_mm=54)

## B UMAP NeFRS feature
pB<-ggplot(PTs[order(NeFRS)],aes(UMAP_1,UMAP_2,colour=NeFRS))+geom_point(size=0.12,stroke=0)+
  scale_colour_viridis_c(option="magma",name="NeFRS",limits=c(0,quantile(PTs$NeFRS,0.90,na.rm=TRUE)),oob=scales::squish)+labs(title="NeFRS (UCell)")+
  uth()+theme(legend.position="right")
save_panel(pB,file.path(PFIG,"MF5_B_umap_nefrs"),w_mm=58,h_mm=54)

## C NeFRS by cell type (all 12 types, median + IQR)
F6<-fread(file.path(SC,"w4_fig6a_ucell_by_celltype.tsv"))
F6<-F6[order(median)]; F6[,Level1:=factor(Level1,levels=Level1)]
F6[,mal:=Level1=="Epithelial (malignant)"]
pC<-ggplot(F6,aes(median,Level1))+
  geom_vline(xintercept=0,linetype=2,linewidth=0.4,colour="grey70")+
  geom_segment(aes(x=p25,xend=p75,yend=Level1),linewidth=0.5,colour="grey55")+
  geom_point(aes(colour=mal),size=1.3)+
  scale_colour_manual(values=c("TRUE"=COL_HIGH,"FALSE"="grey35"),guide="none")+
  labs(x="NeFRS UCell (median, IQR)",y=NULL)+theme_pub(base=6.5)+
  theme(axis.text.y=element_text(size=5.6))
save_panel(pC,file.path(PFIG,"MF5_C_celltype"),w_mm=64,h_mm=54)

## D UMAP by pseudotime
pD<-ggplot(PTs[order(pseudotime)],aes(UMAP_1,UMAP_2,colour=pseudotime))+geom_point(size=0.12,stroke=0)+
  scale_colour_viridis_c(option="viridis",name="Pseudotime")+labs(title="Pseudotime trajectory")+
  uth()+theme(legend.position="right")
save_panel(pD,file.path(PFIG,"MF5_D_umap_pt"),w_mm=58,h_mm=52)

## E NeFRS vs pseudotime
pE<-ggplot(PTs,aes(pseudotime,NeFRS))+geom_point(size=0.08,alpha=0.12,colour="grey50",stroke=0)+
  geom_smooth(method="loess",se=TRUE,colour=COL_HIGH,fill=COL_HIGH,linewidth=0.7,alpha=0.2)+
  labs(x="Pseudotime",y="NeFRS")+theme_pub(base=6.5)
save_panel(pE,file.path(PFIG,"MF5_E_nefrs_pt"),w_mm=58,h_mm=52)

## F Malignant high vs low NeFRS DE volcano
DE<-fread(file.path(PS,"malignant_NeFRS_hilo_DE.tsv"))
DE[,neglogp:=-log10(padj)]; DE[is.infinite(neglogp)|neglogp>300,neglogp:=300]
DE[,sig:=fifelse(padj<0.05 & logFC>0.5,"up",fifelse(padj<0.05 & logFC< -0.5,"down","ns"))]
hub<-c("NRP1","EGFR","BMPR1B","SEMA3A","NR2F2","WWTR1","MAPK10")
lab<-DE[feature %in% hub]
pF<-ggplot(DE,aes(logFC,neglogp))+geom_point(aes(colour=sig),size=0.3,alpha=0.4,stroke=0)+
  scale_colour_manual(values=c(up=COL_HIGH,down=COL_LOW,ns="grey78"),guide="none")+
  geom_point(data=lab,colour="black",size=0.9)+
  geom_text_repel(data=lab,aes(label=feature),size=5.6/.pt,max.overlaps=20,min.segment.length=0,segment.size=0.2)+
  labs(x="log2FC (high vs low NeFRS)",y="-log10 padj")+theme_pub(base=6.5)
save_panel(pF,file.path(PFIG,"MF5_F_volcano"),w_mm=58,h_mm=52)

## G NeFRS-member trajectory association (top by Wald)
TA<-fread(file.path(PS,"nefrs_members_trajectory_assoc.tsv"))
TA<-TA[is_nefrs==TRUE][order(-waldStat)][1:12]; TA[,gene:=factor(gene,levels=rev(gene))]
pG<-ggplot(TA,aes(waldStat,gene))+geom_col(fill=COL_HIGH,width=0.7)+
  labs(x="tradeSeq Wald statistic",y=NULL,title="NeFRS members along trajectory")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=6))
save_panel(pG,file.path(PFIG,"MF5_G_traj_genes"),w_mm=88,h_mm=48)

## H Spatial NeFRS by segment
SP<-fread(file.path(SPA,"spatial_nefrs_by_segment.tsv"))
SP[,segment:=factor(segment,levels=c("CAF","Epithelial","Immune"))]
pH<-ggplot(SP,aes(segment,mean,fill=segment))+geom_col(width=0.6)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2,linewidth=0.5)+
  geom_text(aes(y=mean+sd,label=sprintf("%.2f",mean)),vjust=-0.7,size=6.5/.pt)+coord_cartesian(ylim=c(0,0.92))+
  scale_fill_manual(values=c(CAF="#E69F00",Epithelial=COL_HIGH,Immune=COL_LOW),guide="none")+
  labs(x=NULL,y="Spatial NeFRS (mean \u00b1 SD)",title="Spatial segments")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"))
save_panel(pH,file.path(PFIG,"MF5_H_spatial"),w_mm=88,h_mm=48)
cat("MF5 A-H 完成\n")
