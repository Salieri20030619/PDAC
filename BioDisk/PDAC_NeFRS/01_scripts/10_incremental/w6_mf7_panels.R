#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; IMM<-file.path(ROOT,"02_results/14_immune"); GSE<-file.path(ROOT,"02_results/15_gsea"); PFIG<-file.path(ROOT,"03_figures/panels")

## A immune–NeFRS correlation heatmap (cell x cohort, rho)
COR<-fread(file.path(IMM,"immune_nefrs_cor.tsv"))
ord<-COR[,.(m=mean(rho)),by=cell][order(m)]$cell
COR[,cell:=factor(cell,levels=ord)]
COR[,cohort:=factor(cohort,levels=c("PACA_AU","TCGA_PAAD","GSE71729"))]
cat("免疫细胞数:",length(ord),"| 队列:",paste(levels(COR$cohort),collapse=","),"\n")
pA<-ggplot(COR,aes(cohort,cell,fill=rho))+geom_tile(colour="white",linewidth=0.3)+
  scale_fill_gradient2(low=COL_LOW,mid="white",high=COL_HIGH,midpoint=0,name="\u03c1")+
  labs(x=NULL,y=NULL,title="NeFRS vs immune cells (Spearman)")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=4.6),
    axis.text.x=element_text(size=5.4,angle=30,hjust=1),legend.key.width=unit(2,"mm"),legend.key.height=unit(3,"mm"),
    legend.title=element_text(size=6),legend.text=element_text(size=4.8),panel.grid=element_blank())
save_panel(pA,file.path(PFIG,"MF7_A_cor"),w_mm=62,h_mm=78)

## B mean rho lollipop
MR<-COR[,.(rho=mean(rho)),by=cell][order(rho)]; MR[,cell:=factor(cell,levels=cell)]
pB<-ggplot(MR,aes(rho,cell))+geom_vline(xintercept=0,linewidth=0.3,colour="grey70")+
  geom_segment(aes(x=0,xend=rho,y=cell,yend=cell),linewidth=0.4,colour="grey60")+
  geom_point(aes(colour=rho>0),size=1.4)+
  scale_colour_manual(values=c("TRUE"=COL_HIGH,"FALSE"=COL_LOW),guide="none")+
  labs(x="Mean Spearman \u03c1 (3 cohorts)",y=NULL,title="Immune correlation (mean)")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=4.6))
save_panel(pB,file.path(PFIG,"MF7_B_lollipop"),w_mm=58,h_mm=78)

## C immune desert: high vs low (PACA_AU diff), top |diff|
HL<-fread(file.path(IMM,"immune_highlow_diff.tsv"))[cohort=="PACA_AU"]
HL<-HL[order(diff)]; HL[,cell:=factor(cell,levels=cell)]
pC<-ggplot(HL,aes(diff,cell))+geom_vline(xintercept=0,linewidth=0.3,colour="grey70")+
  geom_col(aes(fill=diff>0),width=0.72)+
  scale_fill_manual(values=c("TRUE"=COL_HIGH,"FALSE"=COL_LOW),labels=c("Down in high-NeFRS","Up in high-NeFRS"),name=NULL)+
  labs(x="ssGSEA diff (high - low NeFRS, PACA-AU)",y=NULL,title="Immune desert (high vs low)")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=4.6),
    legend.position="bottom",legend.key.size=unit(2.4,"mm"),legend.text=element_text(size=5))
save_panel(pC,file.path(PFIG,"MF7_C_desert"),w_mm=64,h_mm=78)

## D GSEA GO-BP consistent: top up + top down by mean_NES
GC<-fread(file.path(GSE,"gsea_consistent_3cohorts.tsv"))
if("consistent" %in% colnames(GC)) GC<-GC[consistent==TRUE]
up<-GC[order(-mean_NES)][1:10]; dn<-GC[order(mean_NES)][1:10]
cat("一致通路:",nrow(GC)," | NES范围:",paste(round(range(GC$mean_NES),2),collapse=" ~ "),"\n")
GG<-rbind(up,dn); GG<-GG[order(mean_NES)]; GG[,Description:=factor(Description,levels=Description)]
pD<-ggplot(GG,aes(mean_NES,Description,fill=mean_NES>0))+geom_vline(xintercept=0,linewidth=0.3,colour="grey70")+
  geom_col(width=0.72)+
  scale_fill_manual(values=c("TRUE"=COL_HIGH,"FALSE"=COL_LOW),labels=c("Down in high-NeFRS","Up in high-NeFRS"),name=NULL)+
  labs(x="Mean NES (3 cohorts)",y=NULL,title="GSEA GO-BP (consistent)")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=5),
    legend.position="bottom",legend.key.size=unit(2.4,"mm"),legend.text=element_text(size=5))
save_panel(pD,file.path(PFIG,"MF7_D_gsea"),w_mm=110,h_mm=78)
cat("MF7 A-D 完成\n")
