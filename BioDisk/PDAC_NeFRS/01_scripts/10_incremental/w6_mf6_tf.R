#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; TF<-file.path(ROOT,"02_results/11_tf"); PFIG<-file.path(ROOT,"03_figures/panels")
res<-fread(file.path(TF,"TF_highlow_diff.tsv"))

res25<-res[order(-abs(diff))][1:25]; res25$TF<-factor(res25$TF, levels=res25$TF[order(res25$diff)])
p1<-ggplot(res25,aes(TF,diff,color=diff>0))+
  geom_segment(aes(xend=TF,y=0,yend=diff),linewidth=0.6)+geom_point(size=2.4)+coord_flip()+
  scale_color_manual(values=c("TRUE"=COL_HIGH,"FALSE"=COL_LOW),guide="none")+
  geom_hline(yintercept=0,lty=2,color="gray50")+
  labs(x=NULL,y="TF activity diff (High - Low NeFRS)")+theme_pub()
save_panel(p1, file.path(PFIG,"MF6_TF_top25"), w_mm=82, h_mm=120)

key<-c("NFIB","TRPS1","REST","SOX4","NEUROD1","GLIS3","WWTR1","NR2F1","ZEB1","INSM1")
kdf<-res[TF%in%key]; kdf$TF<-factor(kdf$TF,levels=kdf$TF[order(kdf$diff)])
p2<-ggplot(kdf,aes(TF,diff,fill=diff>0))+geom_col(width=0.7)+coord_flip()+
  scale_fill_manual(values=c("TRUE"=COL_HIGH,"FALSE"=COL_LOW),guide="none")+
  geom_hline(yintercept=0,lty=2,color="gray50")+
  labs(x=NULL,y="TF activity diff (High - Low NeFRS)")+theme_pub()
save_panel(p2, file.path(PFIG,"MF6_TF_keyTF"), w_mm=80, h_mm=70)
cat("MF6 TF 完成 (top25 棒棒糖 + 关键TF)\n")
