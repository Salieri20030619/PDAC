#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; PFIG<-file.path(ROOT,"03_figures/panels")
TFd<-fread(file.path(ROOT,"02_results/11_tf/TF_highlow_diff.tsv"))
TFd[,neglogp:=-log10(padj)]; TFd[is.infinite(neglogp)|neglogp>320,neglogp:=320]
TFd[,sig:=fifelse(padj<0.05 & diff>0.05,"up",fifelse(padj<0.05 & diff< -0.05,"down","ns"))]
cat("padj<0.05:",sum(TFd$padj<0.05)," | up(diff>0.05):",sum(TFd$sig=="up")," down:",sum(TFd$sig=="down"),"\n")
stress<-c("NFE2L3","BACH1","FOXM1","E2F1","MAX","PROX1")   # 应激/增殖 (主信号)
matureN<-c("TRPS1","INSM1","NEUROD1")                       # 成熟神经分化 (强抑制)
progN<-c("NR2F1","GLIS3","NFIB")                            # 神经祖细胞 (弱)
lab<-TFd[TF %in% c(stress,matureN,progN)]
lab[,grp:=fifelse(TF %in% stress,"Stress/proliferation",fifelse(TF %in% matureN,"Mature neuronal","Neural progenitor"))]
pA<-ggplot(TFd,aes(diff,neglogp))+
  geom_point(aes(colour=sig),size=0.4,alpha=0.35,stroke=0)+
  scale_colour_manual(values=c(up=COL_HIGH,down=COL_LOW,ns="grey82"),guide="none")+
  geom_vline(xintercept=0,linewidth=0.3,colour="grey70")+
  geom_hline(yintercept=-log10(0.05),linetype=2,linewidth=0.4,colour="grey55")+
  geom_point(data=lab,aes(fill=grp),colour="black",size=1.2,shape=21,stroke=0.3)+
  scale_fill_manual(values=c("Stress/proliferation"=COL_HIGH,"Mature neuronal"=COL_LOW,"Neural progenitor"="grey55"),name=NULL)+
  geom_text_repel(data=lab,aes(label=TF),size=5/.pt,max.overlaps=30,min.segment.length=0,segment.size=0.2)+
  labs(x="TF activity diff (high - low 1/3 NeFRS)",y="-log10 padj",title="TF activity (decoupleR)")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),
    legend.position="bottom",legend.key.size=unit(2.4,"mm"),legend.text=element_text(size=5),
    legend.margin=margin(0,0,0,0))
save_panel(pA,file.path(PFIG,"MF6_A_tf"),w_mm=58,h_mm=54)
cat("A 重做完成\n")
