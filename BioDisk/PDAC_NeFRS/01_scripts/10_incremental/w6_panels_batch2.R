#!/usr/bin/env Rscript
# W6 panel 批次2: MF4-E/F/G(NeFRS vs 纯度/增殖/Basal) + MF7-J(免疫相关性热图)
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; PFIG<-file.path(ROOT,"03_figures/panels")

## ---- MF4-E/F/G: NeFRS vs 混杂量 散点 ----
ctl<-fread(file.path(ROOT,"02_results/18_controls/controls_scores_persample.tsv"))
scat<-function(xv,xlab,file,w=80,h=72){
  d<-ctl[is.finite(get(xv)) & is.finite(NeFRS)]
  ct<-suppressWarnings(cor.test(d[[xv]],d$NeFRS,method="spearman"))
  ann<-sprintf("Spearman rho=%.2f, %s", as.numeric(ct$estimate), fmt_p(ct$p.value))
  p<-ggplot(d,aes(.data[[xv]],NeFRS,colour=cohort))+
    geom_point(size=1.4,alpha=0.55)+
    geom_smooth(aes(group=1),method="lm",se=FALSE,colour="black",linewidth=0.7)+
    scale_colour_manual(values=PAL_CAT)+
    annotate("text",x=-Inf,y=Inf,hjust=-0.08,vjust=1.5,label=ann,size=13/.pt)+
    labs(x=xlab,y="NeFRS")+theme(legend.position="none")
  save_panel(p,file,w_mm=w,h_mm=h)
}
scat("purity","Tumor purity (ESTIMATE, z)",  file.path(PFIG,"MF4E_NeFRS_vs_purity"))
scat("prolif","Proliferation score (z)",     file.path(PFIG,"MF4F_NeFRS_vs_prolif"))
scat("basal_prob","PurIST Basal probability",file.path(PFIG,"MF4G_NeFRS_vs_basal"))

## ---- MF7-J: 免疫细胞 × 队列 与 NeFRS 相关性热图 ----
ic<-fread(file.path(ROOT,"02_results/14_immune/immune_nefrs_cor.tsv"))
ic[,star:=sig_star(padj)]
ord<-ic[,.(m=mean(rho,na.rm=TRUE)),by=cell][order(m),cell]
ic[,cell:=factor(cell,levels=ord)]
lim<-max(abs(ic$rho),na.rm=TRUE)
p7j<-ggplot(ic,aes(cohort,cell,fill=rho))+
  geom_tile(colour="white",linewidth=0.3)+
  geom_text(aes(label=star),size=13/.pt,vjust=0.8)+
  scale_fill_gradient2(low=COL_LOW,mid="white",high=COL_HIGH,midpoint=0,
                       limits=c(-lim,lim),name="Spearman\nrho")+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x=element_text(angle=40,hjust=1),
        axis.line=element_blank(),axis.ticks=element_blank())
save_panel(p7j,file.path(PFIG,"MF7J_immune_cor_heatmap"),w_mm=90,h_mm=max(85,length(ord)*4.5))
cat("\n批次2 完成: MF4E/F/G, MF7J (各 pdf+png+tiff)\n")
