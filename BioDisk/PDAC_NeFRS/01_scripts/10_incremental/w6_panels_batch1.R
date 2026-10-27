#!/usr/bin/env Rscript
# W6 panel 批次1: MF2-K(8队列C-index) / MF3-H(超越分期C) / MF3-I(IDI/NRI)
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(survival); library(glmnet); library(data.table) })
ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
PFIG <- file.path(ROOT,"03_figures/panels"); dir.create(PFIG,recursive=TRUE,showWarnings=FALSE)

## ---- MF2-K: 8 队列 C-index 条形 ----
W2IN<-file.path(ROOT,"02_results/02_w2_inputs"); MIME<-file.path(ROOT,"02_results/04_mime")
cohort_data<-readRDS(file.path(W2IN,"cohort_data_8.rds"))
mime_input<-readRDS(file.path(MIME,"mime_input_list.rds")); res<-readRDS(file.path(MIME,"mime_results.rds"))
mime_genes<-colnames(mime_input[[1]])[-(1:3)]; ml<-res$ml.res[["Enet[α=0.1]"]]
compute_rs<-function(em){X<-matrix(0,ncol(em),length(mime_genes),dimnames=list(colnames(em),mime_genes))
  for(g in intersect(mime_genes,rownames(em)))X[,g]<-em[g,]; setNames(as.numeric(predict(ml,newx=X,s="lambda.min")),colnames(em))}
ci<-data.frame()
for(cn in names(cohort_data)){co<-cohort_data[[cn]];cl<-co$clinical;rs<-compute_rs(co$expression)
  d<-data.frame(t=cl$os_time_month,s=cl$os_status,z=as.numeric(rs[cl$sample_id]))
  d<-d[!is.na(d$t)&d$t>0&!is.na(d$s)&!is.na(d$z),]
  C<-concordance(coxph(Surv(t,s)~z,data=d))$concordance
  ci<-rbind(ci,data.frame(cohort=cn,C=C,set=ifelse(cn=="PACA_AU","Training","Validation")))}
ci$cohort<-factor(ci$cohort,levels=ci$cohort[order(-ci$C)])
p2k<-ggplot(ci,aes(cohort,C,fill=set))+
  geom_col(width=0.7)+geom_hline(yintercept=0.5,linetype=2,linewidth=0.4,colour="grey50")+
  geom_text(aes(label=sprintf("%.3f",C)),vjust=-0.4,size=13/.pt)+
  scale_fill_manual(values=c(Training=COL_HIGH,Validation=COL_LOW))+
  coord_cartesian(ylim=c(0.45,0.82))+labs(x=NULL,y="C-index",fill=NULL)+
  theme(axis.text.x=element_text(angle=40,hjust=1),legend.position=c(0.85,0.92))
save_panel(p2k,file.path(PFIG,"MF2K_cindex_8cohort"),w_mm=110,h_mm=80)

## ---- MF3-H: 超越分期 C-index ----
bs<-fread(file.path(ROOT,"02_results/21_nri/beyond_stage_cindex.tsv"))
dh<-data.frame(model=factor(c("AJCC stage","Stage + NeFRS"),levels=c("AJCC stage","Stage + NeFRS")),
               C=c(bs$C_stage,bs$C_stage_NeFRS))
p3h<-ggplot(dh,aes(model,C,fill=model))+
  geom_col(width=0.6)+geom_hline(yintercept=0.5,linetype=2,linewidth=0.4,colour="grey50")+
  geom_text(aes(label=sprintf("%.3f",C)),vjust=-0.4,size=13/.pt)+
  annotate("text",x=1.5,y=max(dh$C)+0.05,label=sprintf("ΔC = +%.3f",bs$dC),fontface=2,size=13/.pt)+
  scale_fill_manual(values=c("AJCC stage"="grey70","Stage + NeFRS"=COL_HIGH),guide="none")+
  coord_cartesian(ylim=c(0.45,0.78))+labs(x=NULL,y="C-index")
save_panel(p3h,file.path(PFIG,"MF3H_beyondstage_C"),w_mm=72,h_mm=74)

## ---- MF3-I: IDI / 连续NRI / 中位改善 (12, 24 月) ----
di<-data.frame(
  measure=factor(rep(c("IDI","Continuous NRI","Median improve"),2),levels=c("IDI","Continuous NRI","Median improve")),
  t0=factor(rep(c("12 mo","24 mo"),each=3),levels=c("12 mo","24 mo")),
  est=c(0.152,0.376,0.146, 0.132,0.280,0.160),
  lo =c(0.107,0.291,0.106, 0.094,0.206,0.104),
  hi =c(0.205,0.473,0.215, 0.173,0.366,0.233))
p3i<-ggplot(di,aes(measure,est,colour=t0))+
  geom_hline(yintercept=0,linewidth=0.4,colour="grey60")+
  geom_point(position=position_dodge(0.5),size=2.4)+
  geom_errorbar(aes(ymin=lo,ymax=hi),position=position_dodge(0.5),width=0.18,linewidth=0.6)+
  scale_colour_manual(values=c("12 mo"=COL_LOW,"24 mo"=COL_HIGH))+
  labs(x=NULL,y="Improvement vs AJCC stage",colour=NULL)+
  theme(axis.text.x=element_text(angle=20,hjust=1),legend.position=c(0.85,0.92))
save_panel(p3i,file.path(PFIG,"MF3I_idi_nri"),w_mm=98,h_mm=74)

cat("\n批次1 完成 -> 03_figures/panels/ : MF2K, MF3H, MF3I (各 pdf+png+tiff)\n")
