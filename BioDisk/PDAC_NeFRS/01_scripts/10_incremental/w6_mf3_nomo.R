#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(survival); library(rms) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; NOMD<-file.path(ROOT,"02_results/05_nomogram"); PFIG<-file.path(ROOT,"03_figures/panels")
clin<-readRDS(file.path(NOMD,"clinical_with_nefrs_3cohorts.rds"))
paca<-as.data.frame(clin$PACA_AU)
paca_cc<-paca[complete.cases(paca[,c("NeFRS","age","Sex","Stage_AJCC","Grade","os_time_month","os_status")]),]
dd<-datadist(paca_cc); options(datadist="dd")
cf<-cph(Surv(os_time_month,os_status)~NeFRS+age+Sex+Stage_AJCC+Grade,data=paca_cc,x=TRUE,y=TRUE,surv=TRUE,time.inc=12)
set.seed(5201314); val<-validate(cf,method="boot",B=200,dxy=TRUE); Cc<-(val["Dxy","index.corrected"]+1)/2
so<-Survival(cf)
nom<-nomogram(cf,fun=list(function(x)so(12,lp=x),function(x)so(36,lp=x),function(x)so(60,lp=x)),
  funlabel=c("1-Year OS","3-Year OS","5-Year OS"),fun.at=c(0.1,0.3,0.5,0.7,0.85,0.95),lp=FALSE)
save_base(function(){par(mar=c(1.6,1.2,1.2,1.2));plot(nom,xfrac=0.32,cex.axis=0.55,cex.var=0.62,
  col.grid=gray(c(0.85,0.96)),lmgp=0.13,label.every=1)
  mtext(sprintf("Bootstrap-corrected C-index = %.3f",Cc),3,line=0.0,cex=0.6,col="gray30")},
  file.path(PFIG,"MF3_nomogram"),w_mm=176,h_mm=68)
cat("nomogram h=68 重出, C=",round(Cc,3),"\n")
