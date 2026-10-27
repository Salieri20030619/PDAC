#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(survival); library(data.table) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; NOMD<-file.path(ROOT,"02_results/05_nomogram"); PFIG<-file.path(ROOT,"03_figures/panels")
clin<-readRDS(file.path(NOMD,"clinical_with_nefrs_3cohorts.rds"))
fp<-function(p){if(is.na(p))return("P=NA");if(p<0.001)return("P<0.001");sprintf("P=%.3f",p)}
pool<-rbindlist(lapply(names(clin),function(cn){d<-as.data.table(clin[[cn]])
  data.table(cohort=cn,t=as.numeric(d$os_time_month),s=as.numeric(d$os_status),
    RS=as.numeric(d$NeFRS),stg=as.character(d$Stage_AJCC))}),fill=TRUE)
pool<-pool[!is.na(RS)&!is.na(t)&t>0&!is.na(s)]
pool[,RS_z:=as.numeric(scale(RS)),by=cohort]
pool[,stage:=ifelse(grepl("III|IV",stg),"Late",ifelse(grepl("I",stg),"Early",NA))]
cat("== 归并 stage ==\n"); print(table(pool$stage,useNA="ifany"))
draw_km<-function(tgt,ttl){                       # 参数改名 tgt, 不再撞列名 stg
  d<-pool[stage==tgt & !is.na(RS_z)]
  cat(sprintf("[%s] n=%d events=%d\n",ttl,nrow(d),sum(d$s)))
  if(nrow(d)<10){plot.new();text(0.5,0.5,"n too small");return(invisible())}
  d[,grp:=factor(ifelse(RS_z>median(RS_z,na.rm=TRUE),"High","Low"),levels=c("Low","High"))]
  fit<-survfit(Surv(t,s)~grp,data=d)
  hr<-tryCatch(exp(coef(coxph(Surv(t,s)~grp+strata(cohort),data=d))),error=function(e)NA)
  pv<-tryCatch(1-pchisq(survdiff(Surv(t,s)~grp,data=d)$chisq,1),error=function(e)NA)
  par(mar=c(2.4,2.9,1.4,0.6),mgp=c(1.4,0.4,0),tcl=-0.2,lwd=0.7)
  plot(fit,col=c(COL_LOW,COL_HIGH),lwd=1.5,conf.int=FALSE,xlab="",ylab="",ylim=c(0,1),
    xaxt="n",yaxt="n",mark.time=TRUE,cex=0.4)
  axis(1,cex.axis=0.6,lwd=0.6,padj=-0.6);axis(2,seq(0,1,.5),las=2,cex.axis=0.6,lwd=0.6)
  mtext(sprintf("%s (n=%d)",ttl,nrow(d)),3,line=0.3,font=2,cex=0.62)
  mtext("Months",1,line=1.2,cex=0.55);mtext("Survival",2,line=1.9,cex=0.55,las=0)
  legend("topright",c(if(!is.na(hr))sprintf("HR %.2f",hr) else "HR NA",fp(pv)),bty="n",cex=0.56,y.intersp=0.95)
  legend("bottomleft",c("Low","High"),col=c(COL_LOW,COL_HIGH),lwd=1.5,bty="n",cex=0.56,seg.len=0.8,y.intersp=0.9)
}
save_base(function() draw_km("Early","Early stage (I-II)"),file.path(PFIG,"MF3_km_early"),w_mm=86,h_mm=40)
save_base(function() draw_km("Late","Late stage (III-IV)"),file.path(PFIG,"MF3_km_late"),w_mm=86,h_mm=40)
