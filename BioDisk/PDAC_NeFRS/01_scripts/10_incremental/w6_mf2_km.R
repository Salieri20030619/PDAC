#!/usr/bin/env Rscript
# MF2 KM 重出 (10pt; 统计量移到右上)
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(survival); library(glmnet) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
W2<-file.path(ROOT,"02_results/02_w2_inputs"); MIME<-file.path(ROOT,"02_results/04_mime"); PFIG<-file.path(ROOT,"03_figures/panels")
cohort_data<-readRDS(file.path(W2,"cohort_data_8.rds"))
mi<-readRDS(file.path(MIME,"mime_input_list.rds")); rr<-readRDS(file.path(MIME,"mime_results.rds"))
mg<-colnames(mi[[1]])[-(1:3)]; ml<-rr$ml.res[["Enet[α=0.1]"]]
compute_rs<-function(cn){co<-cohort_data[[cn]]
  X<-matrix(0,ncol(co$expression),length(mg),dimnames=list(colnames(co$expression),mg))
  for(g in intersect(mg,rownames(co$expression)))X[,g]<-co$expression[g,]
  rs<-as.numeric(predict(ml,newx=X,s="lambda.min")); names(rs)<-colnames(co$expression)
  cl<-co$clinical; data.frame(OS_time=cl$os_time_month,OS_status=cl$os_status,RS=rs[match(cl$sample_id,names(rs))])}
co<-list(c("PACA_AU","PACA-AU (n=267)","TRAIN"),c("GSE71729","GSE71729 (n=123)","SEL"),
  c("GSE85916","GSE85916 (n=79)","SEL"),c("GSE62452","GSE62452 (n=66)","SEL"),
  c("GSE57495","GSE57495 (n=63)","SEL"),c("GSE21501","GSE21501 (n=102)","HO"),
  c("TCGA_PAAD","TCGA-PAAD (n=147)","HO"),c("E_MTAB_6134","E-MTAB-6134 (n=288)","HO"))
fmt_p<-function(p){if(is.na(p))return("NA");if(p<1e-15)return("P < 1e-15")
  if(p<1e-4)return(sprintf("P = %.1e",p));if(p<0.001)return(sprintf("P = %.4f",p));sprintf("P = %.3f",p)}
stat1<-function(cn){d<-compute_rs(cn);d<-d[!is.na(d$RS)&!is.na(d$OS_time)&d$OS_time>0,]
  d$Group<-factor(ifelse(d$RS>median(d$RS,na.rm=TRUE),"High","Low"),levels=c("Low","High"))
  cox<-coxph(Surv(OS_time,OS_status)~Group,data=d);ci<-exp(confint(cox))
  list(fit=survfit(Surv(OS_time,OS_status)~Group,data=d),d=d,hr=exp(coef(cox)),lo=ci[1,1],hi=ci[1,2],
    pv=1-pchisq(survdiff(Surv(OS_time,OS_status)~Group,data=d)$chisq,1),
    C=as.numeric(concordance(Surv(OS_time,OS_status)~RS,data=d,reverse=TRUE)$concordance),
    n=nrow(d),nh=sum(d$Group=="High"),nl=sum(d$Group=="Low"))}
cols<-c(COL_LOW,COL_HIGH)
draw_km<-function(s,lab){
  layout(matrix(c(1,2),ncol=1),heights=c(3.0,1.25)); mt<-max(s$d$OS_time);tb<-pretty(c(0,mt),5)
  par(mar=c(0.6,4,2.1,1),mgp=c(2.3,0.5,0),tcl=-0.25)
  plot(s$fit,col=cols,lwd=2,conf.int=FALSE,xlab="",ylab="Overall survival",xlim=c(0,mt),ylim=c(0,1),
       xaxt="n",yaxt="n",mark.time=TRUE,cex=0.4); axis(2,seq(0,1,.2),seq(0,1,.2),las=2)
  mtext(lab,side=3,line=0.5,font=2,cex=1.05)
  legend("topright",
    legend=c(sprintf("Low (n=%d)",s$nl),sprintf("High (n=%d)",s$nh),"",
             paste("Log-rank",fmt_p(s$pv)),sprintf("HR=%.2f (%.2f-%.2f)",s$hr,s$lo,s$hi),sprintf("C-index=%.3f",s$C)),
    col=c(cols,NA,NA,NA,NA),lwd=c(2,2,NA,NA,NA,NA),bty="n",cex=0.8,seg.len=1.1,y.intersp=1.05)
  par(mar=c(2.9,4,0.5,1),mgp=c(1.7,0.5,0)); plot.new();plot.window(xlim=c(0,mt),ylim=c(0,3))
  text(0,2.7,"Number at risk",adj=0,font=2,cex=0.8)
  rk<-function(g)sapply(tb,function(t)sum(s$d$Group==g&s$d$OS_time>=t))
  text(-mt*0.085,1.55,"Low",adj=1,col=cols[1],font=2,cex=0.8,xpd=NA)
  text(-mt*0.085,0.6,"High",adj=1,col=cols[2],font=2,cex=0.8,xpd=NA)
  for(k in seq_along(tb)){text(tb[k],1.55,rk("Low")[k],cex=0.8,col=cols[1]);text(tb[k],0.6,rk("High")[k],cex=0.8,col=cols[2])}
  axis(1,tb,tb);mtext("Time (months)",side=1,line=1.8,cex=0.9)
}
summ<-data.frame()
for(x in co){s<-stat1(x[1]); save_base(function() draw_km(s,x[2]), file.path(PFIG,paste0("MF2_KM_",x[1])),w_mm=80,h_mm=88)
  summ<-rbind(summ,data.frame(cohort=x[1],role=x[3],n=s$n,c_index=round(s$C,3),HR=round(s$hr,2),HR_lo=round(s$lo,2),HR_hi=round(s$hi,2),KM_p=s$pv))}
write.table(summ,file.path(PFIG,"MF2_summary.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
cat("MF2 KM 完成 (8 张, 统计量右上)\n")
