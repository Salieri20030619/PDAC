#!/usr/bin/env Rscript
# MF2 旧图 10pt 重出: 8 KM + 8 ROC + 森林 (base R -> save_base, 三格式, 无大标题)
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(survival); library(glmnet); library(timeROC) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
W2<-file.path(ROOT,"02_results/02_w2_inputs"); MIME<-file.path(ROOT,"02_results/04_mime")
PFIG<-file.path(ROOT,"03_figures/panels"); dir.create(PFIG,recursive=TRUE,showWarnings=FALSE)

cohort_data<-readRDS(file.path(W2,"cohort_data_8.rds"))
mi<-readRDS(file.path(MIME,"mime_input_list.rds")); rr<-readRDS(file.path(MIME,"mime_results.rds"))
mg<-colnames(mi[[1]])[-(1:3)]; ml<-rr$ml.res[["Enet[α=0.1]"]]
compute_rs<-function(cn){co<-cohort_data[[cn]]
  X<-matrix(0,ncol(co$expression),length(mg),dimnames=list(colnames(co$expression),mg))
  for(g in intersect(mg,rownames(co$expression)))X[,g]<-co$expression[g,]
  rs<-as.numeric(predict(ml,newx=X,s="lambda.min")); names(rs)<-colnames(co$expression)
  cl<-co$clinical
  data.frame(sample_id=cl$sample_id,OS_time=cl$os_time_month,OS_status=cl$os_status,RS=rs[match(cl$sample_id,names(rs))])}

co<-list(
  c("PACA_AU","PACA-AU (n=267)","TRAIN"), c("GSE71729","GSE71729 (n=123)","SEL"),
  c("GSE85916","GSE85916 (n=79)","SEL"),  c("GSE62452","GSE62452 (n=66)","SEL"),
  c("GSE57495","GSE57495 (n=63)","SEL"),  c("GSE21501","GSE21501 (n=102)","HO"),
  c("TCGA_PAAD","TCGA-PAAD (n=147)","HO"),c("E_MTAB_6134","E-MTAB-6134 (n=288)","HO"))

fmt_p<-function(p){if(is.na(p))return("NA");if(p<1e-15)return("P < 1e-15")
  if(p<1e-4)return(sprintf("P = %.1e",p));if(p<0.001)return(sprintf("P = %.4f",p));sprintf("P = %.3f",p)}
stat1<-function(cn){d<-compute_rs(cn); d<-d[!is.na(d$RS)&!is.na(d$OS_time)&d$OS_time>0,]
  d$Group<-factor(ifelse(d$RS>median(d$RS,na.rm=TRUE),"High","Low"),levels=c("Low","High"))
  cox<-coxph(Surv(OS_time,OS_status)~Group,data=d); ci<-exp(confint(cox))
  list(fit=survfit(Surv(OS_time,OS_status)~Group,data=d), d=d,
       hr=exp(coef(cox)),lo=ci[1,1],hi=ci[1,2],
       pv=1-pchisq(survdiff(Surv(OS_time,OS_status)~Group,data=d)$chisq,1),
       C=as.numeric(concordance(Surv(OS_time,OS_status)~RS,data=d,reverse=TRUE)$concordance),
       n=nrow(d),nh=sum(d$Group=="High"),nl=sum(d$Group=="Low"))}

cols<-c(COL_LOW,COL_HIGH)
## ---- 8 KM ----
draw_km<-function(s,lab){
  layout(matrix(c(1,2),ncol=1),heights=c(3.4,1.0))
  par(mar=c(2.6,3.4,1.6,0.8),mgp=c(1.9,0.5,0),tcl=-0.25)
  mt<-max(s$d$OS_time);tb<-pretty(c(0,mt),5)
  plot(s$fit,col=cols,lwd=2,conf.int=FALSE,xlab="Time (months)",ylab="Overall survival",
       xlim=c(0,mt),ylim=c(0,1),xaxt="n",yaxt="n",mark.time=TRUE,cex=0.4)
  axis(1,tb,tb);axis(2,seq(0,1,.2),seq(0,1,.2),las=2)
  mtext(lab,side=3,line=0.4,font=2,cex=1.0)
  legend("topright",c(sprintf("Low (n=%d)",s$nl),sprintf("High (n=%d)",s$nh)),col=cols,lwd=2,bty="n",cex=0.85,seg.len=1.1)
  legend("bottomleft",c(paste("Log-rank",fmt_p(s$pv)),sprintf("HR=%.2f (%.2f-%.2f)",s$hr,s$lo,s$hi),
         sprintf("C-index=%.3f",s$C)),bty="n",cex=0.82,y.intersp=1.15)
  par(mar=c(2.6,3.4,0.3,0.8));plot.new();plot.window(xlim=c(0,mt),ylim=c(0,2.4))
  mtext("Number at risk",side=3,line=-0.7,adj=0,cex=0.78,font=2)
  rk<-function(g)sapply(tb,function(t)sum(s$d$Group==g&s$d$OS_time>=t))
  text(-mt*0.02,1.6,"Low",adj=1,col=cols[1],font=2,cex=0.8);text(-mt*0.02,0.6,"High",adj=1,col=cols[2],font=2,cex=0.8)
  for(k in seq_along(tb)){text(tb[k],1.6,rk("Low")[k],cex=0.8,col=cols[1]);text(tb[k],0.6,rk("High")[k],cex=0.8,col=cols[2])}
  axis(1,tb,tb)
}
summ<-data.frame()
for(x in co){s<-stat1(x[1])
  save_base(function() draw_km(s,x[2]), file.path(PFIG,paste0("MF2_KM_",x[1])), w_mm=78, h_mm=84)
  summ<-rbind(summ,data.frame(cohort=x[1],role=x[3],n=s$n,c_index=round(s$C,3),
    HR=round(s$hr,2),HR_lo=round(s$lo,2),HR_hi=round(s$hi,2),KM_p=s$pv))}

## ---- 森林 (MF2-J) ----
fp<-function(p){if(p<1e-15)return("< 1e-15");if(p<1e-4)return(sprintf("%.1e",p));if(p<0.001)return(sprintf("%.4f",p));sprintf("%.3f",p)}
draw_forest<-function(df){
  df$rf<-factor(df$role,levels=c("TRAIN","SEL","HO")); df<-df[order(df$rf,-df$HR),]; df$y<-nrow(df):1
  df$col<-ifelse(df$role=="TRAIN",COL_HIGH,ifelse(df$role=="HO",COL_LOW,"#666666"))
  df$cex<-sqrt(df$n)/sqrt(max(df$n))*1.4+0.5
  layout(matrix(c(1,2),nrow=1),widths=c(6,4))
  par(mar=c(3.6,0.5,1.2,0.2));plot.new();plot.window(xlim=c(0,10),ylim=c(-0.6,nrow(df)+1.0))
  cx<-c(0.1,2.6,3.9,5.2,6.6,9.8)
  hd<-c("Cohort","Role","N","C-index","HR (95% CI)","P")
  for(j in 1:6) text(cx[j],nrow(df)+0.6,hd[j],adj=ifelse(j==6,1,0),font=2,cex=0.9)
  abline(h=nrow(df)+0.3,lwd=0.8,col="gray50")
  for(i in 1:nrow(df)){r<-df[i,];y<-r$y
    text(cx[1],y,r$cohort,adj=0,cex=0.9,col=r$col,font=ifelse(r$role=="HO",2,1))
    text(cx[2],y,r$role,adj=0,cex=0.85,col=r$col);text(cx[3],y,r$n,adj=0,cex=0.88)
    text(cx[4],y,sprintf("%.3f",r$c_index),adj=0,cex=0.88,font=2)
    text(cx[5],y,sprintf("%.2f (%.2f-%.2f)",r$HR,r$HR_lo,r$HR_hi),adj=0,cex=0.88)
    text(cx[6],y,fp(r$KM_p),adj=1,cex=0.88,font=ifelse(r$KM_p<0.05,2,1),col=ifelse(r$KM_p<0.05,COL_HIGH,"#666666"))}
  par(mar=c(3.6,0.5,1.2,1.0));xl<-c(0.4,7);xb<-c(0.5,1,2,4,6)
  plot.new();plot.window(xlim=log(xl),ylim=c(-0.6,nrow(df)+1.0))
  text(log(2),nrow(df)+0.6,"HR (95% CI)",adj=0.5,font=2,cex=0.9);abline(h=nrow(df)+0.3,lwd=0.8,col="gray50")
  abline(v=log(1),lwd=1.2,lty=2,col="gray50");abline(v=log(c(0.5,2,4,6)),lwd=0.4,lty=3,col="gray85")
  for(i in 1:nrow(df)){r<-df[i,];y<-r$y
    segments(log(r$HR_lo),y,log(r$HR_hi),y,lwd=1.8,col=r$col)
    segments(log(r$HR_lo),y-0.13,log(r$HR_lo),y+0.13,lwd=1.8,col=r$col)
    segments(log(r$HR_hi),y-0.13,log(r$HR_hi),y+0.13,lwd=1.8,col=r$col)
    points(log(r$HR),y,pch=15,cex=r$cex,col=r$col)}
  axis(1,log(xb),xb,cex.axis=0.85,padj=-0.6);mtext("Hazard ratio (log scale)",side=1,line=1.9,cex=0.9)
}
save_base(function() draw_forest(summ), file.path(PFIG,"MF2_forest_8cohort"), w_mm=175, h_mm=80)

## ---- 8 ROC ----
rc<-c(COL_LOW,"#009E73",COL_HIGH)
draw_roc<-function(cn,lab){d<-compute_rs(cn);d<-d[!is.na(d$RS)&!is.na(d$OS_time)&d$OS_time>0,]
  tms<-c(12,36,60);tms<-tms[tms<=max(d$OS_time)]
  par(mar=c(3.4,3.6,1.6,0.8),mgp=c(2,0.6,0),tcl=-0.25)
  plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="1 - Specificity",ylab="Sensitivity",las=1)
  mtext(lab,side=3,line=0.4,font=2,cex=1.0);abline(0,1,lty=2,col="gray60")
  roc<-tryCatch(timeROC::timeROC(T=d$OS_time,delta=d$OS_status,marker=d$RS,cause=1,weighting="marginal",times=tms,iid=FALSE),error=f
