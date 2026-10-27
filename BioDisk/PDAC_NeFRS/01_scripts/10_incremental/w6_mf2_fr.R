#!/usr/bin/env Rscript
# MF2 森林(列距修复) + 8 ROC
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(survival); library(glmnet); library(timeROC) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
W2<-file.path(ROOT,"02_results/02_w2_inputs"); MIME<-file.path(ROOT,"02_results/04_mime"); PFIG<-file.path(ROOT,"03_figures/panels")
df<-fread(file.path(PFIG,"MF2_summary.tsv"))
fp<-function(p){if(p<1e-15)return("< 1e-15");if(p<1e-4)return(sprintf("%.1e",p));if(p<0.001)return(sprintf("%.4f",p));sprintf("%.3f",p)}
draw_forest<-function(df){
  df$rf<-factor(df$role,levels=c("TRAIN","SEL","HO")); df<-df[order(df$rf,-df$HR),]; df$y<-nrow(df):1
  df$col<-ifelse(df$role=="TRAIN",COL_HIGH,ifelse(df$role=="HO",COL_LOW,"#666666"))
  df$cz<-sqrt(df$n)/sqrt(max(df$n))*1.4+0.5
  layout(matrix(c(1,2),nrow=1),widths=c(6.3,3.7))
  par(mar=c(3.6,0.5,1.2,0.2));plot.new();plot.window(xlim=c(0,11.9),ylim=c(-0.5,nrow(df)+1.0))
  cx<-c(0.1,2.9,4.1,5.25,6.9,11.6); hd<-c("Cohort","Role","N","C-index","HR (95% CI)","P")
  for(j in 1:6) text(cx[j],nrow(df)+0.6,hd[j],adj=ifelse(j==6,1,0),font=2,cex=0.9)
  abline(h=nrow(df)+0.3,lwd=0.8,col="gray50")
  for(i in 1:nrow(df)){r<-df[i,];y<-r$y
    text(cx[1],y,r$cohort,adj=0,cex=0.9,col=r$col,font=ifelse(r$role=="HO",2,1))
    text(cx[2],y,r$role,adj=0,cex=0.85,col=r$col);text(cx[3],y,r$n,adj=0,cex=0.88)
    text(cx[4],y,sprintf("%.3f",r$c_index),adj=0,cex=0.88,font=2)
    text(cx[5],y,sprintf("%.2f (%.2f-%.2f)",r$HR,r$HR_lo,r$HR_hi),adj=0,cex=0.88)
    text(cx[6],y,fp(r$KM_p),adj=1,cex=0.88,font=ifelse(r$KM_p<0.05,2,1),col=ifelse(r$KM_p<0.05,COL_HIGH,"#666666"))}
  par(mar=c(3.6,0.5,1.2,1.0));xb<-c(0.5,1,2,4,6)
  plot.new();plot.window(xlim=log(c(0.4,7)),ylim=c(-0.5,nrow(df)+1.0))
  text(log(2),nrow(df)+0.6,"HR (95% CI)",adj=0.5,font=2,cex=0.9);abline(h=nrow(df)+0.3,lwd=0.8,col="gray50")
  abline(v=log(1),lwd=1.2,lty=2,col="gray50");abline(v=log(c(0.5,2,4,6)),lwd=0.4,lty=3,col="gray85")
  for(i in 1:nrow(df)){r<-df[i,];y<-r$y
    segments(log(r$HR_lo),y,log(r$HR_hi),y,lwd=1.8,col=r$col)
    segments(log(r$HR_lo),y-0.13,log(r$HR_lo),y+0.13,lwd=1.8,col=r$col)
    segments(log(r$HR_hi),y-0.13,log(r$HR_hi),y+0.13,lwd=1.8,col=r$col)
    points(log(r$HR),y,pch=15,cex=r$cz,col=r$col)}
  axis(1,log(xb),xb,cex.axis=0.85,padj=-0.6);mtext("Hazard ratio (log scale)",side=1,line=1.9,cex=0.9)
}
save_base(function() draw_forest(df), file.path(PFIG,"MF2_forest_8cohort"), w_mm=178, h_mm=80)

cohort_data<-readRDS(file.path(W2,"cohort_data_8.rds"))
mi<-readRDS(file.path(MIME,"mime_input_list.rds")); rr<-readRDS(file.path(MIME,"mime_results.rds"))
mg<-colnames(mi[[1]])[-(1:3)]; ml<-rr$ml.res[["Enet[α=0.1]"]]
compute_rs<-function(cn){co<-cohort_data[[cn]]
  X<-matrix(0,ncol(co$expression),length(mg),dimnames=list(colnames(co$expression),mg))
  for(g in intersect(mg,rownames(co$expression)))X[,g]<-co$expression[g,]
  rs<-as.numeric(predict(ml,newx=X,s="lambda.min")); names(rs)<-colnames(co$expression)
  cl<-co$clinical; data.frame(OS_time=cl$os_time_month,OS_status=cl$os_status,RS=rs[match(cl$sample_id,names(rs))])}
co<-list(c("PACA_AU","PACA-AU (n=267)"),c("GSE71729","GSE71729 (n=123)"),c("GSE85916","GSE85916 (n=79)"),
  c("GSE62452","GSE62452 (n=66)"),c("GSE57495","GSE57495 (n=63)"),c("GSE21501","GSE21501 (n=102)"),
  c("TCGA_PAAD","TCGA-PAAD (n=147)"),c("E_MTAB_6134","E-MTAB-6134 (n=288)"))
rc<-c(COL_LOW,"#009E73",COL_HIGH)
draw_roc<-function(cn,lab){d<-compute_rs(cn);d<-d[!is.na(d$RS)&!is.na(d$OS_time)&d$OS_time>0,]
  tms<-c(12,36,60);tms<-tms[tms<=max(d$OS_time)]
  par(mar=c(3.5,3.7,2.0,0.8),mgp=c(2.1,0.6,0),tcl=-0.25)
  plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="1 - Specificity",ylab="Sensitivity",las=1)
  mtext(lab,side=3,line=0.5,font=2,cex=1.05);abline(0,1,lty=2,col="gray60")
  roc<-tryCatch(timeROC::timeROC(T=d$OS_time,delta=d$OS_status,marker=d$RS,cause=1,weighting="marginal",times=tms,iid=FALSE),error=function(e)NULL)
  if(is.null(roc)){text(0.5,0.5,"insufficient");return(invisible())}
  ll<-c();for(j in seq_along(tms)){lines(roc$FP[,j],roc$TP[,j],col=rc[j],lwd=2.2)
    yl<-c("12"="1-year","36"="3-year","60"="5-year")[as.character(tms[j])]
    ll<-c(ll,sprintf("%s AUC=%.3f",yl,roc$AUC[j]))}
  legend("bottomright",ll,col=rc[seq_along(tms)],lwd=2.2,bty="n",cex=0.85)}
for(x in co) save_base(function() draw_roc(x[1],x[2]), file.path(PFIG,paste0("MF2_ROC_",x[1])), w_mm=72, h_mm=72)
cat("MF2 森林 + 8 ROC 完成\n")
