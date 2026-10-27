#!/usr/bin/env Rscript
# FigS2: 8 队列 1/3/5 年 time-dependent ROC (2x4, WYSIWYG 43x40mm)
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(survival); library(glmnet); library(timeROC) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
W2<-file.path(ROOT,"02_results/02_w2_inputs"); MIME<-file.path(ROOT,"02_results/04_mime")
PFIG<-file.path(ROOT,"03_figures/panels")
cohort_data<-readRDS(file.path(W2,"cohort_data_8.rds"))
mi<-readRDS(file.path(MIME,"mime_input_list.rds")); rr<-readRDS(file.path(MIME,"mime_results.rds"))
mg<-colnames(mi[[1]])[-(1:3)]; ml<-rr$ml.res[["Enet[α=0.1]"]]
compute_rs<-function(cn){co<-cohort_data[[cn]]
  X<-matrix(0,ncol(co$expression),length(mg),dimnames=list(colnames(co$expression),mg))
  for(g in intersect(mg,rownames(co$expression)))X[,g]<-co$expression[g,]
  rs<-as.numeric(predict(ml,newx=X,s="lambda.min")); names(rs)<-colnames(co$expression)
  cl<-co$clinical; data.frame(OS_time=cl$os_time_month,OS_status=cl$os_status,RS=rs[match(cl$sample_id,names(rs))])}
co<-list(c("PACA_AU","PACA-AU (n=267)"),c("TCGA_PAAD","TCGA-PAAD (n=147)"),
  c("E_MTAB_6134","E-MTAB-6134 (n=288)"),c("GSE71729","GSE71729 (n=123)"),
  c("GSE85916","GSE85916 (n=79)"),c("GSE62452","GSE62452 (n=66)"),
  c("GSE57495","GSE57495 (n=63)"),c("GSE21501","GSE21501 (n=102)"))
rc<-c("#2166AC","#D6604D","#009E73"); yr<-c("1y","3y","5y"); tt<-c(12,36,60)

draw_roc<-function(d,lab){
  d<-d[!is.na(d$RS)&!is.na(d$OS_time)&d$OS_time>0&!is.na(d$OS_status),]
  tr<-tryCatch(timeROC(T=d$OS_time,delta=d$OS_status,marker=d$RS,cause=1,times=tt,iid=FALSE),
               error=function(e) NULL)
  par(mar=c(2.3,2.9,1.4,0.6),mgp=c(1.4,0.4,0),tcl=-0.2,lwd=0.7)
  plot(NA,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  abline(0,1,lty=2,lwd=0.6,col="gray60")
  leg<-character(0)
  for(j in 1:3){
    if(is.null(tr)||all(is.na(tr$TP[,j]))){leg<-c(leg,sprintf("%s  NA",yr[j]));next}
    fp<-tr$FP[,j]; tp<-tr$TP[,j]; o<-order(fp)
    lines(fp[o],tp[o],col=rc[j],lwd=1.4)
    leg<-c(leg,sprintf("%s  %.2f",yr[j],tr$AUC[j]))
  }
  axis(1,at=c(0,.5,1),cex.axis=0.6,lwd=0.6,padj=-0.6)
  axis(2,at=c(0,.5,1),las=2,cex.axis=0.6,lwd=0.6)
  mtext(lab,3,line=0.3,font=2,cex=0.62)
  mtext("1 - Specificity",1,line=1.25,cex=0.55)
  mtext("Sensitivity",2,line=1.9,cex=0.55,las=0)
  legend("bottomright",legend=leg,col=rc,lwd=1.4,bty="n",cex=0.52,seg.len=0.9,
         y.intersp=0.95,title="AUC",title.adj=0)
}
for(x in co){ d<-compute_rs(x[1])
  save_base(function() draw_roc(d,x[2]),file.path(PFIG,paste0("FigS2_ROC_",x[1])),w_mm=43,h_mm=40) }
cat("FigS2 panels 完成: 8 timeROC\n")
