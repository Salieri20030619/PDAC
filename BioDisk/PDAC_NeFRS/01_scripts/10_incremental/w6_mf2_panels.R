#!/usr/bin/env Rscript
# MF2 主图 panel：8 紧凑KM + 森林(7验证,无框) + C-index（WYSIWYG,紧凑版）
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(survival); library(glmnet) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
W2<-file.path(ROOT,"02_results/02_w2_inputs"); MIME<-file.path(ROOT,"02_results/04_mime")
META<-file.path(ROOT,"02_results/20_meta"); PFIG<-file.path(ROOT,"03_figures/panels")
cohort_data<-readRDS(file.path(W2,"cohort_data_8.rds"))
mi<-readRDS(file.path(MIME,"mime_input_list.rds")); rr<-readRDS(file.path(MIME,"mime_results.rds"))
mg<-colnames(mi[[1]])[-(1:3)]; ml<-rr$ml.res[["Enet[α=0.1]"]]
compute_rs<-function(cn){co<-cohort_data[[cn]]
  X<-matrix(0,ncol(co$expression),length(mg),dimnames=list(colnames(co$expression),mg))
  for(g in intersect(mg,rownames(co$expression)))X[,g]<-co$expression[g,]
  rs<-as.numeric(predict(ml,newx=X,s="lambda.min")); names(rs)<-colnames(co$expression)
  cl<-co$clinical; data.frame(OS_time=cl$os_time_month,OS_status=cl$os_status,RS=rs[match(cl$sample_id,names(rs))])}
stat1<-function(cn){d<-compute_rs(cn);d<-d[!is.na(d$RS)&!is.na(d$OS_time)&d$OS_time>0,]
  d$Group<-factor(ifelse(d$RS>median(d$RS,na.rm=TRUE),"High","Low"),levels=c("Low","High"))
  cox<-coxph(Surv(OS_time,OS_status)~Group,data=d);ci<-exp(confint(cox))
  list(fit=survfit(Surv(OS_time,OS_status)~Group,data=d),d=d,hr=exp(coef(cox)),lo=ci[1,1],hi=ci[1,2],
    pv=1-pchisq(survdiff(Surv(OS_time,OS_status)~Group,data=d)$chisq,1),
    C=as.numeric(concordance(Surv(OS_time,OS_status)~RS,data=d,reverse=TRUE)$concordance),
    nh=sum(d$Group=="High"),nl=sum(d$Group=="Low"))}
cols<-c(COL_LOW,COL_HIGH)
fmt_p2<-function(p){if(is.na(p))return("P=NA");if(p<0.001)return("P<0.001");sprintf("P=%.3f",p)}
co<-list(c("PACA_AU","PACA-AU (n=267)"),c("TCGA_PAAD","TCGA-PAAD (n=147)"),
  c("E_MTAB_6134","E-MTAB-6134 (n=288)"),c("GSE71729","GSE71729 (n=123)"),
  c("GSE85916","GSE85916 (n=79)"),c("GSE62452","GSE62452 (n=66)"),
  c("GSE57495","GSE57495 (n=63)"),c("GSE21501","GSE21501 (n=102)"))

draw_kmc<-function(s,lab){mt<-max(s$d$OS_time)
  par(mar=c(2.0,2.9,1.2,0.4),mgp=c(1.4,0.4,0),tcl=-0.2,lwd=0.7)
  plot(s$fit,col=cols,lwd=1.5,conf.int=FALSE,xlab="",ylab="",xlim=c(0,mt),ylim=c(0,1),
       xaxt="n",yaxt="n",mark.time=TRUE,cex=0.22)
  axis(1,cex.axis=0.6,lwd=0.6,padj=-0.6);axis(2,seq(0,1,.5),las=2,cex.axis=0.6,lwd=0.6)
  mtext(lab,3,line=0.25,font=2,cex=0.62);mtext("Months",1,line=1.2,cex=0.55)
  mtext("Survival",2,line=1.9,cex=0.55,las=0)
  legend("topright",legend=c(sprintf("HR %.2f",s$hr),fmt_p2(s$pv),sprintf("C %.2f",s$C)),
    bty="n",cex=0.56,y.intersp=0.95)
  legend("bottomleft",legend=c("Low","High"),col=cols,lwd=1.5,bty="n",cex=0.56,seg.len=0.8,y.intersp=0.9)}
for(x in co){s<-stat1(x[1]);save_base(function() draw_kmc(s,x[2]),file.path(PFIG,paste0("MF2_KMc_",x[1])),w_mm=43,h_mm=38)}

draw_forest<-function(){
  pc<-read.table(file.path(META,"meta_per_cohort.tsv"),header=TRUE,sep="\t",stringsAsFactors=FALSE)
  pl<-read.table(file.path(META,"meta_pooled.tsv"),header=TRUE,sep="\t",stringsAsFactors=FALSE)
  pc<-pc[!grepl("PACA",pc$cohort,ignore.case=TRUE),]
  pc<-pc[order(pc$HR,decreasing=TRUE),]
  labs<-gsub("_","-",pc$cohort); pool<-pl[pl$set=="val7",]
  n<-nrow(pc); ys<-n:1; xlim<-c(0.9,2.5)
  par(mar=c(2.2,6.0,0.9,6.0),mgp=c(1.4,0.4,0),tcl=-0.2,lwd=0.7)
  plot(NA,xlim=xlim,ylim=c(-1.2,n+0.65),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
  abline(v=1,lty=2,lwd=0.6,col="gray55")
  for(i in 1:n){y<-ys[i]
    segments(pc$lo[i],y,pc$hi[i],y,lwd=1.0)
    points(pc$HR[i],y,pch=15,cex=0.85,col="black")
    text(2.60,y,sprintf("%.2f (%.2f-%.2f)",pc$HR[i],pc$lo[i],pc$hi[i]),adj=0,xpd=NA,cex=0.55)}
  yc<--0.7;d<-pool
  polygon(c(d$lo,d$HR,d$hi,d$HR),c(yc,yc+0.32,yc,yc-0.32),col=COL_LOW,border=NA)
  text(2.60,yc,sprintf("%.2f (%.2f-%.2f)",d$HR,d$lo,d$hi),adj=0,xpd=NA,cex=0.55,font=2)
  axis(2,at=ys,labels=labs,las=1,tick=FALSE,cex.axis=0.6,hadj=1,line=-0.3)
  axis(2,at=yc,labels="Pooled (7 val)",las=1,tick=FALSE,cex.axis=0.6,font.axis=2,hadj=1,line=-0.3)
  axis(1,at=seq(1.0,2.5,0.5),cex.axis=0.6,lwd=0.6,padj=-0.6)
  mtext("HR per SD (95% CI)",1,line=1.4,cex=0.6,font=2)
  text(1.95,n+0.58,sprintf("Random-effects (7 val): I\u00b2 = %.0f%%",d$I2),cex=0.58,font=2,xpd=NA)}
save_base(draw_forest,file.path(PFIG,"MF2_forest"),w_mm=88,h_mm=46)

draw_cindex<-function(){
  sm<-read.table(file.path(PFIG,"MF2_summary.tsv"),header=TRUE,sep="\t",stringsAsFactors=FALSE)
  sm<-sm[order(sm$c_index),];n<-nrow(sm)
  par(mar=c(2.2,5.6,0.9,2.4),mgp=c(1.4,0.4,0),tcl=-0.2,lwd=0.7)
  plot(NA,xlim=c(0.5,0.80),ylim=c(0.5,n+0.5),xaxt="n",yaxt="n",xlab="",ylab="")
  abline(v=seq(0.5,0.8,0.05),lty=3,lwd=0.4,col="gray85")
  for(i in 1:n){y<-i;col<-ifelse(sm$role[i]=="TRAIN",COL_HIGH,"#4C72B0")
    segments(0.5,y,sm$c_index[i],y,col="gray75",lwd=0.9)
    points(sm$c_index[i],y,pch=19,cex=1.0,col=col)
    text(sm$c_index[i]+0.006,y,sprintf("%.3f",sm$c_index[i]),adj=0,cex=0.56,xpd=NA)}
  axis(2,at=1:n,labels=gsub("_","-",sm$cohort),las=1,tick=FALSE,cex.axis=0.62,hadj=1,line=-0.3)
  axis(1,at=seq(0.5,0.8,0.1),cex.axis=0.6,lwd=0.6,padj=-0.6)
  mtext("C-index",1,line=1.4,cex=0.6,font=2)}
save_base(draw_cindex,file.path(PFIG,"MF2_cindex"),w_mm=88,h_mm=46)
cat("MF2 panels 完成\n")
