#!/usr/bin/env Rscript
# MF6J 微调: 类别标签移出圆外, 其余不变
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; PFIG<-file.path(ROOT,"03_figures/panels")
KOD<-file.path(ROOT,"02_results/10_knockout/Results_final")
gs<-function(f) unique(fread(file.path(KOD,f))$gene)
A<-gs("NRP1_sig_adjP005.csv");B<-gs("EGFR_sig_adjP005.csv");C<-gs("BMPR1B_sig_adjP005.csv")
nABC<-length(Reduce(intersect,list(A,B,C)));nAB<-length(intersect(A,B));nAC<-length(intersect(A,C));nBC<-length(intersect(B,C))
oA<-length(A)-nAB-nAC+nABC;oB<-length(B)-nAB-nBC+nABC;oC<-length(C)-nAC-nBC+nABC
oAB<-nAB-nABC;oAC<-nAC-nABC;oBC<-nBC-nABC
dc<-function(x,y,r,col){a<-seq(0,2*pi,length.out=200);polygon(x+r*cos(a),y+r*sin(a),col=adjustcolor(col,0.30),border=col,lwd=1.5)}
draw<-function(){par(mar=c(0.3,0.3,0.3,0.3));plot(NA,xlim=c(-2.1,2.1),ylim=c(-2.1,2.05),axes=FALSE,xlab="",ylab="",asp=1)
  dc(-0.5,0.35,0.95,COL_HIGH);dc(0.5,0.35,0.95,COL_LOW);dc(0,-0.55,0.95,"#55A868")
  tt<-function(x,y,n,cx=1)text(x,y,n,cex=cx)
  tt(-1,0.72,oA);tt(1,0.72,oB);tt(0,-1.25,oC);tt(0,0.78,oAB);tt(-0.64,-0.3,oAC);tt(0.64,-0.3,oBC);tt(0,0.03,nABC,1.25)
  text(-1.25,1.62,"NRP1 KO",font=2,col=COL_HIGH)
  text( 1.25,1.62,"EGFR KO",font=2,col=COL_LOW)
  text(0,-1.72,"BMPR1B KO",font=2,col="#55A868")}
save_base(draw,file.path(PFIG,"MF6J_KO_convergence_venn"),w_mm=95,h_mm=95)
cat("MF6J 标签已移出圆外\n")
