#!/usr/bin/env Rscript
# W6 panel 批次3: MF6-J 三 KO 下游程序趋同 Venn (NRP1/EGFR/BMPR1B 显著基因)
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; PFIG<-file.path(ROOT,"03_figures/panels")
KOD<-file.path(ROOT,"02_results/10_knockout/Results_final")
gset<-function(f) unique(fread(file.path(KOD,f))$gene)
A<-gset("NRP1_sig_adjP005.csv"); B<-gset("EGFR_sig_adjP005.csv"); C<-gset("BMPR1B_sig_adjP005.csv")
nABC<-length(Reduce(intersect,list(A,B,C)))
nAB<-length(intersect(A,B)); nAC<-length(intersect(A,C)); nBC<-length(intersect(B,C))
rA<-length(A); rB<-length(B); rC<-length(C)
oA<-rA-nAB-nAC+nABC; oB<-rB-nAB-nBC+nABC; oC<-rC-nAC-nBC+nABC
oAB<-nAB-nABC; oAC<-nAC-nABC; oBC<-nBC-nABC
cat(sprintf("NRP1=%d EGFR=%d BMPR1B=%d | AB=%d AC=%d BC=%d ABC=%d\n",rA,rB,rC,nAB,nAC,nBC,nABC))

drawcirc<-function(x,y,r,col){a<-seq(0,2*pi,length.out=200)
  polygon(x+r*cos(a),y+r*sin(a),col=adjustcolor(col,0.30),border=col,lwd=1.5)}
draw<-function(){
  par(mar=c(0.4,0.4,0.4,0.4))
  plot(NA,xlim=c(-1.9,1.9),ylim=c(-2.0,1.7),axes=FALSE,xlab="",ylab="",asp=1)
  drawcirc(-0.5,0.35,0.95,COL_HIGH); drawcirc(0.5,0.35,0.95,COL_LOW); drawcirc(0,-0.55,0.95,"#55A868")
  t<-function(x,y,n,cx=1) text(x,y,n,cex=cx)
  t(-1.0,0.72,oA); t(1.0,0.72,oB); t(0,-1.25,oC)
  t(0,0.78,oAB); t(-0.64,-0.30,oAC); t(0.64,-0.30,oBC); t(0,0.03,nABC,1.2)
  text(-1.2,1.3,"NRP1 KO",font=2,col=COL_HIGH); text(1.2,1.3,"EGFR KO",font=2,col=COL_LOW)
  text(0,-1.65,"BMPR1B KO",font=2,col="#55A868")
}
save_base(draw, file.path(PFIG,"MF6J_KO_convergence_venn"), w_mm=92, h_mm=92)
cat("\n批次3 完成: MF6J\n")
