#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table) })
MIME<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime"; PFIG<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures/panels"
nz<-fread(file.path(MIME,"nefrs_v3_enet_coefficients.tsv"))
P3COL<-"#8C8C8C"
pcol<-function(pr) ifelse(pr=="P1_TopPriority",COL_HIGH,ifelse(pr=="P2_HighPriority","#FF8C00",P3COL))
manual_leg<-function(){
  u<-par("usr"); xL<-u[1]+0.05*(u[2]-u[1]); yT<-u[3]+0.17*(u[4]-u[3])
  sw<-0.035*(u[2]-u[1]); sh<-0.008*(u[4]-u[3]); dy<-0.042*(u[4]-u[3])
  cols<-c(COL_HIGH,"#FF8C00",P3COL); labs<-c("P1 (3 src)","P2 (2 src)","P3 (1 src)")
  for(k in 1:3){ yy<-yT-(k-1)*dy
    rect(xL,yy-sh,xL+sw,yy+sh,col=cols[k],border=ifelse(k<3,"black",NA),lwd=0.5)
    text(xL+sw*1.7,yy,labs[k],adj=0,cex=0.55) }
}
draw_bar<-function(d,legmode,lab_subset=NULL){
  d$col<-pcol(d$priority); d$hub<-d$priority%in%c("P1_TopPriority","P2_HighPriority"); xr<-range(d$coef)*1.18
  showlab<- if(is.null(lab_subset)) rep(TRUE,nrow(d)) else d$gene %in% lab_subset
  par(mar=c(3.5,5.0,1.2,0.8),mgp=c(1.9,0.5,0),tcl=-0.25,font.lab=2,cex.lab=0.78)
  plot(NA,xlim=xr,ylim=c(0.5,nrow(d)+1.2),xlab="NeFRS coef. (beta)",ylab="",las=1,yaxt="n",cex.axis=0.72,font.axis=2)
  abline(v=0,lwd=1.0,col="gray40");abline(v=seq(-0.15,0.15,0.05),lty=3,lwd=0.4,col="gray85")
  for(i in 1:nrow(d)){r<-d[i,];bh<-ifelse(r$hub,0.42,0.34)
    rect(min(0,r$coef),r$y-bh,max(0,r$coef),r$y+bh,col=r$col,border=ifelse(r$hub,"black","white"),lwd=ifelse(r$hub,0.7,0.3))}
  for(i in 1:nrow(d)){ if(!showlab[i]) next; r<-d[i,]
    axis(2,at=r$y,labels=r$gene,las=2,cex.axis=ifelse(r$hub,0.72,0.62),font=ifelse(r$hub,2,1),
      col.axis=ifelse(r$hub,"black","gray25"),lwd=0,lwd.ticks=0)}
  text(xr[1]*0.95,nrow(d)+0.7,"<- Protective",cex=0.62,col=COL_LOW,adj=0,font=2)
  text(xr[2]*0.95,nrow(d)+0.7,"Risk ->",cex=0.62,col=COL_HIGH,adj=1,font=2)
  if(identical(legmode,"BL")) manual_leg()
  else legend(legmode,c("P1 (3 sources)","P2 (2 sources)","P3 (1 source)"),
    fill=c(COL_HIGH,"#FF8C00",P3COL),border=c("black","black",NA),bty="n",cex=0.7)
}
da<-nz[order(coef)]; da$y<-seq_len(nrow(da))
hub_g <- da$gene[da$priority %in% c("P1_TopPriority","P2_HighPriority")]
top10_g <- da$gene[order(-abs(da$coef))][1:10]; lab_sub <- union(hub_g, top10_g)
save_base(function() draw_bar(da,"topleft",lab_subset=lab_sub), file.path(PFIG,"MF1_coef_all"), w_mm=88, h_mm=max(120,nrow(da)*2.6))
dt<-nz[order(-abs(coef))][1:30][order(abs(coef))]; dt$y<-seq_len(nrow(dt))
save_base(function() draw_bar(dt,"BL"), file.path(PFIG,"MF1_coef_top30"), w_mm=58, h_mm=90)
cat("MF1 系数条形 完成(深色加大+增高)\n")
