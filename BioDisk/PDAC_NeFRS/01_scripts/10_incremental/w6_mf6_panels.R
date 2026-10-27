#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
TFD<-file.path(ROOT,"02_results/11_tf"); KOD<-file.path(ROOT,"02_results/10_knockout/Results_final")
CCD<-file.path(ROOT,"02_results/12_cellchat"); PFIG<-file.path(ROOT,"03_figures/panels")

## A TF activity (nominal; FDR n.s.)
TFd<-fread(file.path(TFD,"TF_activity_NeFRS_hilo_diff.csv"))
TFd[,neglogp:=-log10(pval)]
TFd[,sig:=fifelse(pval<0.05 & diff>0,"up",fifelse(pval<0.05 & diff<0,"down","ns"))]
hi<-c("WWTR1","INSM1","NR2F1","GLIS3","NFIB","ZEB1","NEUROD1")
lab<-TFd[TF %in% hi]
pA<-ggplot(TFd,aes(diff,neglogp))+geom_point(aes(colour=sig),size=0.4,alpha=0.5,stroke=0)+
  scale_colour_manual(values=c(up=COL_HIGH,down=COL_LOW,ns="grey80"),guide="none")+
  geom_vline(xintercept=0,linewidth=0.3,colour="grey70")+
  geom_hline(yintercept=-log10(0.05),linetype=2,linewidth=0.4,colour="grey55")+
  geom_point(data=lab,colour="black",size=0.9)+
  geom_text_repel(data=lab,aes(label=TF),size=5.4/.pt,max.overlaps=20,min.segment.length=0,segment.size=0.2)+
  labs(x="TF activity diff (high - low NeFRS)",y="-log10 nominal P",title="TF activity (FDR n.s.)")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"))
save_panel(pA,file.path(PFIG,"MF6_A_tf"),w_mm=58,h_mm=50)

## B-D receptor KO differential regulation
mkko<-function(g,stub){d<-fread(file.path(KOD,sprintf("%s_diffRegulation.csv",g)))
  d<-d[gene!=g]; d[,logp:=-log10(p.adj)]; d[is.infinite(logp)|logp>300,logp:=300]
  d[,sigf:=p.adj<0.05]; top<-d[order(-Z)][1:8]
  p<-ggplot(d,aes(Z,logp))+geom_point(aes(colour=sigf),size=0.3,alpha=0.4,stroke=0)+
    scale_colour_manual(values=c("TRUE"=COL_HIGH,"FALSE"="grey80"),guide="none")+
    geom_text_repel(data=top,aes(label=gene),size=4.8/.pt,max.overlaps=15,min.segment.length=0,segment.size=0.2)+
    labs(x="Differential regulation Z",y="-log10 padj",title=sprintf("%s knockout",g))+
    theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"))
  save_panel(p,file.path(PFIG,stub),w_mm=58,h_mm=50)}
mkko("NRP1","MF6_B_NRP1"); mkko("EGFR","MF6_C_EGFR"); mkko("BMPR1B","MF6_D_BMPR1B")

## E shared-downstream Venn (hand-drawn from overlap counts)
draw_venn<-function(){par(mar=c(0.3,0.3,1.6,0.3))
  plot(0,0,type="n",xlim=c(0,10),ylim=c(0,9.5),axes=FALSE,xlab="",ylab="",asp=1)
  symbols(c(3.9,6.1,5),c(6.2,6.2,4.1),circles=rep(2.3,3),inches=FALSE,add=TRUE,
    fg=c(COL_HIGH,COL_LOW,"#009E73"),lwd=2)
  text(2.7,6.6,"28",cex=0.8);text(7.3,6.6,"22",cex=0.8);text(5,2.5,"18",cex=0.8)
  text(5,7.1,"2",cex=0.75);text(3.5,4.4,"11",cex=0.75);text(6.5,4.4,"7",cex=0.75)
  text(5,5.2,"23",font=2,cex=1.0)
  text(2.6,8.7,"NRP1",col=COL_HIGH,font=2,cex=0.8);text(7.4,8.7,"EGFR",col=COL_LOW,font=2,cex=0.8)
  text(5,0.9,"BMPR1B",col="#009E73",font=2,cex=0.8)
  mtext("Shared KO downstream",3,line=0.2,font=2,cex=0.72)}
save_base(draw_venn,file.path(PFIG,"MF6_E_venn"),w_mm=58,h_mm=50)

## F shared 23 downstream genes (Z in NRP1 KO)
readsig<-function(g){d<-fread(file.path(KOD,sprintf("%s_sig_adjP005.csv",g)))
  gc<-intersect(c("gene","Gene","feature","x"),colnames(d)); if(length(gc)==0) gc<-colnames(d)[1]; d[[gc[1]]]}
shared<-Reduce(intersect,list(readsig("NRP1"),readsig("EGFR"),readsig("BMPR1B")))
cat("shared genes n=",length(shared),"\n")
nd<-fread(file.path(KOD,"NRP1_diffRegulation.csv"))[gene %in% shared][order(-Z)]
topF<-nd[1:min(15,.N)]; topF[,gene:=factor(gene,levels=rev(gene))]
pF<-ggplot(topF,aes(Z,gene))+geom_col(fill="#009E73",width=0.7)+
  labs(x="Differential regulation Z (NRP1 KO)",y=NULL,
    title=sprintf("Shared downstream (%d genes)",length(shared)))+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=5.6))
save_panel(pF,file.path(PFIG,"MF6_F_shared"),w_mm=58,h_mm=50)

## G CellChat pathway differential (high vs low NeFRS)
CC<-fread(file.path(CCD,"cellchat_B_pathway_diff.tsv"))
CC<-CC[order(-abs(diff))][1:15]; CC[,pathway:=factor(pathway,levels=rev(pathway))]
pG<-ggplot(CC,aes(diff,pathway,fill=diff>0))+geom_col(width=0.72)+
  geom_vline(xintercept=0,linewidth=0.3,colour="grey60")+
  scale_fill_manual(values=c("TRUE"=COL_HIGH,"FALSE"=COL_LOW),
    labels=c("Up in low-NeFRS","Up in high-NeFRS"),name=NULL)+
  labs(x="Interaction-strength diff (high - low NeFRS)",y=NULL,title="CellChat signaling pathways")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=5.8),
    legend.position="bottom",legend.key.size=unit(2.6,"mm"),legend.text=element_text(size=5.6))
save_panel(pG,file.path(PFIG,"MF6_G_cellchat"),w_mm=120,h_mm=52)
cat("MF6 A-G 完成\n")
