#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
TFD<-file.path(ROOT,"02_results/11_tf"); KOD<-file.path(ROOT,"02_results/10_knockout/Results_final")
CCD<-file.path(ROOT,"02_results/12_cellchat"); PFIG<-file.path(ROOT,"03_figures/panels")

## A TF activity volcano (correct table; honest grouping)
TFd<-fread(file.path(TFD,"TF_highlow_diff.tsv"))
TFd[,neglogp:=-log10(padj)]; TFd[is.infinite(neglogp)|neglogp>320,neglogp:=320]
TFd[,sig:=fifelse(padj<0.05 & diff>0.05,"up",fifelse(padj<0.05 & diff< -0.05,"down","ns"))]
stress<-c("NFE2L3","BACH1","FOXM1","E2F1","MAX","PROX1"); matureN<-c("TRPS1","INSM1","NEUROD1"); progN<-c("NR2F1","GLIS3","NFIB")
lab<-TFd[TF %in% c(stress,matureN,progN)]
lab[,grp:=fifelse(TF %in% stress,"Stress / proliferation",fifelse(TF %in% matureN,"Mature neuronal","Neural progenitor"))]
pA<-ggplot(TFd,aes(diff,neglogp))+geom_point(aes(colour=sig),size=0.4,alpha=0.32,stroke=0)+
  scale_colour_manual(values=c(up=COL_HIGH,down=COL_LOW,ns="grey82"),guide="none")+
  geom_vline(xintercept=0,linewidth=0.3,colour="grey70")+geom_hline(yintercept=-log10(0.05),linetype=2,linewidth=0.4,colour="grey55")+
  geom_point(data=lab,aes(fill=grp),colour="black",size=1.2,shape=21,stroke=0.3)+
  scale_fill_manual(values=c("Stress / proliferation"=COL_HIGH,"Mature neuronal"=COL_LOW,"Neural progenitor"="grey55"),name=NULL)+
  geom_text_repel(data=lab,aes(label=TF),size=5/.pt,max.overlaps=30,min.segment.length=0,segment.size=0.2)+
  labs(x="TF activity diff (high - low 1/3 NeFRS)",y="-log10 padj",title="TF activity (decoupleR)")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),legend.position="bottom",
    legend.key.size=unit(2.4,"mm"),legend.text=element_text(size=5),legend.margin=margin(0,0,0,0))
save_panel(pA,file.path(PFIG,"MF6_A_tf"),w_mm=78,h_mm=56)

## shared 23 downstream
readsig<-function(g){d<-fread(file.path(KOD,sprintf("%s_sig_adjP005.csv",g)))
  gc<-intersect(c("gene","Gene","feature","x"),colnames(d)); if(length(gc)==0) gc<-colnames(d)[1]; d[[gc[1]]]}
shared<-Reduce(intersect,list(readsig("NRP1"),readsig("EGFR"),readsig("BMPR1B")))
cat("shared n=",length(shared),"\n")
getZ<-function(g){d<-fread(file.path(KOD,sprintf("%s_diffRegulation.csv",g)))[gene %in% shared,.(gene,Z)]; d[,recep:=g]; d}
L<-rbindlist(list(getZ("NRP1"),getZ("EGFR"),getZ("BMPR1B")))
ord<-L[,.(m=mean(Z)),by=gene][order(m)]$gene
L[,gene:=factor(gene,levels=ord)]; L[,recep:=factor(recep,levels=c("NRP1","EGFR","BMPR1B"))]

## B convergence heatmap (23 shared x 3 receptors)
pB<-ggplot(L,aes(recep,gene,fill=Z))+geom_tile(colour="white",linewidth=0.5)+
  scale_fill_gradient(low="#FDEAE0",high=COL_HIGH,name="Z")+
  labs(x=NULL,y=NULL,title="Shared downstream (Z)")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),
    axis.text.y=element_text(size=5),axis.text.x=element_text(size=5.6),
    legend.key.width=unit(2.2,"mm"),legend.key.height=unit(3,"mm"),
    legend.title=element_text(size=5.5),legend.text=element_text(size=4.8),panel.grid=element_blank())
save_panel(pB,file.path(PFIG,"MF6_B_heat"),w_mm=46,h_mm=66)

## C Venn (hand-drawn)
draw_venn<-function(){par(mar=c(0.3,0.3,1.6,0.3))
  plot(0,0,type="n",xlim=c(0,10),ylim=c(0,9.5),axes=FALSE,xlab="",ylab="",asp=1)
  symbols(c(3.9,6.1,5),c(6.2,6.2,4.1),circles=rep(2.3,3),inches=FALSE,add=TRUE,fg=c(COL_HIGH,COL_LOW,"#009E73"),lwd=2)
  text(2.7,6.6,"28",cex=0.8);text(7.3,6.6,"22",cex=0.8);text(5,2.5,"18",cex=0.8)
  text(5,7.1,"2",cex=0.75);text(3.5,4.4,"11",cex=0.75);text(6.5,4.4,"7",cex=0.75);text(5,5.2,"23",font=2,cex=1.05)
  text(2.6,8.7,"NRP1",col=COL_HIGH,font=2,cex=0.82);text(7.4,8.7,"EGFR",col=COL_LOW,font=2,cex=0.82);text(5,0.9,"BMPR1B",col="#009E73",font=2,cex=0.82)
  mtext("Shared KO downstream",3,line=0.2,font=2,cex=0.72)}
save_base(draw_venn,file.path(PFIG,"MF6_C_venn"),w_mm=56,h_mm=56)

## D shared genes ranked by mean Z, coloured by functional class
nd<-L[,.(Z=mean(Z)),by=gene][order(-Z)][1:min(16,.N)]
neural<-c("LINGO2","NR2F1-AS1","FOXP2","NFIB","ROBO1","EPHA6"); emt<-c("CDH11","MTSS1","WWTR1","FNDC3B","COL4A5")
metab<-c("AKT3","PRKD1","PPP4R4","CP","GMDS")
nd[,cls:=fifelse(gene %in% neural,"Neural/axon",fifelse(gene %in% emt,"EMT/matrix",fifelse(gene %in% metab,"PI3K/metab","Other")))]
nd[,gene:=factor(gene,levels=rev(gene))]
pD<-ggplot(nd,aes(Z,gene,fill=cls))+geom_col(width=0.72)+
  scale_fill_manual(values=c("Neural/axon"="#009E73","EMT/matrix"=COL_HIGH,"PI3K/metab"="#E69F00","Other"="grey70"),name=NULL)+
  labs(x="Mean differential-regulation Z",y=NULL,title="Shared 23 downstream")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=5.4),
    legend.position="bottom",legend.key.size=unit(2.2,"mm"),legend.text=element_text(size=4.8),legend.margin=margin(0,0,0,0))
save_panel(pD,file.path(PFIG,"MF6_D_shared"),w_mm=64,h_mm=66)

## E CellChat pathways
CC<-fread(file.path(CCD,"cellchat_B_pathway_diff.tsv"))[order(-abs(diff))][1:15]
CC[,pathway:=factor(pathway,levels=rev(pathway))]
pE<-ggplot(CC,aes(diff,pathway,fill=diff>0))+geom_col(width=0.72)+geom_vline(xintercept=0,linewidth=0.3,colour="grey60")+
  scale_fill_manual(values=c("TRUE"=COL_HIGH,"FALSE"=COL_LOW),labels=c("Up in low-NeFRS","Up in high-NeFRS"),name=NULL)+
  labs(x="Interaction-strength diff (high - low NeFRS)",y=NULL,title="CellChat signaling pathways")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=5.8),
    legend.position="bottom",legend.key.size=unit(2.6,"mm"),legend.text=element_text(size=5.4))
save_panel(pE,file.path(PFIG,"MF6_E_cellchat"),w_mm=112,h_mm=66)
cat("v2 A-E 完成\n")
