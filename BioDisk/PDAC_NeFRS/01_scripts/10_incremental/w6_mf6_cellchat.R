#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(CellChat); library(data.table); library(ggplot2) })
options(bitmapType="cairo")
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; CC<-file.path(ROOT,"02_results/12_cellchat"); PFIG<-file.path(ROOT,"03_figures/panels")
cc_A<-readRDS(file.path(CC,"cellchat_A_allcelltypes.rds"))
vw<-as.numeric(table(cc_A@meta$labels)); names(vw)<-names(table(cc_A@meta$labels))
save_base(function(){ par(mar=c(2,2,2,2))
  netVisual_circle(cc_A@net$count, vertex.weight=vw, weight.scale=TRUE, label.edge=FALSE, vertex.label.cex=0.65) },
  file.path(PFIG,"MF6_cellchat_count"), w_mm=96, h_mm=96)
save_base(function(){ par(mar=c(2,2,2,2))
  netVisual_circle(cc_A@net$weight, vertex.weight=vw, weight.scale=TRUE, label.edge=FALSE, vertex.label.cex=0.65) },
  file.path(PFIG,"MF6_cellchat_weight"), w_mm=96, h_mm=96)

neuro<-grep("^(NRXN|NLGN|NEGR|NCAM|SEMA|EPH|SLIT|GDNF|NGF|NTN|BDNF|CADM|PSAP)", cc_A@netP$pathways, value=TRUE)
cat("神经通路(已剔除WNT等误匹配):",paste(neuro,collapse=", "),"\n")
if(length(neuro)>0){
  ctl<-levels(cc_A@meta$labels); mal_i<-which(ctl=="Epithelial (malignant)")
  pb<-tryCatch(
    netVisual_bubble(cc_A, sources.use=1:length(ctl), targets.use=mal_i, signaling=neuro, remove.isolate=TRUE)+
      scale_x_discrete(labels=function(x) sub("\\s*->.*","",x))+
      theme_pub()+labs(title=NULL, x=NULL, caption="Target: Epithelial (malignant)")+
      theme(axis.text.x=element_text(angle=45,hjust=1)),
    error=function(e){cat("bubble fail:",conditionMessage(e),"\n");NULL})
  if(!is.null(pb)) save_panel(pb, file.path(PFIG,"MF6_cellchat_neural_bubble"), w_mm=110, h_mm=88)
}
cat("MF6 CellChat 完成\n")
