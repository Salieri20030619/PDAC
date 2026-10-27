#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(VennDiagram); library(grid); library(UpSetR) })
MIME<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime"; PFIG<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures/panels"
nz<-fread(file.path(MIME,"nefrs_v3_enet_coefficients.tsv"))
psrc<-function(s){if(is.na(s))return(character(0));unlist(strsplit(s,"\\|"))}

smap_v<-c("FerrDb_v2"="Ferroptosis","Hwang2022_NatGenet_SupTable2"="Hwang_NRP",
  "GO_BP_high_specificity"="Neural_GO","Literature_review_PDAC_PNI"="PDAC_PNI")
sets<-list(); for(sr in names(smap_v)) sets[[smap_v[[sr]]]]<-nz$gene[sapply(nz$Sources,function(s) sr%in%psrc(s))]
vc<-c("#E63946","#F4A261","#2A9D8F","#264653")
vp<-venn.diagram(x=sets, category.names=c("Ferr.","NRP","Neural","PNI"),
  filename=NULL, lwd=2, lty=1, fill=vc, alpha=0.45,
  cat.cex=0.62, cat.fontface="bold", cat.col=vc, cat.default.pos="outer",
  cat.dist=c(0.22,0.22,0.12,0.12), cat.pos=c(-15,15,-15,15),
  cex=0.78, fontface="bold", margin=0.10, disable.logging=TRUE)
save_base(function(){ grid::grid.newpage(); grid::grid.draw(vp) }, file.path(PFIG,"MF1_venn"), w_mm=54, h_mm=48)

smap_u<-c("FerrDb_v2"="Ferroptosis","Hwang2022_NatGenet"="Hwang_NRP",
  "GO_BP_high_specificity"="Neural_GO","Literature_review_PDAC_PNI"="PDAC_PNI")
ud<-data.frame(gene=nz$gene)
for(sr in names(smap_u)) ud[[smap_u[[sr]]]]<-as.integer(sapply(nz$Sources,function(s) sr%in%psrc(s)))
save_base(function(){ print(upset(ud, sets=c("Ferroptosis","Hwang_NRP","Neural_GO","PDAC_PNI"),
  order.by="freq", keep.order=TRUE, point.size=3.0, line.size=1.0,
  mainbar.y.label="NeFRS genes", sets.x.label="Genes per source",
  text.scale=c(1.6,1.3,1.5,1.2,1.5,1.3), mb.ratio=c(0.62,0.38),
  matrix.color="#1F4E79", main.bar.color=COL_LOW, sets.bar.color=COL_HIGH,
  shade.color="gray85", shade.alpha=0.4, matrix.dot.alpha=0.6)) },
  file.path(PFIG,"MF1_upset"), w_mm=155, h_mm=125)
cat("MF1 Venn + UpSet 完成 (UpSet 字号适中+加高)\n")
