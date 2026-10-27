#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(clusterProfiler); library(org.Hs.eg.db) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; KOD<-file.path(ROOT,"02_results/10_knockout/Results_final"); PFIG<-file.path(ROOT,"03_figures/panels")
rs<-function(g){d<-fread(file.path(KOD,sprintf("%s_sig_adjP005.csv",g))); d[[1]]}
genes<-unique(c(rs("NRP1"),rs("EGFR"),rs("BMPR1B")))
shared<-Reduce(intersect,list(rs("NRP1"),rs("EGFR"),rs("BMPR1B")))
bg<-unique(fread(file.path(KOD,"NRP1_diffRegulation.csv"))$gene)
cat("并集:",length(genes)," 共享:",length(shared)," 背景:",length(bg),"\n")
eg<-bitr(genes,"SYMBOL","ENTREZID",org.Hs.eg.db); bge<-bitr(bg,"SYMBOL","ENTREZID",org.Hs.eg.db)
cat("映射成功 并集:",nrow(eg),"/",length(genes)," 背景:",nrow(bge),"/",length(bg),"\n")
ego<-enrichGO(eg$ENTREZID,OrgDb=org.Hs.eg.db,ont="BP",universe=bge$ENTREZID,
  pAdjustMethod="BH",pvalueCutoff=0.1,qvalueCutoff=0.2,readable=TRUE,minGSSize=5,maxGSSize=500)
res<-as.data.table(ego@result)
cat("=== padj<0.1 条目数:",sum(res$p.adjust<0.1)," | pval<0.05:",sum(res$pvalue<0.05),"===\n")
print(res[p.adjust<0.1][order(p.adjust)][1:min(15,.N),.(Description,Count,GeneRatio,p.adjust=signif(p.adjust,2))])
fwrite(res,file.path(ROOT,"02_results/10_knockout/Results_final/shared_KO_GO_BP.tsv"),sep="\t")
## panel
dd<-res[p.adjust<0.1][order(p.adjust)][1:min(12,.N)]
if(nrow(dd)>0){
  dd[,gr:=as.numeric(sub("/.*","",GeneRatio))/as.numeric(sub(".*/","",GeneRatio))]
  dd[,Description:=factor(Description,levels=rev(Description))]
  p<-ggplot(dd,aes(gr,Description))+geom_point(aes(size=Count,colour=p.adjust))+
    scale_colour_gradient(low=COL_HIGH,high="grey75",name="padj")+scale_size(range=c(1.2,4),name="Count")+
    labs(x="Gene ratio",y=NULL,title="Shared KO downstream — GO-BP")+
    theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=5.4),
      legend.key.size=unit(2.6,"mm"),legend.title=element_text(size=5.5),legend.text=element_text(size=5))
  save_panel(p,file.path(PFIG,"MF6_E_enrich"),w_mm=116,h_mm=66)
  cat("富集 panel 完成\n")
} else cat("!! 无 padj<0.1 条目, 需换策略(KEGG/各KO共识/放宽阈值)\n")
