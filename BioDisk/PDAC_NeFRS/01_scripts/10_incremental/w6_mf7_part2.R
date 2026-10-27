#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(enrichplot) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; GSE<-file.path(ROOT,"02_results/15_gsea"); PFIG<-file.path(ROOT,"03_figures/panels")

## D GO-BP full width re-render
GC<-fread(file.path(GSE,"gsea_consistent_3cohorts.tsv")); if("consistent" %in% colnames(GC)) GC<-GC[consistent==TRUE]
up<-GC[order(-mean_NES)][1:10]; dn<-GC[order(mean_NES)][1:10]
GG<-rbind(up,dn)[order(mean_NES)]; GG[,Description:=factor(Description,levels=Description)]
pD<-ggplot(GG,aes(mean_NES,Description,fill=mean_NES>0))+geom_vline(xintercept=0,linewidth=0.3,colour="grey70")+
  geom_col(width=0.72)+scale_fill_manual(values=c("TRUE"=COL_HIGH,"FALSE"=COL_LOW),labels=c("Down in high-NeFRS","Up in high-NeFRS"),name=NULL)+
  labs(x="Mean NES (3 cohorts)",y=NULL,title="GSEA GO-BP (consistent across 3 cohorts)")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=5.4),
    legend.position="bottom",legend.key.size=unit(2.6,"mm"),legend.text=element_text(size=5.4))
save_panel(pD,file.path(PFIG,"MF7_D_gsea"),w_mm=176,h_mm=58)

## E/F GSEA enrichment curves from rds
go<-readRDS(file.path(GSE,"gse_go_PACA_obj.rds"))
res<-as.data.table(go@result); cat("GSEA对象通路:",nrow(res),"\n")
prolif<-res[order(-NES)][1]
imm<-res[grepl("B cell receptor|antigen process|MHC|immune response|lymphocyte|leukocyte",Description,ignore.case=TRUE) & NES<0][order(NES)][1]
cat("增殖曲线:",prolif$Description,"NES",round(prolif$NES,2),"\n免疫曲线:",imm$Description,"NES",round(imm$NES,2),"\n")
shorten<-function(s) ifelse(nchar(s)>34,paste0(substr(s,1,32),"\u2026"),s)
pE<-gseaplot2(go,geneSetID=prolif$ID,title=shorten(prolif$Description),base_size=6,rel_heights=c(1.3,0.25,0.6),subplots=1:3)
ggsave(file.path(PFIG,"MF7_E_curveUp.png"),pE,width=58/25.4,height=54/25.4,dpi=600)
pF<-gseaplot2(go,geneSetID=imm$ID,title=shorten(imm$Description),base_size=6,rel_heights=c(1.3,0.25,0.6),subplots=1:3)
ggsave(file.path(PFIG,"MF7_F_curveDn.png"),pF,width=58/25.4,height=54/25.4,dpi=600)

## G KEGG
KG<-fread(file.path(GSE,"gsea_KEGG_PACA_AU.tsv")); cat("KEGG列:",paste(colnames(KG),collapse=","),"\n")
nescol<-intersect(c("NES","nes"),colnames(KG))[1]; descol<-intersect(c("Description","ID","pathway","Name"),colnames(KG))[1]
KG<-KG[order(-abs(get(nescol)))][1:14]; setnames(KG,c(nescol,descol),c("NES","Desc")); KG[,Desc:=factor(Desc,levels=rev(Desc))]
pG<-ggplot(KG,aes(NES,Desc,fill=NES>0))+geom_vline(xintercept=0,linewidth=0.3,colour="grey70")+geom_col(width=0.72)+
  scale_fill_manual(values=c("TRUE"=COL_HIGH,"FALSE"=COL_LOW),guide="none")+
  labs(x="NES (PACA-AU)",y=NULL,title="GSEA KEGG")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=5))
save_panel(pG,file.path(PFIG,"MF7_G_kegg"),w_mm=64,h_mm=54)
cat("MF7 D/E/F/G 完成\n")
