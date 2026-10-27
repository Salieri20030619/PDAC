#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; DRG<-file.path(ROOT,"02_results/16_drug"); MUT<-file.path(ROOT,"02_results/17_mutation"); SIG<-file.path(ROOT,"02_results/13_sig_compare"); PFIG<-file.path(ROOT,"03_figures/panels")
clean<-function(x) gsub("_[0-9]+$","",x)

## A drug sensitivity consistent (top sensitive + resistant)
DC<-fread(file.path(DRG,"drug_consistent.tsv")); cat("drug_consistent列:",paste(colnames(DC),collapse=","),"\n")
fc<-intersect(c("mean_log2fc","log2fc","mean_lfc","lfc"),colnames(DC))[1]; dn<-intersect(c("drug","Drug","compound"),colnames(DC))[1]
setnames(DC,c(dn,fc),c("drug","lfc")); DC[,drug:=clean(drug)]
sens<-DC[order(lfc)][1:10]; resi<-DC[order(-lfc)][1:10]; DD<-rbind(sens,resi)[order(lfc)]
DD[,drug:=factor(drug,levels=drug)]
pA<-ggplot(DD,aes(lfc,drug,fill=lfc<0))+geom_vline(xintercept=0,linewidth=0.3,colour="grey70")+geom_col(width=0.72)+
  scale_fill_manual(values=c("TRUE"=COL_HIGH,"FALSE"=COL_LOW),labels=c("Resistant in high","Sensitive in high"),name=NULL)+
  labs(x="log2FC drug response (high/low NeFRS)",y=NULL,title="Drug sensitivity (3-cohort consistent)")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=5),
    legend.position="bottom",legend.key.size=unit(2.4,"mm"),legend.text=element_text(size=5))
save_panel(pA,file.path(PFIG,"MF8_A_drug"),w_mm=72,h_mm=72)

## B mutation multivariate Cox forest
MC<-fread(file.path(MUT,"mutation_multicox.tsv")); cat("multicox:",paste(colnames(MC),collapse=","),"| var:",paste(MC$var,collapse=","),"\n")
MC[,lab:=ifelse(!is.na(label)&label!="",label,var)]; MC[,lab:=factor(lab,levels=rev(lab))]
MC[,hl:=ifelse(grepl("NeFRS",var,ignore.case=TRUE),"NeFRS","other")]
pB<-ggplot(MC,aes(HR,lab,colour=hl))+geom_vline(xintercept=1,linetype=2,linewidth=0.4,colour="grey55")+
  geom_errorbarh(aes(xmin=low,xmax=high),height=0.25,linewidth=0.5)+geom_point(size=1.8)+
  scale_colour_manual(values=c(NeFRS=COL_HIGH,other="grey35"),guide="none")+scale_x_log10()+
  labs(x="Hazard ratio (95% CI, log scale)",y=NULL,title="Multivariate Cox (vs driver mutations)")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=5.6))
save_panel(pB,file.path(PFIG,"MF8_B_forest"),w_mm=64,h_mm=50)

## C model C-index comparison
MD<-fread(file.path(MUT,"model_comparison.tsv")); cat("model_comparison:",paste(MD$Model,collapse=" | "),"\n")
MD[,Model:=factor(Model,levels=Model[order(C_index)])]; MD[,hl:=grepl("NeFRS",Model)]
pC<-ggplot(MD,aes(C_index,Model,fill=hl))+geom_col(width=0.66)+
  geom_text(aes(label=sprintf("%.3f",C_index)),hjust=-0.15,size=5/.pt)+
  scale_fill_manual(values=c("TRUE"=COL_HIGH,"FALSE"="grey60"),guide="none")+
  coord_cartesian(xlim=c(0.5,max(MD$C_index)*1.08))+
  labs(x="C-index",y=NULL,title="Model comparison")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=5.6))
save_panel(pC,file.path(PFIG,"MF8_C_model"),w_mm=64,h_mm=50)

## D NeFRS rank percentile across cohorts
NR<-fread(file.path(SIG,"nefrs_rank_summary.tsv")); cat("rank队列:",nrow(NR),"\n")
NR[,cohort:=factor(cohort,levels=cohort[order(percentile)])]
pD<-ggplot(NR,aes(percentile,cohort))+geom_col(fill=COL_HIGH,width=0.66)+
  geom_text(aes(label=sprintf("%d/%d",rank,total)),hjust=-0.1,size=4.6/.pt)+
  coord_cartesian(xlim=c(0,108))+
  labs(x="Percentile vs published signatures (%)",y=NULL,title="NeFRS rank (per cohort)")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"),axis.text.y=element_text(size=5.6))
save_panel(pD,file.path(PFIG,"MF8_D_rank"),w_mm=72,h_mm=50)

## E signature C-index distribution (PACA_AU), NeFRS highlighted
SC<-fread(file.path(SIG,"sig_compare_all_cohorts.tsv"))[cohort=="PACA_AU"]
nefrs_c<-NR[cohort=="PACA_AU"]$NeFRS_cindex
cat("PACA signature数:",nrow(SC)," NeFRS C:",nefrs_c,"\n")
pE<-ggplot(SC,aes(cindex))+geom_histogram(bins=24,fill="grey70",colour="white",linewidth=0.2)+
  geom_vline(xintercept=nefrs_c,colour=COL_HIGH,linewidth=0.8)+
  annotate("text",x=nefrs_c,y=Inf,label=sprintf(" NeFRS %.3f",nefrs_c),colour=COL_HIGH,hjust=0,vjust=1.4,size=5.4/.pt)+
  labs(x="C-index (PACA-AU)",y="Published signatures (n)",title="NeFRS vs 84 published signatures")+
  theme_pub(base=6.5)+theme(plot.title=element_text(size=6.5,face="bold"))
save_panel(pE,file.path(PFIG,"MF8_E_dist"),w_mm=72,h_mm=50)
cat("MF8 A-E 完成\n")
