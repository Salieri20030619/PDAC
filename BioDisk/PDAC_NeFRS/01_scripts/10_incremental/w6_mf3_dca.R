#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(survival); library(dcurves); library(ggplot2) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; OUT<-file.path(ROOT,"02_results/05_nomogram"); PFIG<-file.path(ROOT,"03_figures/panels")
clin<-readRDS(file.path(OUT,"clinical_with_nefrs_3cohorts.rds"))
tcga<-as.data.frame(clin$TCGA_PAAD)
tcga_cc<-tcga[complete.cases(tcga[,c("NeFRS","age","Sex","Stage_AJCC","os_time_month","os_status")]),]
sp<-function(f){sf<-survfit(coxph(f,data=tcga_cc),newdata=tcga_cc);i<-which.min(abs(sf$time-36));1-sf$surv[i,]}
tcga_cc$risk_3y<-sp(Surv(os_time_month,os_status)~NeFRS+age+Sex+Stage_AJCC)
tcga_cc$risk_nefrs<-sp(Surv(os_time_month,os_status)~NeFRS)
tcga_cc$risk_clin<-sp(Surv(os_time_month,os_status)~age+Sex+Stage_AJCC)
d<-dca(Surv(os_time_month,os_status)~risk_3y+risk_nefrs+risk_clin,data=tcga_cc,time=36,
  label=list(risk_3y="NeFRS + Clinical",risk_nefrs="NeFRS only",risk_clin="Clinical only"))
p<-plot(d,smooth=TRUE)+
  scale_color_manual(values=c("NeFRS + Clinical"=COL_HIGH,"NeFRS only"=COL_LOW,"Clinical only"="#009E73",
    "Treat All"="gray60","Treat None"="gray30","all"="gray60","none"="gray30"))+
  theme_pub()+theme(legend.position="right",legend.title=element_blank())+labs(title=NULL)
save_panel(p, file.path(PFIG,"MF3_DCA_TCGA"), w_mm=120, h_mm=78)
cat("MF3 DCA 完成\n")
