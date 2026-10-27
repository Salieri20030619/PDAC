#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(survival); library(data.table); library(ggplot2) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; PFIG<-file.path(ROOT,"03_figures/panels")
cd8<-readRDS(file.path(ROOT,"02_results/02_w2_inputs/cohort_data_8.rds"))
P<-cd8[["PACA_AU"]]
expr<-P$expression                     # gene x sample (20741 x 267)
clin<-as.data.frame(P$clinical)        # sample x 3
stopifnot(ncol(expr)==nrow(clin))      # 假设列序==行序
cat("clinical cols:",paste(colnames(clin),collapse=","),"\n"); print(head(clin,2))
tc<-grep("time",colnames(clin),ignore.case=TRUE,value=TRUE)[1]
sc<-setdiff(grep("status|event|^OS$|death|vital",colnames(clin),ignore.case=TRUE,value=TRUE),tc)[1]
cat("time col:",tc," | status col:",sc,"\n")
sv<-Surv(as.numeric(clin[[tc]]),as.numeric(clin[[sc]]))

G<-t(expr)                             # sample x gene
keep<-apply(G,2,function(x) sd(x,na.rm=TRUE)>0); G<-G[,keep]; Gz<-scale(G)
cat("usable genes:",ncol(Gz),"\n")

cf<-fread(file.path(ROOT,"02_results/04_mime/nefrs_v3_enet_coefficients.tsv"))
cf<-cf[gene %in% colnames(Gz)]; k<-nrow(cf)
cat("NeFRS genes matched:",k,"\n")
sc_true<-as.numeric(Gz[,cf$gene] %*% cf$coef)
c_true<-summary(coxph(sv~sc_true))$concordance[1]; if(c_true<0.5) c_true<-1-c_true
cat("true NeFRS C-index:",round(c_true,3)," (应≈0.757)\n")

cat("预计算全基因单变量 Cox 符号 (约 1-2 分钟)...\n")
signs<-apply(Gz,2,function(x) sign(coef(coxph(sv~x))))
set.seed(20260614); np<-1000; nc<-numeric(np)
for(i in 1:np){idx<-sample(ncol(Gz),k); s<-as.numeric(Gz[,idx]%*%signs[idx])
  nc[i]<-summary(coxph(sv~s))$concordance[1]}
nc<-ifelse(nc<0.5,1-nc,nc)
p_emp<-(sum(nc>=c_true)+1)/(np+1)
cat(sprintf("random panel null: mean=%.3f max=%.3f trueNeFRS=%.3f P_emp=%.3g\n",mean(nc),max(nc),c_true,p_emp))

pH<-ggplot(data.frame(nc=nc),aes(nc))+
  geom_histogram(bins=30,fill="grey75",colour="white",linewidth=0.1)+
  geom_vline(xintercept=c_true,colour=COL_HIGH,linewidth=0.8)+
  annotate("text",x=c_true,y=Inf,vjust=1.6,hjust=1.08,
    label=sprintf("NeFRS %.3f\nP %s",c_true,ifelse(p_emp<0.001,"< 0.001",sprintf("= %.3f",p_emp))),
    size=6/.pt,colour=COL_HIGH,fontface=2,lineheight=0.9)+
  labs(x=sprintf("C-index of random %d-gene panels",k),y="Count")+theme_pub(base=6.5)
save_panel(pH,file.path(PFIG,"MF4_H_randomnull"),w_mm=60,h_mm=46)
cat("MF4 H 完成\n")
