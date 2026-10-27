#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(survival); library(data.table); library(ggplot2) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
CTL<-file.path(ROOT,"02_results/18_controls"); SUB<-file.path(ROOT,"02_results/19_subgroup")
PAN<-file.path(ROOT,"02_results/08_pancancer"); NOMD<-file.path(ROOT,"02_results/05_nomogram")
PFIG<-file.path(ROOT,"03_figures/panels")
fp<-function(p){if(is.na(p))return("NA");if(p<0.001)return("P<0.001");sprintf("P=%.3f",p)}

## ---- A 控制森林: 校正纯度/增殖/basal 后各队列 adj_HR ----
A<-fread(file.path(CTL,"controls_summary.tsv"))
A<-A[order(-adj_HR)]
draw_A<-function(){n<-nrow(A);ys<-n:1
  par(mar=c(2.2,0.3,1.4,0.3),mgp=c(1.3,0.4,0),tcl=-0.2,lwd=0.7)
  plot(NA,xlim=c(0,1),ylim=c(0.3,n+0.9),axes=FALSE,xlab="",ylab="")
  text(0.003,n+0.6,"Cohort",adj=0,font=2,cex=0.6)
  text(0.50,n+0.6,"Adjusted HR (95% CI)",adj=0.5,font=2,cex=0.6)
  text(0.995,n+0.6,"P",adj=1,font=2,cex=0.6)
  xl<-0.30;xr<-0.68;hmin<-0.8;hmax<-5
  fx<-function(h) xl+(log(h)-log(hmin))/(log(hmax)-log(hmin))*(xr-xl)
  segments(fx(1),0.5,fx(1),n+0.3,lty=2,col="gray55",lwd=0.6)
  for(i in 1:n){y<-ys[i];r<-A[i];col<-ifelse(r$adj_p<0.05,COL_HIGH,"#666666")
    text(0.003,y,r$cohort,adj=0,cex=0.54,col=col,font=ifelse(r$adj_p<0.05,2,1))
    segments(fx(max(r$HR_lo,hmin)),y,fx(min(r$HR_hi,hmax)),y,lwd=1.1,col=col)
    points(fx(r$adj_HR),y,pch=15,cex=0.85,col=col)
    text(0.72,y,sprintf("%.2f (%.2f-%.2f)",r$adj_HR,r$HR_lo,r$HR_hi),adj=0,cex=0.5)
    text(0.995,y,fp(r$adj_p),adj=1,cex=0.5,col=col,font=ifelse(r$adj_p<0.05,2,1))}
  xb<-c(1,2,3,5);axis(1,at=fx(xb),labels=xb,cex.axis=0.55,lwd=0.6,padj=-0.7)
  mtext("Adj. HR (purity+prolif+basal corrected)",1,line=1.2,cex=0.55)}
save_base(draw_A,file.path(PFIG,"MF4_A_controls"),w_mm=88,h_mm=52)

## ---- B 亚组森林 + 交互P ----
B<-fread(file.path(SUB,"subgroup_summary.tsv"))
clin<-readRDS(file.path(NOMD,"clinical_with_nefrs_3cohorts.rds"))
pool<-rbindlist(lapply(names(clin),function(cn){d<-as.data.table(clin[[cn]])
  data.table(cohort=cn,t=as.numeric(d$os_time_month),s=as.numeric(d$os_status),
    RS=as.numeric(d$NeFRS),age=as.numeric(d$age),sex=as.character(d$Sex),
    stg=as.character(d$Stage_AJCC),grade=as.character(d$Grade))}),fill=TRUE)
pool<-pool[!is.na(RS)&!is.na(t)&t>0&!is.na(s)]
pool[,RSz:=as.numeric(scale(RS)),by=cohort]
pool[,stage2:=ifelse(grepl("III|IV",stg),"III-IV",ifelse(grepl("I",stg),"I-II",NA))]
pool[,age2:=ifelse(age>=median(age,na.rm=TRUE),">=med","<med")]
pool[,grade2:=ifelse(grepl("3|4",grade),"G3-4",ifelse(grepl("1|2",grade),"G1-2",NA))]
intP<-function(v){d<-pool[!is.na(get(v))]
  m1<-coxph(Surv(t,s)~RSz*factor(get(v))+strata(cohort),data=d)
  m0<-coxph(Surv(t,s)~RSz+factor(get(v))+strata(cohort),data=d)
  pchisq(2*(m1$loglik[2]-m0$loglik[2]),df=length(coef(m1))-length(coef(m0)),lower.tail=FALSE)}
ipv<-c(Stage=intP("stage2"),Age=intP("age2"),Sex=intP("sex"),Grade=intP("grade2"))
cat("交互P:\n");print(round(ipv,3))
draw_B<-function(){n<-nrow(B);ys<-n:1
  par(mar=c(2.4,0.3,1.4,0.3),mgp=c(1.3,0.4,0),tcl=-0.2,lwd=0.7)
  plot(NA,xlim=c(0,1),ylim=c(0.3,n+0.9),axes=FALSE,xlab="",ylab="")
  text(0.003,n+0.6,"Subgroup",adj=0,font=2,cex=0.6)
  text(0.50,n+0.6,"HR (95% CI)",adj=0.5,font=2,cex=0.6)
  text(0.995,n+0.6,"P",adj=1,font=2,cex=0.6)
  xl<-0.34;xr<-0.70;hmin<-0.7;hmax<-4
  fx<-function(h) xl+(log(h)-log(hmin))/(log(hmax)-log(hmin))*(xr-xl)
  segments(fx(1),0.5,fx(1),n+0.3,lty=2,col="gray55",lwd=0.6)
  for(i in 1:n){y<-ys[i];r<-B[i];col<-ifelse(r$p<0.05,COL_HIGH,"#666666")
    text(0.003,y,r$subgroup,adj=0,cex=0.5,col=col,font=ifelse(r$p<0.05,2,1))
    segments(fx(max(r$HR_lo,hmin)),y,fx(min(r$HR_hi,hmax)),y,lwd=1.1,col=col)
    points(fx(r$HR),y,pch=15,cex=0.8,col=col)
    text(0.72,y,sprintf("%.2f (%.2f-%.2f)",r$HR,r$HR_lo,r$HR_hi),adj=0,cex=0.48)
    text(0.995,y,fp(r$p),adj=1,cex=0.48,col=col,font=ifelse(r$p<0.05,2,1))}
  xb<-c(1,2,3);axis(1,at=fx(xb),labels=xb,cex.axis=0.55,lwd=0.6,padj=-0.7)
  mtext(sprintf("Interaction P: Stage %.2f, Age %.2f, Sex %.2f, Grade %.2f",
    ipv["Stage"],ipv["Age"],ipv["Sex"],ipv["Grade"]),1,line=1.15,cex=0.46)}
save_base(draw_B,file.path(PFIG,"MF4_B_subgroup"),w_mm=88,h_mm=52)

## ---- C/D/E 散点: NeFRS vs purity/prolif/basal ----
S<-fread(file.path(CTL,"controls_scores_persample.tsv"))
mk<-function(yv,ylab,stub){rho<-cor(S$NeFRS,S[[yv]],method="spearman",use="complete.obs")
  p<-ggplot(S,aes(NeFRS,.data[[yv]]))+
    geom_point(size=0.25,alpha=0.28,colour="grey45")+
    geom_smooth(method="lm",se=FALSE,colour=COL_HIGH,linewidth=0.6)+
    annotate("text",x=-Inf,y=Inf,hjust=-0.08,vjust=1.4,
      label=sprintf("Spearman \u03c1 = %.2f",rho),size=6.5/.pt)+
    labs(x="NeFRS",y=ylab)+theme_pub(base=6.5)
  save_panel(p,file.path(PFIG,stub),w_mm=58,h_mm=46)}
mk("purity","Tumor purity","MF4_C_purity")
mk("prolif","Proliferation score","MF4_D_prolif")
mk("basal_prob","Basal-like probability","MF4_E_basal")

## ---- F 泛癌外推 C-index: PDAC vs HCC vs CRC ----
PC<-fread(file.path(PAN,"pancancer_summary.tsv"))
dF<-data.frame(can=factor(c("PDAC","HCC","CRC"),levels=c("PDAC","HCC","CRC")),
  C=c(0.757,PC[cohort=="HCC"]$c_index,PC[cohort=="CRC"]$c_index))
pF<-ggplot(dF,aes(can,C,fill=can))+geom_col(width=0.62)+
  geom_hline(yintercept=0.5,linetype=2,linewidth=0.4,colour="grey50")+
  geom_text(aes(label=sprintf("%.3f",C)),vjust=-0.4,size=6.5/.pt)+
  scale_fill_manual(values=c("PDAC"=COL_HIGH,"HCC"="grey70","CRC"="grey70"),guide="none")+
  coord_cartesian(ylim=c(0.4,0.82))+labs(x=NULL,y="C-index (locked weights)")+theme_pub(base=6.5)
save_panel(pF,file.path(PFIG,"MF4_F_pancancer"),w_mm=56,h_mm=46)

## ---- G refit 零分布: 锁权重 vs 重训 null (HCC/CRC) ----
nH<-fread(file.path(PAN,"refit_null_HCC.tsv"))$null_cindex
nC<-fread(file.path(PAN,"refit_null_CRC.tsv"))$null_cindex
trueH<-PC[cohort=="HCC"]$c_index; trueC<-PC[cohort=="CRC"]$c_index
dG<-rbind(data.frame(can="HCC",null=nH),data.frame(can="CRC",null=nC))
vl<-data.frame(can=c("HCC","CRC"),x=c(trueH,trueC))
pG<-ggplot(dG,aes(null))+geom_histogram(bins=28,fill="grey75",colour="white",linewidth=0.1)+
  geom_vline(data=vl,aes(xintercept=x),colour=COL_HIGH,linewidth=0.7)+
  facet_wrap(~can,ncol=1,scales="free_y")+
  labs(x="Refit C-index (null)",y="Count")+theme_pub(base=6.5)+
  theme(strip.text=element_text(face="bold",size=6))
save_panel(pG,file.path(PFIG,"MF4_G_refitnull"),w_mm=60,h_mm=46)
cat("MF4 A-G 完成\n")
