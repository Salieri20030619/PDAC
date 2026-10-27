#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(survival); library(rms) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; OUT<-file.path(ROOT,"02_results/05_nomogram"); PFIG<-file.path(ROOT,"03_figures/panels")
clin<-readRDS(file.path(OUT,"clinical_with_nefrs_3cohorts.rds"))
run_cal<-function(df,f,u,m){ct<-cph(f,data=df,x=TRUE,y=TRUE,surv=TRUE,time.inc=u)
  tryCatch({capture.output({r<-calibrate(ct,cmethod="KM",method="boot",u=u,m=m,B=200)});r},
    error=function(e){cat(sprintf("  cal t=%d fail: %s\n",u,conditionMessage(e)));NULL})}
dl<-function(cal,col){if(is.null(cal)||!"mean.predicted"%in%colnames(cal))return(FALSE)
  x<-cal[,"mean.predicted"];y<-cal[,"KM"];ok<-!is.na(x)&!is.na(y);if(sum(ok)<2)return(FALSE)
  lines(x[ok],y[ok],col=col,lwd=2,type="b",pch=16,cex=1.0);TRUE}
cc3<-c(COL_LOW,"#009E73",COL_HIGH)
basec<-function(){plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="Predicted OS probability",ylab="Observed OS (KM)",las=1)
  abline(0,1,lty=2,col="gray60");abline(h=seq(0,1,.2),v=seq(0,1,.2),col="gray92",lty=3,lwd=0.5)}

paca<-as.data.frame(clin$PACA_AU)
paca_cc<-paca[complete.cases(paca[,c("NeFRS","age","Sex","Stage_AJCC","Grade","os_time_month","os_status")]),]
dd<-datadist(paca_cc);options(datadist="dd")
fP<-Surv(os_time_month,os_status)~NeFRS+age+Sex+Stage_AJCC+Grade
p12<-run_cal(paca_cc,fP,12,65);p36<-run_cal(paca_cc,fP,36,65);p60<-run_cal(paca_cc,fP,60,65)
labP<-sprintf("PACA-AU (n=%d)",nrow(paca_cc))
save_base(function(){par(mar=c(3.4,3.6,1.5,1),mgp=c(2.1,0.6,0),tcl=-0.25);basec();L<-c();C<-c()
  if(dl(p12,cc3[1])){L<-c(L,"1-year");C<-c(C,cc3[1])};if(dl(p36,cc3[2])){L<-c(L,"3-year");C<-c(C,cc3[2])}
  if(dl(p60,cc3[3])){L<-c(L,"5-year");C<-c(C,cc3[3])}
  mtext(labP,side=3,line=0.4,font=2,cex=1.0);legend("topleft",L,col=C,lwd=2,pch=16,bty="n",cex=0.85)},
  file.path(PFIG,"MF3_calib_PACA"),w_mm=74,h_mm=74)

tcga<-as.data.frame(clin$TCGA_PAAD)
tcga_cc<-tcga[complete.cases(tcga[,c("NeFRS","age","Sex","Stage_AJCC","os_time_month","os_status")]),
  c("sample_id","os_time_month","os_status","NeFRS","age","Sex","Stage_AJCC")]
tcga_cc$Stage_AJCC<-droplevels(tcga_cc$Stage_AJCC);dd_t<-datadist(tcga_cc);options(datadist="dd_t")
fT<-Surv(os_time_month,os_status)~NeFRS+age+Sex+Stage_AJCC
t12<-run_cal(tcga_cc,fT,12,36);t36<-run_cal(tcga_cc,fT,36,36)
labT<-sprintf("TCGA-PAAD (n=%d)",nrow(tcga_cc))
save_base(function(){par(mar=c(3.4,3.6,1.5,1),mgp=c(2.1,0.6,0),tcl=-0.25);basec();L<-c();C<-c()
  if(dl(t12,cc3[1])){L<-c(L,"1-year");C<-c(C,cc3[1])};if(dl(t36,cc3[2])){L<-c(L,"3-year");C<-c(C,cc3[2])}
  mtext(labT,side=3,line=0.4,font=2,cex=1.0);legend("topleft",L,col=C,lwd=2,pch=16,bty="n",cex=0.85)},
  file.path(PFIG,"MF3_calib_TCGA"),w_mm=74,h_mm=74)

emtab<-as.data.frame(clin$E_MTAB_6134)
emtab_cc<-emtab[complete.cases(emtab[,c("NeFRS","Sex","Stage_AJCC","Grade","os_time_month","os_status")]),
  c("sample_id","os_time_month","os_status","NeFRS","Sex","Stage_AJCC","Grade")]
emtab_cc$Stage_AJCC<-droplevels(emtab_cc$Stage_AJCC);emtab_cc$Grade<-droplevels(emtab_cc$Grade)
dd_e<-datadist(emtab_cc);options(datadist="dd_e")
fE<-Surv(os_time_month,os_status)~NeFRS+Sex+Stage_AJCC+Grade
e12<-run_cal(emtab_cc,fE,12,70);e36<-run_cal(emtab_cc,fE,36,70)
labE<-sprintf("E-MTAB-6134 (n=%d)",nrow(emtab_cc))
save_base(function(){par(mar=c(3.4,3.6,1.5,1),mgp=c(2.1,0.6,0),tcl=-0.25);basec();L<-c();C<-c()
  if(dl(e12,cc3[1])){L<-c(L,"1-year");C<-c(C,cc3[1])};if(dl(e36,cc3[2])){L<-c(L,"3-year");C<-c(C,cc3[2])}
  mtext(labE,side=3,line=0.4,font=2,cex=1.0);legend("topleft",L,col=C,lwd=2,pch=16,bty="n",cex=0.85)},
  file.path(PFIG,"MF3_calib_EMTAB"),w_mm=74,h_mm=74)
cat("MF3 校准 完成 (PACA/TCGA/E-MTAB)\n")
