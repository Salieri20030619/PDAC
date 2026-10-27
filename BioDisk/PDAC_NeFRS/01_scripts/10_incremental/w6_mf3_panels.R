#!/usr/bin/env Rscript
# MF3 全面版 8panel (WYSIWYG): A列线图 B多变量森林 C超越分期 D IDI/NRI E校准 F DCA G/H分期内KM
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(survival); library(rms); library(data.table); library(ggplot2) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; NOMD<-file.path(ROOT,"02_results/05_nomogram")
NRID<-file.path(ROOT,"02_results/21_nri"); PFIG<-file.path(ROOT,"03_figures/panels")
clin<-readRDS(file.path(NOMD,"clinical_with_nefrs_3cohorts.rds"))
fp<-function(p){if(is.na(p))return("NA");if(p<0.001)return("P<0.001");sprintf("P=%.3f",p)}
ok<-function(nm,expr) tryCatch({expr; cat("[OK]",nm,"\n")}, error=function(e) cat("[FAIL]",nm,":",conditionMessage(e),"\n"))

paca<-as.data.frame(clin$PACA_AU)
paca_cc<-paca[complete.cases(paca[,c("NeFRS","age","Sex","Stage_AJCC","Grade","os_time_month","os_status")]),]
dd<-datadist(paca_cc); options(datadist="dd")

## A 列线图
ok("A_nomogram",{
  cf<-cph(Surv(os_time_month,os_status)~NeFRS+age+Sex+Stage_AJCC+Grade,data=paca_cc,x=TRUE,y=TRUE,surv=TRUE,time.inc=12)
  set.seed(5201314); val<-validate(cf,method="boot",B=200,dxy=TRUE); Cc<-(val["Dxy","index.corrected"]+1)/2
  so<-Survival(cf)
  nom<-nomogram(cf,fun=list(function(x)so(12,lp=x),function(x)so(36,lp=x),function(x)so(60,lp=x)),
    funlabel=c("1-Year OS","3-Year OS","5-Year OS"),fun.at=c(0.1,0.3,0.5,0.7,0.85,0.95),lp=FALSE)
  save_base(function(){par(mar=c(1.6,1.2,1.2,1.2));plot(nom,xfrac=0.32,cex.axis=0.55,cex.var=0.62,
    col.grid=gray(c(0.85,0.96)),lmgp=0.13,label.every=1)
    mtext(sprintf("Bootstrap-corrected C-index = %.3f",Cc),3,line=0.0,cex=0.6,col="gray30")},
    file.path(PFIG,"MF3_nomogram"),w_mm=176,h_mm=68)
})

## B 多变量 Cox 森林
ok("B_forest",{
  cs<-summary(coxph(Surv(os_time_month,os_status)~NeFRS+age+Sex+Stage_AJCC+Grade,data=paca_cc))
  sh<-c("NeFRS"="NeFRS (per unit)","age"="Age (per year)","SexMale"="Sex (M vs F)",
    "Stage_AJCCIIA"="Stage IIA (vs I)","Stage_AJCCIIB"="Stage IIB (vs I)","Stage_AJCCIII-IV"="Stage III-IV (vs I)",
    "GradeG2"="Grade G2 (vs G1)","GradeG3-G4"="Grade G3-G4 (vs G1)")
  ht<-data.frame(short=unname(sh[rownames(cs$coefficients)]),HR=cs$coefficients[,"exp(coef)"],
    lo=cs$conf.int[,"lower .95"],hi=cs$conf.int[,"upper .95"],P=cs$coefficients[,"Pr(>|z|)"],stringsAsFactors=FALSE)
  draw_f<-function(){n<-nrow(ht);ys<-n:1
    par(mar=c(2.2,0.3,1.0,0.3),mgp=c(1.3,0.4,0),tcl=-0.2,lwd=0.7)
    plot(NA,xlim=c(0,1),ylim=c(0.3,n+0.9),axes=FALSE,xlab="",ylab="")
    text(0.003,n+0.65,"Variable",adj=0,font=2,cex=0.6)
    text(0.55,n+0.65,"Hazard ratio (95% CI)",adj=0.5,font=2,cex=0.6)
    text(0.995,n+0.65,"P",adj=1,font=2,cex=0.6)
    xl<-0.27;xr<-0.66;hmin<-0.3;hmax<-12
    fx<-function(hr) xl+(log(hr)-log(hmin))/(log(hmax)-log(hmin))*(xr-xl)
    segments(fx(1),0.5,fx(1),n+0.35,lty=2,col="gray55",lwd=0.6)
    for(i in 1:n){y<-ys[i];r<-ht[i,];col<-ifelse(r$P<0.05,COL_HIGH,"#666666")
      text(0.003,y,r$short,adj=0,cex=0.55,col=col,font=ifelse(r$P<0.05,2,1))
      segments(fx(max(r$lo,hmin)),y,fx(min(r$hi,hmax)),y,lwd=1.1,col=col)
      points(fx(r$HR),y,pch=15,cex=0.9,col=col)
      text(0.70,y,sprintf("%.2f (%.2f-%.2f)",r$HR,r$lo,r$hi),adj=0,cex=0.53)
      text(0.995,y,fp(r$P),adj=1,cex=0.53,col=col,font=ifelse(r$P<0.05,2,1))}
    xb<-c(0.5,1,2,5,10);axis(1,at=fx(xb),labels=xb,cex.axis=0.55,lwd=0.6,padj=-0.7)
    mtext("Hazard ratio (log scale)",1,line=1.2,cex=0.58)}
  save_base(draw_f,file.path(PFIG,"MF3_forest"),w_mm=124,h_mm=44)
})

## C 超越分期 ΔC
ok("C_beyond",{
  bs<-fread(file.path(NRID,"beyond_stage_cindex.tsv"))
  dh<-data.frame(model=factor(c("AJCC\nstage","Stage +\nNeFRS"),levels=c("AJCC\nstage","Stage +\nNeFRS")),
    C=c(bs$C_stage,bs$C_stage_NeFRS))
  p<-ggplot(dh,aes(model,C,fill=model))+geom_col(width=0.62)+
    geom_hline(yintercept=0.5,linetype=2,linewidth=0.4,colour="grey50")+
    geom_text(aes(label=sprintf("%.3f",C)),vjust=-0.4,size=6.5/.pt)+
    annotate("text",x=1.5,y=max(dh$C)+0.045,label=sprintf("\u0394C = +%.3f",bs$dC),fontface=2,size=6.5/.pt)+
    scale_fill_manual(values=c("AJCC\nstage"="grey70","Stage +\nNeFRS"=COL_HIGH),guide="none")+
    coord_cartesian(ylim=c(0.45,max(dh$C)+0.09))+labs(x=NULL,y="C-index")+theme_pub(base=6.5)
  save_panel(p,file.path(PFIG,"MF3_beyond"),w_mm=48,h_mm=44)
})

## D IDI / NRI (沿用 w5_nri survIDINRI 值)
ok("D_idinri",{
  di<-data.frame(measure=factor(rep(c("IDI","Cont. NRI","Median impr."),2),levels=c("IDI","Cont. NRI","Median impr.")),
    t0=factor(rep(c("12 mo","24 mo"),each=3),levels=c("12 mo","24 mo")),
    est=c(0.152,0.376,0.146,0.132,0.280,0.160),lo=c(0.107,0.291,0.106,0.094,0.206,0.104),hi=c(0.205,0.473,0.215,0.173,0.366,0.233))
  p<-ggplot(di,aes(measure,est,colour=t0))+geom_hline(yintercept=0,linewidth=0.4,colour="grey60")+
    geom_point(position=position_dodge(0.5),size=1.5)+
    geom_errorbar(aes(ymin=lo,ymax=hi),position=position_dodge(0.5),width=0.2,linewidth=0.5)+
    scale_colour_manual(values=c("12 mo"=COL_LOW,"24 mo"=COL_HIGH))+scale_y_continuous(limits=c(0,NA))+
    theme_pub(base=6.5)+labs(x=NULL,y="Improvement vs stage",colour=NULL)+
    theme(axis.text.x=element_text(angle=20,hjust=1),legend.position=c(0.8,0.9),
      legend.key.size=unit(2.6,"mm"),legend.text=element_text(size=5.6))
  save_panel(p,file.path(PFIG,"MF3_idinri"),w_mm=58,h_mm=48)
})

## E 校准 PACA
ok("E_calib",{
  cc3<-c(COL_LOW,"#009E73",COL_HIGH)
  rc<-function(u,m){ct<-cph(Surv(os_time_month,os_status)~NeFRS+age+Sex+Stage_AJCC+Grade,data=paca_cc,x=TRUE,y=TRUE,surv=TRUE,time.inc=u)
    tryCatch({capture.output({r<-calibrate(ct,cmethod="KM",method="boot",u=u,m=m,B=200)});r},error=function(e)NULL)}
  p12<-rc(12,65);p36<-rc(36,65);p60<-rc(60,65)
  dl<-function(cal,col){if(is.null(cal)||!"mean.predicted"%in%colnames(cal))return(FALSE)
    x<-cal[,"mean.predicted"];y<-cal[,"KM"];v<-!is.na(x)&!is.na(y);if(sum(v)<2)return(FALSE)
    lines(x[v],y[v],col=col,lwd=1.6,type="b",pch=16,cex=0.7);TRUE}
  save_base(function(){par(mar=c(2.6,2.8,1.2,0.6),mgp=c(1.5,0.4,0),tcl=-0.2)
    plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="Predicted OS prob.",ylab="Observed OS (KM)",las=1,cex.axis=0.6,cex.lab=0.62)
    abline(0,1,lty=2,col="gray60",lwd=0.7)
    L<-c();C<-c();if(dl(p12,cc3[1])){L<-c(L,"1-y");C<-c(C,cc3[1])};if(dl(p36,cc3[2])){L<-c(L,"3-y");C<-c(C,cc3[2])}
    if(dl(p60,cc3[3])){L<-c(L,"5-y");C<-c(C,cc3[3])}
    legend("topleft",L,col=C,lwd=1.6,pch=16,bty="n",cex=0.6)},
    file.path(PFIG,"MF3_calib"),w_mm=52,h_mm=48)
})

## F DCA TCGA
ok("F_dca",{
  suppressPackageStartupMessages(library(dcurves))
  tcga<-as.data.frame(clin$TCGA_PAAD)
  tc<-tcga[complete.cases(tcga[,c("NeFRS","age","Sex","Stage_AJCC","os_time_month","os_status")]),]
  sp<-function(f){sf<-survfit(coxph(f,data=tc),newdata=tc);i<-which.min(abs(sf$time-36));1-sf$surv[i,]}
  tc$r_full<-sp(Surv(os_time_month,os_status)~NeFRS+age+Sex+Stage_AJCC)
  tc$r_nefrs<-sp(Surv(os_time_month,os_status)~NeFRS)
  tc$r_clin<-sp(Surv(os_time_month,os_status)~age+Sex+Stage_AJCC)
  d<-dca(Surv(os_time_month,os_status)~r_full+r_nefrs+r_clin,data=tc,time=36,
    label=list(r_full="NeFRS+Clinical",r_nefrs="NeFRS only",r_clin="Clinical only"))
  p<-plot(d,smooth=TRUE)+theme_pub(base=6.5)+
    scale_colour_manual(values=c("NeFRS+Clinical"=COL_HIGH,"NeFRS only"=COL_LOW,"Clinical only"="#009E73",
      "Treat all"="gray60","Treat none"="gray30","all"="gray60","none"="gray30"))+
    theme(legend.position=c(0.70,0.80),legend.title=element_blank(),legend.key.size=unit(2.6,"mm"),
      legend.text=element_text(size=5.4))+labs(title=NULL,x="Threshold probability",y="Net benefit")
  save_panel(p,file.path(PFIG,"MF3_dca"),w_mm=62,h_mm=48)
})

## G/H 分期内分层 KM (pool 3 cohorts)
ok("GH_stageKM",{
  pool<-rbindlist(lapply(names(clin),function(cn){d<-as.data.table(clin[[cn]])
    data.table(cohort=cn,t=d$os_time_month,s=d$os_status,RS=as.numeric(d$NeFRS),stg=as.character(d$Stage_AJCC))}),fill=TRUE)
  pool<-pool[!is.na(RS)&!is.na(t)&t>0&!is.na(s)]
  pool[,RS_z:=as.numeric(scale(RS)),by=cohort]
  pool[,stage:=ifelse(grepl("III|IV",stg),"Late",ifelse(grepl("I",stg),"Early",NA))]
  draw_km<-function(tgt,ttl){d<-pool[stage==tgt&!is.na(RS_z)]
    d[,grp:=factor(ifelse(RS_z>median(RS_z),"High","Low"),levels=c("Low","High"))]
    fit<-survfit(Surv(t,s)~grp,data=d)
    hr<-exp(coef(coxph(Surv(t,s)~grp+strata(cohort),data=d)))
    pv<-1-pchisq(survdiff(Surv(t,s)~grp,data=d)$chisq,1)
    par(mar=c(2.4,2.9,1.4,0.6),mgp=c(1.4,0.4,0),tcl=-0.2,lwd=0.7)
    plot(fit,col=c(COL_LOW,COL_HIGH),lwd=1.5,conf.int=FALSE,xlab="",ylab="",ylim=c(0,1),
      xaxt="n",yaxt="n",mark.time=TRUE,cex=0.4)
    axis(1,cex.axis=0.6,lwd=0.6,padj=-0.6);axis(2,seq(0,1,.5),las=2,cex.axis=0.6,lwd=0.6)
    mtext(sprintf("%s (n=%d)",ttl,nrow(d)),3,line=0.3,font=2,cex=0.62)
    mtext("Months",1,line=1.2,cex=0.55);mtext("Survival",2,line=1.9,cex=0.55,las=0)
    legend("topright",c(sprintf("HR %.2f",hr),fp(pv)),bty="n",cex=0.56,y.intersp=0.95)
    legend("bottomleft",c("Low","High"),col=c(COL_LOW,COL_HIGH),lwd=1.5,bty="n",cex=0.56,seg.len=0.8,y.intersp=0.9)}
  save_base(function() draw_km("Early","Early stage (I-II)"),file.path(PFIG,"MF3_km_early"),w_mm=86,h_mm=40)
  save_base(function() draw_km("Late","Late stage (III-IV)"),file.path(PFIG,"MF3_km_late"),w_mm=86,h_mm=40)
})
cat("MF3 panels 完成\n")
