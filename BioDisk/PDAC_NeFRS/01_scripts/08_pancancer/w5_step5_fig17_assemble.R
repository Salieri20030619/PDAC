#!/usr/bin/env Rscript
# W5 Step 5: Fig17 拼图 (2×2): KM(HCC/CRC 阴性) + 随机零分布(HCC/CRC 特异性)
suppressPackageStartupMessages({ library(data.table); library(survival) })
ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
OUT  <- file.path(ROOT, "02_results/08_pancancer")
FIG  <- file.path(ROOT, "03_figures/Fig17_pancancer")

rs   <- fread(file.path(OUT, "pancancer_nefrs_rs.tsv"))
obs  <- list(HCC=0.656, CRC=0.496); pval <- list(HCC=0.057, CRC=0.726)  # 来自 Step4 log

km_panel <- function(tag, title) {
  d <- rs[cohort==tag]; med <- median(d$RS)
  d[, Group := factor(ifelse(RS>med,"High","Low"), levels=c("Low","High"))]
  fit <- survfit(Surv(OS_time,OS_status)~Group, d)
  cox <- coxph(Surv(OS_time,OS_status)~Group, d); hr <- exp(coef(cox)); ci <- exp(confint(cox))
  p   <- 1-pchisq(survdiff(Surv(OS_time,OS_status)~Group,d)$chisq,1)
  cidx<- as.numeric(concordance(Surv(OS_time,OS_status)~RS,d,reverse=TRUE)$concordance)
  cols<- c("#0072B2","#C00000"); mt <- max(d$OS_time); tb <- pretty(c(0,mt),5)
  plot(fit,col=cols,lwd=2.5,conf.int=FALSE,xlab="Time (months)",
       ylab="Overall Survival",xlim=c(0,mt),ylim=c(0,1),xaxt="n",yaxt="n",mark.time=TRUE,cex=0.5)
  axis(1,at=tb); axis(2,at=seq(0,1,.2),las=2); title(main=title,cex.main=1.1,font.main=2)
  legend("topright",c(sprintf("Low (n=%d)",sum(d$Group=="Low")),
         sprintf("High (n=%d)",sum(d$Group=="High"))),col=cols,lwd=2.5,bty="n",cex=0.85)
  legend("bottomleft",c(sprintf("Log-rank P = %.3f",p),
         sprintf("HR = %.2f (%.2f-%.2f)",hr,ci[1],ci[2]),
         sprintf("C-index = %.3f",cidx)),bty="n",cex=0.85,y.intersp=1.15)
}

null_panel <- function(tag, title) {
  nl <- fread(file.path(OUT, sprintf("refit_null_%s.tsv",tag)))$null_cindex
  oc <- obs[[tag]]; q <- quantile(nl,c(.025,.975))
  hist(nl,breaks=40,col="grey85",border="white",main=title,cex.main=1.1,font.main=2,
       xlab="Refit 5-fold CV C-index",xlim=range(c(nl,oc))+c(-.02,.02))
  abline(v=oc,col="#C00000",lwd=3); abline(v=q,col="#0072B2",lwd=1.5,lty=2)
  legend("topright",c(sprintf("NeFRS = %.3f",oc),sprintf("empirical p = %.3f",pval[[tag]]),
         "random 95%"),col=c("#C00000",NA,"#0072B2"),lwd=c(3,NA,1.5),lty=c(1,NA,2),bty="n",cex=0.8)
}

pdf(file.path(FIG,"Fig17_pancancer_combined.pdf"), width=11, height=9, pointsize=11)
par(mfrow=c(2,2), mar=c(4,4.5,3,1.5), mgp=c(2.3,0.7,0), oma=c(0,0,2.5,0))
km_panel("HCC","A  HCC (TCGA-LIHC): fixed-weight NeFRS")
km_panel("CRC","B  CRC (TCGA-COAD+READ): fixed-weight NeFRS")
null_panel("HCC","C  HCC: NeFRS-69 vs random panels (refit)")
null_panel("CRC","D  CRC: NeFRS-69 vs random panels (refit)")
mtext("Fig 17  NeFRS is PDAC-specific: no transfer to HCC/CRC",
      side=3,line=0.3,outer=TRUE,cex=1.3,font=2)
dev.off()
cat("✓ 03_figures/Fig17_pancancer/Fig17_pancancer_combined.pdf\n")
