#!/usr/bin/env Rscript
# 增量-1: 8队列 NeFRS per-SD HR 随机效应 meta (DerSimonian-Laird, 无外部依赖)
suppressPackageStartupMessages({ library(survival); library(glmnet) })
ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
W2IN <- file.path(ROOT,"02_results/02_w2_inputs"); MIME <- file.path(ROOT,"02_results/04_mime")
OUT  <- file.path(ROOT,"02_results/20_meta"); FIG <- file.path(ROOT,"03_figures/Fig20_meta")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE); dir.create(FIG,recursive=TRUE,showWarnings=FALSE)
cat("=== 增量-1: 8队列 meta (per-SD HR) ===\n")

cohort_data <- readRDS(file.path(W2IN,"cohort_data_8.rds"))
mime_input  <- readRDS(file.path(MIME,"mime_input_list.rds"))
res <- readRDS(file.path(MIME,"mime_results.rds"))
mime_genes <- colnames(mime_input[[1]])[-(1:3)]; ml <- res$ml.res[["Enet[α=0.1]"]]
compute_rs <- function(em){ X <- matrix(0,ncol(em),length(mime_genes),dimnames=list(colnames(em),mime_genes))
  for(g in intersect(mime_genes,rownames(em))) X[,g] <- em[g,]
  setNames(as.numeric(predict(ml,newx=X,s="lambda.min")),colnames(em)) }

rows <- data.frame()
for(cn in names(cohort_data)){
  co <- cohort_data[[cn]]; cl <- co$clinical; rs <- compute_rs(co$expression)
  d <- data.frame(t=cl$os_time_month, s=cl$os_status, z=as.numeric(scale(rs[cl$sample_id])))
  d <- d[!is.na(d$t) & d$t>0 & !is.na(d$s) & !is.na(d$z),]
  f <- coxph(Surv(t,s)~z, data=d); b <- coef(f)["z"]; se <- sqrt(vcov(f)["z","z"])
  rows <- rbind(rows, data.frame(cohort=cn,n=nrow(d),ev=sum(d$s),logHR=as.numeric(b),se=se,
    HR=exp(b),lo=exp(b-1.96*se),hi=exp(b+1.96*se)))
}
rownames(rows) <- NULL

dl <- function(yi,vi){ k<-length(yi); wf<-1/vi; muf<-sum(wf*yi)/sum(wf)
  Q<-sum(wf*(yi-muf)^2); C<-sum(wf)-sum(wf^2)/sum(wf); tau2<-max(0,(Q-(k-1))/C)
  wr<-1/(vi+tau2); mu<-sum(wr*yi)/sum(wr); se<-sqrt(1/sum(wr))
  list(HR=exp(mu),lo=exp(mu-1.96*se),hi=exp(mu+1.96*se),p=2*pnorm(-abs(mu/se)),
       I2=max(0,(Q-(k-1))/Q)*100, tau2=tau2, Q=Q, Qp=pchisq(Q,k-1,lower.tail=FALSE)) }

P  <- dl(rows$logHR, rows$se^2)                       # 全部 8 队列
val<- rows$cohort!="PACA_AU"; V <- dl(rows$logHR[val], rows$se[val]^2)  # 7 验证队列

cat("\n各队列:\n"); print(rows[,c("cohort","n","ev","HR","lo","hi")], row.names=FALSE, digits=3)
cat(sprintf("\n[全 8 队列] 随机效应合并 HR = %.2f (%.2f-%.2f), p=%.2e | I2=%.1f%%, Q p=%.3f\n",
  P$HR,P$lo,P$hi,P$p,P$I2,P$Qp))
cat(sprintf("[7 验证队列, 排除训练 PACA_AU] 合并 HR = %.2f (%.2f-%.2f), p=%.2e | I2=%.1f%%\n",
  V$HR,V$lo,V$hi,V$p,V$I2))
write.table(rows, file.path(OUT,"meta_per_cohort.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(set=c("all8","val7"),
  HR=c(P$HR,V$HR),lo=c(P$lo,V$lo),hi=c(P$hi,V$hi),p=c(P$p,V$p),I2=c(P$I2,V$I2)),
  file.path(OUT,"meta_pooled.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

pdf(file.path(FIG,"Fig20_meta_forest.pdf"), width=8, height=5.5)
k1<-nrow(rows); par(mar=c(4,11,2,7)); xr<-range(c(rows$lo,rows$hi,P$lo,P$hi))
plot(NA,xlim=c(min(0.8,xr[1]*0.9),xr[2]*1.1),ylim=c(0,k1+2),log="x",yaxt="n",xlab="NeFRS HR per SD (95% CI)",ylab="")
abline(v=1,lty=2,col="grey50")
for(i in 1:k1){ y<-k1-i+2; segments(rows$lo[i],y,rows$hi[i],y,lwd=2)
  points(rows$HR[i],y,pch=15,cex=0.9+rows$n[i]/max(rows$n)*1.1,col="#0072B2")
  axis(2,at=y,labels=sprintf("%s (n=%d)",rows$cohort[i],rows$n[i]),las=2,tick=FALSE,cex.axis=0.78)
  text(xr[2]*1.13,y,sprintf("%.2f (%.2f-%.2f)",rows$HR[i],rows$lo[i],rows$hi[i]),adj=0,cex=0.7,xpd=NA) }
polygon(c(P$lo,P$HR,P$hi,P$HR),c(1,1.3,1,0.7),col="#C00000",border="#C00000")
axis(2,at=1,labels="Random-effects (all 8)",las=2,tick=FALSE,cex.axis=0.82,font=2)
text(xr[2]*1.13,1,sprintf("%.2f (%.2f-%.2f)",P$HR,P$lo,P$hi),adj=0,cex=0.76,font=2,xpd=NA)
title(sprintf("8-cohort meta-analysis (I2=%.0f%%)",P$I2)); dev.off()
cat("\n✓ Fig20_meta_forest.pdf\nDone.\n")
