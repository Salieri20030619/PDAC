#!/usr/bin/env Rscript
# 增量-2: NeFRS 超越 AJCC 分期 (C-index/LRT + IDI/NRI), 3 全分期队列
suppressPackageStartupMessages({ library(data.table); library(survival) })
ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
NOM  <- file.path(ROOT,"02_results/05_nomogram/clinical_with_nefrs_3cohorts.rds")
OUT  <- file.path(ROOT,"02_results/21_nri"); dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
cat("=== 增量-2: NeFRS beyond AJCC stage ===\n")

cl <- readRDS(NOM)
stage_num <- function(x){ x<-as.character(x)
  ifelse(grepl("III|IV",x),4, ifelse(grepl("IIB",x),3, ifelse(grepl("IIA",x),2, ifelse(grepl("II",x),2, ifelse(grepl("I",x),1,NA))))) }
pool <- rbindlist(lapply(names(cl), function(cn){ d<-as.data.table(cl[[cn]])
  data.table(cohort=cn, t=d$os_time_month, s=d$os_status, RS=as.numeric(d$NeFRS), stg=stage_num(d$Stage_AJCC)) }),fill=TRUE)
pool <- pool[!is.na(t)&t>0&!is.na(s)&!is.na(RS)&!is.na(stg)]
pool[, z := as.numeric(scale(RS)), by=cohort]
cat(sprintf("n=%d, 事件=%d, 队列=%d; 分期(1=I,2=II,3=IIB,4=III/IV):\n", nrow(pool),sum(pool$s),uniqueN(pool$cohort))); print(table(pool$stg))

# (A) C-index + 似然比 (无外部依赖)
cb <- coxph(Surv(t,s)~stg+strata(cohort), data=pool)
cf <- coxph(Surv(t,s)~stg+z+strata(cohort), data=pool)
Cb <- concordance(cb)$concordance; Cf <- concordance(cf)$concordance
lr <- 2*(as.numeric(logLik(cf))-as.numeric(logLik(cb))); pl <- pchisq(lr,1,lower.tail=FALSE)
cat(sprintf("\n[C-index] 分期模型 C=%.3f  ->  分期+NeFRS C=%.3f  (ΔC=+%.3f)\n", Cb,Cf,Cf-Cb))
cat(sprintf("[似然比] 加入 NeFRS: LR=%.1f, df=1, p=%.2e\n", lr,pl))
cat(sprintf("[NeFRS|stage] HR per SD=%.2f, p=%.2e\n", exp(coef(cf)["z"]), summary(cf)$coefficients["z","Pr(>|z|)"]))
fwrite(data.frame(C_stage=Cb,C_stage_NeFRS=Cf,dC=Cf-Cb,LR=lr,LR_p=pl), file.path(OUT,"beyond_stage_cindex.tsv"), sep="\t")

# (B) IDI/NRI (survIDINRI)
if (requireNamespace("survIDINRI", quietly=TRUE)) {
  suppressPackageStartupMessages(library(survIDINRI))
  D <- as.matrix(pool[, .(t,s)])
  cohD <- model.matrix(~cohort, data=as.data.frame(pool))[,-1,drop=FALSE]
  covs0 <- as.matrix(cbind(stg=pool$stg, cohD))
  covs1 <- as.matrix(cbind(stg=pool$stg, cohD, z=pool$z))
  for (t0 in c(12,24)) {
    x <- IDI.INF(D, covs0, covs1, t0, npert=500, seed1=1)
    cat(sprintf("\n--- t0=%d 月 (M1=IDI, M2=连续NRI, M3=中位风险改善) ---\n", t0)); print(round(IDI.INF.OUT(x),4))
  }
} else cat("\n[提示] 未装 survIDINRI,已完成 C-index/LRT 部分。装后重跑获取 IDI/NRI。\n")
cat("\nDone.\n")
