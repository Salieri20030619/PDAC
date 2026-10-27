#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(survival) })
ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
cl <- readRDS(file.path(ROOT,"02_results/05_nomogram/clinical_with_nefrs_3cohorts.rds"))
pool <- rbindlist(lapply(names(cl), function(cn){ d<-as.data.table(cl[[cn]])
  data.table(cohort=cn, OS_time=d$os_time_month, OS_status=d$os_status, RS=as.numeric(d$NeFRS),
    age=suppressWarnings(as.numeric(d$age)), sex=d$Sex, stage_raw=d$Stage_AJCC, grade_raw=d$Grade) }),fill=TRUE)
pool <- pool[!is.na(RS)&!is.na(OS_time)&OS_time>0&!is.na(OS_status)]
pool[, RS_z:=as.numeric(scale(RS)), by=cohort]
norm_sex<-function(x){x<-tolower(trimws(as.character(x)));ifelse(grepl("^f|female",x),"Female",ifelse(grepl("^m|male",x),"Male",NA))}
stage_grp<-function(x){x<-as.character(x);ifelse(grepl("III|IV",x),"Late",ifelse(grepl("I",x),"Early",NA))}
grade_grp<-function(x){x<-as.character(x);ifelse(grepl("3|4|[Pp]oor",x),"High",ifelse(grepl("1|2|[Ww]ell|[Mm]oder",x),"Low",NA))}
pool[,sex:=norm_sex(sex)][,stage:=stage_grp(stage_raw)][,grade:=grade_grp(grade_raw)]
pool[,age_grp:=ifelse(is.na(age),NA,ifelse(age>65,">65","<=65"))]
ip <- function(var){ d<-pool[!is.na(get(var))&!is.na(RS_z)]
  if(uniqueN(d[[var]])<2 || sum(d$OS_status)<10) return(NA_real_)
  f0<-coxph(as.formula(sprintf("Surv(OS_time,OS_status)~RS_z+%s+strata(cohort)",var)),data=d)
  f1<-coxph(as.formula(sprintf("Surv(OS_time,OS_status)~RS_z*%s+strata(cohort)",var)),data=d)
  lr<-2*(as.numeric(logLik(f1))-as.numeric(logLik(f0))); dfd<-length(coef(f1))-length(coef(f0))
  as.numeric(pchisq(lr,dfd,lower.tail=FALSE)) }
cat("=== 交互检验 (似然比, RS_z x 亚组) ===\n")
for(v in c("stage","age_grp","sex","grade")) cat(sprintf("  %-10s  p_interaction = %.4f\n", v, ip(v)))
