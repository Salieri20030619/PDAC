#!/usr/bin/env Rscript
# 补充B v5: 亚组分析 (单行append + 加固fit_sub + 逐亚组打印HR)
suppressPackageStartupMessages({ library(data.table); library(survival) })
ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
NOM  <- file.path(ROOT,"02_results/05_nomogram/clinical_with_nefrs_3cohorts.rds")
OUT  <- file.path(ROOT,"02_results/19_subgroup"); FIG <- file.path(ROOT,"03_figures/Fig19_subgroup")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE); dir.create(FIG,recursive=TRUE,showWarnings=FALSE)
cat("=== 补充B: 亚组分析 ===\n")

cl <- readRDS(NOM)
pool <- rbindlist(lapply(names(cl), function(cn){ d <- as.data.table(cl[[cn]])
  data.table(cohort=cn, sample_id=d$sample_id, OS_time=d$os_time_month, OS_status=d$os_status,
    RS=as.numeric(d$NeFRS), age=suppressWarnings(as.numeric(d$age)),
    sex=d$Sex, stage_raw=d$Stage_AJCC, grade_raw=d$Grade) }), fill=TRUE)
pool <- pool[!is.na(RS) & !is.na(OS_time) & OS_time>0 & !is.na(OS_status)]
pool[, RS_z := as.numeric(scale(RS)), by=cohort]
norm_sex  <- function(x){x<-tolower(trimws(as.character(x))); ifelse(grepl("^f|female",x),"Female",ifelse(grepl("^m|male",x),"Male",NA))}
stage_grp <- function(x){x<-as.character(x); ifelse(grepl("III|IV",x),"Late (III-IV)",ifelse(grepl("I",x),"Early (I-II)",NA))}
grade_grp <- function(x){x<-as.character(x); ifelse(grepl("3|4|[Pp]oor|[Uu]ndiff",x),"High (G3-4)",ifelse(grepl("1|2|[Ww]ell|[Mm]oder",x),"Low (G1-2)",NA))}
pool[, sex:=norm_sex(sex)][, stage:=stage_grp(stage_raw)][, grade:=grade_grp(grade_raw)]
pool[, age_grp:=ifelse(is.na(age),NA,ifelse(age>65,">65","<=65"))]
cat(sprintf("池化: %d 样本, %d 事件\n", nrow(pool), sum(pool$OS_status)))

fit_sub <- function(d){ d <- d[!is.na(RS_z)]
  if (nrow(d)<20 || sum(d$OS_status)<10) return(list(ok=FALSE, why="n/ev不足"))
  tryCatch({
    f <- coxph(Surv(OS_time,OS_status)~RS_z+strata(cohort), data=d)
    s <- summary(f)$coefficients; ci <- suppressWarnings(confint(f))
    cir <- if (is.matrix(ci)) ci["RS_z",] else ci
    list(ok=TRUE, n=nrow(d), ev=sum(d$OS_status), hr=as.numeric(s["RS_z","exp(coef)"]),
         lo=as.numeric(exp(cir[1])), hi=as.numeric(exp(cir[2])), p=as.numeric(s["RS_z","Pr(>|z|)"]))
  }, error=function(e) list(ok=FALSE, why=conditionMessage(e))) }

interaction_p <- function(var){ d<-pool[!is.na(get(var)) & !is.na(RS_z)]
  if(uniqueN(d[[var]])<2 || sum(d$OS_status)<10) return(NA_real_)
  f0<-tryCatch(coxph(as.formula(sprintf("Surv(OS_time,OS_status)~RS_z+%s+strata(cohort)",var)),data=d),error=function(e)NULL)
  f1<-tryCatch(coxph(as.formula(sprintf("Surv(OS_time,OS_status)~RS_z*%s+strata(cohort)",var)),data=d),error=function(e)NULL)
  if(is.null(f0)||is.null(f1)) return(NA_real_); tryCatch(as.numeric(anova(f0,f1)[2,"P(>|Chi|)"]),error=function(e)NA_real_) }
pint_map <- sapply(c("stage","age_grp","sex","grade"), interaction_p)

sgdef <- list(c("All patients","",""),
  c("Stage Early (I-II)","stage","Early (I-II)"), c("Stage Late (III-IV)","stage","Late (III-IV)"),
  c("Age <=65","age_grp","<=65"), c("Age >65","age_grp",">65"),
  c("Sex Female","sex","Female"), c("Sex Male","sex","Male"),
  c("Grade Low (G1-2)","grade","Low (G1-2)"), c("Grade High (G3-4)","grade","High (G3-4)"))
res_rows <- list()
for (ii in seq_along(sgdef)) {
  sg<-sgdef[[ii]]; lab<-sg[1]; v<-sg[2]; lv<-sg[3]
  if (v=="") d <- pool else { keep <- pool[[v]]==lv; keep[is.na(keep)] <- FALSE; d <- pool[keep] }
  r <- fit_sub(d)
  if (!isTRUE(r$ok)) { cat(sprintf("  [%d/9] %-20s n=%d ev=%d 跳过: %s\n", ii,lab,nrow(d),sum(d$OS_status),r$why)); next }
  pint <- if (v=="") NA_real_ else as.numeric(pint_map[[v]])
  cat(sprintf("  [%d/9] %-20s n=%d ev=%d HR=%.2f (%.2f-%.2f) p=%.1e\n", ii,lab,r$n,r$ev,r$hr,r$lo,r$hi,r$p))
  res_rows[[length(res_rows)+1]] <- data.table(subgroup=lab,n=r$n,events=r$ev,HR=round(r$hr,2),HR_lo=round(r$lo,2),HR_hi=round(r$hi,2),p=signif(r$p,3))
}
cat(sprintf("res_rows 累计 %d 行\n", length(res_rows)))
summ <- rbindlist(res_rows, fill=TRUE)
fwrite(summ, file.path(OUT,"subgroup_summary.tsv"), sep="\t")
cat("\n=== 亚组汇总 ===\n"); print(summ)

pdf(file.path(FIG,"Fig19A_subgroup_forest.pdf"),width=9,height=6); k<-nrow(summ); par(mar=c(4,11,3,7))
plot(NA,xlim=c(0.8,max(summ$HR_hi,na.rm=TRUE)*1.1),ylim=c(0.5,k+0.5),log="x",yaxt="n",xlab="NeFRS HR per SD (95% CI)",ylab="")
abline(v=1,lty=2,col="grey50")
for(i in 1:k){ y<-k-i+1; segments(summ$HR_lo[i],y,summ$HR_hi[i],y,lwd=2)
  points(summ$HR[i],y,pch=18,cex=1.5,col=ifelse(summ$p[i]<0.05,"#C00000","grey50"))
  axis(2,at=y,labels=summ$subgroup[i],las=2,tick=FALSE,cex.axis=0.8)
  t<-sprintf("%.2f (%.2f-%.2f)",summ$HR[i],summ$HR_lo[i],summ$HR_hi[i])
  text(max(summ$HR_hi,na.rm=TRUE)*1.12,y,t,adj=0,cex=0.7,xpd=NA) }
title("NeFRS prognostic across clinical subgroups (3 cohorts)"); dev.off()

km_strat <- function(stg,ttl){ d<-pool[stage==stg & !is.na(RS_z)]; if(nrow(d)<20) return(invisible())
  d[, grp:=factor(ifelse(RS_z>median(RS_z),"High","Low"),levels=c("Low","High"))]
  fit<-survfit(Surv(OS_time,OS_status)~grp,data=d); hr<-exp(coef(coxph(Surv(OS_time,OS_status)~grp+strata(cohort),data=d)))
  pv<-1-pchisq(survdiff(Surv(OS_time,OS_status)~grp,data=d)$chisq,1); cols<-c("#0072B2","#C00000")
  plot(fit,col=cols,lwd=2.5,conf.int=FALSE,xlab="Time (months)",ylab="OS",ylim=c(0,1),mark.time=TRUE,cex=0.5)
  title(sprintf("%s (n=%d)",ttl,nrow(d)),cex.main=1.05)
  legend("topright",c(sprintf("Low (n=%d)",sum(d$grp=="Low")),sprintf("High (n=%d)",sum(d$grp=="High"))),col=cols,lwd=2.5,bty="n",cex=0.85)
  legend("bottomleft",c(sprintf("Log-rank P=%.1e",pv),sprintf("HR=%.2f",hr)),bty="n",cex=0.85) }
pdf(file.path(FIG,"Fig19B_stage_stratified_KM.pdf"),width=11,height=5); par(mfrow=c(1,2),mar=c(4,4.5,3,1))
km_strat("Early (I-II)","Early stage I-II"); km_strat("Late (III-IV)","Late stage III-IV"); dev.off()
cat("\n✓ Fig19A + Fig19B\nDone.\n")
cat("\n=== 交互检验 p (RS_z x 亚组) ===\n"); print(pint_map)
