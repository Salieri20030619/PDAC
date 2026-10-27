#!/usr/bin/env Rscript
# 补充: 三大对照 — NeFRS 是否独立于 纯度/增殖/分子亚型  (v2: 修TCGA样本名)
suppressPackageStartupMessages({
  library(data.table); library(survival); library(glmnet); library(estimate)
})
ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
W2IN <- file.path(ROOT, "02_results/02_w2_inputs"); MIME <- file.path(ROOT, "02_results/04_mime")
OUT  <- file.path(ROOT, "02_results/18_controls"); FIG <- file.path(ROOT, "03_figures/Fig18_controls")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE); dir.create(FIG, recursive=TRUE, showWarnings=FALSE)
cat("=== 补充对照: 纯度 / 增殖 / 亚型 ===\n")

cohort_data <- readRDS(file.path(W2IN, "cohort_data_8.rds"))
mime_input  <- readRDS(file.path(MIME, "mime_input_list.rds"))
res         <- readRDS(file.path(MIME, "mime_results.rds"))
mime_genes  <- colnames(mime_input[[1]])[-(1:3)]; ml <- res$ml.res[["Enet[α=0.1]"]]
compute_rs <- function(em) {
  X <- matrix(0, ncol(em), length(mime_genes), dimnames=list(colnames(em), mime_genes))
  for (g in intersect(mime_genes, rownames(em))) X[, g] <- em[g, ]
  setNames(as.numeric(predict(ml, newx=X, s="lambda.min")), colnames(em))
}

# 名字稳健匹配 (修复 ESTIMATE/PurIST 把 TCGA 连字符改成点的问题)
getv <- function(v, ids) {
  if (is.null(v) || is.null(names(v))) return(rep(NA_real_, length(ids)))
  o <- v[ids]; m <- is.na(o)
  if (any(m)) { nm <- setNames(as.numeric(v), make.names(names(v))); o[m] <- nm[make.names(ids[m])] }
  as.numeric(o)
}
safecor <- function(a,b) tryCatch(suppressWarnings(cor(a,b,method="spearman",use="complete.obs")),
                                  error=function(e) NA_real_)

cc <- unique(c(Seurat::cc.genes.updated.2019$s.genes, Seurat::cc.genes.updated.2019$g2m.genes))

load(file.path(ROOT, "tools/PurIST/fitteds_public_2019-02-12.Rdata"))
source(file.path(ROOT, "tools/PurIST/functions.R"))
purist_clf <- classifs[[1]]

run_estimate <- function(em) {
  ti <- tempfile(fileext=".tsv"); tg <- tempfile(fileext=".gct"); to <- tempfile(fileext=".gct")
  write.table(data.frame(GeneSymbol=rownames(em), em, check.names=FALSE),
              ti, sep="\t", row.names=FALSE, quote=FALSE)
  filterCommonGenes(input.f=ti, output.f=tg, id="GeneSymbol")
  estimateScore(input.ds=tg, output.ds=to, platform="affymetrix")
  es <- read.delim(to, skip=2, row.names=1, check.names=FALSE)[,-1, drop=FALSE]
  setNames(as.numeric(es["ESTIMATEScore", ]), colnames(es))
}

summ <- data.frame(); scores_all <- list()
for (cn in names(cohort_data)) {
  co <- cohort_data[[cn]]; em <- co$expression; clin <- co$clinical
  cat(sprintf("\n[%s] expr %d genes x %d samples\n", cn, nrow(em), ncol(em)))
  if (nrow(em) < 5000) { cat("  warn: 基因过少, 跳过\n"); next }

  rs   <- compute_rs(em)
  est  <- tryCatch(run_estimate(em), error=function(e){cat("  ESTIMATE err:",conditionMessage(e),"\n"); NULL})
  pro  <- { sub <- em[intersect(cc, rownames(em)), , drop=FALSE]; colMeans(t(scale(t(sub))), na.rm=TRUE) }
  pur  <- tryCatch({ p <- apply_classifier(data=as.data.frame(em), classifier=purist_clf)
                     setNames(p$Pred_prob_basal, rownames(p)) },
                   error=function(e){cat("  PurIST err:",conditionMessage(e),"\n"); NULL})

  ids <- clin$sample_id
  d <- data.frame(sample_id=ids, OS_time=clin$os_time_month, OS_status=clin$os_status,
                  NeFRS=as.numeric(scale(getv(rs, ids))),
                  purity=as.numeric(scale(getv(est, ids))),
                  prolif=as.numeric(scale(getv(pro, ids))),
                  basal_prob=getv(pur, ids))
  d$subtype <- ifelse(is.na(d$basal_prob), NA, ifelse(d$basal_prob>0.5,"basal","classical"))
  d <- d[!is.na(d$OS_time) & d$OS_time>0 & !is.na(d$OS_status), ]
  scores_all[[cn]] <- data.frame(cohort=cn, d)

  c_pur <- safecor(d$NeFRS, d$purity); c_pro <- safecor(d$NeFRS, d$prolif); c_bas <- safecor(d$NeFRS, d$basal_prob)

  vars <- "NeFRS"
  if (!all(is.na(d$purity))) vars <- c(vars, "purity")
  vars <- c(vars, "prolif")
  if (length(unique(na.omit(d$subtype)))==2) vars <- c(vars, "subtype")
  f <- tryCatch(coxph(as.formula(paste("Surv(OS_time,OS_status) ~", paste(vars, collapse="+"))), data=d),
                error=function(e) NULL)
  if (!is.null(f)) { sm <- summary(f)$coefficients
    hr <- exp(sm["NeFRS",1]); ci <- exp(confint(f)["NeFRS",]); pv <- sm["NeFRS",5]
  } else { hr<-ci<-pv<-NA }
  cat(sprintf("  cor purity=%.2f prolif=%.2f basalProb=%.2f | 调整后 NeFRS HR=%.2f (%.2f-%.2f) p=%.2g\n",
              c_pur,c_pro,c_bas,hr,ci[1],ci[2],pv))
  summ <- rbind(summ, data.frame(cohort=cn, n=nrow(d),
    cor_purity=round(c_pur,2), cor_prolif=round(c_pro,2), cor_basalprob=round(c_bas,2),
    pct_basal=round(mean(d$subtype=="basal",na.rm=TRUE),2),
    adj_HR=round(hr,2), HR_lo=round(ci[1],2), HR_hi=round(ci[2],2), adj_p=signif(pv,3)))
}
fwrite(summ, file.path(OUT, "controls_summary.tsv"), sep="\t")
fwrite(rbindlist(scores_all), file.path(OUT, "controls_scores_persample.tsv"), sep="\t")
cat("\n=== 汇总 (adj_HR = 调整纯度+增殖+亚型后 NeFRS HR) ===\n"); print(summ)
cat(sprintf("\n调整三者后仍显著(p<0.05): %d/%d 队列\n", sum(summ$adj_p<0.05, na.rm=TRUE), nrow(summ)))

pdf(file.path(FIG, "Fig18_NeFRS_adjusted_forest.pdf"), width=8, height=5.5)
s <- summ[!is.na(summ$adj_HR), ]; k <- nrow(s); par(mar=c(4,8,3,2))
plot(NA, xlim=range(c(0.8, s$HR_hi, 4)), ylim=c(0.5,k+0.5), log="x", yaxt="n",
     xlab="NeFRS HR (adjusted for purity + proliferation + subtype)", ylab="")
abline(v=1, lty=2, col="grey50")
for (i in 1:k) { y<-k-i+1
  segments(s$HR_lo[i], y, s$HR_hi[i], y, lwd=2)
  points(s$adj_HR[i], y, pch=18, cex=1.6, col=ifelse(s$adj_p[i]<0.05,"#C00000","grey50"))
  axis(2, at=y, labels=s$cohort[i], las=2, tick=FALSE)
  text(max(s$HR_hi)*1.05, y, sprintf("%.2f (%.2f-%.2f)", s$adj_HR[i],s$HR_lo[i],s$HR_hi[i]), adj=0, cex=0.8, xpd=NA) }
title("NeFRS prognostic after adjusting purity, proliferation & subtype")
dev.off()
cat("OK Fig18_NeFRS_adjusted_forest.pdf\nDone.\n")
