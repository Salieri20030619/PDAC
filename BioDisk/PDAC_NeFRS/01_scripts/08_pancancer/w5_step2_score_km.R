#!/usr/bin/env Rscript
# ============================================================
# W5 Step 2: NeFRS 泛癌验证 (HCC + CRC)
# 照搬 w3_step1 口径: predict(ml, lambda.min) over 170 mime_genes (缺失填0,
#   不z-score) → median 分组 (Low=ref) → KM/Cox/concordance(reverse=TRUE)
# 输入: 8队列模型 + pancancer/{HCC,CRC}/{HCC,CRC}.gene.tsv.gz + .clinical_std.tsv
# 输出: Fig17 泛癌 KM + pancancer_summary.tsv + pancancer_nefrs_rs.tsv
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(survival); library(glmnet) })

ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
MIME <- file.path(ROOT, "02_results/04_mime")
PAN  <- file.path(ROOT, "00_raw_data/pancancer")
OUT  <- file.path(ROOT, "02_results/08_pancancer")
FIG  <- file.path(ROOT, "03_figures/Fig17_pancancer")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

cat("=== W5 Step 2: NeFRS pancancer validation (HCC + CRC) ===\n")

# ---- 1) 加载 NeFRS 模型 (与 w3_step1 完全一致) ----
mime_input <- readRDS(file.path(MIME, "mime_input_list.rds"))
res        <- readRDS(file.path(MIME, "mime_results.rds"))
mime_genes <- colnames(mime_input[[1]])[-(1:3)]            # 170 funnel 基因
ml         <- res$ml.res[["Enet[α=0.1]"]]
nz <- fread(file.path(MIME, "nefrs_v3_enet_coefficients.tsv"))$gene   # 69 非零基因
cat(sprintf("  模型: Enet[α=0.1], %d funnel 基因, %d 非零(NeFRS)\n", length(mime_genes), length(nz)))

# ---- 2) 评分函数 (复刻 compute_rs: 缺失填0, 不z-score) ----
compute_rs <- function(expr_mat) {   # expr_mat: gene(行) × sample(列), log2
  X <- matrix(0, nrow=ncol(expr_mat), ncol=length(mime_genes),
              dimnames=list(colnames(expr_mat), mime_genes))
  avail <- intersect(mime_genes, rownames(expr_mat))
  for (g in avail) X[, g] <- expr_mat[g, ]
  rs <- as.numeric(predict(ml, newx=X, s="lambda.min")); names(rs) <- rownames(X)
  list(rs=rs, cov170=length(avail), cov69=length(intersect(nz, rownames(expr_mat))))
}

load_cohort <- function(tag) {
  gm  <- fread(file.path(PAN, tag, paste0(tag, ".gene.tsv.gz")))
  cl  <- fread(file.path(PAN, tag, paste0(tag, ".clinical_std.tsv")))
  em  <- as.matrix(gm[, -1]); rownames(em) <- gm$gene
  list(expression=em, clinical=cl)
}

format_pval <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-15) return("P < 1e-15")
  if (p < 1e-4)  return(sprintf("P = %.1e", p))
  if (p < 0.001) return(sprintf("P = %.4f", p))
  sprintf("P = %.3f", p)
}

compute_stats <- function(co, label) {
  sc <- compute_rs(co$expression)
  cl <- co$clinical
  m  <- match(cl$sample_id, names(sc$rs))
  df <- data.frame(sample_id=cl$sample_id, OS_time=cl$os_time_month,
                   OS_status=cl$os_status, RS=sc$rs[m],
                   age=cl$age, sex=cl$sex, stage=cl$stage)
  df <- df[!is.na(df$RS) & !is.na(df$OS_time) & df$OS_time > 0 & !is.na(df$OS_status), ]
  med <- median(df$RS, na.rm=TRUE)
  df$Group <- factor(ifelse(df$RS > med, "High", "Low"), levels=c("Low","High"))

  fit <- survfit(Surv(OS_time, OS_status) ~ Group, data=df)
  cox <- coxph(Surv(OS_time, OS_status) ~ Group, data=df)
  hr  <- exp(coef(cox)); ci <- exp(confint(cox))
  sd_obj <- survdiff(Surv(OS_time, OS_status) ~ Group, data=df)
  pv  <- 1 - pchisq(sd_obj$chisq, df=1)
  cidx <- as.numeric(concordance(Surv(OS_time, OS_status) ~ RS, data=df, reverse=TRUE)$concordance)

  # 多变量 Cox (调整可用临床, RS per SD) — 看独立性
  mv_hr <- mv_p <- NA
  df$RSz <- as.numeric(scale(df$RS))
  mv <- tryCatch({
    d <- df; d$stage <- ifelse(grepl("III|IV", d$stage), "late",
                        ifelse(grepl("I", d$stage), "early", NA))
    f <- coxph(Surv(OS_time, OS_status) ~ RSz + age + sex + stage, data=d)
    s <- summary(f)$coefficients["RSz", ]; c(hr=exp(s[1]), p=s[5])
  }, error=function(e) c(hr=NA, p=NA))
  mv_hr <- mv[1]; mv_p <- mv[2]

  cat(sprintf("  [%s] n=%d events=%d | 170cov=%d 69cov=%d | C=%.3f HR=%.2f %s | mvHR(perSD)=%.2f p=%.2g\n",
              label, nrow(df), sum(df$OS_status), sc$cov170, sc$cov69,
              cidx, hr, format_pval(pv), mv_hr, mv_p))
  list(fit=fit, df=df, hr=hr, lo=ci[1], hi=ci[2], pv=pv, cidx=cidx,
       cov170=sc$cov170, cov69=sc$cov69, mv_hr=mv_hr, mv_p=mv_p, label=label,
       n=nrow(df), nh=sum(df$Group=="High"), nl=sum(df$Group=="Low"))
}

# ---- 3) KM 绘图 (复刻 w3_step1 样式: 蓝低红高 + risk table) ----
plot_km <- function(s) {
  cols <- c("#0072B2", "#C00000")
  layout(matrix(c(1,2), ncol=1), heights=c(3.5, 1.2))
  par(mar=c(3.5,5,3.5,1.5), mgp=c(2.2,0.7,0))
  mt <- max(s$df$OS_time, na.rm=TRUE); tb <- pretty(c(0, mt), n=5)
  plot(s$fit, col=cols, lwd=3, conf.int=FALSE, xlab="Time (months)",
       ylab="Overall Survival Probability", xlim=c(0,mt), ylim=c(0,1),
       xaxt="n", yaxt="n", mark.time=TRUE, cex=0.6, cex.lab=1.15)
  axis(1, at=tb); axis(2, at=seq(0,1,0.2), las=2)
  title(main=s$label, line=1.5, cex.main=1.2, font.main=2)
  abline(h=seq(0,1,0.2), col="gray90", lty=3, lwd=0.5)
  legend("topright", legend=c(sprintf("Low Risk (n=%d)", s$nl),
         sprintf("High Risk (n=%d)", s$nh)), col=cols, lwd=3, bty="n", cex=1.05)
  legend("bottomleft", legend=c(paste("Log-rank", format_pval(s$pv)),
         sprintf("HR = %.2f (%.2f-%.2f)", s$hr, s$lo, s$hi),
         sprintf("C-index = %.3f", s$cidx)), bty="n", cex=1.0, y.intersp=1.2)
  par(mar=c(3.5,5,0,1.5)); plot.new(); plot.window(xlim=c(0,mt), ylim=c(0,3))
  text(0, 2.5, "Number at risk", adj=0, font=2)
  rar <- function(grp) sapply(tb, function(t) sum(s$df$Group==grp & s$df$OS_time>=t))
  text(-mt*0.02, 1.5, "Low Risk",  adj=1, col=cols[1], font=2, cex=0.95)
  text(-mt*0.02, 0.6, "High Risk", adj=1, col=cols[2], font=2, cex=0.95)
  lr <- rar("Low"); hrk <- rar("High")
  for (i in seq_along(tb)) { text(tb[i],1.5,lr[i],col=cols[1],cex=0.95)
                             text(tb[i],0.6,hrk[i],col=cols[2],cex=0.95) }
  axis(1, at=tb); mtext("Time (months)", side=1, line=2.2, cex=0.9)
}

# ---- 4) 跑 HCC + CRC ----
defs <- list(HCC=list(tag="HCC", label="TCGA-LIHC / HCC"),
             CRC=list(tag="CRC", label="TCGA-COAD+READ / CRC"))
stats <- list(); summ <- data.frame(); rs_all <- data.frame()
for (k in names(defs)) {
  co <- load_cohort(defs[[k]]$tag)
  s  <- compute_stats(co, defs[[k]]$label); stats[[k]] <- s
  summ <- rbind(summ, data.frame(cohort=k, n=s$n, events=sum(s$df$OS_status),
    cov_170=s$cov170, cov_69=s$cov69, c_index=round(s$cidx,3),
    HR=round(s$hr,2), HR_lo=round(s$lo,2), HR_hi=round(s$hi,2),
    KM_p=signif(s$pv,3), mv_HR_perSD=round(s$mv_hr,2), mv_p=signif(s$mv_p,3)))
  rs_all <- rbind(rs_all, data.frame(cohort=k, s$df[,c("sample_id","RS","OS_time","OS_status","Group")]))
}
fwrite(summ,   file.path(OUT, "pancancer_summary.tsv"), sep="\t")
fwrite(rs_all, file.path(OUT, "pancancer_nefrs_rs.tsv"), sep="\t")
cat("\n=== 泛癌验证汇总 ===\n"); print(summ)

# ---- 5) 出图 ----
for (k in names(defs)) {
  pdf(file.path(FIG, sprintf("Fig17_KM_%s.pdf", k)), width=7, height=7.5, pointsize=10)
  plot_km(stats[[k]]); dev.off()
}
pdf(file.path(FIG, "Fig17_pancancer_KM_combined.pdf"), width=14, height=8, pointsize=10)
par(mfrow=c(1,2), oma=c(0,0,3,0))
for (k in names(defs)) plot_km(stats[[k]])
mtext("NeFRS Pan-cancer Validation (HCC, CRC)", side=3, line=0.5, outer=TRUE, cex=1.4, font=2)
dev.off()
cat("\n✓ 图: 03_figures/Fig17_pancancer/  数据: 02_results/08_pancancer/\nDone.\n")
