# W3 Day 3 Step 1: 8-cohort KM curves (publication quality)
suppressPackageStartupMessages({
  library(data.table); library(survival); library(glmnet)
})

W2_IN <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/02_w2_inputs"
MIME  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime"
FIG   <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"
dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

cat("=== W3 Step 1: 8 cohort KM curves (publication quality) ===\n")

cohort_data <- readRDS(file.path(W2_IN, "cohort_data_8.rds"))
mime_input  <- readRDS(file.path(MIME, "mime_input_list.rds"))
res         <- readRDS(file.path(MIME, "mime_results.rds"))
mime_genes  <- colnames(mime_input[[1]])[-(1:3)]
ml <- res$ml.res[["Enet[α=0.1]"]]

compute_rs <- function(cn) {
  co <- cohort_data[[cn]]
  X <- matrix(0, nrow=ncol(co$expression), ncol=length(mime_genes),
              dimnames=list(colnames(co$expression), mime_genes))
  avail <- intersect(mime_genes, rownames(co$expression))
  for (g in avail) X[, g] <- co$expression[g, ]
  rs <- as.numeric(predict(ml, newx=X, s="lambda.min"))
  names(rs) <- rownames(X)
  clin <- co$clinical
  m <- match(clin$sample_id, names(rs))
  data.frame(sample_id=clin$sample_id, OS_time=clin$os_time_month,
             OS_status=clin$os_status, RS=rs[m])
}

cohorts <- list(
  list(name="PACA_AU",     label="ICGC PACA-AU (Training, n=267)",      role="TRAIN"),
  list(name="GSE71729",    label="GSE71729 (Sel-Val, n=123)",            role="SEL"),
  list(name="GSE85916",    label="GSE85916 (Sel-Val, n=79)",             role="SEL"),
  list(name="GSE62452",    label="GSE62452 (Sel-Val, n=66)",             role="SEL"),
  list(name="GSE57495",    label="GSE57495 (Sel-Val, n=63)",             role="SEL"),
  list(name="GSE21501",    label="GSE21501 (Hold-out, n=102)",           role="HO"),
  list(name="TCGA_PAAD",   label="TCGA-PAAD (Hold-out, n=147)",          role="HO"),
  list(name="E_MTAB_6134", label="E-MTAB-6134 (Hold-out, n=288)",        role="HO")
)

# 修复 1: factor 用 c("Low","High") -> Low 作 reference, "High vs Low" HR > 1
compute_stats <- function(co_info) {
  rs_df <- compute_rs(co_info$name)
  rs_df <- rs_df[!is.na(rs_df$RS) & !is.na(rs_df$OS_time) & rs_df$OS_time > 0, ]
  
  med <- median(rs_df$RS, na.rm=TRUE)
  rs_df$Group <- factor(ifelse(rs_df$RS > med, "High", "Low"),
                        levels=c("Low","High"))  # Low = reference
  
  fit_km <- survfit(Surv(OS_time, OS_status) ~ Group, data=rs_df)
  cox <- coxph(Surv(OS_time, OS_status) ~ Group, data=rs_df)
  hr <- as.numeric(exp(coef(cox)))   # 系数对应 "GroupHigh", 是 High vs Low
  hr_ci <- exp(confint(cox))
  sd_obj <- survdiff(Surv(OS_time, OS_status) ~ Group, data=rs_df)
  pv <- 1 - pchisq(sd_obj$chisq, df=1)
  ci_val <- as.numeric(concordance(Surv(OS_time, OS_status) ~ RS,
                                    data=rs_df, reverse=TRUE)$concordance)
  
  list(fit_km=fit_km, rs_df=rs_df, hr=hr, hr_lo=hr_ci[1,1], hr_hi=hr_ci[1,2],
       pv=pv, ci_val=ci_val, n=nrow(rs_df),
       n_high=sum(rs_df$Group=="High"), n_low=sum(rs_df$Group=="Low"))
}

# 修复 2: P 值科学计数,显式处理边界值
format_pval <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-15) return("P < 1e-15")
  if (p < 1e-4) return(sprintf("P = %.1e", p))
  if (p < 0.001) return(sprintf("P = %.4f", p))
  return(sprintf("P = %.3f", p))
}

# 修复 3: 高质量 KM (加 risk table + 优化字体)
plot_km <- function(stats, co_info) {
  pv_text <- paste("Log-rank", format_pval(stats$pv))
  hr_text <- sprintf("HR = %.2f (%.2f-%.2f)", stats$hr, stats$hr_lo, stats$hr_hi)
  ci_text <- sprintf("C-index = %.3f", stats$ci_val)
  
  cols <- c("#0072B2", "#C00000")  # Low blue, High red
  
  # 主图区 (上 75%) + risk table 区 (下 25%)
  layout(matrix(c(1,2), ncol=1), heights=c(3.5, 1.2))
  
  # 主 KM 图
  par(mar=c(3.5, 5, 3.5, 1.5), mgp=c(2.2, 0.7, 0))
  max_time <- max(stats$rs_df$OS_time, na.rm=TRUE)
  time_breaks <- pretty(c(0, max_time), n=5)
  
  plot(stats$fit_km, col=cols, lwd=3, conf.int=FALSE,
       xlab="Time (months)", ylab="Overall Survival Probability",
       xlim=c(0, max_time), ylim=c(0,1),
       xaxt="n", yaxt="n",
       mark.time=TRUE, cex=0.6,
       cex.lab=1.15, font.lab=1)
  
  # 自定义坐标轴
  axis(1, at=time_breaks, labels=time_breaks, cex.axis=1.0)
  axis(2, at=seq(0,1,0.2), labels=seq(0,1,0.2), las=2, cex.axis=1.0)
  
  # 主标题
  title(main=co_info$label, line=1.5, cex.main=1.2, font.main=2)
  
  # 网格(可选,科学论文喜欢)
  abline(h=seq(0,1,0.2), col="gray90", lty=3, lwd=0.5)
  
  # Legend 右上 (risk 分组 + n)
  legend("topright",
         legend=c(sprintf("Low Risk (n=%d)", stats$n_low),
                  sprintf("High Risk (n=%d)", stats$n_high)),
         col=cols, lwd=3, bty="n", cex=1.05, x.intersp=0.5)
  
  # Legend 左下 (统计量)
  legend("bottomleft",
         legend=c(pv_text, hr_text, ci_text),
         bty="n", cex=1.0, text.col="black", y.intersp=1.2)
  
  # === Risk table ===
  par(mar=c(3.5, 5, 0, 1.5))
  plot.new()
  plot.window(xlim=c(0, max_time), ylim=c(0, 3))
  
  # Risk table 标题
  text(0, 2.5, "Number at risk", adj=0, font=2, cex=1.0)
  
  # 计算每个时间点剩余样本数
  risk_at_times <- function(group_name) {
    sapply(time_breaks, function(t) {
      sum(stats$rs_df$Group == group_name & stats$rs_df$OS_time >= t)
    })
  }
  
  low_at_risk <- risk_at_times("Low")
  high_at_risk <- risk_at_times("High")
  
  # 行标签
  text(-max_time*0.02, 1.5, "Low Risk", adj=1, col=cols[1], font=2, cex=0.95)
  text(-max_time*0.02, 0.6, "High Risk", adj=1, col=cols[2], font=2, cex=0.95)
  
  # 数字
  for (i in seq_along(time_breaks)) {
    text(time_breaks[i], 1.5, low_at_risk[i], cex=0.95, col=cols[1])
    text(time_breaks[i], 0.6, high_at_risk[i], cex=0.95, col=cols[2])
  }
  
  # X 轴用同样的刻度
  axis(1, at=time_breaks, labels=time_breaks, cex.axis=1.0)
  mtext("Time (months)", side=1, line=2.2, cex=0.9)
}

# === 计算所有 cohort ===
cat("\nComputing statistics for 8 cohorts...\n")
all_stats <- list()
summary_df <- data.frame()
for (i in seq_along(cohorts)) {
  ci <- cohorts[[i]]
  cat(sprintf("  [%d/8] %s ...\n", i, ci$name))
  all_stats[[i]] <- compute_stats(ci)
  summary_df <- rbind(summary_df, data.frame(
    cohort=ci$name, role=ci$role, n=all_stats[[i]]$n,
    c_index=round(all_stats[[i]]$ci_val, 3),
    HR=round(all_stats[[i]]$hr, 2),
    HR_lo=round(all_stats[[i]]$hr_lo, 2),
    HR_hi=round(all_stats[[i]]$hr_hi, 2),
    KM_p=signif(all_stats[[i]]$pv, 3)
  ))
}

fwrite(summary_df, file.path(FIG, "Fig3_KM_summary.tsv"), sep="\t")

cat("\n=== Summary (HR should all be > 1 now) ===\n")
print(summary_df)

# === 单 cohort PDF (高分辨率,放大尺寸) ===
cat("\nGenerating individual KM PDFs...\n")
for (i in seq_along(cohorts)) {
  outpath <- file.path(FIG, sprintf("Fig3_KM_%02d_%s.pdf", i, cohorts[[i]]$name))
  pdf(outpath, width=7, height=7.5, pointsize=10)
  plot_km(all_stats[[i]], cohorts[[i]])
  dev.off()
  sz <- file.info(outpath)$size
  cat(sprintf("  saved %-35s %d bytes\n", basename(outpath), sz))
}

# === 8-panel combined PDF ===
cat("\nGenerating combined 8-panel PDF...\n")
pdf(file.path(FIG, "Fig3_KM_8cohort_combined.pdf"), width=18, height=24, pointsize=10)
par(mfrow=c(4, 2), oma=c(0, 0, 3, 0))
for (i in seq_along(cohorts)) {
  plot_km(all_stats[[i]], cohorts[[i]])
}
mtext("NeFRS v3 Risk Stratification across 8 PDAC Cohorts (Enet 69 genes)",
      side=3, line=0.5, outer=TRUE, cex=1.5, font=2)
dev.off()
combined_sz <- file.info(file.path(FIG, "Fig3_KM_8cohort_combined.pdf"))$size
cat(sprintf("  Combined PDF: %d bytes\n", combined_sz))

cat("\nDone.\n")
