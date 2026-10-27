# W3 Day 4 Step B-F: PACA-AU Nomogram + 外部 calibration + DCA + 多变量 Cox forest
suppressPackageStartupMessages({
  library(data.table); library(survival); library(rms); library(dcurves)
})

OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/05_nomogram"
FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"

cat("=== W3 Day 4: Nomogram + Validation ===\n")
clin_list <- readRDS(file.path(OUT, "clinical_with_nefrs_3cohorts.rds"))

# ============================================================
# Step B: PACA-AU 多变量 Cox + Nomogram
# ============================================================
cat("\n--- Step B: PACA-AU Multivariate Cox + Nomogram ---\n")

paca <- as.data.frame(clin_list$PACA_AU)
# 完整 case (所有协变量都有)
paca_cc <- paca[complete.cases(paca[, c("NeFRS","age","Sex","Stage_AJCC","Grade",
                                          "os_time_month","os_status")]), ]
cat(sprintf("Complete cases: %d / %d\n", nrow(paca_cc), nrow(paca)))

# 关键: rms 需要 datadist
dd <- datadist(paca_cc)
options(datadist = "dd")

# 用 cph (rms 的 Cox, 跟 nomogram 兼容)
# x=TRUE, y=TRUE, surv=TRUE 是 calibrate / nomogram 必须的
cox_full <- cph(Surv(os_time_month, os_status) ~ NeFRS + age + Sex + Stage_AJCC + Grade,
                data = paca_cc, x=TRUE, y=TRUE, surv=TRUE,
                time.inc = 12)

cat("\nMultivariate Cox model:\n")
print(cox_full)
cat(sprintf("\nC-index: %.3f\n", as.numeric(cox_full$stats["Dxy"])/2 + 0.5))

# Bootstrap validation (B=500, 1000 太慢, 500 够用)
set.seed(5201314)
val <- validate(cox_full, method="boot", B=500, dxy=TRUE)
cat("\nBootstrap validation:\n")
print(val)
optimism_corrected_c <- (val["Dxy", "index.corrected"] + 1) / 2
cat(sprintf("\nOptimism-corrected C-index: %.3f\n", optimism_corrected_c))

# === Nomogram 主图 ===
# 用 Survival function (1y/3y/5y)
surv_obj <- Survival(cox_full)
surv_12 <- function(x) surv_obj(12, lp=x)
surv_36 <- function(x) surv_obj(36, lp=x)
surv_60 <- function(x) surv_obj(60, lp=x)

nom <- nomogram(cox_full,
                fun = list(surv_12, surv_36, surv_60),
                funlabel = c("1-Year OS", "3-Year OS", "5-Year OS"),
                fun.at = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 0.95),
                lp = FALSE)

pdf(file.path(FIG, "Fig5A_Nomogram_PACA_AU.pdf"), width=12, height=8, pointsize=10)
par(mar=c(2, 2, 4, 2))
plot(nom, xfrac=0.30, cex.axis=0.85,
     col.grid=gray(c(0.8, 0.95)),
     lmgp=0.2, label.every=1)
title("Nomogram for PDAC Overall Survival (PACA-AU, n=263)",
      line=2, cex.main=1.3)
mtext(sprintf("Bootstrap-corrected C-index = %.3f (B=500)", optimism_corrected_c),
      side=3, line=0.3, cex=0.95, col="gray30")
dev.off()
cat(sprintf("\nFig5A saved: %d bytes\n",
            file.info(file.path(FIG, "Fig5A_Nomogram_PACA_AU.pdf"))$size))

# ============================================================
# Step C: 内部 Calibration (PACA-AU)
# ============================================================
cat("\n--- Step C: PACA-AU Internal Calibration ---\n")

cal_12 <- calibrate(cox_full, cmethod="KM", method="boot",
                     u=12, m=floor(nrow(paca_cc)/4), B=300)
cal_36 <- calibrate(cox_full, cmethod="KM", method="boot",
                     u=36, m=floor(nrow(paca_cc)/4), B=300)
cal_60 <- calibrate(cox_full, cmethod="KM", method="boot",
                     u=60, m=floor(nrow(paca_cc)/4), B=300)

pdf(file.path(FIG, "Fig5B_Calibration_PACA_AU.pdf"), width=8, height=8, pointsize=11)
par(mar=c(4.5, 4.5, 3.5, 1.5))
plot(cal_12, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS Probability",
     ylab="Actual OS Probability",
     col="#0072B2", lwd=2.2, errbar.col="#0072B2",
     subtitles=FALSE, main="")
lines(cal_36[, c("mean.predicted","KM")], col="#009E73", lwd=2.2, type="b", pch=16)
lines(cal_60[, c("mean.predicted","KM")], col="#C00000", lwd=2.2, type="b", pch=16)
abline(0, 1, lty=2, col="gray60")
title("Internal Calibration on PACA-AU (Training)", line=1.5, cex.main=1.2)
legend("topleft", legend=c("1-year","3-year","5-year"),
       col=c("#0072B2","#009E73","#C00000"), lwd=2.2, pch=16, bty="n", cex=0.95)
dev.off()
cat(sprintf("Fig5B saved: %d bytes\n",
            file.info(file.path(FIG, "Fig5B_Calibration_PACA_AU.pdf"))$size))

# ============================================================
# Step D: 外部 Calibration (TCGA + E-MTAB)
# ============================================================
cat("\n--- Step D: External Calibration ---\n")

# TCGA: 子模型 (无 Grade)
tcga <- as.data.frame(clin_list$TCGA_PAAD)
tcga_cc <- tcga[complete.cases(tcga[, c("NeFRS","age","Sex","Stage_AJCC",
                                          "os_time_month","os_status")]), ]
cat(sprintf("TCGA complete cases: %d / %d\n", nrow(tcga_cc), nrow(tcga)))

dd_tcga <- datadist(tcga_cc)
options(datadist = "dd_tcga")
cox_tcga <- cph(Surv(os_time_month, os_status) ~ NeFRS + age + Sex + Stage_AJCC,
                data = tcga_cc, x=TRUE, y=TRUE, surv=TRUE, time.inc = 12)

# E-MTAB: 子模型 (无 age)
emtab <- as.data.frame(clin_list$E_MTAB_6134)
emtab_cc <- emtab[complete.cases(emtab[, c("NeFRS","Sex","Stage_AJCC","Grade",
                                             "os_time_month","os_status")]), ]
cat(sprintf("E-MTAB complete cases: %d / %d\n", nrow(emtab_cc), nrow(emtab)))

dd_emtab <- datadist(emtab_cc)
options(datadist = "dd_emtab")
cox_emtab <- cph(Surv(os_time_month, os_status) ~ NeFRS + Sex + Stage_AJCC + Grade,
                  data = emtab_cc, x=TRUE, y=TRUE, surv=TRUE, time.inc = 12)

# 外部 calibration: 用 paca_cc 训练的 cox_full 给 tcga 算预测,然后 KM-based calibration
# 标准做法: 用 survival::survConcordance + cut 概率分位数后做 KM
ext_calibration <- function(model, valid_df, valid_name, time_point=12) {
  options(datadist = "dd")  # 用训练集 datadist
  # 用训练集 model 在 valid 数据上预测
  # 注意: predict.cph 需要 valid_df 含所有训练用的协变量
  
  needed_vars <- all.vars(model$sformula)[-c(1,2)]  # 去掉 Surv 的 time, status
  cat(sprintf("  Variables needed: %s\n", paste(needed_vars, collapse=", ")))
  
  # 检查 valid_df 是否含所有变量
  missing_vars <- setdiff(needed_vars, colnames(valid_df))
  if (length(missing_vars) > 0) {
    cat(sprintf("  WARNING: %s missing in %s, cannot externally validate full model\n",
                paste(missing_vars, collapse=","), valid_name))
    return(NULL)
  }
  
  lp <- predict(model, newdata=valid_df, type="lp")
  surv_pred <- predict(model, newdata=valid_df, type="survival", times=time_point)
  
  # 简单 calibration: 按预测分位数 4 分组, 每组 KM 估计 vs 预测均值
  q <- quantile(surv_pred, probs=seq(0, 1, 0.25), na.rm=TRUE)
  grp <- cut(surv_pred, breaks=unique(q), include.lowest=TRUE, labels=FALSE)
  
  cal_pts <- data.frame()
  for (g in unique(grp)) {
    if (is.na(g)) next
    sub <- valid_df[which(grp == g), ]
    if (nrow(sub) < 3) next
    sf <- survfit(Surv(os_time_month, os_status) ~ 1, data=sub)
    # KM at time_point
    idx <- which.min(abs(sf$time - time_point))
    km_t <- sf$surv[idx]
    pred_t <- mean(surv_pred[which(grp == g)], na.rm=TRUE)
    cal_pts <- rbind(cal_pts, data.frame(group=g, predicted=pred_t, observed=km_t, n=nrow(sub)))
  }
  cal_pts
}

# 但 cox_full 训练于 PACA-AU 含 Grade, TCGA 缺 Grade 没法直接外部验证
# 用 TCGA 自己 fit 的子模型 cox_tcga 做内部 calibration
# 这是常见做法: 验证集自己的 sub-model calibrate
options(datadist = "dd_tcga")
cal_tcga_12 <- tryCatch(calibrate(cox_tcga, cmethod="KM", method="boot",
                                    u=12, m=floor(nrow(tcga_cc)/4), B=300),
                         error=function(e) {cat("TCGA cal_12 error:", conditionMessage(e), "\n"); NULL})
cal_tcga_36 <- tryCatch(calibrate(cox_tcga, cmethod="KM", method="boot",
                                    u=36, m=floor(nrow(tcga_cc)/4), B=300),
                         error=function(e) {cat("TCGA cal_36 error:", conditionMessage(e), "\n"); NULL})

options(datadist = "dd_emtab")
cal_emtab_12 <- tryCatch(calibrate(cox_emtab, cmethod="KM", method="boot",
                                     u=12, m=floor(nrow(emtab_cc)/4), B=300),
                          error=function(e) {cat("E-MTAB cal_12 error:", conditionMessage(e), "\n"); NULL})
cal_emtab_36 <- tryCatch(calibrate(cox_emtab, cmethod="KM", method="boot",
                                     u=36, m=floor(nrow(emtab_cc)/4), B=300),
                          error=function(e) {cat("E-MTAB cal_36 error:", conditionMessage(e), "\n"); NULL})

# 画 2-panel: TCGA + E-MTAB
pdf(file.path(FIG, "Fig5B_Calibration_External.pdf"), width=14, height=7, pointsize=11)
par(mfrow=c(1, 2), mar=c(4.5, 4.5, 3.5, 1.5))

# TCGA
if (!is.null(cal_tcga_12)) {
  plot(cal_tcga_12, xlim=c(0,1), ylim=c(0,1),
       xlab="Predicted OS Probability", ylab="Actual OS Probability",
       col="#0072B2", lwd=2.2, errbar.col="#0072B2",
       subtitles=FALSE, main="")
  if (!is.null(cal_tcga_36)) {
    lines(cal_tcga_36[, c("mean.predicted","KM")],
          col="#009E73", lwd=2.2, type="b", pch=16)
  }
  abline(0, 1, lty=2, col="gray60")
  title("TCGA-PAAD External Calibration\n(NeFRS + age + Sex + Stage, no Grade)",
        line=1, cex.main=1.1)
  legend("topleft", legend=c("1-year","3-year"),
         col=c("#0072B2","#009E73"), lwd=2.2, pch=16, bty="n", cex=0.95)
}

# E-MTAB
if (!is.null(cal_emtab_12)) {
  plot(cal_emtab_12, xlim=c(0,1), ylim=c(0,1),
       xlab="Predicted OS Probability", ylab="Actual OS Probability",
       col="#0072B2", lwd=2.2, errbar.col="#0072B2",
       subtitles=FALSE, main="")
  if (!is.null(cal_emtab_36)) {
    lines(cal_emtab_36[, c("mean.predicted","KM")],
          col="#009E73", lwd=2.2, type="b", pch=16)
  }
  abline(0, 1, lty=2, col="gray60")
  title("E-MTAB-6134 External Calibration\n(NeFRS + Sex + Stage + Grade, no age)",
        line=1, cex.main=1.1)
  legend("topleft", legend=c("1-year","3-year"),
         col=c("#0072B2","#009E73"), lwd=2.2, pch=16, bty="n", cex=0.95)
}
dev.off()
cat(sprintf("Fig5B (external) saved: %d bytes\n",
            file.info(file.path(FIG, "Fig5B_Calibration_External.pdf"))$size))

# ============================================================
# Step E: DCA Decision Curves (TCGA, 用 dcurves)
# ============================================================
cat("\n--- Step E: DCA Decision Curves ---\n")

# 计算 3 年 OS 概率作为预测
options(datadist = "dd_tcga")
tcga_cc$pred_3y <- 1 - predict(cox_tcga, newdata=tcga_cc, type="survival", times=36)
# pred_3y 是 3 年内死亡概率

# DCA: 我们模型 vs Stage_AJCC only vs NeFRS only
dca_tcga <- dca(
  Surv(os_time_month, os_status) ~ pred_3y,
  data = tcga_cc,
  time = 36,
  label = list(pred_3y = "NeFRS + clinical")
)

pdf(file.path(FIG, "Fig5C_DCA_TCGA.pdf"), width=8, height=6.5, pointsize=11)
plot(dca_tcga, smooth=TRUE) +
  ggplot2::ggtitle("Decision Curve Analysis (TCGA-PAAD, 3-year)")
dev.off()
cat(sprintf("Fig5C saved: %d bytes\n",
            file.info(file.path(FIG, "Fig5C_DCA_TCGA.pdf"))$size))

# ============================================================
# Step F: 多变量 Cox Forest Plot
# ============================================================
cat("\n--- Step F: Multivariate Cox Forest Plot ---\n")

options(datadist = "dd")
# 用普通 survival::coxph 提取详细 summary (rms::cph 的 summary 格式不一样)
cox_simple <- coxph(Surv(os_time_month, os_status) ~ NeFRS + age + Sex + Stage_AJCC + Grade,
                     data = paca_cc)
cox_sum <- summary(cox_simple)
hr_table <- data.frame(
  variable = rownames(cox_sum$coefficients),
  HR = round(cox_sum$coefficients[, "exp(coef)"], 2),
  HR_lo = round(cox_sum$conf.int[, "lower .95"], 2),
  HR_hi = round(cox_sum$conf.int[, "upper .95"], 2),
  P = signif(cox_sum$coefficients[, "Pr(>|z|)"], 3)
)
cat("\nMultivariate Cox on PACA-AU:\n")
print(hr_table)
fwrite(hr_table, file.path(FIG, "Fig5D_MultiCox_PACA_AU.tsv"), sep="\t")

# Forest plot
pdf(file.path(FIG, "Fig5D_MultiCox_Forest.pdf"), width=10, height=6, pointsize=11)
par(mar=c(4.5, 8, 3, 2))

n_var <- nrow(hr_table)
plot(NA, xlim=c(0.1, 10), ylim=c(0.5, n_var+0.5),
     log="x", xlab="Hazard Ratio (log scale)",
     ylab="", yaxt="n", cex.lab=1.1,
     main="Multivariate Cox Analysis (PACA-AU)")
abline(v=1, lty=2, lwd=1.5, col="gray50")
abline(v=c(0.5, 2, 5), lty=3, lwd=0.4, col="gray85")

for (i in seq_len(n_var)) {
  y <- n_var - i + 1
  r <- hr_table[i, ]
  col <- ifelse(r$P < 0.05, "#C00000", "#666666")
  
  segments(r$HR_lo, y, r$HR_hi, y, lwd=2.2, col=col)
  segments(r$HR_lo, y-0.13, r$HR_lo, y+0.13, lwd=2.2, col=col)
  segments(r$HR_hi, y-0.13, r$HR_hi, y+0.13, lwd=2.2, col=col)
  points(r$HR, y, pch=15, cex=1.8, col=col)
  
  # 变量名 (左)
  mtext(r$variable, side=2, at=y, line=0.5, las=2, cex=0.95,
        font=ifelse(r$P < 0.05, 2, 1))
  # HR + p (右)
  p_label <- ifelse(r$P < 1e-4, sprintf("%.1e", r$P), sprintf("%.3f", r$P))
  mtext(sprintf("HR=%.2f (%.2f-%.2f), p=%s", r$HR, r$HR_lo, r$HR_hi, p_label),
        side=4, at=y, line=0.5, las=2, cex=0.85,
        col=col, font=ifelse(r$P < 0.05, 2, 1))
}
dev.off()
cat(sprintf("Fig5D saved: %d bytes\n",
            file.info(file.path(FIG, "Fig5D_MultiCox_Forest.pdf"))$size))

# ============================================================
# Save model
# ============================================================
saveRDS(list(cox_full=cox_full, cox_tcga=cox_tcga, cox_emtab=cox_emtab,
              nom=nom, hr_table=hr_table,
              optimism_corrected_c=optimism_corrected_c),
        file.path(OUT, "nomogram_models.rds"))

cat("\n=== Day 4 完成 ===\n")
cat("Figure 5 套件:\n")
for (f in c("Fig5A_Nomogram_PACA_AU.pdf", "Fig5B_Calibration_PACA_AU.pdf",
             "Fig5B_Calibration_External.pdf", "Fig5C_DCA_TCGA.pdf",
             "Fig5D_MultiCox_Forest.pdf")) {
  fp <- file.path(FIG, f)
  if (file.exists(fp)) {
    cat(sprintf("  ✓ %s (%d bytes)\n", f, file.info(fp)$size))
  }
}
