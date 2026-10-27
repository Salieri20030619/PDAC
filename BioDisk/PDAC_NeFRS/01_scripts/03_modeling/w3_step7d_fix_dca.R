# W3 Day 4 - 修复 Step E (DCA) + 跑 Step F
suppressPackageStartupMessages({
  library(data.table); library(survival); library(rms); library(dcurves); library(ggplot2)
})

OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/05_nomogram"
FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"

clin_list <- readRDS(file.path(OUT, "clinical_with_nefrs_3cohorts.rds"))

# === Reload 数据 ===
paca <- as.data.frame(clin_list$PACA_AU)
paca_cc <- paca[complete.cases(paca[, c("NeFRS","age","Sex","Stage_AJCC","Grade",
                                          "os_time_month","os_status")]), ]

tcga <- as.data.frame(clin_list$TCGA_PAAD)
tcga_cc <- tcga[complete.cases(tcga[, c("NeFRS","age","Sex","Stage_AJCC",
                                          "os_time_month","os_status")]),
                 c("sample_id","os_time_month","os_status","NeFRS","age","Sex","Stage_AJCC")]
cat(sprintf("TCGA complete: %d\n", nrow(tcga_cc)))

# ============================================================
# Step E: DCA on TCGA — 用 survival::coxph 而不是 rms::cph
# ============================================================
cat("\n--- Step E: DCA (using survival::coxph) ---\n")

# 用 base survival::coxph (支持 type="expected" 来算 cumulative hazard)
cox_tcga_surv <- coxph(Surv(os_time_month, os_status) ~ NeFRS + age + Sex + Stage_AJCC,
                        data = tcga_cc)

# Survival probability at t=36 months
# 方法: 用 survfit() 给每个新样本预测 survival curve
# 然后在 t=36 取值
sf_pred <- survfit(cox_tcga_surv, newdata = tcga_cc)
# sf_pred$time 是时间向量, sf_pred$surv 是矩阵 (n_times × n_samples)

# 找最接近 36 月的时间点
t_idx <- which.min(abs(sf_pred$time - 36))
cat(sprintf("Closest time point to 36 months: %.2f\n", sf_pred$time[t_idx]))

# 提取每个样本在该时点的 survival prob
surv_36 <- sf_pred$surv[t_idx, ]
tcga_cc$risk_3y <- 1 - surv_36  # 3 年内死亡概率

cat(sprintf("Risk_3y range: [%.3f, %.3f], median=%.3f\n",
            min(tcga_cc$risk_3y), max(tcga_cc$risk_3y), median(tcga_cc$risk_3y)))

# === DCA ===
# 同时对比: NeFRS only / Stage only / NeFRS + Clinical (full model)
# 这样可以显示 NeFRS 加入临床变量的增量价值

# Sub-model 1: NeFRS only
cox_nefrs_only <- coxph(Surv(os_time_month, os_status) ~ NeFRS, data = tcga_cc)
sf_nefrs <- survfit(cox_nefrs_only, newdata = tcga_cc)
t_idx <- which.min(abs(sf_nefrs$time - 36))
tcga_cc$risk_nefrs <- 1 - sf_nefrs$surv[t_idx, ]

# Sub-model 2: Clinical only (no NeFRS)
cox_clin_only <- coxph(Surv(os_time_month, os_status) ~ age + Sex + Stage_AJCC, data = tcga_cc)
sf_clin <- survfit(cox_clin_only, newdata = tcga_cc)
t_idx <- which.min(abs(sf_clin$time - 36))
tcga_cc$risk_clin <- 1 - sf_clin$surv[t_idx, ]

cat(sprintf("risk_nefrs range: [%.3f, %.3f]\n",
            min(tcga_cc$risk_nefrs), max(tcga_cc$risk_nefrs)))
cat(sprintf("risk_clin  range: [%.3f, %.3f]\n",
            min(tcga_cc$risk_clin), max(tcga_cc$risk_clin)))

# DCA with 3 models
dca_tcga <- tryCatch({
  dca(Surv(os_time_month, os_status) ~ risk_3y + risk_nefrs + risk_clin,
      data = tcga_cc, time = 36,
      label = list(
        risk_3y    = "NeFRS + Clinical (full)",
        risk_nefrs = "NeFRS only",
        risk_clin  = "Clinical only"
      ))
}, error = function(e) {
  cat(sprintf("DCA error: %s\n", conditionMessage(e))); NULL
})

if (!is.null(dca_tcga)) {
  pdf(file.path(FIG, "Fig5C_DCA_TCGA.pdf"), width=9, height=6.5, pointsize=11)
  p_dca <- plot(dca_tcga, smooth=TRUE) +
    ggplot2::scale_color_manual(
      values = c(
        "NeFRS + Clinical (full)" = "#C00000",
        "NeFRS only" = "#0072B2",
        "Clinical only" = "#009E73",
        "Treat All" = "gray60",
        "Treat None" = "gray40"
      )
    ) +
    ggplot2::ggtitle("Decision Curve Analysis (TCGA-PAAD, 3-year mortality)") +
    ggplot2::theme(plot.title = ggplot2::element_text(size=12, face="bold"),
                   legend.position = "right")
  print(p_dca)
  dev.off()
  cat(sprintf("Fig5C saved: %d bytes\n",
              file.info(file.path(FIG, "Fig5C_DCA_TCGA.pdf"))$size))
}

# ============================================================
# Step F: 多变量 Cox Forest (复制之前代码, 确认运行)
# ============================================================
cat("\n--- Step F: Multivariate Cox Forest ---\n")

cox_simple <- coxph(Surv(os_time_month, os_status) ~ NeFRS + age + Sex + Stage_AJCC + Grade,
                     data = paca_cc)
cox_sum <- summary(cox_simple)

var_names_pretty <- c(
  "NeFRS" = "NeFRS (per unit)",
  "age" = "Age (per year)",
  "SexMale" = "Sex (Male vs Female)",
  "Stage_AJCCIIA" = "Stage IIA (vs I)",
  "Stage_AJCCIIB" = "Stage IIB (vs I)",
  "Stage_AJCCIII-IV" = "Stage III-IV (vs I)",
  "GradeG2" = "Grade G2 (vs G1)",
  "GradeG3-G4" = "Grade G3-G4 (vs G1)"
)
hr_table <- data.frame(
  variable_raw = rownames(cox_sum$coefficients),
  variable = unname(var_names_pretty[rownames(cox_sum$coefficients)]),
  HR = round(cox_sum$coefficients[, "exp(coef)"], 2),
  HR_lo = round(cox_sum$conf.int[, "lower .95"], 2),
  HR_hi = round(cox_sum$conf.int[, "upper .95"], 2),
  P = signif(cox_sum$coefficients[, "Pr(>|z|)"], 3)
)
print(hr_table)
fwrite(hr_table, file.path(FIG, "Fig5D_MultiCox_PACA_AU.tsv"), sep="\t")

pdf(file.path(FIG, "Fig5D_MultiCox_Forest.pdf"), width=12, height=6, pointsize=11)
par(mar=c(4.5, 11, 3.5, 13))
n_var <- nrow(hr_table)
plot(NA, xlim=c(0.3, 8), ylim=c(0.5, n_var+0.5),
     log="x", xlab="Hazard Ratio (log scale)", ylab="", yaxt="n",
     cex.lab=1.05,
     main=sprintf("Multivariate Cox Analysis (PACA-AU, n=%d)", nrow(paca_cc)))
abline(v=1, lty=2, lwd=1.5, col="gray50")
abline(v=c(0.5, 2, 5), lty=3, lwd=0.4, col="gray85")

for (i in seq_len(n_var)) {
  y <- n_var - i + 1
  r <- hr_table[i, ]
  col <- ifelse(r$P < 0.05, "#C00000", "#666666")
  segments(r$HR_lo, y, r$HR_hi, y, lwd=2.2, col=col)
  segments(r$HR_lo, y-0.13, r$HR_lo, y+0.13, lwd=2.2, col=col)
  segments(r$HR_hi, y-0.13, r$HR_hi, y+0.13, lwd=2.2, col=col)
  points(r$HR, y, pch=15, cex=1.6, col=col)
  
  mtext(r$variable, side=2, at=y, line=0.5, las=2, cex=0.92,
        font=ifelse(r$P < 0.05, 2, 1), col=col)
  p_label <- ifelse(r$P < 1e-4, sprintf("%.1e", r$P), sprintf("%.3f", r$P))
  mtext(sprintf("%.2f (%.2f-%.2f), p=%s", r$HR, r$HR_lo, r$HR_hi, p_label),
        side=4, at=y, line=0.3, las=2, cex=0.82,
        col=col, font=ifelse(r$P < 0.05, 2, 1))
}
dev.off()
cat(sprintf("Fig5D saved: %d bytes\n",
            file.info(file.path(FIG, "Fig5D_MultiCox_Forest.pdf"))$size))

cat("\n=== 全部完成 ===\n")
for (f in c("Fig5A_Nomogram_PACA_AU.pdf", "Fig5B_Calibration_PACA_AU.pdf",
             "Fig5B_Calibration_External.pdf", "Fig5C_DCA_TCGA.pdf",
             "Fig5D_MultiCox_Forest.pdf")) {
  fp <- file.path(FIG, f)
  if (file.exists(fp)) {
    cat(sprintf("  ✓ %s (%d bytes)\n", f, file.info(fp)$size))
  } else {
    cat(sprintf("  ✗ %s MISSING\n", f))
  }
}
