# W3 Day 4 - 修复 Step C/D/E/F
suppressPackageStartupMessages({
  library(data.table); library(survival); library(rms); library(dcurves); library(ggplot2)
})

OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/05_nomogram"
FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"

cat("=== Fix: Step C/D/E/F ===\n")
clin_list <- readRDS(file.path(OUT, "clinical_with_nefrs_3cohorts.rds"))

# ============================================================
# PACA-AU 主模型
# ============================================================
paca <- as.data.frame(clin_list$PACA_AU)
paca_cc <- paca[complete.cases(paca[, c("NeFRS","age","Sex","Stage_AJCC","Grade",
                                          "os_time_month","os_status")]), ]
cat(sprintf("PACA-AU complete: %d\n", nrow(paca_cc)))

dd <- datadist(paca_cc)
options(datadist = "dd")

# ============================================================
# Step C: PACA-AU 内部 Calibration (修复 m + time.inc)
# ============================================================
cat("\n--- Step C: PACA-AU Internal Calibration ---\n")

# 关键: 每个时间点独立 cph + 大 m
make_cal <- function(df, u_month, label) {
  m_use <- floor(nrow(df) / 4)
  cph_t <- cph(Surv(os_time_month, os_status) ~ NeFRS + age + Sex + Stage_AJCC + Grade,
               data = df, x=TRUE, y=TRUE, surv=TRUE, time.inc = u_month)
  cal <- tryCatch({
    capture.output({
      result <- calibrate(cph_t, cmethod="KM", method="boot",
                          u = u_month, m = m_use, B = 200)
    })
    result
  }, error = function(e) {
    cat(sprintf("  %s t=%d FAILED: %s\n", label, u_month, conditionMessage(e)))
    NULL
  })
  if (!is.null(cal)) cat(sprintf("  %s t=%d OK (m=%d)\n", label, u_month, m_use))
  cal
}

cal_12 <- make_cal(paca_cc, 12, "PACA")
cal_36 <- make_cal(paca_cc, 36, "PACA")
cal_60 <- make_cal(paca_cc, 60, "PACA")

# 画 PACA cal
pdf(file.path(FIG, "Fig5B_Calibration_PACA_AU.pdf"), width=8, height=8, pointsize=11)
par(mar=c(4.5, 4.5, 3.5, 1.5))

cols_cal <- c("12"="#0072B2", "36"="#009E73", "60"="#C00000")
valid_cals <- list()
if (!is.null(cal_12)) valid_cals[["12"]] <- cal_12
if (!is.null(cal_36)) valid_cals[["36"]] <- cal_36
if (!is.null(cal_60)) valid_cals[["60"]] <- cal_60

if (length(valid_cals) > 0) {
  first_t <- names(valid_cals)[1]
  plot(valid_cals[[first_t]], xlim=c(0,1), ylim=c(0,1),
       xlab="Predicted OS Probability", ylab="Observed OS Probability (KM)",
       col=cols_cal[first_t], lwd=2.2, errbar.col=cols_cal[first_t],
       subtitles=FALSE, main="")
  
  for (t in names(valid_cals)[-1]) {
    cal_data <- valid_cals[[t]]
    cn <- colnames(cal_data)
    pred_col <- if ("mean.predicted" %in% cn) "mean.predicted" else cn[1]
    obs_col <- if ("KM" %in% cn) "KM" else cn[2]
    lines(cal_data[, pred_col], cal_data[, obs_col],
          col=cols_cal[t], lwd=2.2, type="b", pch=16)
  }
  abline(0, 1, lty=2, col="gray60")
  title(sprintf("Internal Calibration on PACA-AU (n=%d)", nrow(paca_cc)),
        line=1.5, cex.main=1.2)
  
  legend_labels <- sapply(names(valid_cals), function(t) {
    sprintf("%d-year", as.numeric(t)/12)
  })
  legend("topleft", legend=legend_labels,
         col=cols_cal[names(valid_cals)], lwd=2.2, pch=16, bty="n", cex=0.95)
} else {
  plot.new(); text(0.5, 0.5, "All cal failed", cex=1.2)
}
dev.off()
cat(sprintf("Fig5B PACA: %d bytes\n", file.info(file.path(FIG, "Fig5B_Calibration_PACA_AU.pdf"))$size))

# ============================================================
# Step D: TCGA + E-MTAB 外部 Calibration
# ============================================================
cat("\n--- Step D: External Calibration ---\n")

# TCGA (有 age, 无 Grade)
tcga <- as.data.frame(clin_list$TCGA_PAAD)
tcga_cc <- tcga[complete.cases(tcga[, c("NeFRS","age","Sex","Stage_AJCC",
                                          "os_time_month","os_status")]),
                 c("sample_id","os_time_month","os_status","NeFRS","age","Sex","Stage_AJCC")]
cat(sprintf("TCGA complete: %d\n", nrow(tcga_cc)))
dd_tcga <- datadist(tcga_cc)
options(datadist = "dd_tcga")

make_cal_tcga <- function(u_month) {
  m_use <- floor(nrow(tcga_cc) / 4)
  cph_t <- cph(Surv(os_time_month, os_status) ~ NeFRS + age + Sex + Stage_AJCC,
               data = tcga_cc, x=TRUE, y=TRUE, surv=TRUE, time.inc = u_month)
  tryCatch({
    capture.output({
      result <- calibrate(cph_t, cmethod="KM", method="boot",
                          u = u_month, m = m_use, B = 200)
    })
    cat(sprintf("  TCGA t=%d OK\n", u_month))
    result
  }, error = function(e) {
    cat(sprintf("  TCGA t=%d FAILED: %s\n", u_month, conditionMessage(e)))
    NULL
  })
}
cal_tcga_12 <- make_cal_tcga(12)
cal_tcga_36 <- make_cal_tcga(36)

# E-MTAB (无 age, 有 Grade) — 关键修复: 删除 age 列!
emtab <- as.data.frame(clin_list$E_MTAB_6134)
# 用 emtab 子模型变量, 不含 age
emtab_cc <- emtab[complete.cases(emtab[, c("NeFRS","Sex","Stage_AJCC","Grade",
                                             "os_time_month","os_status")]),
                   c("sample_id","os_time_month","os_status","NeFRS","Sex","Stage_AJCC","Grade")]
# 注意上面 c(...) 不含 age, 所以 emtab_cc 没有 age 列, datadist 不会崩
cat(sprintf("E-MTAB complete: %d\n", nrow(emtab_cc)))
cat(sprintf("  Columns: %s\n", paste(colnames(emtab_cc), collapse=", ")))

dd_emtab <- datadist(emtab_cc)
options(datadist = "dd_emtab")

make_cal_emtab <- function(u_month) {
  m_use <- floor(nrow(emtab_cc) / 4)
  cph_t <- cph(Surv(os_time_month, os_status) ~ NeFRS + Sex + Stage_AJCC + Grade,
               data = emtab_cc, x=TRUE, y=TRUE, surv=TRUE, time.inc = u_month)
  tryCatch({
    capture.output({
      result <- calibrate(cph_t, cmethod="KM", method="boot",
                          u = u_month, m = m_use, B = 200)
    })
    cat(sprintf("  E-MTAB t=%d OK\n", u_month))
    result
  }, error = function(e) {
    cat(sprintf("  E-MTAB t=%d FAILED: %s\n", u_month, conditionMessage(e)))
    NULL
  })
}
cal_emtab_12 <- make_cal_emtab(12)
cal_emtab_36 <- make_cal_emtab(36)

# 画外部 2-panel
plot_external_cal <- function(c12, c36, title_str) {
  if (is.null(c12) && is.null(c36)) {
    plot.new(); text(0.5, 0.5, "Cal failed", cex=1.2); return()
  }
  first_cal <- if (!is.null(c12)) c12 else c36
  plot(first_cal, xlim=c(0,1), ylim=c(0,1),
       xlab="Predicted OS Probability", ylab="Observed OS Probability",
       col="#0072B2", lwd=2.2, errbar.col="#0072B2",
       subtitles=FALSE, main="")
  if (!is.null(c12) && !is.null(c36)) {
    cn <- colnames(c36)
    pred_col <- if ("mean.predicted" %in% cn) "mean.predicted" else cn[1]
    obs_col <- if ("KM" %in% cn) "KM" else cn[2]
    lines(c36[, pred_col], c36[, obs_col], col="#009E73", lwd=2.2, type="b", pch=16)
  }
  abline(0, 1, lty=2, col="gray60")
  title(title_str, line=0.5, cex.main=1.05)
  
  leg_labs <- c()
  leg_cols <- c()
  if (!is.null(c12)) { leg_labs <- c(leg_labs, "1-year"); leg_cols <- c(leg_cols, "#0072B2") }
  if (!is.null(c36)) { leg_labs <- c(leg_labs, "3-year"); leg_cols <- c(leg_cols, "#009E73") }
  legend("topleft", legend=leg_labs, col=leg_cols, lwd=2.2, pch=16, bty="n", cex=0.95)
}

pdf(file.path(FIG, "Fig5B_Calibration_External.pdf"), width=14, height=7, pointsize=11)
par(mfrow=c(1, 2), mar=c(4.5, 4.5, 3.5, 1.5))
plot_external_cal(cal_tcga_12, cal_tcga_36,
                   sprintf("TCGA-PAAD External Cal (n=%d)\nNeFRS + age + Sex + Stage", nrow(tcga_cc)))
plot_external_cal(cal_emtab_12, cal_emtab_36,
                   sprintf("E-MTAB-6134 External Cal (n=%d)\nNeFRS + Sex + Stage + Grade", nrow(emtab_cc)))
dev.off()
cat(sprintf("Fig5B External: %d bytes\n", file.info(file.path(FIG, "Fig5B_Calibration_External.pdf"))$size))

# ============================================================
# Step E: DCA on TCGA
# ============================================================
cat("\n--- Step E: DCA ---\n")
options(datadist = "dd_tcga")
cph_tcga_36 <- cph(Surv(os_time_month, os_status) ~ NeFRS + age + Sex + Stage_AJCC,
                    data = tcga_cc, x=TRUE, y=TRUE, surv=TRUE, time.inc = 36)
tcga_cc$risk_3y <- 1 - predict(cph_tcga_36, newdata=tcga_cc, type="survival", times=36)

dca_tcga <- tryCatch({
  dca(Surv(os_time_month, os_status) ~ risk_3y, data = tcga_cc, time = 36,
      label = list(risk_3y = "NeFRS + Clinical"))
}, error = function(e) {
  cat(sprintf("DCA error: %s\n", conditionMessage(e))); NULL
})

if (!is.null(dca_tcga)) {
  pdf(file.path(FIG, "Fig5C_DCA_TCGA.pdf"), width=8, height=6.5, pointsize=11)
  p_dca <- plot(dca_tcga, smooth=TRUE) +
    ggplot2::ggtitle("Decision Curve Analysis (TCGA-PAAD, 3-year mortality)") +
    ggplot2::theme(plot.title = ggplot2::element_text(size=12, face="bold"))
  print(p_dca)
  dev.off()
  cat(sprintf("Fig5C: %d bytes\n", file.info(file.path(FIG, "Fig5C_DCA_TCGA.pdf"))$size))
}

# ============================================================
# Step F: 多变量 Cox Forest
# ============================================================
cat("\n--- Step F: Multivariate Cox Forest ---\n")
options(datadist = "dd")
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

pdf(file.path(FIG, "Fig5D_MultiCox_Forest.pdf"), width=11, height=6, pointsize=11)
par(mar=c(4.5, 11, 3.5, 12))
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
cat(sprintf("Fig5D: %d bytes\n", file.info(file.path(FIG, "Fig5D_MultiCox_Forest.pdf"))$size))

cat("\n=== 全部完成 ===\n")
for (f in c("Fig5A_Nomogram_PACA_AU.pdf", "Fig5B_Calibration_PACA_AU.pdf",
             "Fig5B_Calibration_External.pdf", "Fig5C_DCA_TCGA.pdf",
             "Fig5D_MultiCox_Forest.pdf")) {
  fp <- file.path(FIG, f)
  if (file.exists(fp)) {
    cat(sprintf("  ✓ %s (%d bytes)\n", f, file.info(fp)$size))
  } else {
    cat(sprintf("  ✗ %s NOT GENERATED\n", f))
  }
}
