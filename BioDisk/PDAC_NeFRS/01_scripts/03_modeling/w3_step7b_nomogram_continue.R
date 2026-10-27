# W3 Day 4 Step C-F (修复): calibration + DCA + forest
# 假设 Step B 的 Nomogram 已生成, 继续后续步骤
suppressPackageStartupMessages({
  library(data.table); library(survival); library(rms); library(dcurves); library(ggplot2)
})

OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/05_nomogram"
FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"

cat("=== W3 Day 4 Step C-F: Calibration + DCA + Forest ===\n")
clin_list <- readRDS(file.path(OUT, "clinical_with_nefrs_3cohorts.rds"))

# === 重建 cox_full (省得依赖之前 session) ===
paca <- as.data.frame(clin_list$PACA_AU)
paca_cc <- paca[complete.cases(paca[, c("NeFRS","age","Sex","Stage_AJCC","Grade",
                                          "os_time_month","os_status")]), ]
cat(sprintf("PACA-AU complete cases: %d\n", nrow(paca_cc)))

dd <- datadist(paca_cc)
options(datadist = "dd")

# 关键修复: cph 的 time.inc 应该和 calibrate 的 u 一致(都用月份)
# 但 calibrate 报错主要是 m 太小 — 让我们对每个时间点分别用合适的 m
cox_full <- cph(Surv(os_time_month, os_status) ~ NeFRS + age + Sex + Stage_AJCC + Grade,
                data = paca_cc, x=TRUE, y=TRUE, surv=TRUE,
                time.inc = 12)  # 1 年, 用月数

# ============================================================
# Step C 修复: PACA-AU 内部 Calibration
# ============================================================
cat("\n--- Step C: Internal Calibration (修复版) ---\n")

# 计算各时间点的实际可用样本量
# m 应该够大让每组都有至少 ~30 样本
n_total <- nrow(paca_cc)
# 用 4-5 个 bin, 每个 bin ~50-65 样本
m_size <- floor(n_total / 4)
cat(sprintf("Using m=%d (n_total/4) per calibration bin\n", m_size))

# 关键: 检查每个时间点有多少事件 (events>0 是 calibrate 必需的)
for (t in c(12, 36, 60)) {
  n_alive_at_t <- sum(paca_cc$os_time_month >= t & paca_cc$os_status == 0) +
                  sum(paca_cc$os_time_month >= t & paca_cc$os_status == 1)
  n_event_before <- sum(paca_cc$os_time_month < t & paca_cc$os_status == 1)
  cat(sprintf("  At t=%d: %d events before t, %d still at risk at t\n",
              t, n_event_before, n_alive_at_t))
}

# 关键修复: 重新拟合 cph, 每个时间点单独一个 cph (避免 surv.summary slot bug)
cat("\nRefitting cph with time.inc matching u...\n")

run_cal <- function(t_month) {
  # 每个时间点单独 fit cph (time.inc = u)
  cph_t <- cph(Surv(os_time_month, os_status) ~ NeFRS + age + Sex + Stage_AJCC + Grade,
               data = paca_cc, x=TRUE, y=TRUE, surv=TRUE,
               time.inc = t_month)
  cal <- tryCatch({
    calibrate(cph_t, cmethod="KM", method="boot",
              u = t_month, m = m_size, B = 200)
  }, error = function(e) {
    cat(sprintf("  Cal at t=%d FAILED: %s\n", t_month, conditionMessage(e)))
    return(NULL)
  })
  if (!is.null(cal)) {
    cat(sprintf("  Cal at t=%d month OK\n", t_month))
  }
  cal
}

cal_12 <- run_cal(12)
cal_36 <- run_cal(36)
cal_60 <- run_cal(60)

# 画图 (即使有的 cal 失败也能画其他的)
pdf(file.path(FIG, "Fig5B_Calibration_PACA_AU.pdf"), width=8, height=8, pointsize=11)
par(mar=c(4.5, 4.5, 3.5, 1.5))

valid_cals <- list()
if (!is.null(cal_12)) valid_cals[["12"]] <- cal_12
if (!is.null(cal_36)) valid_cals[["36"]] <- cal_36
if (!is.null(cal_60)) valid_cals[["60"]] <- cal_60

if (length(valid_cals) > 0) {
  # 第一条 cal 画 base
  first_t <- names(valid_cals)[1]
  cols <- c("12"="#0072B2", "36"="#009E73", "60"="#C00000")
  
  plot(valid_cals[[first_t]], xlim=c(0,1), ylim=c(0,1),
       xlab="Nomogram-predicted OS Probability",
       ylab="Observed OS Probability (KM)",
       col=cols[first_t], lwd=2.2, errbar.col=cols[first_t],
       subtitles=FALSE, main="")
  
  # 添加其他时间点
  for (t in names(valid_cals)[-1]) {
    cal_data <- valid_cals[[t]]
    # cal 对象是矩阵, 列名通常含 "mean.predicted" 和 "KM"
    cn <- colnames(cal_data)
    pred_col <- if ("mean.predicted" %in% cn) "mean.predicted" else cn[1]
    obs_col <- if ("KM" %in% cn) "KM" else if ("KM.corrected" %in% cn) "KM.corrected" else cn[2]
    lines(cal_data[, pred_col], cal_data[, obs_col],
          col=cols[t], lwd=2.2, type="b", pch=16)
  }
  
  abline(0, 1, lty=2, col="gray60")
  title("Internal Calibration on PACA-AU (Training, n=260)",
        line=1.5, cex.main=1.2)
  
  legend_labels <- sapply(names(valid_cals), function(t) {
    sprintf("%d-%s", as.numeric(t)/12, if(as.numeric(t)==12) "year" else "years")
  })
  legend("topleft", legend=legend_labels,
         col=cols[names(valid_cals)], lwd=2.2, pch=16, bty="n", cex=0.95)
} else {
  plot.new()
  text(0.5, 0.5, "Calibration failed at all time points", cex=1.2)
}
dev.off()
cat(sprintf("\nFig5B saved: %d bytes\n",
            file.info(file.path(FIG, "Fig5B_Calibration_PACA_AU.pdf"))$size))

# ============================================================
# Step D: 外部 Calibration (TCGA + E-MTAB)
# ============================================================
cat("\n--- Step D: External Calibration ---\n")

# TCGA 子模型
tcga <- as.data.frame(clin_list$TCGA_PAAD)
tcga_cc <- tcga[complete.cases(tcga[, c("NeFRS","age","Sex","Stage_AJCC",
                                          "os_time_month","os_status")]), ]
cat(sprintf("TCGA complete cases: %d\n", nrow(tcga_cc)))

dd_tcga <- datadist(tcga_cc)
options(datadist = "dd_tcga")
m_tcga <- floor(nrow(tcga_cc) / 4)

run_cal_tcga <- function(t) {
  cph_t <- cph(Surv(os_time_month, os_status) ~ NeFRS + age + Sex + Stage_AJCC,
               data = tcga_cc, x=TRUE, y=TRUE, surv=TRUE, time.inc = t)
  tryCatch(calibrate(cph_t, cmethod="KM", method="boot", u=t, m=m_tcga, B=200),
           error=function(e) {
             cat(sprintf("  TCGA cal t=%d failed: %s\n", t, conditionMessage(e))); NULL })
}
cal_tcga_12 <- run_cal_tcga(12)
cal_tcga_36 <- run_cal_tcga(36)

# E-MTAB 子模型
emtab <- as.data.frame(clin_list$E_MTAB_6134)
emtab_cc <- emtab[complete.cases(emtab[, c("NeFRS","Sex","Stage_AJCC","Grade",
                                             "os_time_month","os_status")]), ]
cat(sprintf("E-MTAB complete cases: %d\n", nrow(emtab_cc)))

dd_emtab <- datadist(emtab_cc)
options(datadist = "dd_emtab")
m_emtab <- floor(nrow(emtab_cc) / 4)

run_cal_emtab <- function(t) {
  cph_t <- cph(Surv(os_time_month, os_status) ~ NeFRS + Sex + Stage_AJCC + Grade,
               data = emtab_cc, x=TRUE, y=TRUE, surv=TRUE, time.inc = t)
  tryCatch(calibrate(cph_t, cmethod="KM", method="boot", u=t, m=m_emtab, B=200),
           error=function(e) {
             cat(sprintf("  E-MTAB cal t=%d failed: %s\n", t, conditionMessage(e))); NULL })
}
cal_emtab_12 <- run_cal_emtab(12)
cal_emtab_36 <- run_cal_emtab(36)

# 画外部 2-panel
pdf(file.path(FIG, "Fig5B_Calibration_External.pdf"), width=14, height=7, pointsize=11)
par(mfrow=c(1, 2), mar=c(4.5, 4.5, 3.5, 1.5))

# TCGA
if (!is.null(cal_tcga_12)) {
  plot(cal_tcga_12, xlim=c(0,1), ylim=c(0,1),
       xlab="Predicted OS Probability", ylab="Observed OS Probability",
       col="#0072B2", lwd=2.2, errbar.col="#0072B2",
       subtitles=FALSE, main="")
  if (!is.null(cal_tcga_36)) {
    cn <- colnames(cal_tcga_36)
    pred_col <- if ("mean.predicted" %in% cn) "mean.predicted" else cn[1]
    obs_col <- if ("KM" %in% cn) "KM" else cn[2]
    lines(cal_tcga_36[, pred_col], cal_tcga_36[, obs_col],
          col="#009E73", lwd=2.2, type="b", pch=16)
  }
  abline(0, 1, lty=2, col="gray60")
  title(sprintf("TCGA-PAAD External Calibration (n=%d)\nNeFRS + age + Sex + Stage (no Grade)",
                nrow(tcga_cc)),
        line=0.5, cex.main=1.05)
  legend("topleft", legend=c("1-year","3-year"),
         col=c("#0072B2","#009E73"), lwd=2.2, pch=16, bty="n", cex=0.95)
} else {
  plot.new(); text(0.5, 0.5, "TCGA Cal failed", cex=1.2)
}

# E-MTAB
if (!is.null(cal_emtab_12)) {
  plot(cal_emtab_12, xlim=c(0,1), ylim=c(0,1),
       xlab="Predicted OS Probability", ylab="Observed OS Probability",
       col="#0072B2", lwd=2.2, errbar.col="#0072B2",
       subtitles=FALSE, main="")
  if (!is.null(cal_emtab_36)) {
    cn <- colnames(cal_emtab_36)
    pred_col <- if ("mean.predicted" %in% cn) "mean.predicted" else cn[1]
    obs_col <- if ("KM" %in% cn) "KM" else cn[2]
    lines(cal_emtab_36[, pred_col], cal_emtab_36[, obs_col],
          col="#009E73", lwd=2.2, type="b", pch=16)
  }
  abline(0, 1, lty=2, col="gray60")
  title(sprintf("E-MTAB-6134 External Calibration (n=%d)\nNeFRS + Sex + Stage + Grade (no age)",
                nrow(emtab_cc)),
        line=0.5, cex.main=1.05)
  legend("topleft", legend=c("1-year","3-year"),
         col=c("#0072B2","#009E73"), lwd=2.2, pch=16, bty="n", cex=0.95)
} else {
  plot.new(); text(0.5, 0.5, "E-MTAB Cal failed", cex=1.2)
}
dev.off()
cat(sprintf("Fig5B_External saved: %d bytes\n",
            file.info(file.path(FIG, "Fig5B_Calibration_External.pdf"))$size))

# ============================================================
# Step E: DCA on TCGA
# ============================================================
cat("\n--- Step E: DCA ---\n")

options(datadist = "dd_tcga")
cph_tcga_36 <- cph(Surv(os_time_month, os_status) ~ NeFRS + age + Sex + Stage_AJCC,
                    data = tcga_cc, x=TRUE, y=TRUE, surv=TRUE, time.inc = 36)
# 3 年内死亡概率
tcga_cc$risk_3y <- 1 - predict(cph_tcga_36, newdata=tcga_cc, type="survival", times=36)

dca_tcga <- tryCatch({
  dca(Surv(os_time_month, os_status) ~ risk_3y,
      data = tcga_cc, time = 36,
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
  cat(sprintf("Fig5C saved: %d bytes\n",
              file.info(file.path(FIG, "Fig5C_DCA_TCGA.pdf"))$size))
}

# ============================================================
# Step F: 多变量 Cox Forest
# ============================================================
cat("\n--- Step F: Multivariate Cox Forest ---\n")

options(datadist = "dd")
cox_simple <- coxph(Surv(os_time_month, os_status) ~ NeFRS + age + Sex + Stage_AJCC + Grade,
                     data = paca_cc)
cox_sum <- summary(cox_simple)

# 重命名变量名让更可读
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
cat("\n"); print(hr_table)
fwrite(hr_table, file.path(FIG, "Fig5D_MultiCox_PACA_AU.tsv"), sep="\t")

# Forest plot
pdf(file.path(FIG, "Fig5D_MultiCox_Forest.pdf"), width=11, height=6, pointsize=11)
par(mar=c(4.5, 11, 3.5, 2))

n_var <- nrow(hr_table)
plot(NA, xlim=c(0.3, 8), ylim=c(0.5, n_var+0.5),
     log="x", xlab="Hazard Ratio (log scale)", ylab="", yaxt="n",
     cex.lab=1.05,
     main="Multivariate Cox Analysis (PACA-AU, n=260)")
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
mtext("HR (95% CI), P", side=4, line=4.5, at=n_var+0.7, cex=0.85, font=2)
dev.off()
cat(sprintf("Fig5D saved: %d bytes\n",
            file.info(file.path(FIG, "Fig5D_MultiCox_Forest.pdf"))$size))

# === Save models ===
saveRDS(list(cox_full=cox_full, hr_table=hr_table,
              cal_paca = list(c12=cal_12, c36=cal_36, c60=cal_60),
              cal_tcga = list(c12=cal_tcga_12, c36=cal_tcga_36),
              cal_emtab = list(c12=cal_emtab_12, c36=cal_emtab_36)),
        file.path(OUT, "nomogram_calibration.rds"))

cat("\n=== Day 4 Step C-F 完成 ===\n")
for (f in c("Fig5A_Nomogram_PACA_AU.pdf", "Fig5B_Calibration_PACA_AU.pdf",
             "Fig5B_Calibration_External.pdf", "Fig5C_DCA_TCGA.pdf",
             "Fig5D_MultiCox_Forest.pdf")) {
  fp <- file.path(FIG, f)
  if (file.exists(fp)) {
    cat(sprintf("  ✓ %s (%d bytes)\n", f, file.info(fp)$size))
  } else {
    cat(sprintf("  ✗ %s (NOT GENERATED)\n", f))
  }
}
