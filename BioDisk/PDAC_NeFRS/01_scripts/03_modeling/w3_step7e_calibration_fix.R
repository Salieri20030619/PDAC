# W3 Day 4 - 修复 calibration 图 (手动画, 不用 plot.calibrate)
suppressPackageStartupMessages({
  library(data.table); library(survival); library(rms)
})

OUT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/05_nomogram"
FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"

cat("=== Fix: Manual calibration plots ===\n")
clin_list <- readRDS(file.path(OUT, "clinical_with_nefrs_3cohorts.rds"))

# === PACA-AU ===
paca <- as.data.frame(clin_list$PACA_AU)
paca_cc <- paca[complete.cases(paca[, c("NeFRS","age","Sex","Stage_AJCC","Grade",
                                          "os_time_month","os_status")]), ]
dd <- datadist(paca_cc); options(datadist="dd")

run_cal <- function(df, formula, u, m_size) {
  cph_t <- cph(formula, data = df, x=TRUE, y=TRUE, surv=TRUE, time.inc = u)
  tryCatch({
    capture.output({
      result <- calibrate(cph_t, cmethod="KM", method="boot",
                          u=u, m=m_size, B=200)
    })
    result
  }, error = function(e) {
    cat(sprintf("  Cal failed t=%d: %s\n", u, conditionMessage(e))); NULL
  })
}

# 关键: 手动从 cal 矩阵提取 mean.predicted vs KM 画图
draw_cal_line <- function(cal, col, pch=16, lwd=2.2) {
  if (is.null(cal)) return(FALSE)
  if (!"mean.predicted" %in% colnames(cal)) {
    cat("  WARN: cal missing 'mean.predicted' column\n"); return(FALSE)
  }
  x <- cal[, "mean.predicted"]
  y <- cal[, "KM"]
  # 过滤 NA
  valid <- !is.na(x) & !is.na(y)
  if (sum(valid) < 2) return(FALSE)
  
  # 画点 + 连线
  lines(x[valid], y[valid], col=col, lwd=lwd, type="b", pch=pch, cex=1.3)
  
  # 误差棒 (如果有)
  if ("Lower" %in% colnames(cal) && "Upper" %in% colnames(cal)) {
    # Lower/Upper 是 KM-corrected 的偏差,要加到 mean.predicted 上
    # 不画了, 避免混乱
  }
  return(TRUE)
}

# PACA-AU calibration
cat("\n--- PACA-AU calibration ---\n")
m_paca <- 65
formula_paca <- Surv(os_time_month, os_status) ~ NeFRS + age + Sex + Stage_AJCC + Grade

cal_12 <- run_cal(paca_cc, formula_paca, 12, m_paca)
cal_36 <- run_cal(paca_cc, formula_paca, 36, m_paca)
cal_60 <- run_cal(paca_cc, formula_paca, 60, m_paca)

pdf(file.path(FIG, "Fig5B_Calibration_PACA_AU.pdf"), width=8, height=8, pointsize=11)
par(mar=c(4.5, 4.5, 3.5, 1.5))

# 空 canvas
plot(NA, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS Probability",
     ylab="Observed OS Probability (KM)",
     cex.lab=1.1, las=1,
     main="")

abline(0, 1, lty=2, col="gray60", lwd=1.5)
abline(h=seq(0,1,0.2), v=seq(0,1,0.2), col="gray92", lty=3, lwd=0.5)

cols_cal <- c("#0072B2", "#009E73", "#C00000")
labels <- c()
cols_used <- c()
if (draw_cal_line(cal_12, cols_cal[1])) { labels <- c(labels, "1-year"); cols_used <- c(cols_used, cols_cal[1]) }
if (draw_cal_line(cal_36, cols_cal[2])) { labels <- c(labels, "3-year"); cols_used <- c(cols_used, cols_cal[2]) }
if (draw_cal_line(cal_60, cols_cal[3])) { labels <- c(labels, "5-year"); cols_used <- c(cols_used, cols_cal[3]) }

title(sprintf("Internal Calibration on PACA-AU (n=%d, Bootstrap B=200)", nrow(paca_cc)),
      line=1.5, cex.main=1.15)

if (length(labels) > 0) {
  legend("topleft", legend=labels, col=cols_used, lwd=2.2, pch=16, bty="n", cex=1.0)
}
text(0.95, 0.05, "Diagonal = perfect calibration", cex=0.85, col="gray50", adj=1)
dev.off()
cat(sprintf("Fig5B PACA: %d bytes\n",
            file.info(file.path(FIG, "Fig5B_Calibration_PACA_AU.pdf"))$size))

# === TCGA ===
cat("\n--- TCGA calibration ---\n")
tcga <- as.data.frame(clin_list$TCGA_PAAD)
tcga_cc <- tcga[complete.cases(tcga[, c("NeFRS","age","Sex","Stage_AJCC",
                                          "os_time_month","os_status")]),
                 c("sample_id","os_time_month","os_status","NeFRS","age","Sex","Stage_AJCC")]
tcga_cc$Stage_AJCC <- droplevels(tcga_cc$Stage_AJCC)
dd_tcga <- datadist(tcga_cc); options(datadist="dd_tcga")

m_tcga <- 36
formula_tcga <- Surv(os_time_month, os_status) ~ NeFRS + age + Sex + Stage_AJCC

cal_tcga_12 <- run_cal(tcga_cc, formula_tcga, 12, m_tcga)
cal_tcga_36 <- run_cal(tcga_cc, formula_tcga, 36, m_tcga)

# === E-MTAB ===
cat("\n--- E-MTAB calibration ---\n")
emtab <- as.data.frame(clin_list$E_MTAB_6134)
emtab_cc <- emtab[complete.cases(emtab[, c("NeFRS","Sex","Stage_AJCC","Grade",
                                             "os_time_month","os_status")]),
                   c("sample_id","os_time_month","os_status","NeFRS","Sex","Stage_AJCC","Grade")]
emtab_cc$Stage_AJCC <- droplevels(emtab_cc$Stage_AJCC)
emtab_cc$Grade <- droplevels(emtab_cc$Grade)
dd_emtab <- datadist(emtab_cc); options(datadist="dd_emtab")

m_emtab <- 70
formula_emtab <- Surv(os_time_month, os_status) ~ NeFRS + Sex + Stage_AJCC + Grade

cal_emtab_12 <- run_cal(emtab_cc, formula_emtab, 12, m_emtab)
cal_emtab_36 <- run_cal(emtab_cc, formula_emtab, 36, m_emtab)

# === 外部 calibration 2-panel PDF ===
pdf(file.path(FIG, "Fig5B_Calibration_External.pdf"), width=14, height=7, pointsize=11)
par(mfrow=c(1, 2), mar=c(4.5, 4.5, 3.5, 1.5))

# TCGA panel
plot(NA, xlim=c(0,1), ylim=c(0,1),
     xlab="Predicted OS Probability", ylab="Observed OS Probability (KM)",
     cex.lab=1.05, las=1, main="")
abline(0, 1, lty=2, col="gray60", lwd=1.5)
abline(h=seq(0,1,0.2), v=seq(0,1,0.2), col="gray92", lty=3, lwd=0.5)
tcga_labels <- c(); tcga_cols <- c()
if (draw_cal_line(cal_tcga_12, "#0072B2")) {tcga_labels <- c(tcga_labels, "1-year"); tcga_cols <- c(tcga_cols, "#0072B2")}
if (draw_cal_line(cal_tcga_36, "#009E73")) {tcga_labels <- c(tcga_labels, "3-year"); tcga_cols <- c(tcga_cols, "#009E73")}
title(sprintf("TCGA-PAAD External (n=%d)\nNeFRS + age + Sex + Stage", nrow(tcga_cc)),
      line=0.5, cex.main=1.05)
if (length(tcga_labels) > 0) {
  legend("topleft", legend=tcga_labels, col=tcga_cols, lwd=2.2, pch=16, bty="n", cex=0.95)
}

# E-MTAB panel
plot(NA, xlim=c(0,1), ylim=c(0,1),
     xlab="Predicted OS Probability", ylab="Observed OS Probability (KM)",
     cex.lab=1.05, las=1, main="")
abline(0, 1, lty=2, col="gray60", lwd=1.5)
abline(h=seq(0,1,0.2), v=seq(0,1,0.2), col="gray92", lty=3, lwd=0.5)
emtab_labels <- c(); emtab_cols <- c()
if (draw_cal_line(cal_emtab_12, "#0072B2")) {emtab_labels <- c(emtab_labels, "1-year"); emtab_cols <- c(emtab_cols, "#0072B2")}
if (draw_cal_line(cal_emtab_36, "#009E73")) {emtab_labels <- c(emtab_labels, "3-year"); emtab_cols <- c(emtab_cols, "#009E73")}
title(sprintf("E-MTAB-6134 External (n=%d)\nNeFRS + Sex + Stage + Grade", nrow(emtab_cc)),
      line=0.5, cex.main=1.05)
if (length(emtab_labels) > 0) {
  legend("topleft", legend=emtab_labels, col=emtab_cols, lwd=2.2, pch=16, bty="n", cex=0.95)
}

dev.off()
cat(sprintf("Fig5B External: %d bytes\n",
            file.info(file.path(FIG, "Fig5B_Calibration_External.pdf"))$size))

cat("\n=== Done ===\n")
