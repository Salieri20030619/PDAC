# W3 Day 3 Step 3: 1/3/5 year time-dependent AUC (fixed rbind)
suppressPackageStartupMessages({
  library(data.table); library(survival); library(glmnet); library(timeROC)
})

W2_IN <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/02_w2_inputs"
MIME  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime"
FIG   <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures"

cat("=== W3 Step 3: Time-dependent AUC ===\n")

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
  list(name="PACA_AU",     label="ICGC PACA-AU (Training, n=267)",       role="TRAIN"),
  list(name="GSE71729",    label="GSE71729 (Sel-Val, n=123)",             role="SEL"),
  list(name="GSE85916",    label="GSE85916 (Sel-Val, n=79)",              role="SEL"),
  list(name="GSE62452",    label="GSE62452 (Sel-Val, n=66)",              role="SEL"),
  list(name="GSE57495",    label="GSE57495 (Sel-Val, n=63)",              role="SEL"),
  list(name="GSE21501",    label="GSE21501 (Hold-out, n=102)",            role="HO"),
  list(name="TCGA_PAAD",   label="TCGA-PAAD (Hold-out, n=147)",           role="HO"),
  list(name="E_MTAB_6134", label="E-MTAB-6134 (Hold-out, n=288)",         role="HO")
)

# 计算 timeROC at 12/36/60 months
compute_roc <- function(co_info) {
  rs_df <- compute_rs(co_info$name)
  rs_df <- rs_df[!is.na(rs_df$RS) & !is.na(rs_df$OS_time) & rs_df$OS_time > 0, ]
  max_t <- max(rs_df$OS_time, na.rm=TRUE)
  times <- c(12, 36, 60)
  valid_times <- times[times <= max_t]
  if (length(valid_times) == 0) return(NULL)
  
  roc <- tryCatch({
    timeROC::timeROC(
      T = rs_df$OS_time, delta = rs_df$OS_status,
      marker = rs_df$RS, cause = 1, weighting = "marginal",
      times = valid_times, iid = FALSE
    )
  }, error = function(e) {
    cat(sprintf("  Warning %s: %s\n", co_info$name, conditionMessage(e)))
    return(NULL)
  })
  
  list(roc=roc, times=valid_times, n=nrow(rs_df), max_t=max_t)
}

plot_roc <- function(roc_obj, co_info) {
  if (is.null(roc_obj) || is.null(roc_obj$roc)) {
    plot.new(); text(0.5, 0.5, sprintf("%s\n(Insufficient data)", co_info$name))
    return(invisible(NULL))
  }
  roc <- roc_obj$roc; times <- roc_obj$times
  cols <- c("#0072B2", "#009E73", "#C00000")
  
  plot(NA, xlim=c(0,1), ylim=c(0,1),
       xlab="1 - Specificity", ylab="Sensitivity",
       main=co_info$label,
       cex.main=1.05, cex.lab=1.1, cex.axis=0.95, las=1)
  abline(0, 1, lty=2, col="gray60", lwd=1)
  abline(h=seq(0,1,0.2), v=seq(0,1,0.2), col="gray92", lty=3, lwd=0.5)
  
  legend_labels <- c()
  for (j in seq_along(times)) {
    t <- times[j]; fpr <- roc$FP[, j]; tpr <- roc$TP[, j]; auc_val <- roc$AUC[j]
    lines(fpr, tpr, col=cols[j], lwd=2.8)
    yr_label <- ifelse(t == 12, "1-year",
                ifelse(t == 36, "3-year",
                ifelse(t == 60, "5-year", sprintf("%d-month", t))))
    legend_labels <- c(legend_labels, sprintf("%s AUC = %.3f", yr_label, auc_val))
  }
  legend("bottomright", legend=legend_labels, col=cols[seq_along(times)],
         lwd=2.8, bty="n", cex=0.92)
}

# === 计算 ===
cat("\nComputing time-dependent ROC for 8 cohorts...\n")
all_roc <- list()
# 修复: 使用 list 收集每个 cohort 的 row,最后用 rbindlist(fill=TRUE) 合并
row_list <- list()

for (i in seq_along(cohorts)) {
  ci <- cohorts[[i]]
  cat(sprintf("  [%d/8] %s ...\n", i, ci$name))
  all_roc[[i]] <- compute_roc(ci)
  
  # 始终生成包含 AUC_1yr/3yr/5yr 三列的 row, 缺失填 NA
  row <- list(cohort=ci$name, role=ci$role, n=NA_integer_,
              AUC_1yr=NA_real_, AUC_3yr=NA_real_, AUC_5yr=NA_real_)
  
  if (!is.null(all_roc[[i]]) && !is.null(all_roc[[i]]$roc)) {
    row$n <- all_roc[[i]]$n
    times <- all_roc[[i]]$times
    aucs <- all_roc[[i]]$roc$AUC
    for (j in seq_along(times)) {
      if (times[j] == 12) row$AUC_1yr <- round(aucs[j], 3)
      else if (times[j] == 36) row$AUC_3yr <- round(aucs[j], 3)
      else if (times[j] == 60) row$AUC_5yr <- round(aucs[j], 3)
    }
  }
  row_list[[i]] <- row
}

# 用 rbindlist 合并 (fill=TRUE 在 data.table 里是合法的)
summary_df <- rbindlist(row_list, fill=TRUE)

cat("\n=== Time-dependent AUC ===\n")
print(summary_df)
fwrite(summary_df, file.path(FIG, "Fig3_AUC_summary.tsv"), sep="\t")

# === 单 cohort PDF ===
cat("\nGenerating individual ROC PDFs...\n")
for (i in seq_along(cohorts)) {
  outpath <- file.path(FIG, sprintf("Fig3_ROC_%02d_%s.pdf", i, cohorts[[i]]$name))
  pdf(outpath, width=6, height=6, pointsize=10)
  par(mar=c(4.5, 4.5, 3.5, 1))
  plot_roc(all_roc[[i]], cohorts[[i]])
  dev.off()
  cat(sprintf("  saved %-35s %d bytes\n",
              basename(outpath), file.info(outpath)$size))
}

# === 合并 8 panel ===
cat("\nGenerating combined 8-panel ROC PDF...\n")
pdf(file.path(FIG, "Fig3_ROC_8cohort_combined.pdf"),
    width=16, height=20, pointsize=10)
par(mfrow=c(4, 2), mar=c(4.5, 4.5, 3.5, 1), oma=c(0, 0, 3, 0))
for (i in seq_along(cohorts)) plot_roc(all_roc[[i]], cohorts[[i]])
mtext("NeFRS v2 Time-dependent ROC (1/3/5 years) across 8 PDAC Cohorts",
      side=3, line=0.5, outer=TRUE, cex=1.4, font=2)
dev.off()
cat(sprintf("  Combined PDF: %d bytes\n",
            file.info(file.path(FIG, "Fig3_ROC_8cohort_combined.pdf"))$size))

cat("\nDone.\n")
