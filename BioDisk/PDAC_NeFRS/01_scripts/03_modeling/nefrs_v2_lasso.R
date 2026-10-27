suppressPackageStartupMessages({
  library(data.table); library(survival); library(glmnet)
})

W2_IN <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/02_w2_inputs"
MIME  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime"
GP_IN <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/01_gene_pool"

cohort_data <- readRDS(file.path(W2_IN, "cohort_data_8.rds"))
mime_input  <- readRDS(file.path(MIME, "mime_input_list.rds"))
mime_genes  <- colnames(mime_input[[1]])[-(1:3)]
gp <- fread(file.path(GP_IN, "final_in_all_cohorts.tsv"))

# Train Lasso on PACA-AU with 132 genes
train_df <- mime_input[["PACA_AU"]]
X_train <- as.matrix(train_df[, mime_genes])
y_train <- Surv(train_df$OS.time, train_df$OS)

cat("=== Training Lasso on PACA-AU ===\n")
set.seed(5201314)
cv_lasso <- cv.glmnet(X_train, y_train, family="cox", alpha=1, nfolds=10)
cat(sprintf("Lasso lambda.min: %.5f\n", cv_lasso$lambda.min))
cat(sprintf("Lasso lambda.1se: %.5f\n", cv_lasso$lambda.1se))

# Show non-zero coefs at lambda.min and lambda.1se
for (s_choice in c("lambda.min", "lambda.1se")) {
  coefs <- as.numeric(coef(cv_lasso, s=s_choice))
  names(coefs) <- mime_genes
  nz <- coefs[abs(coefs) > 1e-8]
  cat(sprintf("\n=== %s: %d non-zero genes ===\n", s_choice, length(nz)))
  nz_df <- data.table(gene=names(nz), coef=nz, abs_coef=abs(nz))
  nz_df <- merge(nz_df, gp[, .(GeneSymbol, priority, n_sources)],
                 by.x="gene", by.y="GeneSymbol", all.x=TRUE)
  setorder(nz_df, -abs_coef)
  print(nz_df)
}

# Evaluate on 8 cohorts
cat("\n=== 8-cohort evaluation (using lambda.min) ===\n")
holdout <- c("GSE21501","TCGA_PAAD","E_MTAB_6134")
selval  <- c("GSE71729","GSE85916","GSE62452","GSE57495")
all_cohorts <- c("PACA_AU", selval, holdout)

results <- data.frame()
for (cn in all_cohorts) {
  co <- cohort_data[[cn]]
  X <- matrix(0, nrow=ncol(co$expression), ncol=length(mime_genes),
              dimnames=list(colnames(co$expression), mime_genes))
  avail <- intersect(mime_genes, rownames(co$expression))
  for (g in avail) X[, g] <- co$expression[g, ]
  
  rs <- as.numeric(predict(cv_lasso, newx=X, s="lambda.min"))
  names(rs) <- rownames(X)
  
  clin <- co$clinical
  m <- match(clin$sample_id, names(rs))
  rs_v <- rs[m]
  
  surv <- Surv(clin$os_time_month, clin$os_status)
  ci <- concordance(surv ~ rs_v, reverse=TRUE)
  med <- median(rs_v, na.rm=TRUE)
  grp <- factor(ifelse(rs_v > med, "High","Low"), levels=c("Low","High"))
  sd_obj <- survdiff(surv ~ grp)
  km_p <- 1 - pchisq(sd_obj$chisq, df=1)
  cox <- coxph(surv ~ grp)
  hr <- exp(coef(cox))
  hr_ci <- exp(confint(cox))
  
  role <- ifelse(cn=="PACA_AU","TRAIN",
                 ifelse(cn %in% holdout,"HOLDOUT","SEL-VAL"))
  cat(sprintf("  %-12s [%s] n=%3d C=%.3f KM_p=%.3g HR=%.2f (%.2f-%.2f)\n",
              cn, role, nrow(co$clinical), as.numeric(ci$concordance),
              km_p, hr, hr_ci[1], hr_ci[2]))
  results <- rbind(results, data.frame(
    cohort=cn, role=role, n=nrow(co$clinical),
    c_index=round(as.numeric(ci$concordance),3),
    km_p=signif(km_p,3),
    HR=round(hr,2),
    HR_CI=sprintf("%.2f-%.2f", hr_ci[1], hr_ci[2])
  ))
}

cat("\n=== Summary ===\n")
print(results)
cat(sprintf("\nHold-out mean C-index: %.3f (vs StepCox+Ridge 0.626)\n",
            mean(results$c_index[results$role=="HOLDOUT"])))
cat(sprintf("Sel-val mean C-index:  %.3f (vs StepCox+Ridge 0.621)\n",
            mean(results$c_index[results$role=="SEL-VAL"])))

# Save
coefs_min <- as.numeric(coef(cv_lasso, s="lambda.min"))
names(coefs_min) <- mime_genes
nz_min <- coefs_min[abs(coefs_min) > 1e-8]
nz_min_df <- data.table(gene=names(nz_min), coef=nz_min, abs_coef=abs(nz_min))
nz_min_df <- merge(nz_min_df, gp[, .(GeneSymbol, priority, Sources, n_sources)],
                    by.x="gene", by.y="GeneSymbol", all.x=TRUE)
setorder(nz_min_df, -abs_coef)
nz_min_df[, rank := .I]

fwrite(nz_min_df, file.path(MIME, "nefrs_v2_lasso_coefficients.tsv"), sep="\t")
fwrite(results, file.path(MIME, "nefrs_v2_lasso_8cohort.tsv"), sep="\t")
saveRDS(cv_lasso, file.path(MIME, "nefrs_v2_lasso_model.rds"))
cat("\nSaved: nefrs_v2_lasso_coefficients.tsv, nefrs_v2_lasso_8cohort.tsv, nefrs_v2_lasso_model.rds\n")
