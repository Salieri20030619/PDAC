suppressPackageStartupMessages({
  library(data.table); library(survival); library(glmnet); library(superpc)
})

W2_IN <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/02_w2_inputs"
MIME  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime"

res         <- readRDS(file.path(MIME, "mime_results.rds"))
cohort_data <- readRDS(file.path(W2_IN, "cohort_data_8.rds"))
mime_input  <- readRDS(file.path(MIME, "mime_input_list.rds"))
mime_genes  <- colnames(mime_input[[1]])[-(1:3)]
gp <- fread("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/01_gene_pool/final_in_all_cohorts.tsv")

cat("============================================================\n")
cat("Part 1: NeFRS v1 main formula extraction\n")
cat("============================================================\n\n")

best_alg <- "StepCox[forward] + Ridge"
ml <- res$ml.res[[best_alg]]
cat(sprintf("Algorithm: %s\n", best_alg))
cat(sprintf("Lambda.min: %.4f\n", ml$lambda.min))
cat(sprintf("Total genes: %d\n\n", length(mime_genes)))

coefs <- as.numeric(coef(ml, s="lambda.min"))
names(coefs) <- mime_genes

coef_df <- data.table(gene=mime_genes, coef=coefs, abs_coef=abs(coefs))
coef_df <- merge(coef_df, gp[, .(GeneSymbol, priority, Sources, n_sources)],
                 by.x="gene", by.y="GeneSymbol", all.x=TRUE)
setorder(coef_df, -abs_coef)
coef_df[, rank := .I]

cat("=== Top 30 NeFRS contributors (by |coef|) ===\n")
print(coef_df[1:30, .(rank, gene, coef=round(coef,4), priority, n_sources)])

cat("\n=== P1 genes (3+ source priority) in NeFRS ===\n")
p1_in <- coef_df[priority == "P1_TopPriority"][order(rank)]
print(p1_in[, .(rank, gene, coef=round(coef,4), Sources)])

cat(sprintf("\nRisk-up genes (coef>0): %d (mean=%.3f)\n",
            sum(coefs>0), mean(coefs[coefs>0])))
cat(sprintf("Protective (coef<0):    %d (mean=%.3f)\n",
            sum(coefs<0), mean(coefs[coefs<0])))

fwrite(coef_df, file.path(MIME, "nefrs_v1_main_coefficients.tsv"), sep="\t")
cat("\nSaved: nefrs_v1_main_coefficients.tsv\n\n")

cat("============================================================\n")
cat("Part 2: 10-fold CV on PACA-AU\n")
cat("============================================================\n\n")

train_df <- mime_input[["PACA_AU"]]
X_full <- as.matrix(train_df[, mime_genes])
y_time <- train_df$OS.time
y_stat <- train_df$OS
N <- nrow(train_df)

set.seed(5201314)
folds <- sample(rep(1:10, length.out = N))

run_cv_glmnet <- function(alpha_v) {
  fold_cis <- numeric(10)
  for (f in 1:10) {
    tr <- folds != f; te <- folds == f
    Xtr <- X_full[tr, , drop=FALSE]; Xte <- X_full[te, , drop=FALSE]
    ttr <- y_time[tr]; str_ <- y_stat[tr]
    tte <- y_time[te]; ste <- y_stat[te]
    fold_cis[f] <- tryCatch({
      cv <- cv.glmnet(Xtr, Surv(ttr, str_), family="cox", alpha=alpha_v, nfolds=5)
      beta <- as.numeric(coef(cv, s="lambda.min"))
      lp <- as.numeric(Xte %*% beta)
      if (length(unique(lp)) < 2) return(NA)
      as.numeric(concordance(Surv(tte, ste) ~ lp, reverse=TRUE)$concordance)
    }, error=function(e) NA)
  }
  fold_cis
}

run_cv_stepcox_ridge <- function() {
  fold_cis <- numeric(10)
  for (f in 1:10) {
    tr <- folds != f; te <- folds == f
    Xtr <- X_full[tr, , drop=FALSE]; Xte <- X_full[te, , drop=FALSE]
    ttr <- y_time[tr]; str_ <- y_stat[tr]
    tte <- y_time[te]; ste <- y_stat[te]
    fold_cis[f] <- tryCatch({
      df <- data.frame(time=ttr, status=str_, Xtr, check.names=FALSE)
      null_m <- coxph(Surv(time, status) ~ 1, data=df)
      full_m <- coxph(Surv(time, status) ~ ., data=df)
      step_m <- step(null_m, scope=list(lower=null_m, upper=full_m),
                     direction="forward", trace=0, steps=20)
      sel <- names(coef(step_m))
      sel <- gsub("`", "", sel)
      if (length(sel) < 2) return(NA)
      Xtr_s <- Xtr[, sel, drop=FALSE]
      Xte_s <- Xte[, sel, drop=FALSE]
      cv <- cv.glmnet(Xtr_s, Surv(ttr, str_), family="cox", alpha=0, nfolds=5)
      beta <- as.numeric(coef(cv, s="lambda.min"))
      lp <- as.numeric(Xte_s %*% beta)
      if (length(unique(lp)) < 2) return(NA)
      as.numeric(concordance(Surv(tte, ste) ~ lp, reverse=TRUE)$concordance)
    }, error=function(e) NA)
  }
  fold_cis
}

run_cv_superpc <- function() {
  fold_cis <- numeric(10)
  for (f in 1:10) {
    tr <- folds != f; te <- folds == f
    Xtr <- X_full[tr, , drop=FALSE]; Xte <- X_full[te, , drop=FALSE]
    fold_cis[f] <- tryCatch({
      data_tr <- list(x=t(Xtr), y=y_time[tr], censoring.status=y_stat[tr],
                      featurenames=colnames(Xtr))
      data_te <- list(x=t(Xte), y=y_time[te], censoring.status=y_stat[te],
                      featurenames=colnames(Xte))
      fit <- superpc::superpc.train(data_tr, type="survival")
      pred <- superpc::superpc.predict(object=fit, data=data_tr, newdata=data_te,
                                       threshold=0.5, n.components=1,
                                       prediction.type="continuous")
      lp <- as.numeric(if(is.matrix(pred$v.pred)) pred$v.pred[,1] else pred$v.pred)
      if (length(unique(lp)) < 2) return(NA)
      as.numeric(concordance(Surv(y_time[te], y_stat[te]) ~ lp, reverse=TRUE)$concordance)
    }, error=function(e) NA)
  }
  fold_cis
}

cv_results <- list()
cat(">>> Lasso CV...\n")
cv_results[["Lasso"]] <- run_cv_glmnet(1)
cat(">>> Ridge CV...\n")
cv_results[["Ridge"]] <- run_cv_glmnet(0)
cat(">>> Enet[alpha=0.1] CV...\n")
cv_results[["Enet[alpha=0.1]"]] <- run_cv_glmnet(0.1)
cat(">>> StepCox+Ridge CV (slowest)...\n")
cv_results[["StepCox+Ridge"]] <- run_cv_stepcox_ridge()
cat(">>> SuperPC CV...\n")
cv_results[["SuperPC"]] <- run_cv_superpc()

cv_summary <- data.frame(
  algorithm = names(cv_results),
  cv_mean = sapply(cv_results, function(x) round(mean(x, na.rm=TRUE), 3)),
  cv_sd   = sapply(cv_results, function(x) round(sd(x, na.rm=TRUE), 3)),
  cv_min  = sapply(cv_results, function(x) round(min(x, na.rm=TRUE), 3)),
  cv_max  = sapply(cv_results, function(x) round(max(x, na.rm=TRUE), 3)),
  n_valid = sapply(cv_results, function(x) sum(!is.na(x)))
)
cat("\n=== CV Summary ===\n")
print(cv_summary)

cat("\n=== Per-fold C-index ===\n")
cv_detail <- as.data.frame(do.call(rbind, lapply(cv_results, function(x) round(x, 3))))
colnames(cv_detail) <- paste0("f", 1:10)
cv_detail$algorithm <- names(cv_results)
print(cv_detail[, c("algorithm", paste0("f", 1:10))])

# Composite score
holdout_means <- c("Lasso"=NA, "Ridge"=NA, "Enet[alpha=0.1]"=0.620,
                   "StepCox+Ridge"=0.626, "SuperPC"=0.619)
cv_summary$holdout_mean <- holdout_means[cv_summary$algorithm]
cv_summary$composite <- 0.5*cv_summary$holdout_mean + 0.5*cv_summary$cv_mean - 0.2*cv_summary$cv_sd
cv_summary <- cv_summary[order(-cv_summary$composite), ]
cat("\n=== Composite ranking (0.5*HO + 0.5*CV - 0.2*CV_SD) ===\n")
print(cv_summary[, c("algorithm","cv_mean","cv_sd","holdout_mean","composite")])

fwrite(cv_summary, file.path(MIME, "nefrs_v1_cv_summary.tsv"), sep="\t")
fwrite(cv_detail, file.path(MIME, "nefrs_v1_cv_per_fold.tsv"), sep="\t")
saveRDS(cv_results, file.path(MIME, "nefrs_v1_cv_raw.rds"))
cat("\nAll saved.\n")

cat("\n=== WINNER ===\n")
winner_idx <- which(!is.na(cv_summary$composite))[1]
if (!is.na(winner_idx)) {
  cat(sprintf("Top by composite: %s (composite=%.3f)\n",
              cv_summary$algorithm[winner_idx],
              cv_summary$composite[winner_idx]))
}
