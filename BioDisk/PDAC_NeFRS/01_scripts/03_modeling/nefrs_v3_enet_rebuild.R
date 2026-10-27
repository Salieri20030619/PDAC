# ============================================================================
# NeFRS v3 重建版: Enet[α=0.1] 在新候选池 170 基因上训练
# 跟旧 v3 方法学完全一致, 只是输入候选池从 132 → 170
# 修正: W1 错标 NRP 来源 → 真实 Hwang Sup Table 2 NRP 200
# ============================================================================
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

cat(sprintf("=== Mime input: %d genes ===\n", length(mime_genes)))

# Train Enet alpha=0.1 on PACA-AU
train_df <- mime_input[["PACA_AU"]]
X_train <- as.matrix(train_df[, mime_genes])
y_train <- Surv(train_df$OS.time, train_df$OS)

cat("\n=== Training Enet alpha=0.1 on PACA-AU ===\n")
set.seed(5201314)
cv_enet <- cv.glmnet(X_train, y_train, family="cox", alpha=0.1, nfolds=10)
cat(sprintf("Enet lambda.min: %.5f\n", cv_enet$lambda.min))
cat(sprintf("Enet lambda.1se: %.5f\n", cv_enet$lambda.1se))

# Show non-zero coefs at lambda.min and lambda.1se
for (s_choice in c("lambda.min", "lambda.1se")) {
  coefs <- as.numeric(coef(cv_enet, s=s_choice))
  names(coefs) <- mime_genes
  nz <- coefs[abs(coefs) > 1e-8]
  cat(sprintf("\n=== %s: %d non-zero genes ===\n", s_choice, length(nz)))
  nz_df <- data.table(gene=names(nz), coef=nz, abs_coef=abs(nz))
  nz_df <- merge(nz_df, gp[, .(GeneSymbol, priority, n_sources, Sources)],
                 by.x="gene", by.y="GeneSymbol", all.x=TRUE)
  setorder(nz_df, -abs_coef)
  cat("Top 20 by |coef|:\n")
  print(head(nz_df, 20))
}

# === 关键: 旧版用 lambda.min (53 基因), 新版我们也用 lambda.min ===
S_CHOICE <- "lambda.min"
LAMBDA_USED <- cv_enet$lambda.min
cat(sprintf("\n>>> Using %s = %.5f <<<\n", S_CHOICE, LAMBDA_USED))

# Evaluate on 8 cohorts
cat("\n=== 8-cohort evaluation ===\n")
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

  rs <- as.numeric(predict(cv_enet, newx=X, s=S_CHOICE))
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
cat(sprintf("\nTrain (PACA-AU) C-index: %.3f\n",
            results$c_index[results$role=="TRAIN"]))
cat(sprintf("Sel-val mean C-index:    %.3f (n=4 GEO)\n",
            mean(results$c_index[results$role=="SEL-VAL"])))
cat(sprintf("Hold-out mean C-index:   %.3f (n=3, GSE21501+TCGA+E-MTAB)\n",
            mean(results$c_index[results$role=="HOLDOUT"])))

# Save with full Sources info
coefs <- as.numeric(coef(cv_enet, s=S_CHOICE))
names(coefs) <- mime_genes
nz <- coefs[abs(coefs) > 1e-8]
nz_df <- data.table(gene=names(nz), coef=nz, abs_coef=abs(nz))
nz_df <- merge(nz_df, gp[, .(GeneSymbol, priority, Sources, n_sources, n_cohorts)],
               by.x="gene", by.y="GeneSymbol", all.x=TRUE)
setorder(nz_df, -abs_coef)
nz_df[, rank := .I]

# 标记每个基因是否含 真实 Hwang NRP 200
nz_df[, contains_real_NRP := grepl("Hwang2022_NatGenet_SupTable2", Sources)]

cat(sprintf("\n=== 新 NeFRS v3 (rebuilt): %d 基因 ===\n", nrow(nz_df)))
cat(sprintf("含 真实 Hwang NRP 200 来源: %d\n", sum(nz_df$contains_real_NRP)))
cat("\n=== Top 20 by |coef| (含 Sources) ===\n")
print(head(nz_df, 20))

cat("\n=== Priority 分布 ===\n")
print(table(nz_df$priority, useNA="always"))

# 保存
fwrite(nz_df, file.path(MIME, "nefrs_v3_enet_coefficients.tsv"), sep="\t")
fwrite(results, file.path(MIME, "nefrs_v3_enet_8cohort.tsv"), sep="\t")
saveRDS(cv_enet, file.path(MIME, "nefrs_v3_enet_model.rds"))

cat("\nSaved:\n")
cat("  nefrs_v3_enet_coefficients.tsv\n")
cat("  nefrs_v3_enet_8cohort.tsv\n")
cat("  nefrs_v3_enet_model.rds\n")

# === 跟旧 NeFRS 53 对比 ===
old_nefrs <- fread("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/_archive_pre_pathA_rebuild_20260529/04_mime/nefrs_v3_enet_coefficients.tsv")
cat(sprintf("\n=== 新旧 NeFRS 对比 ===\n"))
cat(sprintf("旧 NeFRS: %d 基因\n", nrow(old_nefrs)))
cat(sprintf("新 NeFRS: %d 基因\n", nrow(nz_df)))
common <- intersect(old_nefrs$gene, nz_df$gene)
cat(sprintf("交集: %d 基因\n", length(common)))
new_only <- setdiff(nz_df$gene, old_nefrs$gene)
old_only <- setdiff(old_nefrs$gene, nz_df$gene)
cat(sprintf("新版独有: %d (%s%s)\n", length(new_only),
            paste(head(new_only, 15), collapse=", "),
            if (length(new_only) > 15) "..." else ""))
cat(sprintf("旧版独有 (丢失): %d (%s)\n", length(old_only),
            paste(old_only, collapse=", ")))
