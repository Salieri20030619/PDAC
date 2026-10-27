# ============================================================================
# W2 Day 1 Step B: PACA-AU univariate Cox screen
# ----------------------------------------------------------------------------
# Input:
#   - 02_w2_inputs/cohort_data_8.rds  (PACA_AU training cohort)
#   - 02_results/01_gene_pool/final_in_all_cohorts.tsv  (1173 candidates)
#   - 02_w2_inputs/common_genes_8cohorts.rds  (8-cohort common)
#
# Process:
#   1. Candidate pool = 1173 INTERSECT 8-cohort common
#   2. PACA-AU univariate Cox on each gene
#   3. Filter p < 0.05 -> ~150-200 genes
#
# Output:
#   - 03_funnel/step1_uniCox_paca_au.tsv (all results)
#   - 03_funnel/step1_sig_genes.tsv (p<0.05, passed to Step 2)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

W2_IN  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/02_w2_inputs"
GP_IN  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/01_gene_pool"
OUT_DIR <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/03_funnel"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(OUT_DIR, sprintf("step1_uniCox_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
sink(LOG, split = TRUE)

cat("===========================================================\n")
cat("W2 Day 1 Step B: PACA-AU univariate Cox\n")
cat("Start:", date(), "\n")
cat("===========================================================\n\n")

# === Load ===
cohort_data    <- readRDS(file.path(W2_IN, "cohort_data_8.rds"))
common_genes_8 <- readRDS(file.path(W2_IN, "common_genes_8cohorts.rds"))
gene_pool      <- fread(file.path(GP_IN, "final_in_all_cohorts.tsv"))

cat(sprintf("Gene pool (1173 candidates): %d\n", nrow(gene_pool)))
cat(sprintf("8-cohort common genes: %d\n", length(common_genes_8)))

# === Step 1.1: Candidate pool = 1173 INTERSECT 8-cohort common ===
candidates_in_8 <- intersect(gene_pool$GeneSymbol, common_genes_8)
cat(sprintf("\nCandidates available in 8 cohorts: %d / %d\n",
            length(candidates_in_8), nrow(gene_pool)))

# Subset gene_pool to keep priority/categories info
gene_pool_8 <- gene_pool[GeneSymbol %in% candidates_in_8]
cat(sprintf("Priority distribution after 8-cohort filter:\n"))
print(table(gene_pool_8$priority))

# === Step 1.2: PACA-AU univariate Cox ===
paca <- cohort_data$PACA_AU
exp_train <- paca$expression
clin_train <- paca$clinical
cat(sprintf("\n--- PACA-AU training cohort ---\n"))
cat(sprintf("  Samples: %d\n", ncol(exp_train)))
cat(sprintf("  Genes (all): %d\n", nrow(exp_train)))
cat(sprintf("  Death events: %d\n", sum(clin_train$os_status == 1)))

surv_obj <- Surv(clin_train$os_time_month, clin_train$os_status)

cat(sprintf("\nRunning univariate Cox on %d candidates...\n", length(candidates_in_8)))
t0 <- Sys.time()

cox_res <- data.table(GeneSymbol=character(), HR=numeric(),
                      HR_lower=numeric(), HR_upper=numeric(),
                      P_value=numeric(), Risk=character())

# Vectorized loop (suppress warnings about ties)
suppressWarnings({
  for (i in seq_along(candidates_in_8)) {
    g <- candidates_in_8[i]
    expr_g <- as.numeric(exp_train[g, ])
    
    # Skip if zero variance
    if (length(unique(expr_g)) < 2 || sd(expr_g, na.rm=TRUE) < 1e-6) next
    
    cox_obj <- tryCatch(coxph(surv_obj ~ expr_g), error=function(e) NULL)
    if (is.null(cox_obj)) next
    cs <- summary(cox_obj)
    hr <- exp(coef(cox_obj))
    ci <- tryCatch(exp(confint(cox_obj)), error=function(e) c(NA, NA))
    p  <- cs$coefficients[5]
    
    cox_res <- rbind(cox_res, data.table(
      GeneSymbol = g,
      HR         = as.numeric(hr),
      HR_lower   = as.numeric(ci[1]),
      HR_upper   = as.numeric(ci[2]),
      P_value    = as.numeric(p),
      Risk       = ifelse(hr > 1, "High", "Low")
    ))
    
    if (i %% 200 == 0) {
      cat(sprintf("  Progress: %d / %d\n", i, length(candidates_in_8)))
    }
  }
})

cat(sprintf("\nDone in %.1f sec | %d genes with Cox result\n",
            as.numeric(difftime(Sys.time(), t0, units="secs")), nrow(cox_res)))

# === Step 1.3: merge with priority info + filter p < 0.05 ===
cox_res <- merge(cox_res,
                 gene_pool_8[, .(GeneSymbol, priority, GeneCategories, Sources, n_sources)],
                 by="GeneSymbol", all.x=TRUE)

cox_res <- cox_res[order(P_value)]

sig <- cox_res[P_value < 0.05]
cat(sprintf("\n--- Filter p < 0.05 ---\n"))
cat(sprintf("  Significant: %d / %d (%.1f%%)\n",
            nrow(sig), nrow(cox_res), 100*nrow(sig)/nrow(cox_res)))
cat(sprintf("  High risk (HR>1): %d\n", sum(sig$Risk == "High")))
cat(sprintf("  Low risk  (HR<1): %d\n", sum(sig$Risk == "Low")))

cat("\nPriority distribution in significant set:\n")
print(table(sig$priority, useNA="always"))

cat("\n--- Top 20 by p-value ---\n")
print(sig[1:20, .(GeneSymbol, HR=round(HR,3), P_value=signif(P_value,3),
                  Risk, priority, Sources)])

# === Special: P1 (20 genes) status in this filter ===
p1_genes <- gene_pool[priority == "P1_TopPriority", GeneSymbol]
cat(sprintf("\n--- P1 (20 top-priority genes) status ---\n"))
p1_in_cox <- cox_res[GeneSymbol %in% p1_genes]
cat(sprintf("  P1 in 8-cohort + Cox-able: %d / %d\n", nrow(p1_in_cox), length(p1_genes)))
p1_sig <- sig[GeneSymbol %in% p1_genes]
cat(sprintf("  P1 passing p<0.05: %d / %d\n", nrow(p1_sig), length(p1_genes)))
cat("  P1 genes passing:\n")
print(p1_sig[, .(GeneSymbol, HR=round(HR,3), P_value=signif(P_value,3), Risk, Sources)])

# === Save ===
fwrite(cox_res, file.path(OUT_DIR, "step1_uniCox_paca_au.tsv"), sep="\t")
fwrite(sig,     file.path(OUT_DIR, "step1_sig_genes.tsv"), sep="\t")
saveRDS(list(cox_all=cox_res, sig=sig), file.path(OUT_DIR, "step1_uniCox.rds"))

cat("\nSaved:\n")
cat("  step1_uniCox_paca_au.tsv  (all", nrow(cox_res), "genes)\n")
cat("  step1_sig_genes.tsv       (sig", nrow(sig), "genes -> Step 2)\n")
cat("  step1_uniCox.rds\n")

cat("\n===========================================================\n")
cat("STEP B COMPLETE\n")
cat("End:", date(), "\n")
cat("===========================================================\n")
sink()
