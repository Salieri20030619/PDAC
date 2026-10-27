# ============================================================================
# W2 Day 1 Step C: 5-GSE cross-cohort consistency screening
# ----------------------------------------------------------------------------
# Input:  step1_sig_genes.tsv (PACA-AU p<0.05 genes)
# Process: For each Step 1 gene, run univariate Cox in each of 5 GEO cohorts;
#          record HR direction (>1 or <1).
#          Filter: STRICT default = direction agreement in all 5 GEO (5/5)
#          Fallback (auto): if STRICT yields <30 genes, relax to 4/5
# Output:
#   step2_gse_consist_matrix.tsv  (per-gene direction per cohort)
#   step2_consist_genes.tsv  (genes passing consistency -> NeFRS candidate pool)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

W2_IN   <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/02_w2_inputs"
FUNNEL  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/03_funnel"

LOG <- file.path(FUNNEL, sprintf("step2_consist_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
sink(LOG, split = TRUE)

cat("===========================================================\n")
cat("W2 Day 1 Step C: 5-GSE direction-consistency screening\n")
cat("Start:", date(), "\n")
cat("===========================================================\n\n")

# === Load ===
cohort_data <- readRDS(file.path(W2_IN, "cohort_data_8.rds"))
sig_step1   <- fread(file.path(FUNNEL, "step1_sig_genes.tsv"))

# 5 GEO validation cohorts (TCGA + E-MTAB held out)
gse_cohorts <- c("GSE71729", "GSE21501", "GSE85916", "GSE62452", "GSE57495")
cat(sprintf("Step 1 input: %d significant genes (PACA-AU p<0.05)\n", nrow(sig_step1)))
cat(sprintf("5 GEO cohorts for consistency: %s\n", paste(gse_cohorts, collapse=", ")))

# PACA-AU direction (reference)
ref_dir <- sig_step1[, .(GeneSymbol, paca_HR=HR, paca_dir=ifelse(HR>1, "High", "Low"))]
candidates <- sig_step1$GeneSymbol

# === Run Cox in each GEO cohort ===
all_results <- copy(ref_dir)

for (cohort_name in gse_cohorts) {
  cat(sprintf("\n--- [%s] ---\n", cohort_name))
  co <- cohort_data[[cohort_name]]
  if (is.null(co)) {
    cat("  SKIP: cohort data missing\n")
    next
  }
  
  exp_mat <- co$expression
  clin    <- co$clinical
  surv    <- Surv(clin$os_time_month, clin$os_status)
  
  # Genes available in this cohort
  avail <- intersect(candidates, rownames(exp_mat))
  cat(sprintf("  Genes available: %d / %d\n", length(avail), length(candidates)))
  
  # Cox per gene
  hr_vec <- rep(NA_real_, length(candidates))
  p_vec  <- rep(NA_real_, length(candidates))
  names(hr_vec) <- candidates
  names(p_vec)  <- candidates
  
  suppressWarnings({
    for (g in avail) {
      expr_g <- as.numeric(exp_mat[g, ])
      if (length(unique(expr_g)) < 2 || sd(expr_g, na.rm=TRUE) < 1e-6) next
      
      cobj <- tryCatch(coxph(surv ~ expr_g), error=function(e) NULL)
      if (is.null(cobj)) next
      
      hr_vec[g] <- exp(coef(cobj))
      p_vec[g]  <- summary(cobj)$coefficients[5]
    }
  })
  
  dir_vec <- ifelse(is.na(hr_vec), NA_character_,
                    ifelse(hr_vec > 1, "High", "Low"))
  
  cohort_dt <- data.table(
    GeneSymbol = candidates,
    HR = hr_vec, dir = dir_vec, P = p_vec
  )
  setnames(cohort_dt, c("HR","dir","P"),
           c(sprintf("%s_HR", cohort_name),
             sprintf("%s_dir", cohort_name),
             sprintf("%s_p", cohort_name)))
  
  all_results <- merge(all_results, cohort_dt, by="GeneSymbol", all.x=TRUE)
  cat(sprintf("  Cox done: %d non-NA HR\n", sum(!is.na(hr_vec))))
}

# === Count direction agreement ===
dir_cols <- sprintf("%s_dir", gse_cohorts)

count_agree <- function(row) {
  paca_d <- row[["paca_dir"]]
  gse_d <- unlist(row[dir_cols])
  gse_d_nona <- gse_d[!is.na(gse_d)]
  if (length(gse_d_nona) == 0) return(c(0L, 0L))
  agree <- sum(gse_d_nona == paca_d)
  total <- length(gse_d_nona)
  c(agree, total)
}

agree_mat <- t(apply(all_results, 1, count_agree))
all_results[, n_agree   := agree_mat[, 1]]
all_results[, n_tested  := agree_mat[, 2]]
all_results[, agree_ratio := ifelse(n_tested>0, n_agree/n_tested, NA)]

cat("\n--- Direction agreement summary ---\n")
cat(sprintf("Genes tested in all 5 GEO: %d\n", sum(all_results$n_tested == 5)))
cat(sprintf("Tested in 4 of 5:          %d\n", sum(all_results$n_tested == 4)))
cat(sprintf("Tested in <4:              %d\n", sum(all_results$n_tested < 4)))

cat("\nAgree count distribution (among 5-tested):\n")
tested5 <- all_results[n_tested == 5]
print(table(tested5$n_agree, useNA="ifany"))

# === Strict filter: 5/5 ===
STRICT_thr <- 5L
strict <- all_results[n_tested >= STRICT_thr & n_agree == STRICT_thr]
cat(sprintf("\n--- STRICT (5/5 agreement) ---\n"))
cat(sprintf("  Genes passing: %d\n", nrow(strict)))

# === Auto fallback if <30 genes ===
MIN_PASS <- 30L
if (nrow(strict) < MIN_PASS) {
  cat(sprintf("\n!!! STRICT (5/5) yields %d < %d, AUTO-RELAX to 4/5 !!!\n", nrow(strict), MIN_PASS))
  RELAX_thr <- 4L
  relax <- all_results[n_tested >= RELAX_thr & n_agree >= RELAX_thr]
  cat(sprintf("  Relaxed (4/5) passing: %d\n", nrow(relax)))
  final_pass <- relax
  used_thr <- "4/5"
} else {
  final_pass <- strict
  used_thr <- "5/5"
}

cat(sprintf("\nFinal threshold used: %s\n", used_thr))
cat(sprintf("Final consistent genes: %d\n", nrow(final_pass)))

# === Sort by mean significance + agreement ===
final_pass[, mean_p := rowMeans(.SD, na.rm=TRUE),
           .SDcols=c(sprintf("%s_p", gse_cohorts))]
final_pass <- final_pass[order(mean_p)]

# === P1 status ===
gene_pool <- fread("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/01_gene_pool/final_in_all_cohorts.tsv")
p1 <- gene_pool[priority == "P1_TopPriority", GeneSymbol]
cat(sprintf("\n--- P1 genes status in final consistent set ---\n"))
cat(sprintf("  P1 passing: %d / %d\n", sum(final_pass$GeneSymbol %in% p1), length(p1)))
cat("  P1 passers:\n")
print(final_pass[GeneSymbol %in% p1, .(GeneSymbol, paca_HR=round(paca_HR,3),
                                        n_agree, n_tested, mean_p=signif(mean_p,3))])

# === Top 20 by mean_p ===
cat("\n--- Top 20 consistent genes (sorted by mean_p across 5 GEO) ---\n")
top_cols <- c("GeneSymbol", "paca_HR", "paca_dir", "n_agree", "n_tested", "mean_p")
print(final_pass[1:min(20, nrow(final_pass)), ..top_cols])

# === Save ===
fwrite(all_results, file.path(FUNNEL, "step2_gse_consist_matrix.tsv"), sep="\t")
fwrite(final_pass, file.path(FUNNEL, "step2_consist_genes.tsv"), sep="\t")
saveRDS(list(all=all_results, pass=final_pass, used_thr=used_thr),
        file.path(FUNNEL, "step2_consist.rds"))

cat("\nSaved:\n")
cat("  step2_gse_consist_matrix.tsv (all", nrow(all_results), "genes x 5 GSE)\n")
cat(sprintf("  step2_consist_genes.tsv (%d passing -> NeFRS candidate pool)\n", nrow(final_pass)))

cat("\n===========================================================\n")
cat("STEP C COMPLETE\n")
cat(sprintf("Funnel so far: 1173 -> Step1 %d -> Step2 %d (%s)\n",
            nrow(sig_step1), nrow(final_pass), used_thr))
cat("End:", date(), "\n")
cat("===========================================================\n")
sink()
