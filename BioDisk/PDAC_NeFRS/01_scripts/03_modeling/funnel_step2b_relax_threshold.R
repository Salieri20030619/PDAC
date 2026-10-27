# ============================================================================
# W2 Day 1 Step C-2 (FIXED): Relax to 4/5 + tabulate all tiers
# Includes mean_p computation before any subsetting
# ============================================================================

suppressPackageStartupMessages({ library(data.table) })

FUNNEL <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/03_funnel"
GP_IN  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/01_gene_pool"

LOG <- file.path(FUNNEL, sprintf("step2b_relax_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
sink(LOG, split = TRUE)

cat("===========================================================\n")
cat("W2 Day 1 Step C-2 (FIXED): Relax to 4/5\n")
cat("Start:", date(), "\n")
cat("===========================================================\n\n")

all_results <- fread(file.path(FUNNEL, "step2_gse_consist_matrix.tsv"))
gene_pool   <- fread(file.path(GP_IN, "final_in_all_cohorts.tsv"))
p1 <- gene_pool[priority == "P1_TopPriority", GeneSymbol]
p2 <- gene_pool[priority == "P2_HighPriority", GeneSymbol]

cat("Matrix dim:", nrow(all_results), "x", ncol(all_results), "\n")

# === Compute mean_p across 5 GEO p columns ===
gse_p_cols <- grep("_p$", colnames(all_results), value = TRUE)
cat("p-value columns found:", paste(gse_p_cols, collapse = ", "), "\n")
all_results[, mean_p := rowMeans(.SD, na.rm = TRUE), .SDcols = gse_p_cols]
cat("mean_p added, range:", round(min(all_results$mean_p, na.rm=TRUE),3), "-",
    round(max(all_results$mean_p, na.rm=TRUE),3), "\n\n")

# === Tabulate every threshold ===
cat("=== Genes passing each threshold ===\n")
for (thr in 1:5) {
  pass <- all_results[n_tested == 5 & n_agree >= thr]
  p1_in <- sum(pass$GeneSymbol %in% p1)
  p2_in <- sum(pass$GeneSymbol %in% p2)
  cat(sprintf("  >= %d/5: %d genes (P1: %d/%d, P2: %d/%d)\n",
              thr, nrow(pass), p1_in, length(p1), p2_in, length(p2)))
}

# === Build 3 tiers (sorted by mean_p) ===
tier_strict <- all_results[n_tested == 5 & n_agree == 5][order(mean_p)]
tier_main   <- all_results[n_tested == 5 & n_agree >= 4][order(mean_p)]
tier_loose  <- all_results[n_tested == 5 & n_agree >= 3][order(mean_p)]

cat(sprintf("\nFinal selection:\n"))
cat(sprintf("  STRICT 5/5 (paper Fig2A highlight): %d genes\n", nrow(tier_strict)))
cat(sprintf("  MAIN   4/5 (Mime input):            %d genes\n", nrow(tier_main)))
cat(sprintf("  LOOSE  3/5 (backup):                %d genes\n", nrow(tier_loose)))

# === Per-cohort outlier diagnosis ===
cat("\n=== Per-cohort outlier rate (in 4/5 main tier) ===\n")
cat("For each GEO cohort, how often it's THE outlier (1 disagree among 5):\n")
dir_cols <- grep("_dir$", colnames(all_results), value = TRUE)
gse_names <- gsub("_dir$", "", dir_cols)
for (i in seq_along(dir_cols)) {
  # Count cases where this cohort is the only one disagreeing with paca_dir
  outlier_count <- 0
  for (j in 1:nrow(tier_main)) {
    if (tier_main$n_agree[j] != 4) next  # only look at 4/5 cases
    paca_d <- tier_main$paca_dir[j]
    other_dirs <- unlist(tier_main[j, ..dir_cols])
    if (is.na(other_dirs[i])) next
    if (other_dirs[i] != paca_d) outlier_count <- outlier_count + 1
  }
  cat(sprintf("  %-10s: outlier in %d / %d (%.1f%%) of 4/5 cases\n",
              gse_names[i], outlier_count, sum(tier_main$n_agree == 4),
              100 * outlier_count / sum(tier_main$n_agree == 4)))
}

# === P1 in MAIN 4/5 ===
cat("\n=== P1 in MAIN 4/5 ===\n")
p1_in_main <- tier_main[GeneSymbol %in% p1]
cat(sprintf("  %d / %d P1 genes pass 4/5\n", nrow(p1_in_main), length(p1)))
if (nrow(p1_in_main) > 0) {
  print(p1_in_main[, .(GeneSymbol, paca_HR=round(paca_HR,3), paca_dir,
                       n_agree, mean_p=signif(mean_p, 3))])
}

# === Key focus genes ===
key_focus <- c("NTRK2","L1CAM","SEMA3A","NRP1","NRP2","BDNF",
               "NGF","NTRK1","TUBB3","NCAM1","GAP43","EZH2","NTN1","ROBO1")
focus_df <- all_results[GeneSymbol %in% key_focus,
                        .(GeneSymbol, paca_HR=round(paca_HR,3),
                          paca_dir, n_agree, mean_p=signif(mean_p,3))]
focus_df <- focus_df[match(intersect(key_focus, focus_df$GeneSymbol), GeneSymbol)]
cat("\n=== Key focus genes (all tiers) ===\n")
print(focus_df)

# === Save 3 tiers + overwrite step2_consist_genes.tsv with MAIN 4/5 ===
fwrite(tier_strict, file.path(FUNNEL, "step2c_strict_5of5.tsv"), sep="\t")
fwrite(tier_main,   file.path(FUNNEL, "step2c_main_4of5.tsv"), sep="\t")
fwrite(tier_loose,  file.path(FUNNEL, "step2c_loose_3of5.tsv"), sep="\t")
fwrite(tier_main,   file.path(FUNNEL, "step2_consist_genes.tsv"), sep="\t")

saveRDS(list(strict=tier_strict, main=tier_main, loose=tier_loose),
        file.path(FUNNEL, "step2c_tiers.rds"))

cat("\nSaved:\n")
cat(sprintf("  step2c_strict_5of5.tsv (%d genes, paper highlight)\n", nrow(tier_strict)))
cat(sprintf("  step2c_main_4of5.tsv   (%d genes, Mime input)\n", nrow(tier_main)))
cat(sprintf("  step2c_loose_3of5.tsv  (%d genes, backup)\n", nrow(tier_loose)))
cat("  step2_consist_genes.tsv  <- MAIN 4/5 (overwritten)\n")

cat("\n===========================================================\n")
cat("STEP C-2 COMPLETE\n")
cat(sprintf("Funnel: 1173 -> 213 -> MAIN 4/5: %d (Mime input)\n", nrow(tier_main)))
cat("End:", date(), "\n")
cat("===========================================================\n")
sink()
