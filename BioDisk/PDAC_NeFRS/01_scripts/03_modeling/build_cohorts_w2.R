# ============================================================================
# W2 Day 1 Step A: 8-cohort entry unification for Mime
# ----------------------------------------------------------------------------
# Output: Mime-compatible list, each element:
#   $name $role $expression (genes x samples) $clinical (sample_id/os_time_month/os_status)
# ============================================================================

suppressPackageStartupMessages({ library(data.table) })

ROOT    <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk"
OUT_DIR <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/02_w2_inputs"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(OUT_DIR, sprintf("build_cohorts_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
sink(LOG, split = TRUE)

cat("===========================================================\n")
cat("W2 Day 1 Step A: 8-cohort unification\n")
cat("Start:", date(), "\n")
cat("===========================================================\n\n")

cohorts <- list(
  list(name="PACA_AU",     dir="ICGC_PACA_AU",
       clin="clinical_std.PACA-AU.tsv",          exp="exp_array.PACA-AU.gene.tsv.gz",
       time_col="os_time_month", status_col="os_status", time_unit="month", role="train"),
  list(name="GSE71729",    dir="GSE71729",
       clin="GSE71729.clinical_std.tsv",         exp="GSE71729.gene.tsv.gz",
       time_col="os_time_month", status_col="os_status", time_unit="month", role="gse_val"),
  list(name="GSE21501",    dir="GSE21501",
       clin="GSE21501.clinical_std.tsv",         exp="GSE21501.gene.tsv.gz",
       time_col="os_time_month", status_col="os_status", time_unit="month", role="gse_val"),
  list(name="GSE85916",    dir="GSE85916",
       clin="GSE85916.clinical_std.tsv",         exp="GSE85916.gene.tsv.gz",
       time_col="os_time_month", status_col="os_status", time_unit="month", role="gse_val"),
  list(name="GSE62452",    dir="GSE62452",
       clin="GSE62452.clinical_std.tsv",         exp="GSE62452.gene.tsv.gz",
       time_col="os_time_month", status_col="os_status", time_unit="month", role="gse_val"),
  list(name="GSE57495",    dir="GSE57495",
       clin="GSE57495.clinical_std.tsv",         exp="GSE57495.gene.tsv.gz",
       time_col="os_time_month", status_col="os_status", time_unit="month", role="gse_val"),
  list(name="TCGA_PAAD",   dir="TCGA_PAAD",
       clin="TCGA-PAAD.curated.clinical_std.tsv", exp="TCGA-PAAD.curated.gene.tsv.gz",
       time_col="os_time_day", status_col="os_status", time_unit="day", role="holdout"),
  list(name="E_MTAB_6134", dir="E_MTAB_6134",
       clin="clinical_std.E-MTAB-6134.valid.tsv", exp="exp_array.E-MTAB-6134.gene.tsv.gz",
       time_col="OS_time", status_col="OS_status", time_unit="month", role="holdout")
)

load_cohort <- function(co) {
  cat(sprintf("\n--- [%s] (%s) ---\n", co$name, co$role))
  clin_path <- file.path(ROOT, co$dir, co$clin)
  exp_path  <- file.path(ROOT, co$dir, co$exp)
  stopifnot(file.exists(clin_path), file.exists(exp_path))
  
  # Clinical
  clin <- fread(clin_path)
  t_vec <- as.numeric(clin[[co$time_col]])
  s_vec <- as.integer(clin[[co$status_col]])
  if (co$time_unit == "day") {
    t_vec <- t_vec / 30.44
    cat("  Converted day->month\n")
  }
  clin_std <- data.table(sample_id=clin$sample_id, os_time_month=t_vec, os_status=s_vec)
  clin_std <- clin_std[!is.na(os_time_month) & !is.na(os_status) &
                       os_time_month > 0 & os_status %in% c(0L, 1L)]
  cat(sprintf("  Clinical valid: %d\n", nrow(clin_std)))
  
  # Expression
  exp_dt <- fread(exp_path)
  symbols <- as.character(exp_dt[[1]])
  exp_mat <- as.matrix(exp_dt[, -1, with=FALSE])
  rownames(exp_mat) <- symbols
  
  # Handle multi-mapped 'A///B///C' -> keep first
  bad <- grepl("///", rownames(exp_mat))
  if (any(bad)) {
    cat(sprintf("  Multi-mapped SYMBOLs: %d -> keep first\n", sum(bad)))
    rownames(exp_mat)[bad] <- sapply(strsplit(rownames(exp_mat)[bad], "///"), `[`, 1)
  }
  
  # Aggregate duplicated SYMBOLs by mean
  if (anyDuplicated(rownames(exp_mat))) {
    n_dup <- sum(duplicated(rownames(exp_mat)))
    cat(sprintf("  Duplicated SYMBOLs: %d -> mean aggregate\n", n_dup))
    df_tmp <- data.table(SYMBOL=rownames(exp_mat), as.data.table(exp_mat))
    df_agg <- df_tmp[, lapply(.SD, mean, na.rm=TRUE), by=SYMBOL]
    syms <- df_agg$SYMBOL
    exp_mat <- as.matrix(df_agg[, -1, with=FALSE])
    rownames(exp_mat) <- syms
  }
  cat(sprintf("  Expression: %d genes x %d samples\n", nrow(exp_mat), ncol(exp_mat)))
  
  # Intersect
  common <- intersect(clin_std$sample_id, colnames(exp_mat))
  if (length(common) < 30) warning(sprintf("[%s] only %d common samples", co$name, length(common)))
  
  clin_std <- clin_std[match(common, clin_std$sample_id)]
  exp_mat  <- exp_mat[, common, drop=FALSE]
  stopifnot(identical(clin_std$sample_id, colnames(exp_mat)))
  
  death_rate <- mean(clin_std$os_status == 1)
  cat(sprintf("  Final: %d samples, deaths %d (%.1f%%), OS median %.1f mo\n",
              ncol(exp_mat), sum(clin_std$os_status==1), death_rate*100,
              median(clin_std$os_time_month)))
  
  list(name=co$name, role=co$role, expression=exp_mat,
       clinical=as.data.frame(clin_std), n_samples=ncol(exp_mat),
       n_genes=nrow(exp_mat), death_rate=death_rate)
}

cohort_data <- list()
for (co in cohorts) {
  cohort_data[[co$name]] <- tryCatch(load_cohort(co), error=function(e) {
    cat(sprintf("!!! [%s] FAILED: %s\n", co$name, conditionMessage(e))); NULL
  })
}

cat("\n===========================================================\n")
cat("=== Summary ===\n")
cat("===========================================================\n")
summary_df <- data.frame(
  cohort     = sapply(cohort_data, function(x) if (is.null(x)) NA else x$name),
  role       = sapply(cohort_data, function(x) if (is.null(x)) NA else x$role),
  n_samples  = sapply(cohort_data, function(x) if (is.null(x)) NA else x$n_samples),
  n_genes    = sapply(cohort_data, function(x) if (is.null(x)) NA else x$n_genes),
  death_pct  = sapply(cohort_data, function(x) if (is.null(x)) NA else round(x$death_rate*100, 1))
)
print(summary_df, row.names=FALSE)
cat(sprintf("\nTotal samples: %d\n", sum(summary_df$n_samples, na.rm=TRUE)))

# Common genes across all 8
gene_lists <- lapply(cohort_data, function(x) if (is.null(x)) NULL else rownames(x$expression))
gene_lists <- gene_lists[!sapply(gene_lists, is.null)]
common_genes_8 <- Reduce(intersect, gene_lists)
cat(sprintf("\n8-cohort common genes: %d\n", length(common_genes_8)))

# Sanity: key genes
key_genes <- c("KRAS","TP53","SMAD4","CDKN2A","NTRK2","NTRK1","GPX4","SLC7A11","NGF","BDNF",
               "ACSL4","GLIS3","SOX4","NEUROD1","INSM1","NFIB","NRP1","NRP2","SEMA3A","TUBB3",
               "L1CAM","GAP43","EZH2","NTN1","REST")
key_in <- intersect(unique(key_genes), common_genes_8)
cat(sprintf("Key gene sanity: %d / %d in 8-cohort common\n",
            length(key_in), length(unique(key_genes))))
key_missing <- setdiff(unique(key_genes), common_genes_8)
if (length(key_missing) > 0) cat("MISSING:", paste(key_missing, collapse=", "), "\n")

# Save
saveRDS(cohort_data, file.path(OUT_DIR, "cohort_data_8.rds"))
saveRDS(common_genes_8, file.path(OUT_DIR, "common_genes_8cohorts.rds"))
fwrite(data.table(gene=common_genes_8), file.path(OUT_DIR, "common_genes_8cohorts.tsv"), sep="\t")
fwrite(summary_df, file.path(OUT_DIR, "cohort_summary.tsv"), sep="\t")
cat("\nSaved:\n")
cat("  cohort_data_8.rds  (Mime input ready)\n")
cat("  common_genes_8cohorts.rds/.tsv\n")
cat("  cohort_summary.tsv\n")

cat("\n===========================================================\n")
cat("STEP A COMPLETE\n")
cat("End:", date(), "\n")
cat("===========================================================\n")
sink()
