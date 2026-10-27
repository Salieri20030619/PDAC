# W2 Day 1 entry diagnostic
# Goal: verify cross-cohort field/format consistency before unification
suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk"

cohorts <- list(
  list(name="PACA_AU",   dir="ICGC_PACA_AU", clin="clinical_std.PACA-AU.tsv",          exp="exp_array.PACA-AU.gene.tsv.gz"),
  list(name="GSE71729",  dir="GSE71729",     clin="GSE71729.clinical_std.tsv",          exp="GSE71729.gene.tsv.gz"),
  list(name="GSE21501",  dir="GSE21501",     clin="GSE21501.clinical_std.tsv",          exp="GSE21501.gene.tsv.gz"),
  list(name="GSE85916",  dir="GSE85916",     clin="GSE85916.clinical_std.tsv",          exp="GSE85916.gene.tsv.gz"),
  list(name="GSE62452",  dir="GSE62452",     clin="GSE62452.clinical_std.tsv",          exp="GSE62452.gene.tsv.gz"),
  list(name="GSE57495",  dir="GSE57495",     clin="GSE57495.clinical_std.tsv",          exp="GSE57495.gene.tsv.gz"),
  list(name="TCGA_basic", dir="TCGA_PAAD",   clin="TCGA-PAAD.clinical_std.tsv",         exp="TCGA-PAAD.curated.gene.tsv.gz"),
  list(name="TCGA_curated", dir="TCGA_PAAD", clin="curated.clinical_std.tsv",            exp="TCGA-PAAD.curated.gene.tsv.gz"),
  list(name="E_MTAB",    dir="E_MTAB_6134",  clin="clinical_std.E-MTAB-6134.valid.tsv", exp="exp_array.E-MTAB-6134.gene.tsv.gz")
)

cat("=== Per-cohort diagnostic ===\n")
report <- data.frame()
for (co in cohorts) {
  clin_path <- file.path(ROOT, co$dir, co$clin)
  exp_path  <- file.path(ROOT, co$dir, co$exp)
  
  cat(sprintf("\n[%s]\n", co$name))
  
  # Clinical
  if (file.exists(clin_path)) {
    clin <- fread(clin_path)
    cat(sprintf("  Clinical: %s (%d rows, %d cols)\n", basename(clin_path), nrow(clin), ncol(clin)))
    cat("    Columns:", paste(colnames(clin), collapse=", "), "\n")
    # Detect time/status columns (case-insensitive)
    cn_low <- tolower(colnames(clin))
    time_col <- colnames(clin)[grepl("^os_?time", cn_low)][1]
    status_col <- colnames(clin)[grepl("^os_?status|^os_?event", cn_low)][1]
    sid_col <- colnames(clin)[grepl("^sample_?id|^sampleid", cn_low)][1]
    cat(sprintf("    -> sample_id col: %s\n", sid_col))
    cat(sprintf("    -> time col: %s\n", time_col))
    cat(sprintf("    -> status col: %s\n", status_col))
    if (!is.na(time_col) && !is.na(status_col)) {
      t_vec <- clin[[time_col]]
      s_vec <- clin[[status_col]]
      cat(sprintf("    -> time NA: %d / range: [%.1f, %.1f]\n",
                  sum(is.na(t_vec)), min(t_vec, na.rm=TRUE), max(t_vec, na.rm=TRUE)))
      cat(sprintf("    -> status table: %s\n",
                  paste(names(table(s_vec, useNA="always")), table(s_vec, useNA="always"),
                        sep="=", collapse=", ")))
      valid <- sum(!is.na(t_vec) & !is.na(s_vec) & t_vec > 0)
      cat(sprintf("    -> valid (time>0 & both present): %d\n", valid))
    }
  } else {
    cat(sprintf("  Clinical: MISSING (%s)\n", clin_path))
  }
  
  # Expression: peek first column name only
  if (file.exists(exp_path)) {
    hdr <- fread(exp_path, nrows = 0)
    first_col <- colnames(hdr)[1]
    # Read first data row to see gene format
    one_row <- fread(exp_path, nrows = 1)
    first_gene <- as.character(one_row[[1]])
    cat(sprintf("  Expression: %s\n", basename(exp_path)))
    cat(sprintf("    First col: '%s' | First gene example: '%s'\n", first_col, first_gene))
    cat(sprintf("    Genes (rows): %s | Samples (cols-1): %d\n",
                format(nrow(fread(exp_path, select=1)), big.mark=","), ncol(hdr)-1))
    # Detect ENSG vs SYMBOL
    is_ensg <- grepl("^ENSG", first_gene)
    cat(sprintf("    Type: %s\n", ifelse(is_ensg, "ENSG (NEEDS conversion to SYMBOL!)", "SYMBOL")))
  }
}

cat("\n\n=== Cross-cohort summary ===\n")
cat("Action items for W2 Day 1 entry unification:\n")
cat("  1. E_MTAB OS_time/OS_status -> lowercase os_time_month/os_status\n")
cat("  2. E_MTAB use .valid.tsv (288), not full 309\n")
cat("  3. TCGA: confirm using curated 147, not basic 183\n")
cat("  4. TCGA gene matrix: confirm SYMBOL vs ENSG\n")
