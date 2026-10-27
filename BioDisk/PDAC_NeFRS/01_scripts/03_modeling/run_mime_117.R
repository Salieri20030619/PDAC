# ============================================================================
# W2 Day 2: Mime 117 algorithm combinations on PDAC NeFRS
# Input: 132 genes (paca p<0.05 + 4-GEO dir_agree >= 3/4)
# 5 selection cohorts: PACA-AU (train) + GSE71729 + GSE85916 + GSE62452 + GSE57495
# 3 hold-out cohorts NOT used here: GSE21501 + TCGA-PAAD + E_MTAB_6134
# ============================================================================

t_start <- Sys.time()
cat("=== Mime 117 PDAC run START ===\n")
cat("Start time:", format(t_start), "\n\n")

# === Load var.select patch BEFORE Mime1 ===
patch_path <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/03_modeling/patch_var_select.R"
if (!file.exists(patch_path)) stop("patch_var_select.R missing at: ", patch_path)
source(patch_path)
if (!exists("var.select", envir = .GlobalEnv))
  stop("var.select patch failed to load to GlobalEnv")
cat("var.select patch loaded OK\n\n")

# === Load packages (strict order from AML proven recipe) ===
cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(randomForestSRC)
  library(survival)
  library(glmnet)
  library(survivalsvm)
  library(gbm)
  library(superpc)
  library(CoxBoost)
  library(plsRcox)
  library(Mime1)
})

# Re-verify patch still in GlobalEnv after Mime1 (plsRcox can mask)
if (!exists("var.select", envir = .GlobalEnv)) {
  cat("WARNING: var.select got masked, re-sourcing patch\n")
  source(patch_path)
}
cat("All packages loaded\n\n")

# === Load data ===
W2_IN  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/02_w2_inputs"
FUNNEL <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/03_funnel"
OUT_DIR <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/04_mime"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Loading cohort_data_8.rds + step2c_loose_3of5.tsv...\n")
cohort_data <- readRDS(file.path(W2_IN, "cohort_data_8.rds"))
chosen_genes <- data.table::fread(file.path(FUNNEL, "step2c_loose_3of5.tsv"))$GeneSymbol
cat(sprintf("Selected candidates: %d genes\n", length(chosen_genes)))

# === Build Mime-compatible list ===
# Mime format: each cohort is data.frame with columns: ID, OS.time, OS, gene1, gene2, ...
build_mime_df <- function(co_name) {
  co <- cohort_data[[co_name]]
  exp_mat <- co$expression
  clin <- co$clinical
  # subset to chosen genes available in this cohort
  g_in <- intersect(chosen_genes, rownames(exp_mat))
  df <- data.frame(
    ID      = clin$sample_id,
    OS.time = clin$os_time_month,
    OS      = clin$os_status,
    t(exp_mat[g_in, , drop = FALSE]),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  df
}

# 5 selection cohorts (PACA-AU as train, 4 GEO as validation in Mime's internal CV)
sel_cohorts <- c("PACA_AU","GSE71729","GSE85916","GSE62452","GSE57495")

list_data <- list()
for (cn in sel_cohorts) {
  df <- build_mime_df(cn)
  list_data[[cn]] <- df
  cat(sprintf("  [%s]: %d samples x %d genes (+ID/OS.time/OS)\n",
              cn, nrow(df), ncol(df) - 3))
}

# Sanity: ensure all cohorts share same gene set (Mime requires it)
gene_sets <- lapply(list_data, function(d) colnames(d)[-(1:3)])
common_g <- Reduce(intersect, gene_sets)
cat(sprintf("\nGenes common to all 5 selection cohorts: %d\n", length(common_g)))

# Subset to common genes (Mime requires identical column structure)
for (cn in names(list_data)) {
  list_data[[cn]] <- list_data[[cn]][, c("ID","OS.time","OS", common_g)]
}

# Save input snapshot
saveRDS(list_data, file.path(OUT_DIR, "mime_input_list.rds"))
cat(sprintf("Saved mime_input_list.rds\n\n"))

# === Run Mime ML.Dev.Prog.Sig() ===
cat("=== Launching ML.Dev.Prog.Sig (117 combinations) ===\n")
cat("Train cohort: PACA_AU (267)\n")
cat("Validation:   GSE71729 + GSE85916 + GSE62452 + GSE57495\n")
cat("Genes:       ", length(common_g), "\n")
cat("Expected time: 30 min - 6 hours\n\n")

t_mime_start <- Sys.time()
res_mime <- tryCatch({
  ML.Dev.Prog.Sig(
    train_data           = list_data[[1]],
    list_train_vali_Data = list_data,
    candidate_genes      = common_g,
    unicox.filter.for.candi = FALSE,
    mode                 = "all",
    nodesize             = 5,
    seed                 = 5201314,
    cores_for_parallel   = 4
  )
}, error = function(e) {
  cat("\n!!! Mime FAILED with error:\n")
  cat(conditionMessage(e), "\n")
  cat("call:", deparse(e$call)[1], "\n")
  NULL
})

elapsed <- as.numeric(difftime(Sys.time(), t_mime_start, units = "mins"))
cat(sprintf("\nMime runtime: %.1f minutes\n", elapsed))

if (is.null(res_mime)) {
  cat("=== Mime aborted ===\n")
  quit(status = 1)
}

# Save full result
saveRDS(res_mime, file.path(OUT_DIR, "mime_results.rds"))
cat(sprintf("Saved mime_results.rds (full Mime output)\n"))

# === Extract C-index summary ===
cat("\n=== Top fields in result ===\n")
cat(paste(head(names(res_mime), 20), collapse=", "), "\n")

ci <- NULL
for (fname in c("Cindex.res","cindex","Cindex")) {
  if (!is.null(res_mime[[fname]])) { ci <- res_mime[[fname]]; break }
}

if (!is.null(ci)) {
  cat(sprintf("\nCindex table: %d rows x %d cols\n", nrow(ci), ncol(ci)))
  cat("Columns:", paste(colnames(ci), collapse=", "), "\n")
  
  # Try to make wide format if long
  wide <- NULL
  if (all(c("Model","Cindex") %in% colnames(ci))) {
    cohort_col <- intersect(c("Cohort","ID"), colnames(ci))[1]
    if (!is.na(cohort_col)) {
      library(reshape2)
      wide <- dcast(ci, as.formula(paste("Model ~", cohort_col)),
                    value.var = "Cindex")
    }
  } else {
    wide <- ci
  }
  
  if (!is.null(wide)) {
    model_col <- intersect(c("Model","method"), colnames(wide))[1]
    cohort_cols <- setdiff(colnames(wide), c(model_col,"ave","mean_cindex"))
    wide$mean_cindex <- rowMeans(wide[, cohort_cols, drop=FALSE], na.rm=TRUE)
    wide <- wide[order(-wide$mean_cindex), ]
    
    data.table::fwrite(wide, file.path(OUT_DIR, "mime_all_cindex.tsv"), sep="\t")
    data.table::fwrite(head(wide, 30), file.path(OUT_DIR, "mime_top30.tsv"), sep="\t")
    
    cat("\n=== Top 10 models by mean C-index ===\n")
    print(head(wide, 10))
    cat("\nSaved: mime_all_cindex.tsv, mime_top30.tsv\n")
  }
}

# Heatmap (best-effort)
tryCatch({
  pdf(file.path(OUT_DIR, "mime_cindex_heatmap.pdf"), width=10, height=16)
  cindex_dis_all(res_mime, validate_set = setdiff(names(list_data), "PACA_AU"))
  dev.off()
  cat("Saved heatmap\n")
}, error = function(e) cat("heatmap skipped:", conditionMessage(e), "\n"))

t_total <- as.numeric(difftime(Sys.time(), t_start, units="mins"))
cat(sprintf("\n=== TOTAL ELAPSED: %.1f minutes ===\n", t_total))
cat("=== Mime 117 PDAC run END ===\n")
