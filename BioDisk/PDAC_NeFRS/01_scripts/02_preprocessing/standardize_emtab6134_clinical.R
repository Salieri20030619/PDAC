# Standardize E-MTAB-6134 clinical info from sdrf
# Output schema (align with PACA-AU clinical_std.PACA-AU.tsv):
#   sample_id, OS_time (month), OS_status (0/1), sex, stage_TNM, grade,
#   DFS_time, DFS_event, subtype_wholetum, kras_mut, tp53_mut, cdkn2a_mut

suppressPackageStartupMessages({
  library(data.table)
})

CEL_DIR  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/E_MTAB_6134"
SDRF     <- file.path(CEL_DIR, "E-MTAB-6134.sdrf.txt")
OUT      <- file.path(CEL_DIR, "clinical_std.E-MTAB-6134.tsv")

cat("=== Read sdrf ===\n")
sdrf <- fread(SDRF, sep = "\t", header = TRUE, quote = "")
cat("sdrf dims:", dim(sdrf), "\n")

cat("\n=== Build clinical_std ===\n")
clin <- sdrf[, .(
  sample_id          = `Source Name`,
  OS_time            = as.numeric(`Characteristics[os.delay]`),
  OS_status          = as.integer(`Characteristics[os.event]`),
  DFS_time           = as.numeric(`Characteristics[dfs.delay]`),
  DFS_event          = as.integer(`Characteristics[dfs.event]`),
  sex                = `Characteristics[sex]`,
  stage_TNM          = `Characteristics[TNM tumour grading]`,
  grade              = `Characteristics[tumor grading]`,
  resection_margin   = `Characteristics[resection margin]`,
  clinical_center    = `Characteristics[clinical center]`,
  subtype_wholetum   = `Characteristics[wholetumclassif]`,
  subtype_hightum    = `Characteristics[hightumcellclassif]`,
  kras_mut           = `Characteristics[krasmut]`,
  tp53_mut           = `Characteristics[tp53mut]`,
  cdkn2a_mut         = `Characteristics[cdkn2amut]`,
  avg_vaf            = as.numeric(`Characteristics[average vaf]`)
)]
clin <- unique(clin)
cat("Clinical rows after dedup:", nrow(clin), "\n")

cat("\n=== Sanity check ===\n")
cat("sample_id unique:", length(unique(clin$sample_id)), "\n")
cat("\nOS_time summary (months):\n"); print(summary(clin$OS_time))
cat("\nOS_status table:\n"); print(table(clin$OS_status, useNA = "always"))
cat("\nDFS_event table:\n"); print(table(clin$DFS_event, useNA = "always"))
cat("\nSex table:\n"); print(table(clin$sex, useNA = "always"))
cat("\nGrade table:\n"); print(table(clin$grade, useNA = "always"))
cat("\nstage_TNM (top 10):\n"); print(head(sort(table(clin$stage_TNM), decreasing = TRUE), 10))
cat("\nsubtype_wholetum table:\n"); print(table(clin$subtype_wholetum, useNA = "always"))
cat("\nKRAS mut table:\n"); print(table(clin$kras_mut, useNA = "always"))

cat("\n=== Filter: Primary + complete OS (alignment with PACA-AU) ===\n")
clin_valid <- clin[!is.na(OS_time) & !is.na(OS_status)]
cat("After OS filter:", nrow(clin_valid), "\n")

cat("\n=== Save ===\n")
fwrite(clin, OUT, sep = "\t")
cat("Saved (all 309):", OUT, "\n")

OUT_VALID <- file.path(CEL_DIR, "clinical_std.E-MTAB-6134.valid.tsv")
fwrite(clin_valid, OUT_VALID, sep = "\t")
cat("Saved (valid OS):", OUT_VALID, "(", nrow(clin_valid), "samples )\n")

cat("\n=== Death events:", sum(clin_valid$OS_status == 1, na.rm=TRUE), 
    "/", nrow(clin_valid), 
    "(", round(100*mean(clin_valid$OS_status == 1, na.rm=TRUE),1), "%)\n")
cat("=== OS_time median:", round(median(clin_valid$OS_time, na.rm=TRUE), 1), "months\n")
cat("\n=== CLINICAL STANDARDIZATION COMPLETE ===\n")
