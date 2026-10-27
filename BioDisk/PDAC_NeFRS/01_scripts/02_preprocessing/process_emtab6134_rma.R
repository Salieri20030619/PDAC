# E-MTAB-6134 RMA normalization + probe-to-gene aggregation
# Platform: Affymetrix HG-U219 (A-GEOD-13667 = GPL13667)
# Output: gene-level expression matrix with Source Name (PANC189) as column names

suppressPackageStartupMessages({
  library(data.table)
  library(affy)          # ReadAffy + rma
  library(hgu219.db)     # probe2gene annotation (already installed for GSE85916)
  library(AnnotationDbi)
})

CEL_DIR    <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/E_MTAB_6134"
OUT_DIR    <- CEL_DIR
SDRF       <- file.path(CEL_DIR, "E-MTAB-6134.sdrf.txt")

setwd(CEL_DIR)
cat("=== Step 1: Read sdrf for CEL <-> Source Name mapping ===\n")
sdrf <- fread(SDRF, sep = "\t", header = TRUE, quote = "")
cat("sdrf dims:", dim(sdrf), "\n")

# Build mapping: CEL filename -> Source Name (patient ID)
map_dt <- sdrf[, .(cel_file = `Array Data File`, source_name = `Source Name`)]
map_dt <- unique(map_dt)
cat("Unique CEL files:", nrow(map_dt), "\n")
cat("Unique Source Names:", length(unique(map_dt$source_name)), "\n")
stopifnot(nrow(map_dt) == 309L)
stopifnot(length(unique(map_dt$source_name)) == 309L)

cat("\n=== Step 2: ReadAffy (309 .CEL) ===\n")
cel_files <- map_dt$cel_file
missing_cel <- cel_files[!file.exists(cel_files)]
if (length(missing_cel) > 0) stop("Missing .CEL files: ", paste(head(missing_cel), collapse=","))
cat("All 309 .CEL present on disk.\n")

t0 <- Sys.time()
af <- ReadAffy(filenames = cel_files)
cat("ReadAffy done in:", round(difftime(Sys.time(), t0, units="mins"), 2), "min\n")
cat("AffyBatch dims:", dim(exprs(af)), "\n")  # probes x samples
cat("CDF detected:", cdfName(af), "\n")

cat("\n=== Step 3: RMA normalization ===\n")
t0 <- Sys.time()
eset <- affy::rma(af)
cat("RMA done in:", round(difftime(Sys.time(), t0, units="mins"), 2), "min\n")
exp_mat <- exprs(eset)
cat("RMA matrix dims:", dim(exp_mat), "(probes x samples)\n")
cat("Value range (log2):", round(range(exp_mat), 2), "\n")
cat("Median:", round(median(exp_mat), 2), "\n")

cat("\n=== Step 4: Rename columns CEL -> Source Name ===\n")
# colnames currently = "H08C5473_E09.CEL", need to map to "PANC189"
m <- match(colnames(exp_mat), map_dt$cel_file)
stopifnot(!anyNA(m))
colnames(exp_mat) <- map_dt$source_name[m]
cat("Sample columns now (first 5):", head(colnames(exp_mat), 5), "\n")

cat("\n=== Step 5: Save probe-level matrix (backup) ===\n")
probe_out <- file.path(OUT_DIR, "exp_array.E-MTAB-6134.probe.tsv.gz")
fwrite(data.table(probe = rownames(exp_mat), exp_mat),
       probe_out, sep = "\t", compress = "gzip")
cat("Saved:", probe_out, "(", nrow(exp_mat), "probes )\n")

cat("\n=== Step 6: probe -> gene via hgu219.db ===\n")
# HG-U219 has ~49293 probes; hgu219.db maps them to Entrez/Symbol
probes <- rownames(exp_mat)
anno <- AnnotationDbi::select(hgu219.db,
                              keys = probes,
                              columns = c("SYMBOL", "ENTREZID"),
                              keytype = "PROBEID")
anno <- as.data.table(anno)
cat("Annotation rows:", nrow(anno), "; with SYMBOL:", sum(!is.na(anno$SYMBOL)), "\n")

# Keep only probes with non-NA SYMBOL
anno_ok <- anno[!is.na(SYMBOL) & SYMBOL != ""]
cat("Probes with valid SYMBOL:", nrow(anno_ok), "\n")
cat("Unique genes:", length(unique(anno_ok$SYMBOL)), "\n")

cat("\n=== Step 7: Aggregate probe -> gene (mean of probes per gene) ===\n")
# Convert matrix to long, merge with anno, aggregate
exp_dt <- as.data.table(exp_mat, keep.rownames = "PROBEID")
exp_dt <- merge(anno_ok[, .(PROBEID, SYMBOL)], exp_dt, by = "PROBEID")
exp_dt[, PROBEID := NULL]

# mean expression per SYMBOL
gene_mat_dt <- exp_dt[, lapply(.SD, mean, na.rm = TRUE), by = SYMBOL]
cat("Gene-level matrix dims:", dim(gene_mat_dt), "(genes x (1+samples))\n")

# Sanity: PDAC marker check
sanity_genes <- c("KRAS","TP53","SMAD4","CDKN2A","NTRK2","GPX4","ACSL4","NGF","NTRK1","BDNF")
cat("\nSanity gene presence:\n")
for (g in sanity_genes) {
  present <- g %in% gene_mat_dt$SYMBOL
  if (present) {
    v <- as.numeric(gene_mat_dt[SYMBOL == g, -1])
    cat(sprintf("  %-8s OK | mean=%.2f, sd=%.2f, range=[%.2f, %.2f]\n",
                g, mean(v), sd(v), min(v), max(v)))
  } else {
    cat(sprintf("  %-8s MISSING\n", g))
  }
}

cat("\n=== Step 8: Save final gene-level matrix ===\n")
gene_out <- file.path(OUT_DIR, "exp_array.E-MTAB-6134.gene.tsv.gz")
fwrite(gene_mat_dt, gene_out, sep = "\t", compress = "gzip")
cat("Saved:", gene_out, "\n")
cat("Final gene matrix:", nrow(gene_mat_dt), "genes x", ncol(gene_mat_dt)-1, "samples\n")

cat("\n=== RMA + probe2gene COMPLETE ===\n")
