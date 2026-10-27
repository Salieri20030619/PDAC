# Rebuild preprocessCore with --disable-threading to bypass
# the pthread_attr_setstacksize EINVAL bug on glibc 2.34+

options(timeout = 1200)
options(Ncpus = 4)

# Use official BioC source (Tsinghua mirror BioC 3.22 is broken)
options(repos = c(
  "BioCsoft"      = "https://bioconductor.org/packages/3.22/bioc",
  "BioCann"       = "https://bioconductor.org/packages/3.22/data/annotation",
  "BioCexp"       = "https://bioconductor.org/packages/3.22/data/experiment",
  "BioCworkflows" = "https://bioconductor.org/packages/3.22/workflows",
  "CRAN"          = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
))

cat("=== Before rebuild ===\n")
cat("preprocessCore version:", as.character(packageVersion("preprocessCore")), "\n")

cat("\n=== Rebuilding preprocessCore from source with --disable-threading ===\n")
BiocManager::install(
  "preprocessCore",
  configure.args = "--disable-threading",
  update = FALSE,
  ask = FALSE,
  force = TRUE,
  type = "source",
  site_repository = "https://bioconductor.org/packages/3.22/bioc"
)

cat("\n=== After rebuild ===\n")
# Force reload (unload old DLL first)
try(unloadNamespace("preprocessCore"), silent = TRUE)
cat("preprocessCore version:", as.character(packageVersion("preprocessCore")), "\n")

cat("\n=== Verify no pthread symbols (should be EMPTY) ===\n")
so_path <- system.file("libs", "preprocessCore.so", package = "preprocessCore")
cat("SO path:", so_path, "\n")
result <- system(paste0("nm -D '", so_path, "' | grep pthread"), intern = TRUE)
if (length(result) == 0) {
  cat("OK: no pthread symbols (threading disabled successfully)\n")
} else {
  cat("WARNING: pthread symbols still present:\n")
  print(result)
}

cat("\n=== Quick smoke test ===\n")
suppressMessages(library(preprocessCore))
# Tiny RMA-like test
set.seed(42)
x <- matrix(rlnorm(1000), 100, 10)
out <- normalize.quantiles(x)
cat("normalize.quantiles output dims:", dim(out), "(expected 100 10)\n")
cat("Output range:", round(range(out), 3), "\n")
cat("\n=== REBUILD COMPLETE ===\n")
