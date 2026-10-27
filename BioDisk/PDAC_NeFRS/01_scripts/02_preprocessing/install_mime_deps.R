# Install missing Mime dependencies
options(timeout = 600)
options(Ncpus = 4)

cat("=== R session info ===\n")
cat("R version:", R.version.string, "\n")
cat("BioC version:", as.character(BiocManager::version()), "\n")
cat(".libPaths():\n"); print(.libPaths())

cat("\n=== Step 1: Bioconductor packages ===\n")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_pkgs <- c("GSVA", "cancerclass", "sparrow", "sva")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("\n>>> Installing %s ...\n", p))
    tryCatch(
      BiocManager::install(p, update = FALSE, ask = FALSE),
      error = function(e) cat(sprintf("!!! %s FAILED: %s\n", p, conditionMessage(e)))
    )
  } else {
    cat(sprintf("    %s already installed\n", p))
  }
}

cat("\n=== Step 2: GitHub packages ===\n")
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", repos = "https://cloud.r-project.org")

if (!requireNamespace("CoxBoost", quietly = TRUE)) {
  cat("\n>>> Installing CoxBoost from github binderh/CoxBoost ...\n")
  tryCatch(
    remotes::install_github("binderh/CoxBoost", upgrade = "never"),
    error = function(e) cat("!!! CoxBoost FAILED:", conditionMessage(e), "\n")
  )
}

if (!requireNamespace("fastAdaboost", quietly = TRUE)) {
  cat("\n>>> Installing fastAdaboost from github souravc83/fastAdaboost ...\n")
  tryCatch(
    remotes::install_github("souravc83/fastAdaboost", upgrade = "never"),
    error = function(e) cat("!!! fastAdaboost FAILED:", conditionMessage(e), "\n")
  )
}

cat("\n=== Step 3: Final verification ===\n")
all_deps <- c("CoxBoost","GSVA","cancerclass","sparrow","sva","fastAdaboost")
status <- sapply(all_deps, function(p) requireNamespace(p, quietly = TRUE))
print(status)

if (all(status)) {
  cat("\nALL 6 DEPS OK. Reloading Mime1:\n")
  suppressMessages(library(Mime1))
  cat("Mime1 v", as.character(packageVersion("Mime1")), "loaded.\n")
  cat("\n=== INSTALL COMPLETE ===\n")
} else {
  cat("\nSTILL FAILING:\n")
  print(names(status)[!status])
  cat("\n=== INSTALL INCOMPLETE ===\n")
}
