# Force install Bioconductor packages from OFFICIAL source
# (Tsinghua mirror BioC 3.22 index is broken)

options(timeout = 600)
options(Ncpus = 4)

# Step 0: override R's repos to force official BioC + CRAN
cat("=== Step 0: Override repos to official BioC ===\n")
local({
  r <- getOption("repos")
  # Keep CRAN at Tsinghua (CRAN mirror works fine)
  r["CRAN"] <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
  options(repos = r)
})

# CRITICAL: override BioC repos via BiocManager::repositories override
options(BioC_mirror = "https://bioconductor.org")
Sys.setenv(R_BIOC_MIRROR = "https://bioconductor.org")

cat("Current repos:\n")
print(getOption("repos"))
cat("\nBiocManager repositories():\n")
print(BiocManager::repositories())

cat("\nR:", R.version.string, "\n")
cat("BioC:", as.character(BiocManager::version()), "\n")

# Step 1: Force BioC packages from official source
cat("\n=== Step 1: BioC packages from OFFICIAL source ===\n")
bioc_pkgs <- c("GSVA", "cancerclass", "sparrow", "sva", "affy", "oligo")

# Use site_repository to explicitly point BiocManager::install at the official source
official_bioc_repos <- c(
  "BioCsoft"       = "https://bioconductor.org/packages/3.22/bioc",
  "BioCann"        = "https://bioconductor.org/packages/3.22/data/annotation",
  "BioCexp"        = "https://bioconductor.org/packages/3.22/data/experiment",
  "BioCworkflows"  = "https://bioconductor.org/packages/3.22/workflows",
  "CRAN"           = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
)
options(repos = official_bioc_repos)
cat("\nRepos after override:\n")
print(getOption("repos"))

for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("\n>>> Installing %s ...\n", p))
    tryCatch(
      BiocManager::install(p, update = FALSE, ask = FALSE,
                           site_repository = "https://bioconductor.org/packages/3.22/bioc"),
      error = function(e) cat(sprintf("!!! %s FAILED: %s\n", p, conditionMessage(e)))
    )
  } else {
    cat(sprintf("    %s already installed\n", p))
  }
}

cat("\n=== Step 2: Final verification ===\n")
all_targets <- c("GSVA","cancerclass","sparrow","sva","affy","oligo",
                 "CoxBoost","fastAdaboost")
status <- sapply(all_targets, function(p) requireNamespace(p, quietly = TRUE))
print(status)

if (all(status)) {
  cat("\nALL 8 TARGETS OK. Reloading Mime1:\n")
  suppressMessages(library(Mime1))
  cat("Mime1 v", as.character(packageVersion("Mime1")), "loaded.\n")
  cat("\n=== INSTALL COMPLETE ===\n")
} else {
  cat("\nSTILL FAILING:\n")
  print(names(status)[!status])
  cat("\n=== INSTALL INCOMPLETE ===\n")
}
