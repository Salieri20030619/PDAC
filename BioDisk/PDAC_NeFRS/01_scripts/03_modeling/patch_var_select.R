# ============================================================================
# patch_var_select.R  (PDAC NeFRS adaptation from AML LCDI project)
# Original: D:/26_01_AML/AML_Lactylation_PCD/scripts/phaseC1_patch_var_select.R
# Adapted: Linux server path, content unchanged
#
# Why: New randomForestSRC removed var.select(); Mime still calls it.
# Fix: Equivalent implementation using max.subtree() minimal depth method.
# ============================================================================

if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
  stop("randomForestSRC not installed")
}
suppressPackageStartupMessages(library(randomForestSRC))

`%||%` <- function(a, b) if (!is.null(a)) a else b

var.select <- function(object, conservative = "high", verbose = FALSE, ...) {
  if (!inherits(object, "rfsrc")) {
    stop("object must be rfsrc result")
  }
  
  tryCatch({
    maxsubtree <- randomForestSRC::max.subtree(object,
                    conservative = (conservative == "high"),
                    sub.order = FALSE)
    md <- maxsubtree$order[, 1]
    threshold <- maxsubtree$threshold
    
    topvars <- names(md)[md <= threshold]
    topvars <- topvars[!is.na(topvars)]
    
    if (verbose) {
      cat(sprintf("  var.select: threshold=%.3f, selected %d\n",
                  threshold, length(topvars)))
    }
    
    if (length(topvars) == 0) {
      sorted_md <- sort(md)
      topvars <- names(sorted_md)[seq_len(min(3, length(sorted_md)))]
    }
    
    return(list(
      topvars = topvars,
      varselect = data.frame(
        depth = md[topvars],
        vimp = object$importance[topvars]
      ),
      md.obj = maxsubtree,
      threshold = threshold
    ))
  }, error = function(e) {
    if (verbose) {
      cat(sprintf("  max.subtree failed (%s), fallback to importance\n", e$message))
    }
    
    imp <- object$importance
    if (is.null(imp)) stop("no importance, set importance=TRUE in rfsrc()")
    
    if (is.matrix(imp)) imp <- imp[, "all"] %||% imp[, 1]
    
    if (conservative == "high") {
      topvars <- names(imp)[imp > 0]
    } else {
      topvars <- names(imp)[imp > mean(imp, na.rm = TRUE)]
    }
    
    if (length(topvars) == 0) {
      sorted <- sort(imp, decreasing = TRUE)
      topvars <- names(sorted)[seq_len(min(3, length(sorted)))]
    }
    
    return(list(
      topvars = topvars,
      varselect = data.frame(vimp = imp[topvars]),
      threshold = NA
    ))
  })
}

cat("-- var.select patch loaded (PDAC NeFRS adaptation) --\n")
cat("  in .GlobalEnv:", exists("var.select", envir = .GlobalEnv), "\n")

# Smoke test
test_patch <- function() {
  cat("\n  Smoke test: small rfsrc + var.select...\n")
  set.seed(1)
  n <- 200
  df <- data.frame(
    time = rexp(n, 0.01),
    status = rbinom(n, 1, 0.6),
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), x4 = rnorm(n), x5 = rnorm(n)
  )
  df$time <- df$time + df$x1 * 100 + df$x2 * 50
  
  fit <- rfsrc(Surv(time, status) ~ ., data = df,
               ntree = 100, nodesize = 15, importance = TRUE)
  res <- var.select(object = fit, conservative = "high", verbose = TRUE)
  cat(sprintf("    topvars: %s\n", paste(res$topvars, collapse = ", ")))
  cat("    PATCH OK\n")
}

tryCatch(test_patch(), error = function(e) cat(sprintf("  TEST FAILED: %s\n", e$message)))

cat("\nNext: source this before library(Mime1) in main analysis scripts.\n")
