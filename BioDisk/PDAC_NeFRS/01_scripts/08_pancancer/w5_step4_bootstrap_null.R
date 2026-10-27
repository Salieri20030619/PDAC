#!/usr/bin/env Rscript
# ============================================================
# W5 Step 4: HCC/CRC 重拟合稳健性 — Bootstrap + 随机面板零分布
# 1) Bootstrap 乐观度校正 (Harrell): NeFRS-69 重拟合 C-index 校正值 + 区间
# 2) 随机 69 基因面板 ×N 的重拟合 CV C-index 零分布 → NeFRS 经验 p (特异性)
#    NeFRS 与随机面板用完全相同的 5折CV重拟合流程 (折内z-score,防泄漏)
# 输出: refit_null_{HCC,CRC}.tsv + Fig17C 零分布图 + 终端汇总
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(survival); library(glmnet); library(parallel)
})

ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
MIME <- file.path(ROOT, "02_results/04_mime")
PAN  <- file.path(ROOT, "00_raw_data/pancancer")
OUT  <- file.path(ROOT, "02_results/08_pancancer")
FIG  <- file.path(ROOT, "03_figures/Fig17_pancancer")
N_NULL <- 1000        # 随机面板数
B_BOOT <- 500         # bootstrap 次数
CORES  <- min(32, detectCores())
cat(sprintf("=== W5 Step 4: Bootstrap + 随机零分布 (cores=%d) ===\n", CORES))

genes69 <- fread(file.path(MIME, "nefrs_v3_enet_coefficients.tsv"))$gene

load_cohort <- function(tag) {
  gm <- fread(file.path(PAN, tag, paste0(tag, ".gene.tsv.gz")))
  cl <- fread(file.path(PAN, tag, paste0(tag, ".clinical_std.tsv")))
  em <- as.matrix(gm[, -1]); rownames(em) <- gm$gene
  cl <- cl[!is.na(os_time_month) & os_time_month > 0 & !is.na(os_status)]
  ids <- intersect(cl$sample_id, colnames(em))
  X <- t(em[, ids, drop=FALSE])                 # 样本 × 全基因 (raw log2)
  list(X = X, time = cl[match(ids, sample_id)]$os_time_month,
       status = cl[match(ids, sample_id)]$os_status)
}

# 折内z-score + LASSO-Cox 拟合 → 在测试集出线性预测
fit_predict <- function(Xtr_raw, ytr, Xte_raw) {
  mu <- colMeans(Xtr_raw); sv <- apply(Xtr_raw, 2, sd); sv[sv==0 | is.na(sv)] <- 1
  Xtr <- scale(Xtr_raw, center=mu, scale=sv)
  Xte <- scale(Xte_raw, center=mu, scale=sv)
  cvf <- tryCatch(suppressWarnings(cv.glmnet(Xtr, ytr, family="cox", nfolds=5)),
                  error=function(e) NULL)
  if (is.null(cvf)) return(rep(NA_real_, nrow(Xte_raw)))
  as.numeric(predict(cvf, newx=Xte, s="lambda.min"))
}

# 一次 5折CV 的 C-index (给定基因子集)
cv_cindex <- function(Xsub, time, status, seed) {
  n <- nrow(Xsub); set.seed(seed)
  folds <- sample(rep(1:5, length.out=n)); lp <- rep(NA_real_, n)
  for (f in 1:5) {
    tr <- folds != f; te <- folds == f
    lp[te] <- fit_predict(Xsub[tr,,drop=FALSE], Surv(time[tr], status[tr]), Xsub[te,,drop=FALSE])
  }
  ok <- !is.na(lp)
  if (sum(ok) < 10 || length(unique(lp[ok])) < 2) return(0.5)
  as.numeric(concordance(Surv(time[ok], status[ok]) ~ lp[ok], reverse=TRUE)$concordance)
}

# Bootstrap 乐观度校正 (Harrell)
boot_optimism <- function(X69, time, status, B, cores) {
  y <- Surv(time, status); n <- nrow(X69)
  lp_app <- fit_predict(X69, y, X69)
  C_app  <- as.numeric(concordance(y ~ lp_app, reverse=TRUE)$concordance)
  out <- mclapply(1:B, function(b) {
    set.seed(10000+b); idx <- sample(n, replace=TRUE)
    yb <- Surv(time[idx], status[idx])
    lp_bb <- fit_predict(X69[idx,,drop=FALSE], yb, X69[idx,,drop=FALSE])
    lp_bo <- fit_predict(X69[idx,,drop=FALSE], yb, X69)
    Cb <- tryCatch(as.numeric(concordance(yb ~ lp_bb, reverse=TRUE)$concordance), error=function(e) NA)
    Co <- tryCatch(as.numeric(concordance(y  ~ lp_bo, reverse=TRUE)$concordance), error=function(e) NA)
    c(opt = Cb - Co, Co = Co)
  }, mc.cores=cores)
  M <- do.call(rbind, out); M <- M[is.finite(M[,"opt"]) & is.finite(M[,"Co"]), , drop=FALSE]
  list(C_app=C_app, optimism=mean(M[,"opt"]), C_corr=C_app - mean(M[,"opt"]),
       Co_med=median(M[,"Co"]), Co_ci=quantile(M[,"Co"], c(.025,.975)))
}

# 随机面板零分布
null_dist <- function(Xall, time, status, cores) {
  pool <- colnames(Xall)[apply(Xall, 2, function(v) sd(v) > 0)]
  C_nef <- cv_cindex(Xall[, genes69, drop=FALSE], time, status, seed=42)
  null <- unlist(mclapply(1:N_NULL, function(i) {
    set.seed(50000+i); gs <- sample(pool, length(genes69))
    cv_cindex(Xall[, gs, drop=FALSE], time, status, seed=50000+i)
  }, mc.cores=cores))
  null <- null[is.finite(null)]
  list(C_nef=C_nef, null=null, p_emp=(sum(null >= C_nef)+1)/(length(null)+1),
       pct=mean(null < C_nef), med=median(null), q=quantile(null, c(.025,.975)))
}

dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
for (tag in c("HCC","CRC")) {
  co <- load_cohort(tag)
  cat(sprintf("\n########## %s (n=%d, events=%d) ##########\n", tag, nrow(co$X), sum(co$status)))

  cat("[1] Bootstrap 乐观度校正 (NeFRS-69 重拟合)...\n")
  bo <- boot_optimism(co$X[, genes69, drop=FALSE], co$time, co$status, B_BOOT, CORES)
  cat(sprintf("    Apparent C = %.3f | optimism = %.3f | 校正后 C = %.3f\n",
              bo$C_app, bo$optimism, bo$C_corr))
  cat(sprintf("    Bootstrap 外样本 C 中位=%.3f, 95%%区间 [%.3f, %.3f]\n",
              bo$Co_med, bo$Co_ci[1], bo$Co_ci[2]))

  cat("[2] 随机 69 基因面板零分布...\n")
  nd <- null_dist(co$X, co$time, co$status, CORES)
  cat(sprintf("    NeFRS CV C = %.3f | 随机面板中位=%.3f, 95%%=[%.3f,%.3f]\n",
              nd$C_nef, nd$med, nd$q[1], nd$q[2]))
  cat(sprintf("    NeFRS 超过 %.1f%% 随机面板 | 经验 p = %.4f\n", 100*nd$pct, nd$p_emp))

  fwrite(data.table(null_cindex=nd$null), file.path(OUT, sprintf("refit_null_%s.tsv", tag)))
  pdf(file.path(FIG, sprintf("Fig17C_refit_null_%s.pdf", tag)), width=7, height=5)
  hist(nd$null, breaks=40, col="grey85", border="white",
       main=sprintf("%s: NeFRS-69 vs random 69-gene panels (refit CV)", tag),
       xlab="Refit 5-fold CV C-index", xlim=range(c(nd$null, nd$C_nef))+c(-.02,.02))
  abline(v=nd$C_nef, col="#C00000", lwd=3)
  abline(v=nd$q, col="#0072B2", lwd=1.5, lty=2)
  legend("topright", c(sprintf("NeFRS = %.3f (p=%.3f)", nd$C_nef, nd$p_emp),
         "random 95%"), col=c("#C00000","#0072B2"), lwd=c(3,1.5), lty=c(1,2), bty="n")
  dev.off()
  cat(sprintf("    ✓ Fig17C_refit_null_%s.pdf\n", tag))
}
cat("\nDone.\n")
