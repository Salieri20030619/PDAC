#!/usr/bin/env Rscript
# ============================================================
# W5 Step 3: 为什么 NeFRS 不转移到 HCC/CRC
# 1) 69 基因单变量 Cox 扫描 (per-SD HR): 几个仍显著? 方向与 PDAC 系数一致?
# 2) 重拟合 CV C-index: 69 基因在各瘤种内重新加权后能否预后 (5折,折内z-score防泄漏)
# 输出: pancancer_gene_univariate.tsv + 终端汇总
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(survival); library(glmnet) })

ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
MIME <- file.path(ROOT, "02_results/04_mime")
PAN  <- file.path(ROOT, "00_raw_data/pancancer")
OUT  <- file.path(ROOT, "02_results/08_pancancer")

cat("=== W5 Step 3: 阴性结果性质诊断 (HCC + CRC) ===\n")

# NeFRS 69 基因 + PDAC 方向 (coef 符号)
nef <- fread(file.path(MIME, "nefrs_v3_enet_coefficients.tsv"))[, .(gene, pdac_coef = coef)]
nef[, pdac_dir := sign(pdac_coef)]
genes69 <- nef$gene

load_cohort <- function(tag) {
  gm <- fread(file.path(PAN, tag, paste0(tag, ".gene.tsv.gz")))
  cl <- fread(file.path(PAN, tag, paste0(tag, ".clinical_std.tsv")))
  em <- as.matrix(gm[, -1]); rownames(em) <- gm$gene
  cl <- cl[!is.na(os_time_month) & os_time_month > 0 & !is.na(os_status)]
  ids <- intersect(cl$sample_id, colnames(em))
  list(expr = em[, ids, drop=FALSE], time = cl[match(ids, sample_id)]$os_time_month,
       status = cl[match(ids, sample_id)]$os_status)
}

# ---- 1) 单变量 Cox 扫描 (per-SD HR) ----
univar_scan <- function(co, tag) {
  res <- rbindlist(lapply(genes69, function(g) {
    if (!(g %in% rownames(co$expr))) return(data.table(gene=g, HR=NA, p=NA, dir=NA))
    x <- as.numeric(co$expr[g, ]); s <- sd(x)
    if (is.na(s) || s == 0) return(data.table(gene=g, HR=NA, p=NA, dir=NA))
    xz <- (x - mean(x)) / s
    f <- tryCatch(coxph(Surv(co$time, co$status) ~ xz), error=function(e) NULL)
    if (is.null(f)) return(data.table(gene=g, HR=NA, p=NA, dir=NA))
    cf <- summary(f)$coefficients
    data.table(gene=g, HR=round(exp(cf[1,1]),3), p=cf[1,5], dir=sign(cf[1,1]))
  }))
  setnames(res, c("HR","p","dir"), paste0(tag, "_", c("HR","p","dir")))
  res
}

# ---- 2) 重拟合 CV C-index (5折, 折内z-score, LASSO-Cox) ----
refit_cv_cindex <- function(co, k=5, seed=42) {
  set.seed(seed)
  X <- t(co$expr[intersect(genes69, rownames(co$expr)), , drop=FALSE])  # 样本 × 基因
  n <- nrow(X); folds <- sample(rep(1:k, length.out=n)); lp <- rep(NA, n)
  for (f in 1:k) {
    tr <- which(folds != f); te <- which(folds == f)
    mu <- colMeans(X[tr,,drop=FALSE]); sv <- apply(X[tr,,drop=FALSE], 2, sd); sv[sv==0] <- 1
    Xtr <- scale(X[tr,,drop=FALSE], center=mu, scale=sv)
    Xte <- scale(X[te,,drop=FALSE], center=mu, scale=sv)
    cvf <- tryCatch(cv.glmnet(Xtr, Surv(co$time[tr], co$status[tr]), family="cox", nfolds=5),
                    error=function(e) NULL)
    if (is.null(cvf)) next
    lp[te] <- as.numeric(predict(cvf, newx=Xte, s="lambda.min"))
  }
  ok <- !is.na(lp)
  as.numeric(concordance(Surv(co$time[ok], co$status[ok]) ~ lp[ok], reverse=TRUE)$concordance)
}

# ---- 跑两瘤种 ----
tbl <- copy(nef)
cv <- list()
for (tag in c("HCC","CRC")) {
  co <- load_cohort(tag)
  cat(sprintf("\n[%s] n=%d events=%d\n", tag, length(co$time), sum(co$status)))
  sc <- univar_scan(co, tag)
  tbl <- merge(tbl, sc, by="gene", all.x=TRUE)
  nsig <- sum(tbl[[paste0(tag,"_p")]] < 0.05, na.rm=TRUE)
  conc <- tbl[[paste0(tag,"_dir")]] == tbl$pdac_dir
  nconc <- sum(conc, na.rm=TRUE); ntot <- sum(!is.na(conc))
  # 显著且方向一致 (真转移的基因)
  sig_conc <- sum(tbl[[paste0(tag,"_p")]] < 0.05 & conc, na.rm=TRUE)
  cat(sprintf("  单变量显著(p<0.05): %d/69\n", nsig))
  cat(sprintf("  方向与PDAC一致: %d/%d (%.0f%%)\n", nconc, ntot, 100*nconc/ntot))
  cat(sprintf("  显著且方向一致: %d/69\n", sig_conc))
  cv[[tag]] <- refit_cv_cindex(co)
  cat(sprintf("  重拟合 5折CV C-index (69基因瘤种内重新加权): %.3f\n", cv[[tag]]))
}

fwrite(tbl[order(-abs(pdac_coef))], file.path(OUT, "pancancer_gene_univariate.tsv"), sep="\t")
cat("\n=== 解读基准 ===\n")
cat("  方向一致率≈50% → 基因在该癌种行为与PDAC无关(随机),阴性=生物学语境不同\n")
cat("  重拟合C-index≈0.5 → 基因本身无预后信息(彻底阴性)\n")
cat("  重拟合C-index明显>0.5(如>0.6) → 程序相关但权重需组织特异校准(叙事可软化)\n")
cat("\n✓ 明细: 02_results/08_pancancer/pancancer_gene_univariate.tsv\nDone.\n")
