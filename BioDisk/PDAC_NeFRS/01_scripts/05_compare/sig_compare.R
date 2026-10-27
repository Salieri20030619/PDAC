# 86 published signature vs NeFRS: 各队列 C-index 对比 (AIDPS 方法)
suppressPackageStartupMessages({
  library(data.table); library(survival); library(survcomp)
})
ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
OUT  <- file.path(ROOT,"02_results/13_sig_compare")
FIG  <- file.path(ROOT,"03_figures/Fig12_sigcompare")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE); dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

cat("=== 86 signature vs NeFRS C-index 对比 ===\n")

# 1. 86 signature
sig <- fread(file.path(ROOT,"00_raw_data/aidps_86_signatures.tsv"))
sig <- sig[!is.na(Coef) & Gene!=""]
models <- unique(sig$Model)
cat(sprintf("86 signature: %d 个 model, %d 基因条目\n", length(models), nrow(sig)))

# 2. NeFRS 系数
nefrs <- fread(file.path(ROOT,"02_results/04_mime/nefrs_v3_enet_coefficients.tsv"))
nefrs_genes <- nefrs$gene; nefrs_coef <- setNames(nefrs$coef, nefrs$gene)
cat(sprintf("NeFRS: %d 基因\n", length(nefrs_genes)))

# 3. 队列清单 (8 队列: 训练 + 验证)
cohorts <- list(
  PACA_AU   = list(expr="00_raw_data/bulk/ICGC_PACA_AU/exp_array.PACA-AU.gene.tsv.gz",
                   clin="00_raw_data/bulk/ICGC_PACA_AU/clinical_std.PACA-AU.tsv"),
  GSE71729  = list(expr="00_raw_data/bulk/GSE71729/GSE71729.gene.tsv.gz",
                   clin="00_raw_data/bulk/GSE71729/GSE71729.clinical_std.tsv"),
  GSE21501  = list(expr="00_raw_data/bulk/GSE21501/GSE21501.gene.tsv.gz",
                   clin="00_raw_data/bulk/GSE21501/GSE21501.clinical_std.tsv"),
  GSE85916  = list(expr="00_raw_data/bulk/GSE85916/GSE85916.gene.tsv.gz",
                   clin="00_raw_data/bulk/GSE85916/GSE85916.clinical_std.tsv"),
  GSE62452  = list(expr="00_raw_data/bulk/GSE62452/GSE62452.gene.tsv.gz",
                   clin="00_raw_data/bulk/GSE62452/GSE62452.clinical_std.tsv"),
  GSE57495  = list(expr="00_raw_data/bulk/GSE57495/GSE57495.gene.tsv.gz",
                   clin="00_raw_data/bulk/GSE57495/GSE57495.clinical_std.tsv"),
  TCGA_PAAD = list(expr="00_raw_data/bulk/TCGA_PAAD/TCGA-PAAD.curated.gene.tsv.gz",
                   clin="00_raw_data/bulk/TCGA_PAAD/TCGA-PAAD.curated.clinical_std.tsv")
)

# C-index 计算函数: signature 基因×系数 → risk → C-index
calc_cindex <- function(expr, clin, genes, coefs) {
  common <- intersect(genes, rownames(expr))
  if (length(common) < 2) return(list(cidx=NA, n_gene=length(common)))
  # risk = Σ z-score(基因) × 系数
  em <- expr[common,,drop=FALSE]
  em_z <- t(scale(t(em)))  # 每基因 z-score
  em_z[is.na(em_z)] <- 0
  risk <- colSums(em_z * coefs[common])
  # 对齐生存 (真实列名: sample_id / os_time_month / os_status)
  s <- clin[match(colnames(expr), clin$sample_id)]
  ok <- !is.na(s$os_time_month) & !is.na(s$os_status) & s$os_time_month>0
  if (sum(ok)<10) return(list(cidx=NA, n_gene=length(common)))
  ci <- tryCatch(
    concordance.index(x=risk[ok], surv.time=s$os_time_month[ok], surv.event=s$os_status[ok],
                      method="noether")$c.index,
    error=function(e) NA)
  list(cidx=ci, n_gene=length(common))
}

# 4. 加载队列 + 算所有 signature C-index
results <- list()
for (cn in names(cohorts)) {
  cat(sprintf("\n--- %s ---\n", cn))
  ef <- file.path(ROOT, cohorts[[cn]]$expr); cf <- file.path(ROOT, cohorts[[cn]]$clin)
  if (!file.exists(ef)||!file.exists(cf)) { cat("文件缺, 跳过\n"); next }
  expr <- as.data.frame(fread(ef)); rownames(expr) <- expr[[1]]; expr[[1]] <- NULL
  expr <- as.matrix(expr)
  clin <- fread(cf)
  cat(sprintf("  %d 基因 x %d 样本, 临床 %d 行\n", nrow(expr), ncol(expr), nrow(clin)))

  res_c <- data.table()
  # 86 signature
  for (m in models) {
    g <- sig[Model==m]$Gene; co <- setNames(sig[Model==m]$Coef, sig[Model==m]$Gene)
    r <- calc_cindex(expr, clin, g, co)
    res_c <- rbind(res_c, data.table(Model=m, cohort=cn, cindex=r$cidx, n_gene_used=r$n_gene,
                                     n_gene_total=length(g),
                                     Author=sig[Model==m]$Author[1], Year=sig[Model==m]$Year[1]))
  }
  # NeFRS
  rn <- calc_cindex(expr, clin, nefrs_genes, nefrs_coef)
  res_c <- rbind(res_c, data.table(Model="NeFRS", cohort=cn, cindex=rn$cidx, n_gene_used=rn$n_gene,
                                   n_gene_total=length(nefrs_genes), Author="Ours", Year=2026))
  # C-index<0.5 翻转 (signature 方向可能相反, AIDPS 也这么处理)
  res_c[, cindex_adj := ifelse(!is.na(cindex) & cindex<0.5, 1-cindex, cindex)]
  results[[cn]] <- res_c
  nf <- res_c[Model=="NeFRS"]$cindex_adj
  rank_nf <- sum(res_c$cindex_adj > nf, na.rm=TRUE)+1
  cat(sprintf("  NeFRS C-index=%.3f, 排名 %d/%d\n", nf, rank_nf, sum(!is.na(res_c$cindex_adj))))
}

all_res <- rbindlist(results)
fwrite(all_res, file.path(OUT,"sig_compare_all_cohorts.tsv"), sep="\t")

# 5. 汇总: NeFRS 各队列排名
cat("\n=== NeFRS 各队列排名汇总 ===\n")
summ <- data.table()
for (cn in names(results)) {
  rc <- results[[cn]]
  nf <- rc[Model=="NeFRS"]$cindex_adj
  valid <- rc[!is.na(cindex_adj)]
  rank_nf <- sum(valid$cindex_adj > nf)+1
  summ <- rbind(summ, data.table(cohort=cn, NeFRS_cindex=round(nf,3),
                                 rank=rank_nf, total=nrow(valid),
                                 percentile=round(100*(1-rank_nf/nrow(valid)),1)))
}
print(summ)
fwrite(summ, file.path(OUT,"nefrs_rank_summary.tsv"), sep="\t")
cat(sprintf("\nNeFRS 平均排名百分位: %.1f%%\n", mean(summ$percentile)))
cat("\n=== 完成 ===\n")
