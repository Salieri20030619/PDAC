# scTenifoldKnk 正式: NRP1 / EGFR / BMPR1B 三个 NeFRS hub, 恶性细胞
suppressPackageStartupMessages({
  library(scTenifoldKnk); library(data.table); library(Matrix); library(ggplot2); library(ggrepel)
})
options(bitmapType="cairo")
SC  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
KO  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/10_knockout"
RES <- file.path(KO,"Results"); FIG <- file.path(KO,"Figures")
dir.create(RES, recursive=TRUE, showWarnings=FALSE); dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

cat("=== scTenifoldKnk 正式 (NRP1/EGFR/BMPR1B) ===\n")

X <- readRDS(file.path(SC,"w4_X_sparse_full.rds"))
annot <- as.data.table(readRDS(file.path(SC,"w4_annot_with_ucell.rds")))
mal <- annot[Level2=="Malignant"]$barcode
idx <- match(mal, rownames(X)); idx <- idx[!is.na(idx)]
Xm <- X[idx,]
cat(sprintf("恶性细胞: %d\n", nrow(Xm)))
rm(X); gc()

# HVG 5000 (方差)
gene_var <- apply(Xm, 2, var)
hvg <- names(sort(gene_var, decreasing=TRUE))[1:5000]

# process_ko_result (改自参考 r_19, 适配)
process_ko <- function(KO_gene, res, RES, FIG) {
  if (is.null(res)) { cat(sprintf("⚠ %s 失败\n", KO_gene)); return(NULL) }
  dr <- res$diffRegulation
  write.csv(dr, file.path(RES, sprintf("%s_diffRegulation.csv",KO_gene)), row.names=FALSE)
  sig_p <- dr[dr$p.value<0.05,]; sig_adj <- dr[dr$p.adj<0.05,]
  cat(sprintf("%s: P<0.05 %d, adjP<0.05 %d\n", KO_gene, nrow(sig_p), nrow(sig_adj)))
  write.csv(sig_adj, file.path(RES, sprintf("%s_significant_adjP005.csv",KO_gene)), row.names=FALSE)

  top20 <- head(dr[order(-dr$FC),],21); top20 <- head(top20[top20$gene!=KO_gene,],20)
  write.csv(top20, file.path(RES, sprintf("%s_top20_byFC.csv",KO_gene)), row.names=FALSE)
  cat(sprintf("%s Top10 扰动 (by FC):\n", KO_gene)); print(head(top20[,c("gene","Z","FC","p.adj")],10))

  # Fig: 火山图 (Z vs -log10 p)
  df <- dr; df$log_pval <- -log10(df$p.value + 1e-300)
  lab <- head(df[order(-abs(df$Z)) & df$p.value<0.05,], 25)
  df$cat <- "ns"
  df$cat[df$p.value<0.05 & abs(df$Z)>2] <- "sig (P<0.05,|Z|>2)"
  df$cat[df$p.adj<0.05] <- "highly sig (adjP<0.05)"
  p <- ggplot(df, aes(Z, log_pval, color=cat)) +
    geom_point(alpha=0.6, size=1) +
    scale_color_manual(values=c("ns"="gray70","sig (P<0.05,|Z|>2)"="#3182bd","highly sig (adjP<0.05)"="#e41a1c")) +
    geom_hline(yintercept=-log10(0.05), lty=2, color="gray50") +
    geom_vline(xintercept=c(-2,2), lty=2, color="gray50") +
    geom_text_repel(data=lab, aes(label=gene), size=2.5, color="black", max.overlaps=30) +
    labs(title=sprintf("Fig9. Perturbed genes after %s knockout (malignant cells)", KO_gene),
         subtitle=sprintf("Total %d | P<0.05 %d | adjP<0.05 %d", nrow(df),nrow(sig_p),nrow(sig_adj)),
         x="Z-score", y="-log10(p)", color="") +
    theme_classic(base_size=11) + theme(legend.position="bottom")
  ggsave(file.path(FIG, sprintf("Fig9_%s_KO_volcano.pdf",KO_gene)), p, width=9, height=7)
  cat(sprintf("Fig9_%s 保存\n", KO_gene))

  data.frame(gene=KO_gene, n_p005=nrow(sig_p), n_adjp005=nrow(sig_adj),
             top1_FC=top20$gene[1], stringsAsFactors=FALSE)
}

genes_ko <- c("NRP1","EGFR","BMPR1B")
summary_list <- list()
for (gene in genes_ko) {
  cat(sprintf("\n========== KO: %s ==========\n", gene))
  if (!gene %in% colnames(Xm)) { cat("基因不在数据\n"); next }
  expr_pct <- 100*sum(Xm[,gene]>0)/nrow(Xm)
  cat(sprintf("%s 恶性细胞表达率: %.1f%%\n", gene, expr_pct))

  genes_use <- unique(c(gene, hvg)); genes_use <- genes_use[genes_use %in% colnames(Xm)]
  cm <- as.matrix(t(round(expm1(Xm[,genes_use]))))  # genes x cells, 整数 counts
  if (which(rownames(cm)==gene)!=1) cm <- cm[c(gene,setdiff(rownames(cm),gene)),]
  cat(sprintf("矩阵: %d x %d\n", nrow(cm), ncol(cm)))

  t0 <- Sys.time()
  res <- tryCatch(scTenifoldKnk(countMatrix=cm, gKO=gene, qc=FALSE),
                  error=function(e){cat("FAIL:",conditionMessage(e),"\n"); NULL})
  if (!is.null(res)) {
    cat(sprintf("%s 完成 %.1f min\n", gene, as.numeric(difftime(Sys.time(),t0,units="mins"))))
    saveRDS(res, file.path(RES, sprintf("scTenifoldKnk_%s.rds",gene)))
  }
  summary_list[[gene]] <- process_ko(gene, res, RES, FIG)
  rm(cm); gc()
}

summary_df <- do.call(rbind, summary_list)
write.csv(summary_df, file.path(RES,"all_KO_summary.csv"), row.names=FALSE)
cat("\n=== 汇总 ===\n"); print(summary_df)
cat("\n=== 全部完成 ===\n")
