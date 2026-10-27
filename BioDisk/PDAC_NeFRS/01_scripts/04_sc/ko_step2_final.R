# 虚拟敲除最终版: 真 counts + 全部 64538 恶性细胞 + 3 hub
suppressPackageStartupMessages({
  library(scTenifoldKnk); library(Matrix); library(data.table); library(ggplot2); library(ggrepel)
})
options(bitmapType="cairo")
KO  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/10_knockout"
RES <- file.path(KO,"Results_final"); FIG <- file.path(KO,"Figures_final")
dir.create(RES, showWarnings=FALSE); dir.create(FIG, showWarnings=FALSE)

counts <- readRDS(file.path(KO,"malignant_real_counts.rds"))
cat(sprintf("真 counts + 全部恶性: %d cells x %d genes\n", nrow(counts), ncol(counts)))

gene_var <- apply(counts, 2, var)
hvg <- names(sort(gene_var, decreasing=TRUE))[1:3000]

process_ko <- function(g, res, RES, FIG) {
  if (is.null(res)) return(NULL)
  dr <- res$diffRegulation
  write.csv(dr, file.path(RES, sprintf("%s_diffRegulation.csv",g)), row.names=FALSE)
  sig_p <- dr[dr$p.value<0.05,]; sig_adj <- dr[dr$p.adj<0.05,]
  write.csv(sig_adj, file.path(RES, sprintf("%s_sig_adjP005.csv",g)), row.names=FALSE)
  cat(sprintf("%s: P<0.05 %d, adjP<0.05 %d\n", g, nrow(sig_p), nrow(sig_adj)))
  top20 <- head(dr[order(-dr$FC),][dr$gene!=g,],20)
  write.csv(top20, file.path(RES, sprintf("%s_top20_byFC.csv",g)), row.names=FALSE)
  cat(sprintf("%s top10:\n", g)); print(head(top20[,c("gene","Z","FC","p.adj")],10))
  df <- dr; df$log_pval <- -log10(df$p.value+1e-300)
  lab <- head(df[order(-abs(df$Z)) & df$p.value<0.05,],25)
  df$cat <- "ns"; df$cat[df$p.value<0.05 & abs(df$Z)>2] <- "sig"; df$cat[df$p.adj<0.05] <- "highly sig"
  p <- ggplot(df, aes(Z,log_pval,color=cat)) + geom_point(alpha=0.6,size=1) +
    scale_color_manual(values=c("ns"="gray70","sig"="#3182bd","highly sig"="#e41a1c")) +
    geom_hline(yintercept=-log10(0.05),lty=2,color="gray50") + geom_vline(xintercept=c(-2,2),lty=2,color="gray50") +
    geom_text_repel(data=lab,aes(label=gene),size=2.5,color="black",max.overlaps=30) +
    labs(title=sprintf("Fig9. %s knockout (real counts, all 64538 malignant cells)",g),
         subtitle=sprintf("Total %d | P<0.05 %d | adjP<0.05 %d",nrow(df),nrow(sig_p),nrow(sig_adj)),
         x="Z-score",y="-log10(p)",color="") + theme_classic(base_size=11)+theme(legend.position="bottom")
  ggsave(file.path(FIG,sprintf("Fig9_%s_KO_volcano.pdf",g)),p,width=9,height=7)
  cat(sprintf("Fig9_%s 保存\n",g))
  data.frame(gene=g,n_p005=nrow(sig_p),n_adjp005=nrow(sig_adj),top1=top20$gene[1])
}

genes_ko <- c("NRP1","EGFR","BMPR1B")
summary_list <- list(); res_all <- list()
for (gene in genes_ko) {
  cat(sprintf("\n===== KO %s =====\n", gene))
  if (!gene %in% colnames(counts)) { cat("不在数据\n"); next }
  expr_pct <- 100*sum(counts[,gene]>0)/nrow(counts)
  cat(sprintf("%s 恶性表达率: %.1f%%\n", gene, expr_pct))
  genes_use <- unique(c(gene,hvg)); genes_use <- genes_use[genes_use %in% colnames(counts)]
  cm <- as.matrix(t(counts[,genes_use]))
  if (which(rownames(cm)==gene)!=1) cm <- cm[c(gene,setdiff(rownames(cm),gene)),]
  cat(sprintf("矩阵: %d x %d\n", nrow(cm), ncol(cm)))
  t0 <- Sys.time()
  res <- tryCatch(scTenifoldKnk(countMatrix=cm,gKO=gene,qc=FALSE),error=function(e){cat("FAIL:",conditionMessage(e),"\n");NULL})
  if (!is.null(res)) { cat(sprintf("%s %.1f min\n",gene,as.numeric(difftime(Sys.time(),t0,units="mins")))); saveRDS(res,file.path(RES,sprintf("scTenifoldKnk_%s.rds",gene))); res_all[[gene]]<-res }
  summary_list[[gene]] <- process_ko(gene,res,RES,FIG); rm(cm); gc()
}
summary_df <- do.call(rbind,summary_list)
write.csv(summary_df,file.path(RES,"all_KO_summary.csv"),row.names=FALSE)

# 共同下游分析 (Top50 distance 重叠)
if (length(res_all)==3) {
  top50 <- lapply(names(res_all), function(gn){d<-res_all[[gn]]$diffRegulation; d<-d[d$gene!=gn,]; head(d[order(-d$distance),]$gene,50)})
  names(top50) <- names(res_all)
  cat("\n=== Top50 重叠 ===\n")
  cat(sprintf("NRP1∩EGFR: %d\n", length(intersect(top50$NRP1,top50$EGFR))))
  cat(sprintf("NRP1∩BMPR1B: %d\n", length(intersect(top50$NRP1,top50$BMPR1B))))
  cat(sprintf("EGFR∩BMPR1B: %d\n", length(intersect(top50$EGFR,top50$BMPR1B))))
  cat(sprintf("三者共有: %d\n", length(Reduce(intersect,top50))))
  cat("NRP1∩EGFR 共同下游:", paste(intersect(top50$NRP1,top50$EGFR),collapse=", "), "\n")
  cat("三者共有:", paste(Reduce(intersect,top50),collapse=", "), "\n")
  ov <- data.table(pair=c("NRP1_EGFR","NRP1_BMPR1B","EGFR_BMPR1B","all3"),
                   overlap=c(length(intersect(top50$NRP1,top50$EGFR)),length(intersect(top50$NRP1,top50$BMPR1B)),
                             length(intersect(top50$EGFR,top50$BMPR1B)),length(Reduce(intersect,top50))))
  fwrite(ov, file.path(RES,"top50_overlap.tsv"), sep="\t")
}
cat("\n=== 汇总 ===\n"); print(summary_df)
cat("\n=== 全部完成 ===\n")
