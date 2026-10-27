# TF 调控网络: decoupleR + CollecTRI, 恶性细胞 TF 活性 + NeFRS 高低对比 + KO 交叉验证
suppressPackageStartupMessages({
  library(decoupleR); library(data.table); library(Matrix)
  library(dplyr); library(tidyr); library(tibble); library(ggplot2)
  library(pheatmap); library(RColorBrewer)
})
options(bitmapType="cairo")
SC  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
KO  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/10_knockout"
TF  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/11_tf"
FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures/Fig10_tf"
dir.create(TF, recursive=TRUE, showWarnings=FALSE); dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

cat("=== TF 调控网络分析 (decoupleR + CollecTRI) ===\n")
t_total <- Sys.time()

# 1. CollecTRI
collectri <- fread("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/collectri_human.csv")
cat(sprintf("CollecTRI: %d 条, %d TFs\n", nrow(collectri), length(unique(collectri$source))))

# 2. 恶性细胞表达 (log-norm, run_ulm 用 log-norm)
X <- readRDS(file.path(SC,"w4_X_sparse_full.rds"))   # cells x genes
annot <- as.data.table(readRDS(file.path(SC,"w4_annot_with_ucell.rds")))
uc <- fread(file.path(SC,"w4_nefrs_ucell_allcells.tsv"))

mal <- annot[Level2=="Malignant"]$barcode
idx <- match(mal, rownames(X)); idx <- idx[!is.na(idx)]
Xm <- X[idx,]
uc_m <- uc[match(rownames(Xm), barcode)]
cat(sprintf("恶性细胞: %d\n", nrow(Xm)))
rm(X); gc()

# run_ulm 要 genes x cells
expr_mat <- t(as.matrix(Xm))   # genes x cells
cat(sprintf("表达矩阵: %d genes x %d cells\n", nrow(expr_mat), ncol(expr_mat)))
rm(Xm); gc()

# 3. run_ulm 计算 TF 活性
cat("\n[run_ulm] 计算 TF 活性 (Holland 2020 benchmark 第一)...\n")
t1 <- Sys.time()
tf_acts <- run_ulm(mat=expr_mat, net=collectri,
                   .source="source", .target="target", .mor="mor", minsize=5)
cat(sprintf("  %.1f min, %d 行\n", as.numeric(difftime(Sys.time(),t1,units="mins")), nrow(tf_acts)))
rm(expr_mat); gc()

# 4. TF x cell 矩阵
tf_mat <- tf_acts %>% filter(statistic=="ulm") %>%
  select(source, condition, score) %>%
  pivot_wider(names_from=condition, values_from=score) %>%
  column_to_rownames("source") %>% as.matrix()
cat(sprintf("TF x cell: %d x %d\n", nrow(tf_mat), ncol(tf_mat)))
fwrite(as.data.frame(tf_mat), file.path(TF,"TF_activity_cell.csv.gz"), row.names=TRUE)

# 5. NeFRS 高低分组 (UCell)
nefrs <- uc_m$nefrs_ucell
q33 <- quantile(nefrs,1/3,na.rm=T); q67 <- quantile(nefrs,2/3,na.rm=T)
grp <- ifelse(nefrs>=q67,"High",ifelse(nefrs<=q33,"Low","Mid"))
names(grp) <- uc_m$barcode

# 对齐 tf_mat 列顺序
common <- intersect(colnames(tf_mat), names(grp))
tf_mat <- tf_mat[, common]; grp_c <- grp[common]
hi_cells <- names(grp_c)[grp_c=="High"]; lo_cells <- names(grp_c)[grp_c=="Low"]
cat(sprintf("\nNeFRS High %d / Low %d\n", length(hi_cells), length(lo_cells)))

# 6. TF 活性 High vs Low 差异
cat("\n=== TF 活性: High vs Low NeFRS 恶性细胞 ===\n")
res <- data.table()
for (tf in rownames(tf_mat)) {
  h <- tf_mat[tf, hi_cells]; l <- tf_mat[tf, lo_cells]
  h <- h[!is.na(h)]; l <- l[!is.na(l)]
  if(length(h)<10||length(l)<10) next
  w <- wilcox.test(h, l)
  res <- rbind(res, data.table(TF=tf, high_mean=mean(h), low_mean=mean(l),
                               diff=mean(h)-mean(l), pval=w$p.value))
}
res[, padj := p.adjust(pval,"BH")]
setorder(res, -diff)
cat("Top 20 高 NeFRS 激活的 TF:\n")
print(head(res[,.(TF,high_mean,low_mean,diff,padj)],20))
cat("\nTop 10 高 NeFRS 抑制的 TF:\n")
print(head(res[order(diff),.(TF,diff,padj)],10))
fwrite(res, file.path(TF,"TF_highlow_diff.tsv"), sep="\t")

# 7. 与虚拟敲除交叉验证: KO hub 的 CollecTRI targets ∩ KO 显著基因
cat("\n=== TF-KO 交叉验证 ===\n")
ko_overlap <- data.table()
for (g in c("NRP1","EGFR","BMPR1B")) {
  ko_f <- file.path(KO,"Results_final",sprintf("%s_sig_adjP005.csv",g))
  if(!file.exists(ko_f)){cat(sprintf("%s KO 文件缺\n",g));next}
  ko_genes <- setdiff(fread(ko_f)$gene, g)
  ct_targets <- setdiff(collectri[source==g]$target, g)
  ov <- intersect(ko_genes, ct_targets)
  cat(sprintf("%s: KO显著 %d, CollecTRI靶 %d, overlap %d %s\n",
              g, length(ko_genes), length(ct_targets), length(ov),
              ifelse(length(ov)>0, paste0("(",paste(ov,collapse=","),")"), "")))
  # 注: NRP1/EGFR/BMPR1B 是受体不是TF, CollecTRI 里可能没有它们作为 source
  ko_overlap <- rbind(ko_overlap, data.table(KO=g, ko_sig=length(ko_genes),
                      ct_targets=length(ct_targets), overlap=length(ov)))
}
fwrite(ko_overlap, file.path(TF,"TF_KO_overlap.tsv"), sep="\t")

# 8. Top TF heatmap (by 高低组均值)
top_tf <- head(res[order(-abs(diff))]$TF, 30)
hm <- cbind(High=rowMeans(tf_mat[top_tf,hi_cells],na.rm=T),
            Low=rowMeans(tf_mat[top_tf,lo_cells],na.rm=T))
hm_s <- t(scale(t(hm))); hm_s[is.na(hm_s)] <- 0
pheatmap(hm_s, color=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
         cluster_cols=FALSE, fontsize_row=9,
         main="Fig10A. Top 30 differential TF activity (High vs Low NeFRS)",
         filename=file.path(FIG,"Fig10A_TF_heatmap.pdf"), width=5, height=9)
cat("Fig10A 保存\n")

# 9. 关键 NeFRS TF 对比图
key_tf <- c("NFIB","TRPS1","REST","SOX4","NEUROD1","GLIS3","WWTR1","NR2F1","ZEB1","INSM1")
key_tf <- intersect(key_tf, res$TF)
kdf <- res[TF %in% key_tf]
kdf[, TF := factor(TF, levels=TF[order(diff)])]
p <- ggplot(kdf, aes(TF, diff, fill=diff>0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values=c("TRUE"="#C0392B","FALSE"="#2980B9"), guide="none") +
  geom_hline(yintercept=0, lty=2) +
  labs(title="Fig10B. Key NeFRS-related TF activity (High - Low)",
       x="", y="TF activity diff (High - Low NeFRS)") +
  theme_bw(base_size=11)
ggsave(file.path(FIG,"Fig10B_keyTF_diff.pdf"), p, width=7, height=5)
cat("Fig10B 保存\n")

cat(sprintf("\n=== 完成, %.1f min ===\n", as.numeric(difftime(Sys.time(),t_total,units="mins"))))
