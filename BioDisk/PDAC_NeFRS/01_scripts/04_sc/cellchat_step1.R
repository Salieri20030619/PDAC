# CellChat 细胞通讯: A 全细胞类型 + B 恶性高/低 NeFRS 对比
suppressPackageStartupMessages({
  library(Seurat); library(CellChat); library(data.table); library(Matrix)
  library(ggplot2); library(dplyr); library(patchwork)
})
options(bitmapType="cairo")
SC  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/06_sc"
CC  <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/12_cellchat"
FIG <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/03_figures/Fig11_cellchat"
dir.create(CC, recursive=TRUE, showWarnings=FALSE); dir.create(FIG, recursive=TRUE, showWarnings=FALSE)
future::plan("multisession", workers=4)
options(future.globals.maxSize=40*1024^3)

cat("=== CellChat 细胞通讯 ===\n"); t0 <- Sys.time()
X <- readRDS(file.path(SC,"w4_X_sparse_full.rds"))   # cells x genes, log-norm
annot <- as.data.table(readRDS(file.path(SC,"w4_annot_with_ucell.rds")))
uc <- fread(file.path(SC,"w4_nefrs_ucell_allcells.tsv"))

# 对齐
annot <- annot[match(rownames(X), barcode)]
data.input <- Matrix::t(X)   # genes x cells, 保持稀疏 dgCMatrix
data.input <- as(data.input, "CsparseMatrix")
# ★ 关键修复: h5ad 提取的基因名是 array 类型, CellChat row.names() 取出当索引会报 S4 错
# 强制 dimnames 为纯 character
rownames(data.input) <- as.character(rownames(data.input))
colnames(data.input) <- as.character(colnames(data.input))
cat(sprintf("表达: %d genes x %d cells, rownames class %s\n", nrow(data.input), ncol(data.input), class(rownames(data.input))))
rm(X); gc()

CellChatDB <- CellChatDB.human
cat(sprintf("CellChatDB.human: %d 配受体对\n", nrow(CellChatDB$interaction)))

run_cellchat <- function(data.input, labels, tag) {
  cat(sprintf("\n===== CellChat: %s =====\n", tag))
  t1 <- Sys.time()
  # 用显式 normalizeData 后的稀疏矩阵 + meta, CellChat 2.2 要 data 是 dgCMatrix
  meta <- data.frame(labels=factor(labels), samples=factor("sample1"),
                     row.names=as.character(colnames(data.input)), stringsAsFactors=FALSE)
  cc <- createCellChat(object=data.input, meta=meta, group.by="labels")
  cc@DB <- CellChatDB
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type="triMean")
  cc <- filterCommunication(cc, min.cells=10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cat(sprintf("  %.1f min, 通讯数 %d, 通路数 %d\n",
              as.numeric(difftime(Sys.time(),t1,units="mins")),
              sum(cc@net$count), length(cc@netP$pathways)))
  cc
}

# ===== A: 全细胞类型通讯 =====
cc_A <- run_cellchat(data.input, annot$Level1, "A_all_celltypes")
saveRDS(cc_A, file.path(CC,"cellchat_A_allcelltypes.rds"))
cat("A 通路:", paste(head(cc_A@netP$pathways,30),collapse=", "), "\n")

# 神经相关通路有哪些
neuro_path <- grep("NRP|SEMA|EPH|NCAM|NGF|NT|NRXN|NLGN|SLIT|NEGR|CADM|GDNF|BDNF",
                   cc_A@netP$pathways, value=TRUE)
cat("神经相关通路:", paste(neuro_path,collapse=", "), "\n")

# ===== B: 恶性拆高/低 NeFRS =====
nefrs <- uc[match(annot$barcode, barcode)]$nefrs_ucell
mal_idx <- which(annot$Level1=="Epithelial (malignant)")
mal_nefrs <- nefrs[mal_idx]
q33 <- quantile(mal_nefrs,1/3,na.rm=T); q67 <- quantile(mal_nefrs,2/3,na.rm=T)
labels_B <- as.character(annot$Level1)
labels_B[mal_idx] <- ifelse(mal_nefrs>=q67, "Malignant_HighNeFRS",
                      ifelse(mal_nefrs<=q33, "Malignant_LowNeFRS", "Malignant_Mid"))
# 去掉 Mid (只留高低对比 + 其他细胞)
keep_B <- labels_B != "Malignant_Mid"
cat(sprintf("\nB: 高 %d, 低 %d\n", sum(labels_B=="Malignant_HighNeFRS"), sum(labels_B=="Malignant_LowNeFRS")))
cc_B <- run_cellchat(data.input[,keep_B], labels_B[keep_B], "B_NeFRS_highlow")
saveRDS(cc_B, file.path(CC,"cellchat_B_NeFRS_highlow.rds"))

# ===== 图 =====
# Fig11A: 全细胞类型通讯网络
vw <- as.numeric(table(cc_A@meta$labels)); names(vw) <- names(table(cc_A@meta$labels))
pdf(file.path(FIG,"Fig11A_network_count.pdf"), width=11, height=11)
par(mar=c(1,1,3,1))
netVisual_circle(cc_A@net$count, vertex.weight=vw, weight.scale=TRUE, label.edge=FALSE,
                 title.name="Fig11A. Cell-cell interactions (count)")
dev.off()
pdf(file.path(FIG,"Fig11A_network_weight.pdf"), width=11, height=11)
par(mar=c(1,1,3,1))
netVisual_circle(cc_A@net$weight, vertex.weight=vw, weight.scale=TRUE, label.edge=FALSE,
                 title.name="Fig11A. Interaction strength")
dev.off()
cat("Fig11A 保存\n")

# Fig11B: 神经通路气泡图 (恶性作为 source/target)
if (length(neuro_path)>0) {
  ct_levels <- levels(cc_A@meta$labels)
  mal_i <- which(ct_levels=="Epithelial (malignant)")
  tryCatch({
    pb <- netVisual_bubble(cc_A, sources.use=1:length(ct_levels), targets.use=mal_i,
                           signaling=neuro_path, remove.isolate=TRUE,
                           title.name="Fig11B. Neural signaling to malignant cells")
    ggsave(file.path(FIG,"Fig11B_neural_bubble.pdf"), pb, width=10, height=8)
    cat("Fig11B 保存\n")
  }, error=function(e) cat("Fig11B 失败:",conditionMessage(e),"\n"))
}

# 汇总
summ <- data.frame(
  analysis=c("A_all","B_highlow"),
  n_celltypes=c(length(unique(cc_A@meta$labels)), length(unique(cc_B@meta$labels))),
  n_interactions=c(sum(cc_A@net$count), sum(cc_B@net$count)),
  n_pathways=c(length(cc_A@netP$pathways), length(cc_B@netP$pathways)))
fwrite(summ, file.path(CC,"cellchat_summary.tsv"), sep="\t")
cat("\n=== 汇总 ===\n"); print(summ)
cat(sprintf("\n=== 完成 %.1f min ===\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
future::plan("sequential")  # 关闭 workers
