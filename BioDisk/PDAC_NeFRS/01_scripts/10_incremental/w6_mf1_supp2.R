suppressMessages({ library(data.table); library(glmnet); library(ComplexHeatmap); library(circlize) })
fs <- "01_scripts/00_style/fig_style.R"; if(!file.exists(fs)) fs <- list.files(".","^fig_style\\.R$",recursive=TRUE,full.names=TRUE)[1]
source(fs); OUT <- "03_figures/panels"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
res <- readRDS("02_results/04_mime/mime_results.rds"); mlnms <- names(res$ml.res)
target <- mlnms[grepl("^Enet",mlnms) & grepl("0\\.1",mlnms) & !grepl("[+]",mlnms)][1]; ml <- res$ml.res[[target]]
fg <- colnames(readRDS("02_results/04_mime/mime_input_list.rds")[[1]])[-(1:3)]
score <- function(e){ X <- matrix(0,ncol(e),length(fg),dimnames=list(colnames(e),fg)); cm<-intersect(rownames(e),fg); X[,cm]<-t(as.matrix(e[cm,,drop=FALSE])); as.numeric(predict(ml,newx=X,s="lambda.min")) }
co <- readRDS("02_results/02_w2_inputs/cohort_data_8.rds"); expr <- as.matrix(co[["PACA_AU"]]$expression)
s <- score(expr); names(s) <- colnames(expr); ord_s <- order(s); risk <- ifelse(s>median(s),"High","Low")[ord_s]
col_fun <- colorRamp2(c(-2,0,2), c(COL_LOW,"white",COL_HIGH))
top_ann <- HeatmapAnnotation(NeFRS=anno_barplot(s[ord_s], gp=gpar(fill=COL_HIGH,col=NA), height=unit(7,"mm")),
  Risk=factor(risk,levels=c("Low","High")), col=list(Risk=c(Low=COL_LOW,High=COL_HIGH)),
  annotation_name_gp=gpar(fontsize=6.5, fontface="bold"), simple_anno_size=unit(3,"mm"))
fgc <- function(dt){ nm<-names(dt); p<-nm[tolower(nm)%in%c("gene","symbol","genes")]; if(length(p)) dt[[p[1]]] else dt[[1]] }

tryCatch({
  cm <- fread("02_results/03_funnel/step2_gse_consist_matrix.tsv")
  hrcol <- names(cm)[grepl("_HR$", names(cm))]
  D <- sign(as.matrix(cm[, hrcol, with=FALSE]) - 1); rownames(D) <- fgc(cm); colnames(D) <- sub("_HR$","",hrcol)
  ar <- if ("agree_ratio" %in% names(cm)) cm$agree_ratio else rowMeans(D==1)
  o <- order(-ar); D <- D[o,,drop=FALSE]; ar <- ar[o]
  ra <- rowAnnotation(agree=anno_barplot(ar, gp=gpar(fill="#009E73",col=NA), width=unit(10,"mm")), annotation_name_gp=gpar(fontsize=6.5, fontface="bold"))
  htD <- Heatmap(D, name="dir", col=colorRamp2(c(-1,0,1), c(COL_LOW,"white",COL_HIGH)),
    cluster_rows=FALSE, cluster_columns=FALSE, show_row_names=FALSE, right_annotation=ra,
    column_names_gp=gpar(fontsize=6, fontface="bold"), column_names_rot=45,
    row_title="candidate genes", row_title_gp=gpar(fontsize=7.5, fontface="bold"),
    heatmap_legend_param=list(at=c(-1,1), labels=c("protective","risk"), labels_gp=gpar(fontsize=6, fontface="bold"), title_gp=gpar(fontsize=6.5, fontface="bold")))
  save_base(function() draw(htD), file.path(OUT,"MF1_D_consist"), w_mm=58, h_mm=74); cat(">>> [D] saved\n")
}, error=function(e) cat("!! [D]", conditionMessage(e), "\n"))

tryCatch({
  ci <- fread("02_results/04_mime/mime_all_cindex.tsv"); modcol <- names(ci)[1]
  coh <- setdiff(names(ci)[sapply(ci,is.numeric)], "mean_cindex"); train <- "PACA_AU"; val <- setdiff(coh, train)
  mat <- as.matrix(ci[, coh, with=FALSE]); rownames(mat) <- ci[[modcol]]
  m_all <- rowMeans(mat, na.rm=TRUE); m_val <- rowMeans(mat[, val, drop=FALSE], na.rm=TRUE)
  o <- order(-m_val)[seq_len(min(12, nrow(mat)))]; mat <- mat[o, c(val, train), drop=FALSE]; m_all <- m_all[o]; m_val <- m_val[o]
  sel <- grepl("^Enet",rownames(mat)) & grepl("0\\.1",rownames(mat)) & !grepl("[+]",rownames(mat))
  body_col <- colorRamp2(quantile(mat,c(0,.5,1),na.rm=TRUE), c("#FEE8C8","#FDBB84","#B30000"))
  cell2 <- function(M) function(j,i,x,y,w,h,f) grid.text(sprintf("%.2f",M[i,j]), x,y, gp=gpar(fontsize=6.5))
  hB <- Heatmap(mat, name="C-index", col=body_col, cluster_rows=FALSE, cluster_columns=FALSE,
    width=unit(9*ncol(mat),"mm"),
    row_names_side="left", row_names_gp=gpar(fontsize=7, col=ifelse(sel,"#C00000","black"), fontface=ifelse(sel,"bold","plain")),
    column_names_gp=gpar(fontsize=6.5, fontface="bold"), column_names_rot=90, cell_fun=cell2(mat),
    heatmap_legend_param=list(labels_gp=gpar(fontsize=6.5, fontface="bold"), title_gp=gpar(fontsize=7, fontface="bold")))
  Mall <- matrix(m_all,ncol=1,dimnames=list(NULL,"Mean all")); Mval <- matrix(m_val,ncol=1,dimnames=list(NULL,"Mean val"))
  hA <- Heatmap(Mall, col=colorRamp2(range(m_all),c("white","#41AB5D")), cluster_rows=FALSE, show_heatmap_legend=FALSE,
    width=unit(8,"mm"), column_names_gp=gpar(fontsize=6.5, fontface="bold"), column_names_rot=90, cell_fun=cell2(Mall))
  hV <- Heatmap(Mval, col=colorRamp2(range(m_val),c("white","#6A51A3")), cluster_rows=FALSE, show_heatmap_legend=FALSE,
    width=unit(8,"mm"), column_names_gp=gpar(fontsize=6.5, fontface="bold"), column_names_rot=90, cell_fun=cell2(Mval))
  save_base(function() draw(hB + hA + hV, merge_legend=TRUE), file.path(OUT,"MF1_E_mime"), w_mm=116, h_mm=90)
  cat(">>> [E] saved (top12, 字号加大)\n")
}, error=function(e) cat("!! [E]", conditionMessage(e), "\n"))

gp <- function(p){ if(file.exists(p)) as.character(fgc(fread(p))) else character(0) }
ferro_g  <- gp("02_results/01_gene_pool/01_ferrdb_v2_genes.tsv")
neural_g <- unique(c(gp("02_results/01_gene_pool/03_neural_dev_genes.tsv"),
                     gp("02_results/01_gene_pool/02_nrp_signature_genes.tsv"),
                     gp("02_results/01_gene_pool/04_pni_literature_genes.tsv")))
cf <- fread("02_results/04_mime/nefrs_v3_enet_coefficients.tsv"); cf <- cf[cf$coef!=0,]; cofv <- setNames(cf$coef, cf$gene)
top_axis <- function(genes){ g <- intersect(intersect(names(cofv),genes), rownames(expr)); g[order(-abs(cofv[g]))][seq_len(min(20,length(g)))] }
I_genes <- top_axis(neural_g); J_genes <- top_axis(ferro_g)
draw_axis <- function(genes, fname, rtitle){
  Mx <- t(scale(t(expr[genes, ord_s, drop=FALSE]))); Mx[Mx>2]<-2; Mx[Mx< -2]<- -2
  ht <- Heatmap(Mx, name="z", col=col_fun, top_annotation=top_ann, cluster_rows=TRUE, cluster_columns=FALSE,
    show_column_names=FALSE, show_row_names=TRUE, row_names_gp=gpar(fontsize=7, fontface="bold"),
    row_title=rtitle, row_title_gp=gpar(fontsize=7.5, fontface="bold"),
    column_split=factor(risk,levels=c("Low","High")), column_title_gp=gpar(fontsize=7.5, fontface="bold"),
    heatmap_legend_param=list(labels_gp=gpar(fontsize=6.5, fontface="bold"), title_gp=gpar(fontsize=7, fontface="bold")))
  save_base(function() draw(ht, heatmap_legend_side="right", annotation_legend_side="right", merge_legend=TRUE),
            file.path(OUT,fname), w_mm=86, h_mm=78)
}
draw_axis(I_genes, "MF1_I_neural_heat", "Neural-axis (top20)")
draw_axis(J_genes, "MF1_J_ferro_heat",  "Ferroptosis-axis (top20)")
cat(">>> D/E/I/J 完成\n")
