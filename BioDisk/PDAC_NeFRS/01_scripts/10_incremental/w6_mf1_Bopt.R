#!/usr/bin/env Rscript
suppressMessages({ library(data.table); library(glmnet); library(ggplot2) })
fs <- "01_scripts/00_style/fig_style.R"; if(!file.exists(fs)) fs <- list.files(".","^fig_style\\.R$",recursive=TRUE,full.names=TRUE)[1]
source(fs); OUT <- "03_figures/panels"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
res <- readRDS("02_results/04_mime/mime_results.rds"); mlnms <- names(res$ml.res)
target <- mlnms[grepl("^Enet",mlnms) & grepl("0\\.1",mlnms) & !grepl("[+]",mlnms)][1]; ml <- res$ml.res[[target]]
fg <- colnames(readRDS("02_results/04_mime/mime_input_list.rds")[[1]])[-(1:3)]
score <- function(e){ X <- matrix(0,ncol(e),length(fg),dimnames=list(colnames(e),fg)); cm<-intersect(rownames(e),fg); X[,cm]<-t(as.matrix(e[cm,,drop=FALSE])); as.numeric(predict(ml,newx=X,s="lambda.min")) }
fgc <- function(dt){ nm<-names(dt); p<-nm[tolower(nm)%in%c("gene","symbol","genes")]; if(length(p)) dt[[p[1]]] else dt[[1]] }
gp <- function(p){ if(file.exists(p)) as.character(fgc(fread(p))) else character(0) }

## ===== 选项1：神经轴 × 铁死亡轴 耦合散点（各取独有基因，诚实 r）=====
tryCatch({
  ferrdb  <- gp("02_results/01_gene_pool/01_ferrdb_v2_genes.tsv")
  neural  <- unique(c(gp("02_results/01_gene_pool/03_neural_dev_genes.tsv"),
                      gp("02_results/01_gene_pool/02_nrp_signature_genes.tsv"),
                      gp("02_results/01_gene_pool/04_pni_literature_genes.tsv")))
  cf <- fread("02_results/04_mime/nefrs_v3_enet_coefficients.tsv"); nefrs <- cf[cf$coef!=0,][[1]]
  co <- readRDS("02_results/02_w2_inputs/cohort_data_8.rds"); expr <- as.matrix(co[["PACA_AU"]]$expression)
  Z <- t(scale(t(expr)))
  nz_neu <- intersect(nefrs[ nefrs %in% neural & !(nefrs %in% ferrdb) ], rownames(Z))
  nz_fer <- intersect(nefrs[ nefrs %in% ferrdb & !(nefrs %in% neural) ], rownames(Z))
  cat(sprintf(">>> 神经轴独有=%d, 铁死亡轴独有=%d\n", length(nz_neu), length(nz_fer)))
  neu_sc <- colMeans(Z[nz_neu,,drop=FALSE]); fer_sc <- colMeans(Z[nz_fer,,drop=FALSE])
  s <- score(expr); risk <- factor(ifelse(s>median(s),"High","Low"), levels=c("Low","High"))
  ct <- cor.test(neu_sc, fer_sc)
  cat(sprintf(">>> 耦合 Pearson r=%.3f, p=%.2e\n", ct$estimate, ct$p.value))
  df <- data.frame(neural=neu_sc, ferro=fer_sc, Risk=risk)
  pSC <- ggplot(df, aes(neural, ferro)) +
    geom_point(aes(color=Risk), size=1.1, alpha=.75) +
    geom_smooth(method="lm", se=TRUE, color="black", linewidth=.5) +
    scale_color_manual(values=c(Low=COL_LOW, High=COL_HIGH)) +
    annotate("text", x=-Inf, y=Inf, hjust=-0.08, vjust=1.4, size=3.6, fontface="bold",
             label=sprintf("r = %.2f, p = %.1e", ct$estimate, ct$p.value)) +
    labs(x="Neural-axis score (z)", y="Ferroptosis-axis score (z)") + theme_pub(base=12) +
    theme(axis.text=element_text(face="bold"), axis.title=element_text(face="bold"),
          legend.position=c(.85,.16), legend.title=element_text(face="bold"))
  save_panel(pSC, file.path(OUT,"MF1_Bopt1_coupling"), w_mm=84, h_mm=72); cat(">>> [Bopt1] saved\n")
}, error=function(e) cat("!! [Bopt1]", conditionMessage(e), "\n"))

## ===== 选项2：Enet 交叉验证曲线 =====
tryCatch({
  cvdf <- data.frame(loglam=log(ml$lambda), cvm=ml$cvm, lo=ml$cvlo, up=ml$cvup, nz=ml$nzero)
  lmin <- log(ml$lambda.min); l1se <- log(ml$lambda.1se)
  ylab_txt <- if(!is.null(ml$name)) as.character(ml$name) else "CV measure"
  pCV <- ggplot(cvdf, aes(loglam, cvm)) +
    geom_ribbon(aes(ymin=lo, ymax=up), fill="grey80", alpha=.6) +
    geom_line(color="grey50", linewidth=.3) +
    geom_point(color=COL_HIGH, size=1.1) +
    geom_vline(xintercept=lmin, linetype=2, color="grey25") +
    geom_vline(xintercept=l1se, linetype=3, color="grey45") +
    annotate("text", x=lmin, y=max(cvdf$up), label="lambda.min", angle=90, vjust=-0.4, hjust=1, size=3.2, fontface="bold") +
    annotate("text", x=l1se, y=max(cvdf$up), label="lambda.1se", angle=90, vjust=-0.4, hjust=1, size=3.2, fontface="bold", color="grey45") +
    labs(x="log(lambda)", y=ylab_txt) + theme_pub(base=12) +
    theme(axis.text=element_text(face="bold"), axis.title=element_text(face="bold"))
  save_panel(pCV, file.path(OUT,"MF1_Bopt2_cvcurve"), w_mm=84, h_mm=68); cat(">>> [Bopt2] saved\n")
}, error=function(e) cat("!! [Bopt2]", conditionMessage(e), "\n"))
cat(">>> 完成 Bopt1 / Bopt2\n")
