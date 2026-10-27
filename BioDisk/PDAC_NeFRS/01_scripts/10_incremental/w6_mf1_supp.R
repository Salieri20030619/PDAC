suppressMessages({ library(ggplot2); library(data.table); library(glmnet) })
fs_cands <- c("01_scripts/00_style/fig_style.R","./01_scripts/00_style/fig_style.R")
fs <- fs_cands[file.exists(fs_cands)][1]
if (is.na(fs)) { hit <- list.files(".", "^fig_style\\.R$", recursive=TRUE, full.names=TRUE); if (length(hit)) fs <- hit[1] }
source(fs); OUT <- "03_figures/panels"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
res   <- readRDS("02_results/04_mime/mime_results.rds"); mlnms <- names(res$ml.res)
target <- mlnms[grepl("^Enet", mlnms) & grepl("0\\.1", mlnms) & !grepl("[+]", mlnms)][1]; ml <- res$ml.res[[target]]
funnel_genes <- colnames(readRDS("02_results/04_mime/mime_input_list.rds")[[1]])[-(1:3)]
nefrs_score <- function(expr){
  X <- matrix(0, nrow=ncol(expr), ncol=length(funnel_genes), dimnames=list(colnames(expr), funnel_genes))
  common <- intersect(rownames(expr), funnel_genes); X[, common] <- t(as.matrix(expr[common, , drop=FALSE]))
  as.numeric(predict(ml, newx=X, s="lambda.min")) }
read_n <- function(p){ if (file.exists(p)) nrow(fread(p)) else NA_integer_ }
coef_tsv <- "02_results/04_mime/nefrs_v3_enet_coefficients.tsv"
n_nefrs  <- if (file.exists(coef_tsv)) { tmp <- fread(coef_tsv); sum(tmp$coef != 0, na.rm=TRUE) } else length(funnel_genes)
stages <- data.frame(
  stage = c("Candidate pool","In 8 cohorts","uniCox sig.","Funnel (enter ML)","NeFRS (Enet a0.1)"),
  n = c(read_n("02_results/01_gene_pool/final_candidate_pool.tsv"),
        read_n("02_results/01_gene_pool/final_in_all_cohorts.tsv"),
        read_n("02_results/03_funnel/step1_sig_genes.tsv"),
        length(funnel_genes), n_nefrs), stringsAsFactors=FALSE)
stages <- stages[!is.na(stages$n), ]
stages$lab <- paste0(stages$stage, " (n=", stages$n, ")"); stages$y <- rev(seq_len(nrow(stages)))
pA <- ggplot(stages) +
  geom_rect(aes(xmin=-n/2, xmax=n/2, ymin=y-0.42, ymax=y+0.42, fill=y)) +
  scale_fill_gradient(low="#9ECAE1", high=COL_HIGH, guide="none") +
  scale_y_continuous(breaks=stages$y, labels=stages$lab, expand=expansion(add=0.6)) +
  labs(x="Number of genes", y=NULL) + theme_pub(base=7) +
  theme(panel.grid=element_blank(), axis.ticks.y=element_blank(),
        axis.text=element_text(face="bold"), axis.title=element_text(face="bold"))
save_panel(pA, file.path(OUT,"MF1_funnel"), w_mm=64, h_mm=44)
co  <- readRDS("02_results/02_w2_inputs/cohort_data_8.rds")
ord <- c("PACA_AU","TCGA_PAAD","E_MTAB_6134","GSE71729","GSE85916","GSE62452","GSE57495","GSE21501")
ord <- ord[ord %in% names(co)]
dd  <- do.call(rbind, lapply(ord, function(k){ s <- nefrs_score(co[[k]]$expression); data.frame(cohort=k, z=as.numeric(scale(s))) }))
dd$cohort <- factor(dd$cohort, levels=ord)
pB <- ggplot(dd, aes(cohort, z)) +
  geom_violin(aes(fill=cohort=="PACA_AU"), colour=NA, alpha=.35, scale="width") +
  geom_boxplot(width=.18, outlier.size=.3, linewidth=.3, fill="white") +
  scale_fill_manual(values=c(`TRUE`=COL_HIGH, `FALSE`="#7F7F7F"), guide="none") +
  geom_hline(yintercept=0, linetype=2, linewidth=.3, colour="grey50") +
  labs(x=NULL, y="NeFRS score (within-cohort z)") + theme_pub(base=10) +
  theme(axis.text.x=element_text(angle=35, hjust=1, face="bold"),
        axis.text.y=element_text(face="bold"), axis.title=element_text(face="bold"))
save_panel(pB, file.path(OUT,"MF1_score_dist"), w_mm=116, h_mm=72)
cat(">>> C/H 完成(字号-2)\n")
