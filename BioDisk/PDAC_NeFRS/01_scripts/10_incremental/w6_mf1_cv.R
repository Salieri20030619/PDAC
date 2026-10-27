#!/usr/bin/env Rscript
suppressMessages({ library(data.table); library(glmnet); library(ggplot2) })
fs <- "01_scripts/00_style/fig_style.R"; if(!file.exists(fs)) fs <- list.files(".","^fig_style\\.R$",recursive=TRUE,full.names=TRUE)[1]
source(fs); OUT <- "03_figures/panels"
res <- readRDS("02_results/04_mime/mime_results.rds"); mlnms <- names(res$ml.res)
target <- mlnms[grepl("^Enet",mlnms) & grepl("0\\.1",mlnms) & !grepl("[+]",mlnms)][1]; ml <- res$ml.res[[target]]
cvdf <- data.frame(loglam=log(ml$lambda), cvm=ml$cvm, lo=ml$cvlo, up=ml$cvup)
lmin <- log(ml$lambda.min); l1se <- log(ml$lambda.1se)
pCV <- ggplot(cvdf, aes(loglam, cvm)) +
  geom_ribbon(aes(ymin=lo, ymax=up), fill="grey80", alpha=.6) +
  geom_line(color="grey50", linewidth=.3) + geom_point(color=COL_HIGH, size=0.7) +
  geom_vline(xintercept=lmin, linetype=2, color="grey25") +
  geom_vline(xintercept=l1se, linetype=3, color="grey45") +
  annotate("text", x=lmin, y=max(cvdf$up), label="lambda.min", angle=90, vjust=-0.4, hjust=1, size=1.9, fontface="bold") +
  annotate("text", x=l1se, y=max(cvdf$up), label="lambda.1se", angle=90, vjust=-0.4, hjust=1, size=1.9, fontface="bold", color="grey45") +
  labs(x="log(lambda)", y="Partial likelihood deviance") + theme_pub(base=7) +
  theme(axis.text=element_text(face="bold", size=6.5), axis.title=element_text(face="bold", size=7.5))
save_panel(pCV, file.path(OUT,"MF1_Bopt2_cvcurve"), w_mm=54, h_mm=45)
cat(">>> CV 重出(最终尺寸字号)\n")
