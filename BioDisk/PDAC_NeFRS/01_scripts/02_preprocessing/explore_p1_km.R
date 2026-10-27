#!/usr/bin/env Rscript
# ============================================================
# P1 候选基因探索性 KM 分析 (PACA-AU 训练集预演)
#
# 输入:
#   exp_array.PACA-AU.gene.tsv.gz (20744 × 269)
#   clinical_std.PACA-AU.tsv (267 含完整生存)
#   final_candidate_pool.tsv  (P1 20 个高优先级基因)
#
# 输出:
#   P1_KM_results.tsv          每基因 log-rank P, HR, median split
#   P1_KM_grid.pdf             20 个 KM 曲线 4×5 网格
#   P1_KM_top.pdf              Top 6 最显著基因单独大图
# ============================================================

suppressPackageStartupMessages({
    library(data.table)
    library(survival)
    library(survminer)
    library(ggplot2)
    library(patchwork)
})

ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
OUT  <- file.path(ROOT, "02_results/01_gene_pool")

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# 1. 加载数据
# ============================================================
log_msg("=== 1. 加载数据 ===")

# 表达
expr <- fread(file.path(ROOT, "00_raw_data/bulk/ICGC_PACA_AU/exp_array.PACA-AU.gene.tsv.gz"))
setnames(expr, 1, "gene")
log_msg(sprintf("  表达: %d 基因 × %d 样本", nrow(expr), ncol(expr) - 1))

# 临床
clin <- fread(file.path(ROOT, "00_raw_data/bulk/ICGC_PACA_AU/clinical_std.PACA-AU.tsv"))
log_msg(sprintf("  临床: %d 样本", nrow(clin)))

# Primary + 完整生存
clin_use <- clin[sample_type_simple == "Primary" & 
                  !is.na(os_time_month) & !is.na(os_status)]
log_msg(sprintf("  Primary + 完整生存: %d", nrow(clin_use)))

# 候选池
pool <- fread(file.path(OUT, "final_candidate_pool.tsv"))
p1_genes <- pool[priority == "P1_TopPriority"]
log_msg(sprintf("  P1 候选基因: %d", nrow(p1_genes)))

# ============================================================
# 2. 提取 P1 基因在 PACA-AU 中的表达
# ============================================================
log_msg("=== 2. 提取 P1 基因表达 ===")

p1_in_paca <- p1_genes[in_PACA_AU == TRUE]
log_msg(sprintf("  P1 in PACA-AU: %d / %d", nrow(p1_in_paca), nrow(p1_genes)))

# 提取这些基因的表达矩阵
p1_expr <- expr[gene %in% p1_in_paca$GeneSymbol]
log_msg(sprintf("  提取的基因数: %d", nrow(p1_expr)))

# 只保留有完整生存的样本
sample_cols <- intersect(colnames(p1_expr), clin_use$sample_id)
log_msg(sprintf("  匹配样本: %d", length(sample_cols)))

p1_expr_mat <- as.matrix(p1_expr[, ..sample_cols])
rownames(p1_expr_mat) <- p1_expr$gene

# 临床对齐
clin_aligned <- clin_use[match(sample_cols, clin_use$sample_id)]
log_msg(sprintf("  对齐后: %d 样本可用", nrow(clin_aligned)))

# ============================================================
# 3. 对每个 P1 基因跑 KM (中位数二分组)
# ============================================================
log_msg("=== 3. 跑 KM 分析 (中位数二分组) ===")

results <- data.table(
    gene = character(),
    median_expr = numeric(),
    n_high = integer(),
    n_low = integer(),
    death_high = integer(),
    death_low = integer(),
    median_OS_high = numeric(),
    median_OS_low = numeric(),
    HR = numeric(),
    HR_lower = numeric(),
    HR_upper = numeric(),
    cox_p = numeric(),
    logrank_p = numeric(),
    direction = character()
)

km_plots <- list()

for (gene in rownames(p1_expr_mat)) {
    expr_vec <- p1_expr_mat[gene, ]
    med <- median(expr_vec, na.rm = TRUE)
    grp <- factor(ifelse(expr_vec > med, "High", "Low"), levels = c("Low", "High"))
    
    # Cox + log-rank
    df <- data.frame(
        time = clin_aligned$os_time_month,
        status = clin_aligned$os_status,
        group = grp,
        expr = expr_vec
    )
    df <- df[complete.cases(df), ]
    
    if (length(unique(df$group)) < 2) {
        log_msg(sprintf("  %-10s ✗ 无法分组 (表达全相同?)", gene))
        next
    }
    
    surv_obj <- Surv(df$time, df$status)
    
    # Cox 模型 (连续变量, 看 HR per 1 单位 log 表达)
    cox_fit <- tryCatch(
        coxph(surv_obj ~ expr, data = df),
        error = function(e) NULL
    )
    
    # KM 分组比较
    fit_km <- survfit(surv_obj ~ group, data = df)
    
    # log-rank
    sd <- survdiff(surv_obj ~ group, data = df)
    lr_p <- 1 - pchisq(sd$chisq, df = 1)
    
    # 提取统计量
    if (!is.null(cox_fit)) {
        cox_summary <- summary(cox_fit)
        hr <- cox_summary$conf.int[1, "exp(coef)"]
        hr_lo <- cox_summary$conf.int[1, "lower .95"]
        hr_hi <- cox_summary$conf.int[1, "upper .95"]
        cox_p <- cox_summary$coefficients[1, "Pr(>|z|)"]
    } else {
        hr <- hr_lo <- hr_hi <- cox_p <- NA
    }
    
    # 各组中位生存
    sf_summary <- summary(fit_km)$table
    med_low <- if ("group=Low"  %in% rownames(sf_summary)) sf_summary["group=Low",  "median"] else NA
    med_high <- if ("group=High" %in% rownames(sf_summary)) sf_summary["group=High", "median"] else NA
    
    direction <- if (!is.na(hr)) {
        if (hr > 1.05) "Risk (high expr = worse)" 
        else if (hr < 0.95) "Protective (high expr = better)" 
        else "Neutral"
    } else "Unknown"
    
    results <- rbind(results, data.table(
        gene = gene,
        median_expr = med,
        n_high = sum(grp == "High"),
        n_low = sum(grp == "Low"),
        death_high = sum(df$status[grp == "High"] == 1),
        death_low = sum(df$status[grp == "Low"] == 1),
        median_OS_high = med_high,
        median_OS_low = med_low,
        HR = hr,
        HR_lower = hr_lo,
        HR_upper = hr_hi,
        cox_p = cox_p,
        logrank_p = lr_p,
        direction = direction
    ))
    
    # KM 图 (简化版, 不用 survminer 避免依赖)
    df_for_plot <- df
    km_data <- broom::tidy(fit_km) 
    # 用基础 R 画 KM
    
    sig_label <- if (!is.na(lr_p)) {
        if (lr_p < 0.001) "p<0.001"
        else sprintf("p=%.3f", lr_p)
    } else "NA"
    
    color_lr_p <- if (!is.na(lr_p) && lr_p < 0.05) "darkred" else "gray40"
    
    p <- ggsurvplot(fit_km, 
                    data = df,
                    pval = TRUE,
                    pval.method = FALSE,
                    pval.size = 3.5,
                    risk.table = FALSE,
                    palette = c("#4575B4", "#D73027"),
                    legend.labs = c("Low", "High"),
                    legend.title = "",
                    legend = c(0.85, 0.85),
                    xlab = "Time (months)",
                    ylab = "OS",
                    title = sprintf("%s (HR=%.2f)", gene, 
                                    ifelse(is.na(hr), NA_real_, round(hr, 2))),
                    font.title = 10,
                    font.x = 8, font.y = 8,
                    font.tickslab = 7,
                    font.legend = 8,
                    ggtheme = theme_classic() +
                        theme(plot.title = element_text(size = 10, face = "bold",
                                                         color = color_lr_p,
                                                         hjust = 0.5),
                              axis.title = element_text(size = 9),
                              legend.text = element_text(size = 7))
    )
    
    km_plots[[gene]] <- p$plot
    
    sig_mark <- if (!is.na(lr_p) && lr_p < 0.001) "***"
                else if (!is.na(lr_p) && lr_p < 0.01) "**"
                else if (!is.na(lr_p) && lr_p < 0.05) "*"
                else "  "
    
    cat(sprintf("  %-10s %s HR=%.3f (P=%.4f) %s\n",
                gene, sig_mark, ifelse(is.na(hr), NA_real_, hr), 
                ifelse(is.na(lr_p), NA_real_, lr_p), direction))
}

# ============================================================
# 4. 排序 + 输出表格
# ============================================================
log_msg("=== 4. 结果排序 + 输出 ===")

results <- results[order(logrank_p)]

cat("\n=== P1 基因 KM 结果排名 (按 log-rank P 升序) ===\n")
print(results[, .(gene, HR = round(HR, 3), 
                  HR_CI = sprintf("%.2f-%.2f", HR_lower, HR_upper),
                  logrank_p = signif(logrank_p, 3),
                  direction)])

cat("\n=== 显著 (P < 0.05) 基因数 ===\n")
n_sig <- sum(results$logrank_p < 0.05, na.rm = TRUE)
n_sig_001 <- sum(results$logrank_p < 0.01, na.rm = TRUE)
n_sig_0001 <- sum(results$logrank_p < 0.001, na.rm = TRUE)
cat(sprintf("  P < 0.05:   %d / %d\n", n_sig, nrow(results)))
cat(sprintf("  P < 0.01:   %d / %d\n", n_sig_001, nrow(results)))
cat(sprintf("  P < 0.001:  %d / %d\n", n_sig_0001, nrow(results)))

# 输出
fwrite(results, file.path(OUT, "P1_KM_results.tsv"), sep = "\t")
log_msg(sprintf("✓ 表格: %s", file.path(OUT, "P1_KM_results.tsv")))

# ============================================================
# 5. 绘图: 4×5 网格
# ============================================================
log_msg("=== 5. 4×5 网格图 ===")

# 按 log-rank P 排序排列
plot_order <- results$gene
km_plots_ordered <- km_plots[plot_order]
km_plots_ordered <- km_plots_ordered[!sapply(km_plots_ordered, is.null)]

grid_plot <- wrap_plots(km_plots_ordered, ncol = 5, nrow = 4)
ggsave(file.path(OUT, "P1_KM_grid.pdf"), grid_plot, 
       width = 18, height = 14, dpi = 150)
log_msg(sprintf("✓ %s", file.path(OUT, "P1_KM_grid.pdf")))

# Top 6 单独图
top6 <- head(plot_order, 6)
top6_plot <- wrap_plots(km_plots[top6], ncol = 3, nrow = 2)
ggsave(file.path(OUT, "P1_KM_top6.pdf"), top6_plot, 
       width = 15, height = 10, dpi = 150)
log_msg(sprintf("✓ %s", file.path(OUT, "P1_KM_top6.pdf")))

cat("\n=== 完成 ===\n")
cat("注意: 这是探索性预筛, W2 用 Mime 117 算法集成才是定稿\n")
