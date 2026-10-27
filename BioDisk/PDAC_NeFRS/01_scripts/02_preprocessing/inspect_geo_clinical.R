#!/usr/bin/env Rscript
# 详细查看 5 个 GEO 的临床字段, 决定标准化策略

suppressPackageStartupMessages({ library(data.table) })

BULK <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk"
GSEs <- c("GSE71729", "GSE21501", "GSE85916", "GSE62452", "GSE57495")

for (gse in GSEs) {
    cat(sprintf("\n========================================\n%s\n========================================\n", gse))
    
    clin <- fread(file.path(BULK, gse, sprintf("%s.clinical.tsv", gse)),
                  sep = "\t", header = TRUE)
    
    # 只看 characteristics_* 字段
    char_cols <- grep("^characteristics", colnames(clin), value = TRUE)
    
    # 对每个 characteristics 字段, 列出所有不同的前缀 (冒号前部分)
    for (cn in char_cols) {
        vals <- clin[[cn]]
        vals <- vals[vals != "" & !is.na(vals)]
        
        # 提取冒号前部分作为字段类型
        prefixes <- sub(":.*$", "", vals)
        prefix_tab <- table(prefixes)
        
        cat(sprintf("\n--- %s ---\n", cn))
        cat(sprintf("  非空: %d / %d\n", length(vals), nrow(clin)))
        if (length(prefix_tab) > 0) {
            for (p in names(prefix_tab)) {
                # 该前缀对应的样本数
                cat(sprintf("    [%d 样本] %s\n", prefix_tab[p], p))
            }
        }
    }
    
    # 总结: 每个样本的所有 characteristics 字段合并起来, 看哪些样本有 "survival" 字眼
    all_char <- apply(clin[, ..char_cols], 1, function(row) paste(row, collapse = " | "))
    has_surv <- grepl("survival|os time|os.year|overall survival", all_char, ignore.case = TRUE)
    has_death <- grepl("death|vital.status|survival status|death_event", all_char, ignore.case = TRUE)
    cat(sprintf("\n>>> 含生存时间字段的样本: %d / %d\n", sum(has_surv), nrow(clin)))
    cat(sprintf(">>> 含生存状态字段的样本: %d / %d\n", sum(has_death), nrow(clin)))
    cat(sprintf(">>> 同时有时间+状态: %d / %d\n", sum(has_surv & has_death), nrow(clin)))
}
