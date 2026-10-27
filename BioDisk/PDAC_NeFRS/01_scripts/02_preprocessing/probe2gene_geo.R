#!/usr/bin/env Rscript
# ============================================================
# 5 个 GEO 队列 probe → gene symbol 映射 + 多 probe 聚合
#
# 输入: {GSE}.expr.tsv.gz (probe × sample) + 对应注释源
# 输出: {GSE}.gene.tsv.gz (gene × sample) + probe2gene_map.tsv
#
# 聚合策略: 同基因多 probe 取最大方差 probe
#
# 队列与注释源:
#   GSE71729: expr 已是 gene level (直接复制)
#   GSE21501: GPL4133.annot.gz       (NCBI GEO 标准, ID + Gene symbol)
#   GSE85916: hgu219.db (R 包)        (Affy HG-U219)
#   GSE62452: GPL6244.annot.gz       (NCBI GEO 标准, ID + Gene symbol)
#   GSE57495: GPL15048_platform.txt.gz (NCBI acc.cgi, ID + GeneSymbol)
# ============================================================

suppressPackageStartupMessages({
    library(data.table)
})

BULK <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk"

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# 通用聚合: 加载表达, 应用 probe→gene 映射, 取最大方差 probe
# ============================================================
aggregate_to_gene <- function(gse, expr_dt, probe2gene) {
    log_msg(sprintf("  表达矩阵: %d probes × %d 样本", nrow(expr_dt), ncol(expr_dt) - 1))
    log_msg(sprintf("  注释表 (probe→gene): %d 条 (有 gene symbol)", nrow(probe2gene)))
    
    # 内连接表达矩阵和注释
    mapped <- merge(expr_dt, probe2gene, by = "probe_id", all.x = FALSE)
    n_total <- nrow(expr_dt)
    n_mapped <- nrow(mapped)
    log_msg(sprintf("  成功映射: %d / %d (%.1f%%)", n_mapped, n_total, 100 * n_mapped / n_total))
    log_msg(sprintf("  映射到的唯一 gene: %d", uniqueN(mapped$gene_symbol)))
    
    # 每基因 probe 数分布
    probe_per_gene <- mapped[, .N, by = gene_symbol]
    log_msg(sprintf("    1 probe: %d | 2: %d | 3: %d | 4+: %d 基因",
                    sum(probe_per_gene$N == 1),
                    sum(probe_per_gene$N == 2),
                    sum(probe_per_gene$N == 3),
                    sum(probe_per_gene$N >= 4)))
    
    # 多 probe 聚合: 最大方差
    sample_cols <- setdiff(colnames(expr_dt), "probe_id")
    log_msg("  计算 probe 方差并选最大方差 probe...")
    mapped[, probe_var := apply(.SD, 1, var, na.rm = TRUE), .SDcols = sample_cols]
    selected <- mapped[, .SD[which.max(probe_var)], by = gene_symbol]
    log_msg(sprintf("  聚合后唯一基因数: %d", nrow(selected)))
    
    # 输出 gene-level expression
    gene_mat <- selected[, c("gene_symbol", sample_cols), with = FALSE]
    setnames(gene_mat, "gene_symbol", "gene")
    
    # probe→gene 映射记录
    map_record <- selected[, .(gene_symbol, selected_probe = probe_id, probe_var)]
    map_record <- merge(map_record, probe_per_gene, by = "gene_symbol")
    setnames(map_record, "N", "n_probes_available")
    
    out_gene <- file.path(BULK, gse, sprintf("%s.gene.tsv.gz", gse))
    out_map  <- file.path(BULK, gse, sprintf("%s.probe2gene_map.tsv", gse))
    fwrite(gene_mat, out_gene, sep = "\t", compress = "gzip")
    fwrite(map_record, out_map, sep = "\t")
    log_msg(sprintf("  ✓ %s (%.1f MB)", basename(out_gene), file.size(out_gene) / 1024^2))
    log_msg(sprintf("  ✓ %s", basename(out_map)))
    
    return(list(gene_mat = gene_mat, n_genes = nrow(gene_mat)))
}

# ============================================================
# 读 NCBI GEO 标准 .annot.gz 格式 (GPL10558 / GPL4133 / GPL6244)
# 结构: 跳过 ! 和 # 开头的注释行, 从 ID<tab>Gene title 开始
# 第 1 列 = ID (probe), 第 3 列 = Gene symbol
# ============================================================
read_geo_annot <- function(annot_file) {
    annot <- fread(annot_file,
                   sep = "\t",
                   skip = "ID\tGene title",
                   header = TRUE,
                   quote = "",
                   fill = TRUE,
                   showProgress = FALSE)
    dt <- data.table(probe_id = annot$ID, gene_symbol = trimws(annot$`Gene symbol`))
    dt <- dt[!is.na(probe_id) & probe_id != ""]
    dt <- dt[!is.na(gene_symbol) & gene_symbol != ""]
    dt
}

# ============================================================
# 队列 1: GSE71729 — 表达矩阵已是 gene level
# ============================================================
process_GSE71729 <- function() {
    log_msg("=== GSE71729 (gene-level 已聚合, 直接复制) ===")
    expr_dt <- fread(file.path(BULK, "GSE71729/GSE71729.expr.tsv.gz"),
                     header = TRUE, showProgress = FALSE)
    setnames(expr_dt, 1, "gene")
    
    log_msg(sprintf("  矩阵: %d genes × %d 样本", nrow(expr_dt), ncol(expr_dt) - 1))
    log_msg(sprintf("  唯一 gene 数: %d", uniqueN(expr_dt$gene)))
    
    # 检查是否有重复 (理论上不应该)
    dups <- expr_dt[duplicated(gene), unique(gene)]
    if (length(dups) > 0) {
        log_msg(sprintf("  ⚠️ 发现 %d 个重复 gene, 取第一个", length(dups)))
        expr_dt <- expr_dt[!duplicated(gene)]
    }
    
    out_gene <- file.path(BULK, "GSE71729/GSE71729.gene.tsv.gz")
    fwrite(expr_dt, out_gene, sep = "\t", compress = "gzip")
    log_msg(sprintf("  ✓ %s (%.1f MB)", basename(out_gene), file.size(out_gene) / 1024^2))
    
    return(list(gene_mat = expr_dt, n_genes = nrow(expr_dt)))
}

# ============================================================
# 队列 2: GSE21501 — GPL4133.annot.gz
# ============================================================
process_GSE21501 <- function() {
    log_msg("=== GSE21501 (Agilent 4×44K, GPL4133) ===")
    expr_dt <- fread(file.path(BULK, "GSE21501/GSE21501.expr.tsv.gz"),
                     header = TRUE, showProgress = FALSE)
    setnames(expr_dt, 1, "probe_id")
    expr_dt[, probe_id := as.character(probe_id)]
    
    annot <- read_geo_annot(file.path(BULK, "GSE21501/GPL4133.annot.gz"))
    annot[, probe_id := as.character(probe_id)]
    
    aggregate_to_gene("GSE21501", expr_dt, annot)
}

# ============================================================
# 队列 3: GSE85916 — hgu219.db R 包
# ============================================================
process_GSE85916 <- function() {
    log_msg("=== GSE85916 (Affy HG-U219, hgu219.db) ===")
    suppressPackageStartupMessages({ library(hgu219.db); library(AnnotationDbi) })
    
    expr_dt <- fread(file.path(BULK, "GSE85916/GSE85916.expr.tsv.gz"),
                     header = TRUE, showProgress = FALSE)
    setnames(expr_dt, 1, "probe_id")
    expr_dt[, probe_id := as.character(probe_id)]
    
    # 从 hgu219.db 取 probe → gene symbol
    log_msg("  从 hgu219.db 提取 probe→symbol 映射...")
    all_probes <- expr_dt$probe_id
    map <- AnnotationDbi::select(hgu219.db,
                                  keys = all_probes,
                                  columns = c("SYMBOL"),
                                  keytype = "PROBEID")
    setDT(map)
    annot <- data.table(probe_id = map$PROBEID, gene_symbol = map$SYMBOL)
    annot <- annot[!is.na(gene_symbol) & gene_symbol != ""]
    annot <- annot[!duplicated(probe_id)]  # 取 1:1
    
    aggregate_to_gene("GSE85916", expr_dt, annot)
}

# ============================================================
# 队列 4: GSE62452 — GPL6244.annot.gz
# ============================================================
process_GSE62452 <- function() {
    log_msg("=== GSE62452 (Affy HuGene 1.0 ST, GPL6244) ===")
    expr_dt <- fread(file.path(BULK, "GSE62452/GSE62452.expr.tsv.gz"),
                     header = TRUE, showProgress = FALSE)
    setnames(expr_dt, 1, "probe_id")
    expr_dt[, probe_id := as.character(probe_id)]
    
    annot <- read_geo_annot(file.path(BULK, "GSE62452/GPL6244.annot.gz"))
    annot[, probe_id := as.character(probe_id)]
    
    aggregate_to_gene("GSE62452", expr_dt, annot)
}

# ============================================================
# 队列 5: GSE57495 — GPL15048_platform.txt.gz (NCBI acc.cgi)
# ============================================================
process_GSE57495 <- function() {
    log_msg("=== GSE57495 (Rosetta/Merck Affy, GPL15048) ===")
    expr_dt <- fread(file.path(BULK, "GSE57495/GSE57495.expr.tsv.gz"),
                     header = TRUE, showProgress = FALSE)
    setnames(expr_dt, 1, "probe_id")
    expr_dt[, probe_id := as.character(probe_id)]
    
    # 读 GPL15048 platform table (从 !platform_table_begin 开始)
    annot_file <- file.path(BULK, "GSE57495/GPL15048_platform.txt.gz")
    annot <- fread(annot_file,
                   sep = "\t",
                   skip = "ID\tGB_LIST",   # 找到表头
                   header = TRUE,
                   quote = "",
                   fill = TRUE,
                   showProgress = FALSE)
    
    # 过滤掉 platform_table_end 行
    annot <- annot[!grepl("^!", ID)]
    
    dt <- data.table(probe_id = as.character(annot$ID),
                     gene_symbol = trimws(annot$GeneSymbol))
    dt <- dt[!is.na(probe_id) & probe_id != ""]
    dt <- dt[!is.na(gene_symbol) & gene_symbol != ""]
    
    aggregate_to_gene("GSE57495", expr_dt, dt)
}

# ============================================================
# 主流程
# ============================================================
results <- list()
for (fn in list(
    list(name = "GSE71729", proc = process_GSE71729),
    list(name = "GSE21501", proc = process_GSE21501),
    list(name = "GSE85916", proc = process_GSE85916),
    list(name = "GSE62452", proc = process_GSE62452),
    list(name = "GSE57495", proc = process_GSE57495)
)) {
    tryCatch({
        results[[fn$name]] <- fn$proc()
    }, error = function(e) {
        log_msg(sprintf("  ✗ %s 失败: %s", fn$name, conditionMessage(e)))
        results[[fn$name]] <<- NULL
    })
    cat("\n")
}

# 总结
log_msg("=== 总结: 5 个 GEO 队列 gene-level 表达矩阵 ===")
for (gse in c("GSE71729", "GSE21501", "GSE85916", "GSE62452", "GSE57495")) {
    r <- results[[gse]]
    if (is.null(r)) {
        log_msg(sprintf("  %s: ✗ 失败", gse))
    } else {
        log_msg(sprintf("  %s: %d genes ✓", gse, r$n_genes))
    }
}

# 关键基因健康检查 (用 GSE85916 作为代表性, 因为它是干净的 80 例)
log_msg("=== 关键基因健康检查 (GSE85916) ===")
check_genes <- c("KRAS", "TP53", "SMAD4", "CDKN2A", "GAPDH", "ACTB",
                 "GLIS3", "NEUROD1", "INSM1", "NFIB", "SOX4",
                 "GPX4", "SLC7A11", "ACSL4", "FTH1")
if (!is.null(results$GSE85916)) {
    gm <- results$GSE85916$gene_mat
    for (g in check_genes) {
        in_mat <- g %in% gm$gene
        cat(sprintf("    %-10s %s\n", g, if (in_mat) "✓" else "✗"))
    }
}
