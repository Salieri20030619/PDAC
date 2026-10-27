#!/usr/bin/env Rscript
# Source 3: 神经发育相关基因 (从 GO 数据库提取)
#
# 策略: 从 org.Hs.eg.db + GO.db 提取以下 GO terms 的所有基因:
#   GO:0007399  nervous system development
#   GO:0030182  neuron differentiation
#   GO:0050769  positive regulation of neurogenesis
#   GO:0007409  axonogenesis
#   GO:0048812  neuron projection morphogenesis
#   GO:0021953  central nervous system neuron differentiation
#   GO:0010975  regulation of neuron projection development
#   GO:0048667  cell morphogenesis involved in neuron differentiation
#   GO:0099536  synaptic signaling
#   GO:0007268  chemical synaptic transmission
#   GO:0050808  synapse organization
# 加上 PDAC PNI 相关:
#   GO:0007411  axon guidance
#   GO:0048846  axon extension involved in axon guidance
#   GO:0050731  positive regulation of peptidyl-tyrosine phosphorylation
#   GO:0099054  presynapse assembly

suppressPackageStartupMessages({
    library(data.table)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(GO.db)
})

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# 1. 目标 GO terms
# ============================================================
log_msg("=== 1. 定义目标 GO terms ===")

target_gos <- list(
    # 神经系统发育主干
    "GO:0007399" = "nervous_system_development",
    "GO:0030182" = "neuron_differentiation",
    "GO:0050769" = "positive_regulation_of_neurogenesis",
    "GO:0021953" = "CNS_neuron_differentiation",
    
    # 神经形态和投射
    "GO:0007409" = "axonogenesis",
    "GO:0048812" = "neuron_projection_morphogenesis",
    "GO:0048667" = "cell_morphogenesis_neuron",
    "GO:0010975" = "regulation_neuron_projection",
    
    # 轴突导向 (PDAC PNI 重要)
    "GO:0007411" = "axon_guidance",
    "GO:0048846" = "axon_extension",
    
    # 突触
    "GO:0099536" = "synaptic_signaling",
    "GO:0007268" = "chemical_synaptic_transmission",
    "GO:0050808" = "synapse_organization",
    "GO:0099054" = "presynapse_assembly"
)

cat(sprintf("目标 GO terms: %d 个\n", length(target_gos)))
for (go_id in names(target_gos)) {
    cat(sprintf("  %s : %s\n", go_id, target_gos[[go_id]]))
}

# ============================================================
# 2. 从 org.Hs.eg.db 提取每个 GO term 的基因
# ============================================================
log_msg("=== 2. 提取每个 GO term 的基因 (含子 term) ===")

# 用 GOALL: 包含所有子 term 的基因 (descendant inclusive)
all_genes_list <- list()
for (go_id in names(target_gos)) {
    tryCatch({
        # 获取 GO term 名称
        go_name_obj <- GOTERM[[go_id]]
        go_name <- if (!is.null(go_name_obj)) Term(go_name_obj) else "?"
        
        # 提取基因 (含所有子 term)
        genes <- AnnotationDbi::select(
            org.Hs.eg.db,
            keys = go_id,
            columns = c("SYMBOL", "ENTREZID"),
            keytype = "GOALL"
        )
        setDT(genes)
        # 去重 + 过滤 NA
        genes <- genes[!is.na(SYMBOL) & SYMBOL != ""]
        genes <- genes[!duplicated(SYMBOL)]
        
        all_genes_list[[go_id]] <- data.table(
            GeneSymbol = genes$SYMBOL,
            EntrezID = genes$ENTREZID,
            GO_ID = go_id,
            GO_Name = target_gos[[go_id]],
            GO_Description = go_name
        )
        log_msg(sprintf("  %s (%s): %d 基因",
                        go_id, target_gos[[go_id]], nrow(genes)))
    }, error = function(e) {
        log_msg(sprintf("  ✗ %s 失败: %s", go_id, conditionMessage(e)))
    })
}

# ============================================================
# 3. 合并所有 GO 来源 + 去重
# ============================================================
log_msg("=== 3. 合并所有 GO 来源 + 去重 ===")

all_genes_dt <- rbindlist(all_genes_list, fill = TRUE)
cat(sprintf("合并前 (按 gene-GO 对): %d 行\n", nrow(all_genes_dt)))
cat(sprintf("唯一基因数: %d\n", uniqueN(all_genes_dt$GeneSymbol)))

# 每基因合并所有归属的 GO terms
genes_by_gene <- all_genes_dt[, .(
    EntrezID = EntrezID[1],
    n_GO_terms = .N,
    GO_IDs = paste(GO_ID, collapse = "|"),
    GO_categories = paste(unique(GO_Name), collapse = "|")
), by = GeneSymbol]

log_msg(sprintf("每基因 GO term 数分布:"))
print(table(genes_by_gene$n_GO_terms))

# ============================================================
# 4. 在 6 个 bulk 队列里的覆盖
# ============================================================
log_msg("=== 4. 在 6 个 bulk 队列里的覆盖 ===")

cohorts <- list(
    "PACA-AU" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/ICGC_PACA_AU/exp_array.PACA-AU.gene.tsv.gz",
    "GSE71729" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/GSE71729/GSE71729.gene.tsv.gz",
    "GSE21501" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/GSE21501/GSE21501.gene.tsv.gz",
    "GSE85916" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/GSE85916/GSE85916.gene.tsv.gz",
    "GSE62452" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/GSE62452/GSE62452.gene.tsv.gz",
    "GSE57495" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/GSE57495/GSE57495.gene.tsv.gz",
    "TCGA_curated" = "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/00_raw_data/bulk/TCGA_PAAD/TCGA-PAAD.curated.gene.tsv.gz"
)

all_cohort_genes <- list()
for (name in names(cohorts)) {
    g <- fread(cohorts[[name]], select = 1, header = TRUE)[[1]]
    all_cohort_genes[[name]] <- g
    in_cohort <- genes_by_gene$GeneSymbol %in% g
    log_msg(sprintf("  %-15s: %d / %d (%.1f%%)",
                    name, sum(in_cohort), nrow(genes_by_gene),
                    100 * mean(in_cohort)))
}

common_genes <- Reduce(intersect, all_cohort_genes)
neural_in_common <- intersect(genes_by_gene$GeneSymbol, common_genes)
log_msg(sprintf("  6 队列共有的神经发育基因: %d / %d (%.1f%%)",
                length(neural_in_common), nrow(genes_by_gene),
                100 * length(neural_in_common) / nrow(genes_by_gene)))

# ============================================================
# 5. 关键基因健康检查
# ============================================================
log_msg("=== 5. 关键基因健康检查 ===")
key_neural <- c(
    # 神经发育主控
    "PAX6", "ASCL1", "NEUROD1", "INSM1", "GLIS3", "REST",
    # 轴突导向
    "ROBO1", "SLIT2", "SEMA3A", "NRP1", "EPHA4",
    # 突触
    "SYN1", "SYP", "SNAP25", "DLG4",
    # 神经骨架/经典 marker
    "TUBB3", "MAP2", "STMN2", "NEFL", "DCX"
)
for (g in key_neural) {
    in_db <- g %in% genes_by_gene$GeneSymbol
    if (in_db) {
        n_gos <- genes_by_gene[GeneSymbol == g, n_GO_terms]
        cat(sprintf("  %-10s ✓ (%d 个 GO term)\n", g, n_gos))
    } else {
        cat(sprintf("  %-10s ✗\n", g))
    }
}

# ============================================================
# 6. 整理为标准格式 + 保存
# ============================================================
log_msg("=== 6. 整理 + 保存 ===")

neural_dev_std <- data.table(
    GeneSymbol = genes_by_gene$GeneSymbol,
    GeneCategory = "NeuralDev",
    SubCategory = genes_by_gene$GO_categories,
    Source = "GO_BP",
    EntrezID = genes_by_gene$EntrezID,
    n_GO_terms = genes_by_gene$n_GO_terms,
    GO_IDs = genes_by_gene$GO_IDs
)

cat(sprintf("最终神经发育基因池: %d 基因\n", nrow(neural_dev_std)))

out_dir <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/01_gene_pool/"
out_file <- file.path(out_dir, "03_neural_dev_genes.tsv")
fwrite(neural_dev_std, out_file, sep = "\t")
log_msg(sprintf("✓ 保存: %s", out_file))

# 也保存"在所有 6 队列里有的"那一部分 (LASSO 用)
neural_in_common_std <- neural_dev_std[GeneSymbol %in% neural_in_common]
out_file2 <- file.path(out_dir, "03_neural_dev_in_all_cohorts.tsv")
fwrite(neural_in_common_std, out_file2, sep = "\t")
log_msg(sprintf("✓ 保存 (6 队列共有): %s (%d 基因)", out_file2, nrow(neural_in_common_std)))
