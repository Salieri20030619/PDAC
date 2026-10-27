#!/usr/bin/env Rscript
# ============================================================
# 最终候选基因池构建
#
# 整合 4 个来源:
#   S1: FerrDb v2          (829 基因)
#   S2: Hwang 2022 NRP     (61 基因)
#   S3: 神经发育 GO ≥ 5    (652 基因, 经聚焦)
#   S4: PDAC PNI 文献      (50 基因)
#
# 输出:
#   final_candidate_pool.tsv       完整候选池
#   final_in_all_cohorts.tsv       6 队列共有的版本 (LASSO 用)
#   gene_pool_summary.tsv          按来源/类别统计
# ============================================================

suppressPackageStartupMessages({
    library(data.table)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
})

POOL <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results/01_gene_pool"
log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# 1. 加载前 3 个来源
# ============================================================
log_msg("=== 1. 加载前 3 个来源 ===")
s1_ferrdb <- fread(file.path(POOL, "01_ferrdb_v2_genes.tsv"))
s2_nrp <- fread(file.path(POOL, "02_nrp_signature_genes.tsv"))
s3_neural_full <- fread(file.path(POOL, "03_neural_dev_genes.tsv"))

# 神经发育 GO 应用 GO term ≥ 5 阈值聚焦
s3_neural <- s3_neural_full[n_GO_terms >= 5]
log_msg(sprintf("  S1 FerrDb v2: %d", nrow(s1_ferrdb)))
log_msg(sprintf("  S2 NRP signature: %d", nrow(s2_nrp)))
log_msg(sprintf("  S3 神经发育 GO (全): %d", nrow(s3_neural_full)))
log_msg(sprintf("  S3 神经发育 GO ≥ 5 (聚焦后): %d", nrow(s3_neural)))

# ============================================================
# 2. 构建 Source 4: PDAC PNI 文献基因
# ============================================================
log_msg("=== 2. 构建 S4: PDAC PNI 文献基因 ===")

# 基于 Chen 2025 Cancer Cell + 经典 PNI 文献综合整理
# 分 4 个子类: 神经营养因子 / 轴突导向 / 神经-肿瘤通讯 / 神经粘附与侵袭
s4_pni <- data.table(
    GeneSymbol = c(
        # ====== 神经营养因子家族 (NTF) ======
        "NGF",        # 神经生长因子
        "NGFR",       # p75NTR, 低亲和力 NGF 受体, PDAC PNI 中关键
        "BDNF",       # 脑源性神经营养因子
        "NTF3",       # 神经营养因子 3
        "NTF4",       # 神经营养因子 4
        "NTRK1",      # TrkA, NGF 高亲和力受体
        "NTRK2",      # TrkB, BDNF 受体
        "NTRK3",      # TrkC
        "GDNF",       # 胶质细胞源性神经营养因子, PDAC PNI 经典
        "GFRA1",      # GDNF family receptor alpha-1
        "GFRA2",
        "RET",        # GDNF 共受体, PDAC 重要
        "ARTN",       # Artemin, GDNF 家族, PDAC PNI 上调
        "NRTN",       # Neurturin
        "PSPN",       # Persephin
        
        # ====== 轴突导向分子 (Axon Guidance) ======
        "SLIT1", "SLIT2", "SLIT3",     # Slit 家族 (排斥)
        "ROBO1", "ROBO2", "ROBO3", "ROBO4",  # Slit 受体
        "SEMA3A", "SEMA3B", "SEMA3C", "SEMA3D", "SEMA3E", "SEMA3F",  # Class 3 semaphorins
        "SEMA4D",     # CD100, PDAC 相关
        "SEMA6A", "SEMA6B",
        "PLXNA1", "PLXNA2", "PLXNA3", "PLXNA4",  # Plexin A
        "PLXNB1", "PLXNB2", "PLXNB3",  # Plexin B
        "NRP1", "NRP2",  # Neuropilins, SEMA 共受体
        "DCC",        # Netrin 受体, 经典轴突导向
        "UNC5A", "UNC5B", "UNC5C", "UNC5D",  # Netrin 排斥受体
        "NTN1", "NTN3", "NTN4", "NTN5",  # Netrins
        "EPHA1", "EPHA2", "EPHA3", "EPHA4", "EPHA5",
        "EPHB1", "EPHB2", "EPHB3", "EPHB4", "EPHB6",
        "EFNA1", "EFNA2", "EFNA3", "EFNA4", "EFNA5",  # Ephrins A
        "EFNB1", "EFNB2", "EFNB3",       # Ephrins B
        
        # ====== 神经-肿瘤通讯 (Neuro-tumor crosstalk) ======
        "CHRNA7",     # nAChR α7
        "ADRB2",      # β2 肾上腺素受体, 慢性应激-PDAC
        "DRD2",       # 多巴胺 D2 受体
        "VIP",        # 血管活性肠肽
        "VIPR1", "VIPR2",
        "CHRM3",      # M3 毒蕈碱受体
        "TAC1",       # 速激肽前体 (P 物质)
        "TACR1",      # 速激肽受体
        "CALCA",      # 降钙素相关
        "GRP",        # 胃泌素释放肽
        "GRPR",
        
        # ====== 神经-肿瘤粘附与侵袭 ======
        "L1CAM",      # 神经粘附分子, PNI 经典
        "NCAM1", "NCAM2",
        "CDH2",       # N-cadherin, 神经粘附 + EMT
        "MUC1",       # 黏蛋白, PDAC 中与 PNI 相关
        "CD24",       # PDAC 干细胞 + PNI
        "MMP2", "MMP9", "MMP14",  # 基质金属蛋白酶, 神经周侵袭
        "PLAU", "PLAUR",  # uPA / uPAR, 侵袭
        
        # ====== PDAC 神经支配相关 (Chen 2025 Cancer Cell) ======
        "S100B",      # 雪旺细胞 marker
        "PMP22",      # 髓鞘蛋白
        "MPZ",        # 髓鞘
        "GFAP",       # 胶质纤维酸性蛋白
        "GAP43",      # 神经再生 marker
        "NES",        # Nestin, 神经干细胞 + PDAC
        "TUBB3",      # βIII-tubulin (神经标志, 也在 NRP)
        
        # ====== 其他 PDAC PNI 关键调控分子 ======
        "MET",        # HGF 受体, PDAC PNI
        "HGF",        # 肝细胞生长因子
        "EGFR",       # PDAC + PNI
        "EREG",       # Epiregulin
        "CXCR4",      # 趋化因子受体
        "CXCL12"      # SDF-1
    ),
    GeneCategory = "PNI",
    SubCategory = c(
        # 对应分类（按顺序）
        rep("NeurotrophicFactor", 15),
        rep("AxonGuidance", 50),
        rep("NeuroTumorComm", 12),
        rep("NeuralAdhesionInvasion", 11),
        rep("PDAC_innervation", 7),
        rep("OtherPNI", 6)
    ),
    Source = "Literature_review_PDAC_PNI"
)

# 去重
s4_pni <- s4_pni[!duplicated(GeneSymbol)]
log_msg(sprintf("  S4 PDAC PNI: %d (去重前 %d, 子类: %d)", 
                nrow(s4_pni), 101, uniqueN(s4_pni$SubCategory)))
log_msg("  子类分布:")
print(table(s4_pni$SubCategory))

# 保存 S4 单独文件
fwrite(s4_pni, file.path(POOL, "04_pni_literature_genes.tsv"), sep = "\t")
log_msg(sprintf("  ✓ 保存 %s", file.path(POOL, "04_pni_literature_genes.tsv")))

# ============================================================
# 3. 标准化 4 个源 (统一格式)
# ============================================================
log_msg("=== 3. 标准化 4 个源 ===")

# S1 标准化
s1_std <- data.table(
    GeneSymbol = s1_ferrdb$GeneSymbol,
    GeneCategory = "Ferroptosis",
    SubCategory = s1_ferrdb$SubCategory,
    Source = "FerrDb_v2"
)

# S2 标准化 (已是标准格式)
s2_std <- s2_nrp[, .(GeneSymbol, GeneCategory, SubCategory, Source)]

# S3 标准化
s3_std <- data.table(
    GeneSymbol = s3_neural$GeneSymbol,
    GeneCategory = "NeuralDev",
    SubCategory = s3_neural$SubCategory,
    Source = "GO_BP_high_specificity"
)

# S4 已经标准格式
s4_std <- s4_pni[, .(GeneSymbol, GeneCategory, SubCategory, Source)]

# ============================================================
# 4. 合并 + 去重 + 多源记录
# ============================================================
log_msg("=== 4. 合并 4 个源 + 去重 ===")

all_pool <- rbind(s1_std, s2_std, s3_std, s4_std)
log_msg(sprintf("  合并前 (含重复): %d 行", nrow(all_pool)))

# 同一基因可能来自多个源, 整合记录
pool_dedup <- all_pool[, .(
    GeneCategories = paste(unique(GeneCategory), collapse = "|"),
    SubCategories = paste(unique(SubCategory), collapse = "|"),
    Sources = paste(unique(Source), collapse = "|"),
    n_sources = uniqueN(Source),
    n_categories = uniqueN(GeneCategory)
), by = GeneSymbol]

log_msg(sprintf("  去重后唯一基因: %d", nrow(pool_dedup)))
log_msg("  按 n_sources (一个基因来自几个源) 分布:")
print(table(pool_dedup$n_sources))
log_msg("  按 n_categories (一个基因属于几个大类) 分布:")
print(table(pool_dedup$n_categories))

# 标记基因优先级 (多源 = 高优先级)
pool_dedup[, priority := fcase(
    n_categories >= 3, "P1_TopPriority",  # 同时是 Ferro + Neural + NRP/PNI
    n_categories == 2, "P2_HighPriority",  # 跨 2 类
    n_categories == 1, "P3_Standard"       # 仅一类
)]
log_msg("  Priority 分布:")
print(table(pool_dedup$priority))

# ============================================================
# 5. 加 EntrezID (org.Hs.eg.db)
# ============================================================
log_msg("=== 5. 加 EntrezID ===")
entrez_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = pool_dedup$GeneSymbol,
    columns = "ENTREZID",
    keytype = "SYMBOL"
)
setDT(entrez_map)
entrez_map <- entrez_map[!duplicated(SYMBOL)]
setnames(entrez_map, c("GeneSymbol", "EntrezID"))
pool_dedup <- merge(pool_dedup, entrez_map, by = "GeneSymbol", all.x = TRUE)
log_msg(sprintf("  有 EntrezID: %d / %d (%.1f%%)",
                sum(!is.na(pool_dedup$EntrezID)), nrow(pool_dedup),
                100 * mean(!is.na(pool_dedup$EntrezID))))

# ============================================================
# 6. 跨 6 个 bulk 队列覆盖率
# ============================================================
log_msg("=== 6. 跨 6 个 bulk 队列覆盖 ===")

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
    in_cohort <- pool_dedup$GeneSymbol %in% g
    log_msg(sprintf("  %-15s : %d / %d (%.1f%%)",
                    name, sum(in_cohort), nrow(pool_dedup),
                    100 * mean(in_cohort)))
}

# 标记每基因的"队列存在记录"
pool_dedup[, in_PACA_AU      := GeneSymbol %in% all_cohort_genes[["PACA-AU"]]]
pool_dedup[, in_GSE71729     := GeneSymbol %in% all_cohort_genes[["GSE71729"]]]
pool_dedup[, in_GSE21501     := GeneSymbol %in% all_cohort_genes[["GSE21501"]]]
pool_dedup[, in_GSE85916     := GeneSymbol %in% all_cohort_genes[["GSE85916"]]]
pool_dedup[, in_GSE62452     := GeneSymbol %in% all_cohort_genes[["GSE62452"]]]
pool_dedup[, in_GSE57495     := GeneSymbol %in% all_cohort_genes[["GSE57495"]]]
pool_dedup[, in_TCGA_curated := GeneSymbol %in% all_cohort_genes[["TCGA_curated"]]]
pool_dedup[, n_cohorts := rowSums(.SD), 
           .SDcols = grep("^in_", colnames(pool_dedup), value = TRUE)]

# 6 队列共有
common_genes <- Reduce(intersect, all_cohort_genes)
pool_in_common <- pool_dedup[GeneSymbol %in% common_genes]
log_msg(sprintf("\n  6 队列全有的候选基因 (LASSO 实际用): %d", nrow(pool_in_common)))

# ============================================================
# 7. 输出
# ============================================================
log_msg("=== 7. 输出 ===")

out_full <- file.path(POOL, "final_candidate_pool.tsv")
fwrite(pool_dedup, out_full, sep = "\t")
log_msg(sprintf("  ✓ %s (%d 基因)", basename(out_full), nrow(pool_dedup)))

out_common <- file.path(POOL, "final_in_all_cohorts.tsv")
fwrite(pool_in_common, out_common, sep = "\t")
log_msg(sprintf("  ✓ %s (%d 基因, 7 队列共有)", basename(out_common), nrow(pool_in_common)))

# Summary 表
summary_dt <- pool_dedup[, .(
    n_genes = .N,
    n_in_common = sum(GeneSymbol %in% common_genes)
), by = GeneCategories][order(-n_genes)]
out_summary <- file.path(POOL, "gene_pool_summary.tsv")
fwrite(summary_dt, out_summary, sep = "\t")
log_msg(sprintf("  ✓ %s", basename(out_summary)))

cat("\n按 GeneCategories 分布 (前 15):\n")
print(head(summary_dt, 15))

# ============================================================
# 8. 关键基因健康检查
# ============================================================
log_msg("=== 8. 关键基因健康检查 ===")
key_genes <- c(
    # 核心 NRP
    "GLIS3", "SOX4", "NEUROD1", "INSM1", "NFIB",
    # 核心铁死亡
    "GPX4", "SLC7A11", "ACSL4", "FTH1", "TFRC",
    # 核心 PNI
    "NGF", "NGFR", "GDNF", "RET", "L1CAM",
    # 神经发育
    "STMN2", "TUBB3", "MAP2", "REST",
    # PDAC 经典
    "KRAS", "TP53", "SMAD4", "CDKN2A",
    # 双重 (FerrDb + Neural)
    "ACSL4", "EGFR", "ZEB1"
)
cat("基因 → 在候选池中 (Priority, Sources):\n")
for (g in unique(key_genes)) {
    if (g %in% pool_dedup$GeneSymbol) {
        info <- pool_dedup[GeneSymbol == g]
        cat(sprintf("  %-10s ✓  [%s] %s\n", g, info$priority, info$Sources))
    } else {
        cat(sprintf("  %-10s ✗ 不在候选池\n", g))
    }
}

log_msg("=== 完成 ===")
