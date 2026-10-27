#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE202051 数据探查脚本
目的: 摸清数据结构, 决定项目故事方向是否走 PNI 路线
作者: PDAC Spatial Project
"""

import os
import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

# ========== 配置 ==========
H5AD_PATH = '/mnt/BioDisk/wangmingyuan/PDAC_spatial_PNI/00_raw_data/GSE202051/GSE202051_totaldata-final-toshare.h5ad'
OUT_DIR = '/mnt/BioDisk/wangmingyuan/PDAC_spatial_PNI/02_results/01_data_exploration'
os.makedirs(OUT_DIR, exist_ok=True)

# 神经/PNI 相关标志物 (用于评估 PNI 故事可行性)
NEURAL_MARKERS = [
    'S100B', 'PLP1', 'MPZ', 'UCHL1', 'GAP43', 'SOX10',  # Schwann/神经
    'NGFR', 'NRP1', 'NRP2',                              # 神经营养因子受体
    'NTRK1', 'NTRK2', 'NTRK3',                           # NGF/BDNF 受体
    'NGF', 'BDNF', 'GDNF', 'NTF3',                       # 神经营养因子配体
    'SEMA3A', 'SEMA3E', 'PLXNA1',                        # 轴突导向
    'TUBB3', 'NEFL', 'NEFM', 'CHGA', 'SYP',              # 神经元/神经内分泌
]

# Hwang 2022 原文报告的关键 program 基因 (验证用)
HWANG_MARKERS = {
    'classical': ['GATA6', 'TFF1', 'TFF2', 'TFF3', 'CLDN18'],
    'basaloid': ['KRT5', 'KRT6A', 'KRT17', 'TP63', 'S100A2'],
    'NRP_neural_like': ['CHGA', 'CHGB', 'SYP', 'INSM1', 'ASCL1'],  # 神经样
}

print("=" * 80)
print("GSE202051 数据探查报告")
print("=" * 80)

# ========== 1. 加载数据 ==========
print("\n[1] 加载 h5ad 文件...")
print(f"路径: {H5AD_PATH}")
print(f"文件大小: {os.path.getsize(H5AD_PATH) / 1e9:.2f} GB")

adata = sc.read_h5ad(H5AD_PATH)
print(f"\n加载完成! AnnData 对象:")
print(adata)

# ========== 2. 基础维度 ==========
print("\n" + "=" * 80)
print("[2] 基础维度")
print("=" * 80)
print(f"obs (细胞/AOI) 总数: {adata.n_obs:,}")
print(f"var (基因) 总数: {adata.n_vars:,}")
print(f"X 矩阵类型: {type(adata.X).__name__}")
print(f"X 矩阵 dtype: {adata.X.dtype}")
if hasattr(adata.X, 'nnz'):
    sparsity = 1 - adata.X.nnz / (adata.n_obs * adata.n_vars)
    print(f"稀疏度: {sparsity * 100:.2f}% (非零比例: {(1-sparsity)*100:.2f}%)")
print(f"内存占用估计: {adata.X.data.nbytes / 1e9:.2f} GB (X 矩阵)" if hasattr(adata.X, 'data') else f"X 大小: {adata.X.nbytes / 1e9:.2f} GB")

# ========== 3. obs 元数据 ==========
print("\n" + "=" * 80)
print("[3] obs 元数据列 (最关键!)")
print("=" * 80)
print(f"obs 总共 {len(adata.obs.columns)} 列:")
for i, col in enumerate(adata.obs.columns):
    dtype = adata.obs[col].dtype
    n_unique = adata.obs[col].nunique()
    n_null = adata.obs[col].isnull().sum()
    print(f"  {i+1:3d}. {col:40s} | dtype: {str(dtype):15s} | unique: {n_unique:6d} | null: {n_null:6d}")

# ========== 4. 每列的具体取值分布 ==========
print("\n" + "=" * 80)
print("[4] obs 每列取值分布 (重点看!)")
print("=" * 80)
for col in adata.obs.columns:
    n_unique = adata.obs[col].nunique()
    if n_unique <= 30:  # 类别型, 显示所有取值
        print(f"\n--- {col} (类别型, {n_unique} 种取值) ---")
        vc = adata.obs[col].value_counts(dropna=False)
        for val, cnt in vc.items():
            print(f"    {str(val)[:50]:50s} : {cnt:8,d} ({cnt/adata.n_obs*100:.1f}%)")
    elif n_unique <= 100:  # 中等基数, 显示前 20
        print(f"\n--- {col} (中基数, {n_unique} 种取值, 显示前 20) ---")
        vc = adata.obs[col].value_counts(dropna=False).head(20)
        for val, cnt in vc.items():
            print(f"    {str(val)[:50]:50s} : {cnt:8,d}")
    else:  # 高基数 (通常是 cell ID 或连续值)
        if adata.obs[col].dtype in ['float32', 'float64', 'int32', 'int64']:
            print(f"\n--- {col} (数值型, {n_unique} 种取值) ---")
            print(f"    min: {adata.obs[col].min():.3f} | max: {adata.obs[col].max():.3f} | mean: {adata.obs[col].mean():.3f} | median: {adata.obs[col].median():.3f}")
        else:
            print(f"\n--- {col} (高基数类别, {n_unique} 种取值, 跳过详细显示) ---")
            print(f"    示例值: {list(adata.obs[col].dropna().unique()[:5])}")

# ========== 5. obsm / uns / layers 内容 ==========
print("\n" + "=" * 80)
print("[5] obsm / varm / uns / layers 内容")
print("=" * 80)

print(f"\nobsm (含 spatial 坐标、降维结果):")
if adata.obsm:
    for key in adata.obsm.keys():
        arr = adata.obsm[key]
        print(f"  - {key}: shape={arr.shape}, dtype={arr.dtype}")
else:
    print("  (空)")

print(f"\nvarm:")
if adata.varm:
    for key in adata.varm.keys():
        print(f"  - {key}: shape={adata.varm[key].shape}")
else:
    print("  (空)")

print(f"\nlayers:")
if adata.layers:
    for key in adata.layers.keys():
        print(f"  - {key}: shape={adata.layers[key].shape}, dtype={adata.layers[key].dtype}")
else:
    print("  (空)")

print(f"\nuns keys:")
if adata.uns:
    for key in adata.uns.keys():
        val = adata.uns[key]
        type_name = type(val).__name__
        if hasattr(val, 'shape'):
            print(f"  - {key}: {type_name}, shape={val.shape}")
        elif isinstance(val, dict):
            print(f"  - {key}: dict with {len(val)} keys: {list(val.keys())[:10]}")
        else:
            print(f"  - {key}: {type_name}")
else:
    print("  (空)")

# ========== 6. var (基因) 信息 ==========
print("\n" + "=" * 80)
print("[6] var (基因) 信息")
print("=" * 80)
print(f"var 列: {list(adata.var.columns)}")
print(f"\n前 10 个基因:")
print(adata.var.head(10))

# ========== 7. 神经标志物检查 (PNI 故事可行性) ==========
print("\n" + "=" * 80)
print("[7] PNI 故事可行性检查 - 神经标志物表达情况")
print("=" * 80)

# 找 var 里的基因名列
gene_name_col = None
for candidate in ['gene_symbols', 'gene_names', 'feature_name', 'symbol', 'gene_symbol']:
    if candidate in adata.var.columns:
        gene_name_col = candidate
        break

if gene_name_col:
    print(f"使用 var['{gene_name_col}'] 作为基因名")
    gene_names = adata.var[gene_name_col].astype(str).values
else:
    print(f"使用 var_names (index) 作为基因名")
    gene_names = adata.var_names.astype(str).values

print(f"\n神经标志物 ({len(NEURAL_MARKERS)} 个) 在数据中的存在情况:")
found_neural = []
missing_neural = []
for marker in NEURAL_MARKERS:
    if marker in gene_names:
        idx = np.where(gene_names == marker)[0][0]
        # 计算该基因的表达细胞比例
        if hasattr(adata.X, 'tocsc'):
            expr = adata.X[:, idx].toarray().flatten()
        else:
            expr = adata.X[:, idx].flatten()
        pct_expressed = (expr > 0).sum() / len(expr) * 100
        mean_expr = expr[expr > 0].mean() if (expr > 0).any() else 0
        found_neural.append((marker, pct_expressed, mean_expr))
        print(f"  ✅ {marker:10s} | 表达细胞比例: {pct_expressed:5.2f}% | 阳性表达均值: {mean_expr:.3f}")
    else:
        missing_neural.append(marker)
        print(f"  ❌ {marker:10s} | 不在数据中")

print(f"\n小结: {len(found_neural)}/{len(NEURAL_MARKERS)} 个神经标志物存在")

# ========== 8. Hwang 原文 marker 验证 ==========
print("\n" + "=" * 80)
print("[8] Hwang 2022 原文 marker 存在情况")
print("=" * 80)
for program, markers in HWANG_MARKERS.items():
    print(f"\n{program}:")
    for m in markers:
        if m in gene_names:
            idx = np.where(gene_names == m)[0][0]
            if hasattr(adata.X, 'tocsc'):
                expr = adata.X[:, idx].toarray().flatten()
            else:
                expr = adata.X[:, idx].flatten()
            pct = (expr > 0).sum() / len(expr) * 100
            print(f"  ✅ {m:10s} ({pct:.1f}% 表达)")
        else:
            print(f"  ❌ {m:10s} 不存在")

# ========== 9. 关键决策信息汇总 ==========
print("\n" + "=" * 80)
print("[9] 项目决策关键信息汇总")
print("=" * 80)

# 检查是否有 spatial 信息
has_spatial = 'spatial' in adata.obsm.keys()
print(f"\n是否有 spatial 坐标 (obsm['spatial']): {has_spatial}")

# 尝试识别 snRNA vs spatial 的区分列
print("\n寻找可能区分 snRNA vs spatial 的列...")
candidates = []
for col in adata.obs.columns:
    col_lower = col.lower()
    if any(kw in col_lower for kw in ['modality', 'tech', 'platform', 'assay', 'method', 'data_type', 'type']):
        candidates.append(col)
        vc = adata.obs[col].value_counts()
        print(f"  候选列 '{col}': {dict(vc.head(10))}")

if not candidates:
    print("  未找到明显的 modality 列, 后续可能需要从样本名/uns 推断")

# 尝试找治疗状态列
print("\n寻找治疗状态列...")
treat_candidates = []
for col in adata.obs.columns:
    col_lower = col.lower()
    if any(kw in col_lower for kw in ['treat', 'therap', 'neoadj', 'naive', 'condition']):
        treat_candidates.append(col)
        vc = adata.obs[col].value_counts()
        print(f"  候选列 '{col}': {dict(vc.head(10))}")

# 尝试找样本/病人 ID 列
print("\n寻找样本/病人 ID 列...")
for col in adata.obs.columns:
    col_lower = col.lower()
    if any(kw in col_lower for kw in ['sample', 'patient', 'donor', 'specimen', 'subject']):
        n_unique = adata.obs[col].nunique()
        print(f"  候选列 '{col}': {n_unique} 个 unique 值")

# 尝试找细胞类型注释列
print("\n寻找细胞类型注释列...")
for col in adata.obs.columns:
    col_lower = col.lower()
    if any(kw in col_lower for kw in ['cell_type', 'celltype', 'cluster', 'annotation', 'label', 'leiden', 'louvain', 'program']):
        n_unique = adata.obs[col].nunique()
        if n_unique < 100:
            vc = adata.obs[col].value_counts().head(10)
            print(f"  候选列 '{col}' ({n_unique} 类): {dict(vc)}")

# ========== 10. 保存关键信息 ==========
print("\n" + "=" * 80)
print("[10] 保存元数据到文件")
print("=" * 80)

# 保存 obs 的前 1000 行预览
obs_preview = adata.obs.head(1000)
obs_preview.to_csv(os.path.join(OUT_DIR, 'obs_preview_1000.csv'))
print(f"✅ obs 前 1000 行已保存: {OUT_DIR}/obs_preview_1000.csv")

# 保存 obs 各列的取值分布
with open(os.path.join(OUT_DIR, 'obs_columns_summary.txt'), 'w') as f:
    f.write(f"GSE202051 obs 列摘要\n")
    f.write(f"总 obs 数: {adata.n_obs}\n")
    f.write(f"总列数: {len(adata.obs.columns)}\n")
    f.write("=" * 80 + "\n\n")
    for col in adata.obs.columns:
        f.write(f"\n--- {col} ---\n")
        f.write(f"dtype: {adata.obs[col].dtype}\n")
        f.write(f"unique: {adata.obs[col].nunique()}\n")
        f.write(f"null: {adata.obs[col].isnull().sum()}\n")
        if adata.obs[col].nunique() <= 50:
            vc = adata.obs[col].value_counts(dropna=False)
            for val, cnt in vc.items():
                f.write(f"  {val}: {cnt}\n")
print(f"✅ obs 列摘要已保存: {OUT_DIR}/obs_columns_summary.txt")

# 保存 var 信息
adata.var.to_csv(os.path.join(OUT_DIR, 'var_info.csv'))
print(f"✅ var 信息已保存: {OUT_DIR}/var_info.csv")

# 保存 neural marker 检查结果
with open(os.path.join(OUT_DIR, 'neural_markers_check.txt'), 'w') as f:
    f.write("PNI 故事可行性 - 神经标志物检查\n")
    f.write("=" * 60 + "\n\n")
    f.write("存在的标志物:\n")
    for marker, pct, mean_expr in found_neural:
        f.write(f"  {marker:10s} | {pct:5.2f}% 表达 | 均值 {mean_expr:.3f}\n")
    f.write(f"\n缺失的标志物: {missing_neural}\n")
print(f"✅ 神经 marker 检查已保存: {OUT_DIR}/neural_markers_check.txt")

print("\n" + "=" * 80)
print("🎉 探查完成! 请将这份输出贴给 Claude, 一起决定下一步方向")
print("=" * 80)
