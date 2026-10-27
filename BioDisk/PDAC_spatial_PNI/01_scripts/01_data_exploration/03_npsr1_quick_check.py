#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
NPSR1 + 同事 marker 快速验证
"""
import os, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, scanpy as sc

H5AD = '/mnt/BioDisk/wangmingyuan/PDAC_spatial_PNI/00_raw_data/GSE202051/GSE202051_totaldata-final-toshare.h5ad'
print("加载数据...")
adata = sc.read_h5ad(H5AD)

# 同事的关键 marker 集合
COLLEAGUE_MARKERS = {
    'High_Risk_Malignant': ['IGF2BP3', 'NPSR1', 'CDH19', 'ZNF536'],
    'Gastric_Malignant': ['ANXA10', 'LGR5', 'GATA4', 'TMEM178B'],
    'Rare_Malignant': ['VSTM2L', 'NPSR1', 'IGF2BP3'],
    'Schwann_colleague': ['TMEM178B', 'CDH19', 'ZNF536'],
    'PNI_axis_重点': ['NPSR1', 'NRP1', 'NRP2', 'NTRK2', 'NTRK3'],
}

genes = adata.var_names.astype(str).values

print("\n" + "=" * 80)
print("同事的 marker 在数据中的总体表达")
print("=" * 80)
all_markers = []
for v in COLLEAGUE_MARKERS.values():
    all_markers.extend(v)
all_markers = list(set([g for g in all_markers if g in genes]))

for g in all_markers:
    idx = np.where(genes == g)[0][0]
    expr = adata.X[:, idx].toarray().flatten() if hasattr(adata.X, 'toarray') else adata.X[:, idx]
    pct = (expr > 0).sum() / len(expr) * 100
    mean_pos = expr[expr > 0].mean() if (expr > 0).any() else 0
    print(f"  {g:10s} | 全局表达 {pct:5.2f}% | 阳性均值 {mean_pos:.3f}")

# 重点: NPSR1 在不同细胞类型的表达
print("\n" + "=" * 80)
print("NPSR1 在各细胞类型 (Level 1) 的表达")
print("=" * 80)
if 'NPSR1' in genes:
    idx = np.where(genes == 'NPSR1')[0][0]
    npsr1 = adata.X[:, idx].toarray().flatten() if hasattr(adata.X, 'toarray') else adata.X[:, idx]
    adata.obs['NPSR1_expr'] = npsr1
    adata.obs['NPSR1_pos'] = (npsr1 > 0).astype(int)
    
    summary = adata.obs.groupby('Level 1 Annotation', observed=True).agg(
        n_cells=('NPSR1_pos', 'size'),
        n_positive=('NPSR1_pos', 'sum'),
        mean_expr=('NPSR1_expr', 'mean'),
    )
    summary['pct_positive'] = summary['n_positive'] / summary['n_cells'] * 100
    summary = summary.sort_values('pct_positive', ascending=False)
    print(summary)
    
    print("\n" + "=" * 80)
    print("NPSR1+ malignant cells 在治疗组的分布")
    print("=" * 80)
    malig = adata.obs[adata.obs['Level 1 Annotation'] == 'Epithelial (malignant)'].copy()
    print(f"Malignant cells 总数: {len(malig)}")
    print(f"NPSR1+ malignant cells: {malig['NPSR1_pos'].sum()} ({malig['NPSR1_pos'].mean()*100:.2f}%)")
    
    ct = pd.crosstab(malig['NPSR1_pos'], malig['new_treatment'])
    print("\nNPSR1+ vs treatment:")
    print(ct)
    print("\nNPSR1+ 比例 by treatment:")
    pct_table = ct.div(ct.sum(axis=0), axis=1) * 100
    print(pct_table.round(2))
    
    print("\n" + "=" * 80)
    print("NPSR1+ malignant cells 在 response 的分布")
    print("=" * 80)
    ct2 = pd.crosstab(malig['NPSR1_pos'], malig['response'])
    pct_table2 = ct2.div(ct2.sum(axis=0), axis=1) * 100
    print(pct_table2.round(2))

# NPSR1 与神经相关基因共表达?
print("\n" + "=" * 80)
print("NPSR1+ malignant cells 是否也表达 NRP1/NTRK?")
print("=" * 80)
malig_mask = adata.obs['Level 1 Annotation'] == 'Epithelial (malignant)'
neural_genes = ['NPSR1', 'NRP1', 'NRP2', 'NTRK2', 'NTRK3', 'GAP43', 'CDH2', 'SEMA3A']
neural_genes = [g for g in neural_genes if g in genes]

malig_expr = pd.DataFrame()
for g in neural_genes:
    idx = np.where(genes == g)[0][0]
    e = adata.X[malig_mask, idx].toarray().flatten() if hasattr(adata.X, 'toarray') else adata.X[malig_mask, idx]
    malig_expr[g] = e

# NPSR1+ 细胞内, 其他基因的共表达比例
npsr1_pos = malig_expr['NPSR1'] > 0
print(f"NPSR1+ malignant cells: {npsr1_pos.sum()}")
print(f"\n在 NPSR1+ cells 中, 其他神经基因的共表达比例:")
for g in neural_genes:
    if g != 'NPSR1':
        coexp = ((malig_expr.loc[npsr1_pos, g] > 0).sum() / npsr1_pos.sum() * 100)
        global_pct = (malig_expr[g] > 0).mean() * 100
        enrichment = coexp / global_pct if global_pct > 0 else 0
        print(f"  {g:8s} | NPSR1+ 内表达 {coexp:5.2f}% | 全局 {global_pct:5.2f}% | 富集 {enrichment:.2f}x")

print("\n" + "=" * 80)
print("✅ NPSR1 快速验证完成")
print("=" * 80)
