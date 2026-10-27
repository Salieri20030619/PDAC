#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE202051 Schwann 细胞聚焦探查
目的: 评估 Schwann 细胞做亚分群和治疗对比的可行性
"""

import os, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

H5AD_PATH = '/mnt/BioDisk/wangmingyuan/PDAC_spatial_PNI/00_raw_data/GSE202051/GSE202051_totaldata-final-toshare.h5ad'
OUT_DIR = '/mnt/BioDisk/wangmingyuan/PDAC_spatial_PNI/02_results/01_data_exploration'
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 80)
print("Schwann 细胞聚焦探查")
print("=" * 80)

print("\n[1] 加载数据...")
adata = sc.read_h5ad(H5AD_PATH)
print(f"  总 obs: {adata.n_obs:,}")

# ========== 1. Schwann 细胞分布 ==========
print("\n[2] Schwann 细胞分布")
print("-" * 60)
schwann_mask = adata.obs['Level 1 Annotation'] == 'Schwann'
schwann_adata = adata[schwann_mask].copy()
print(f"  Schwann 细胞总数: {schwann_adata.n_obs:,}")

# Schwann 在每个病人里的数量
print(f"\n  Schwann 在每个病人 (pid) 的分布:")
pid_schwann = schwann_adata.obs.groupby('pid', observed=True).size().sort_values(ascending=False)
print(pid_schwann)
print(f"\n  有 Schwann 的病人数: {(pid_schwann > 0).sum()} / 43")
print(f"  平均每个病人 Schwann: {pid_schwann.mean():.1f}")
print(f"  最多: {pid_schwann.max()}, 最少: {pid_schwann.min()}")

# Schwann 在治疗组的分布
print(f"\n  Schwann 在不同 treatment_status 的分布:")
print(schwann_adata.obs['treatment_status'].value_counts())

print(f"\n  Schwann 在 new_treatment (4 类) 的分布:")
print(schwann_adata.obs['new_treatment'].value_counts())

print(f"\n  Schwann 在 response 的分布:")
print(schwann_adata.obs['response'].value_counts())

# ========== 2. Schwann 内部异质性 (能不能再分群) ==========
print("\n[3] Schwann 细胞内部异质性 - 关键基因表达")
print("-" * 60)

# Schwann 亚型相关基因
schwann_subtype_markers = {
    'Myelinating': ['MPZ', 'MBP', 'PMP22', 'PRX', 'EGR2', 'MAG'],
    'Non-myelinating': ['NCAM1', 'GAP43', 'GFRA3', 'L1CAM'],
    'Repair_Schwann': ['SHH', 'NGFR', 'GDNF', 'JUN', 'OLIG1'],
    'Activated_Schwann': ['S100B', 'SOX10', 'PLP1', 'NES', 'VIM'],
    'Pro_tumor_PNI': ['CDH2', 'L1CAM', 'NRP1', 'NTRK2', 'BDNF']
}

print(f"  Schwann 细胞数: {schwann_adata.n_obs}")
gene_names = schwann_adata.var_names.astype(str).values
for subtype, markers in schwann_subtype_markers.items():
    print(f"\n  --- {subtype} ---")
    for m in markers:
        if m in gene_names:
            idx = np.where(gene_names == m)[0][0]
            expr = schwann_adata.X[:, idx].toarray().flatten() if hasattr(schwann_adata.X, 'toarray') else schwann_adata.X[:, idx]
            pct = (expr > 0).sum() / len(expr) * 100
            mean_expr = expr[expr > 0].mean() if (expr > 0).any() else 0
            print(f"    {m:8s} | 表达 {pct:5.2f}% | 阳性均值 {mean_expr:.3f}")
        else:
            print(f"    {m:8s} | ❌ 不在数据中")

# ========== 3. Schwann 已有的细致注释 ==========
print("\n[4] Schwann 在 new_celltypes / Level 2/3 里有没有更细的分类")
print("-" * 60)
for col in ['new_celltypes', 'Level 2 Annotation', 'Level 3 Annotation', 'leiden']:
    if col in schwann_adata.obs.columns:
        vc = schwann_adata.obs[col].value_counts()
        vc = vc[vc > 0]
        print(f"\n  {col} 在 Schwann 中的分布:")
        print(vc)

# ========== 4. Neural niche 关键基因在 malignant cells 中的情况 ==========
print("\n[5] PNI 关键基因在 malignant cells 中的表达 (NRP1/NTRK 等)")
print("-" * 60)
malig_mask = adata.obs['Level 1 Annotation'] == 'Epithelial (malignant)'
malig_adata = adata[malig_mask].copy()
print(f"  Malignant cells 总数: {malig_adata.n_obs:,}")

pni_targets = ['NRP1', 'NRP2', 'NTRK1', 'NTRK2', 'NTRK3', 'NGFR', 
               'L1CAM', 'CDH2', 'GAP43', 'NCAM1', 'BDNF', 'NGF',
               'GDNF', 'SEMA3A', 'SEMA3E', 'PLXNA1']

gene_names = malig_adata.var_names.astype(str).values
print(f"\n  在 {malig_adata.n_obs:,} 个 malignant cells 中:")
for g in pni_targets:
    if g in gene_names:
        idx = np.where(gene_names == g)[0][0]
        expr = malig_adata.X[:, idx].toarray().flatten() if hasattr(malig_adata.X, 'toarray') else malig_adata.X[:, idx]
        pct = (expr > 0).sum() / len(expr) * 100
        mean_expr = expr[expr > 0].mean() if (expr > 0).any() else 0
        print(f"    {g:8s} | 表达 {pct:5.2f}% | 阳性均值 {mean_expr:.3f}")

# ========== 5. 每个治疗组的关键细胞数量 ==========
print("\n[6] 每个治疗组的关键细胞数量 (评估统计 power)")
print("-" * 60)
print("\n按 new_treatment (4 类) × Level 1 Annotation 交叉表:")
ct = pd.crosstab(adata.obs['Level 1 Annotation'], adata.obs['new_treatment'])
print(ct)
ct.to_csv(os.path.join(OUT_DIR, 'celltype_x_treatment.csv'))

print("\n按 response × Level 1 Annotation 交叉表:")
ct2 = pd.crosstab(adata.obs['Level 1 Annotation'], adata.obs['response'])
print(ct2)
ct2.to_csv(os.path.join(OUT_DIR, 'celltype_x_response.csv'))

# ========== 6. 病人级别的样本配置 ==========
print("\n[7] 病人级别配置")
print("-" * 60)
patient_info = adata.obs.groupby('pid', observed=True).agg(
    n_cells=('sampleid', 'size'),
    sampleid=('sampleid', 'first'),
    treatment_status=('treatment_status', 'first'),
    new_treatment=('new_treatment', 'first'),
    response=('response', 'first'),
).sort_values('n_cells', ascending=False)
print(patient_info)
patient_info.to_csv(os.path.join(OUT_DIR, 'patient_info.csv'))

print("\n" + "=" * 80)
print("✅ Schwann 聚焦探查完成")
print("=" * 80)
