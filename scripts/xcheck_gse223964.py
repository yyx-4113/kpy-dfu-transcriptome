# -*- coding: utf-8 -*-
"""
xcheck_gse223964.py — GSE223964 独立代码路径交叉验证（任务 #11b 补充）

目的：对 scrna_gse223964.py 的供体级 pseudobulk 结果（32 对比全阴性）做独立复核，
      并单列细胞级描述性统计（明确标注 pseudoreplication 限制）。

路径 A（主脚本）：供体×细胞类型 pseudobulk（log1p 均值）→ Mann-Whitney U (exact)，n=5 vs 3
路径 B（本脚本）：
  B1 细胞级 Mann-Whitney U (asymptotic) + 阳性细胞比例差 —— **描述性/一致性证据**，
     细胞嵌套于供体，非独立观测，p 值存在 pseudoreplication 风险，不作独立推断；
  B2 供体级敏感性分析：细胞类型内中位数 → 供体级 pseudobulk → 置换检验（10000 次）
     —— 与路径 A 不同的统计量与不同的检验，用于确认阴性结论的方向性。

输出：03_结果/单细胞GSE223964/
  xcheck_hubgene_celllevel.csv      路径B1 枢纽基因细胞级描述性统计
  xcheck_pseudobulk_perm.csv        路径B2 供体级置换检验敏感性分析
"""
import warnings

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

BASE = _ROOT
OUT = f"{BASE}/03_结果/单细胞GSE223964"
RNG = np.random.default_rng(42)

a = sc.read_h5ad(f"{BASE}/01_数据/原始数据/GSE223964/GSE223964_allcell.h5ad")
a.var_names_make_unique()
a.obs["cond"] = a.obs["condition"].astype(str).str.replace(" ", "_")
a.obs["sample"] = a.obs["sample"].astype(str)
CT = "Cell Type"
HUB = ["LDHA", "FAHD2B", "ADH1B", "DHX58", "ALDH3A2", "GBP1"]

# ---------- B1 细胞级描述性统计 ----------
cond_v = a.obs["cond"].values
d_mask = cond_v == "diabetic"
n_mask = cond_v == "non_diabetic"
rows = []
for g in HUB:
    x = a[:, g].X
    x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
    x = x.ravel()
    u, p = mannwhitneyu(x[d_mask], x[n_mask], alternative="two-sided", method="asymptotic")
    fd = float((x[d_mask] > 0).mean())
    fn = float((x[n_mask] > 0).mean())
    rows.append({"gene": g, "n_cells_diabetic": int(d_mask.sum()),
                 "n_cells_nondiabetic": int(n_mask.sum()),
                 "frac_pos_diabetic": round(fd, 4), "frac_pos_nondiabetic": round(fn, 4),
                 "frac_diff": round(fd - fn, 4), "celllevel_p": p})
b1 = pd.DataFrame(rows)
b1["celllevel_p_bf"] = multipletests(b1["celllevel_p"], method="fdr_bh")[1]
b1["interpretation"] = ("descriptive only: cells nested within donors -> "
                         "pseudoreplication risk, NOT independent inference")
b1.to_csv(f"{OUT}/xcheck_hubgene_celllevel.csv", index=False)
print("== B1 细胞级（描述性）==")
print(b1.to_string(index=False))

# ---------- B2 供体级置换检验敏感性分析 ----------
K = pd.read_csv(f"{BASE}/01_数据/基因集/KpyRN_modules.csv")
K = K[K["scoring_include"] == "Yes"]
sc_rows = []
for mod in sorted(K["module"].unique()):
    genes = [g for g in K.loc[K["module"] == mod, "gene_symbol"] if g in a.var_names]
    if len(genes) < 3:
        continue
    sc.tl.score_genes(a, gene_list=genes, score_name=f"s_{mod}",
                      ctrl_size=min(len(genes) * 5, 5000), n_bins=25, random_state=42)
    sc_rows.append(f"s_{mod}")
print("\n模块打分完成:", sc_rows)

obs = a.obs[["sample", "cond", CT] + sc_rows].copy()
obs = obs.astype({"sample": str, "cond": str, CT: str})
n_perm = 10000
perm_rows = []
for ct in sorted(obs[CT].unique()):
    sub = obs[obs[CT] == ct]
    if sub["sample"].nunique() < 4:
        continue
    donors = sorted(sub["sample"].unique())
    cond_of = {d: sub.loc[sub["sample"] == d, "cond"].iloc[0] for d in donors}
    for mod in sc_rows:
        # 供体级伪 bulk = 该供体在此细胞类型内的 score 中位数（与路径A的均值不同统计量）
        pb = sub.groupby("sample", observed=True)[mod].median()
        vals = pb.values
        d_idx = [i for i, dd in enumerate(donors) if cond_of[dd] == "diabetic"]
        n_idx = [i for i, dd in enumerate(donors) if cond_of[dd] == "non_diabetic"]
        if len(d_idx) < 2 or len(n_idx) < 2:
            continue
        obs_diff = vals[d_idx].mean() - vals[n_idx].mean()
        # 置换检验：打乱供体标签 10000 次
        cnt = 0
        for _ in range(n_perm):
            perm = RNG.permutation(len(vals))
            pd_ = perm[:len(d_idx)]
            pn_ = perm[len(d_idx):len(d_idx) + len(n_idx)]
            if abs(vals[pd_].mean() - vals[pn_].mean()) >= abs(obs_diff) - 1e-12:
                cnt += 1
        p_perm = (cnt + 1) / (n_perm + 1)
        perm_rows.append({"celltype": ct, "module": mod.replace("s_", ""),
                          "n_donors_diabetic": len(d_idx), "n_donors_nondiabetic": len(n_idx),
                          "observed_diff": obs_diff, "perm_p": p_perm})
b2 = pd.DataFrame(perm_rows)
b2["perm_p_bf"] = multipletests(b2["perm_p"], method="fdr_bh")[1]
b2["note"] = ("donor-level median pseudobulk, 10000-label permutation, "
              "different statistic from path A (mean); sensitivity analysis")
b2.to_csv(f"{OUT}/xcheck_pseudobulk_perm.csv", index=False)
print("\n== B2 供体级置换检验（敏感性分析）==")
print(f"对比数: {len(b2)} | perm_p<0.10: {(b2.perm_p<0.10).sum()} | FDR<0.10: {(b2.perm_p_bf<0.10).sum()}")
print(b2.reindex(b2["observed_diff"].abs().sort_values(ascending=False).index)
      .head(10).round(4).to_string(index=False))
print("\nDONE — xcheck_gse223964.py")
