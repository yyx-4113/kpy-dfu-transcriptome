# -*- coding: utf-8 -*-
"""
scrna_gse223964.py — GSE223964 单细胞验证（任务 #11b）
项目：丙酮酰化（Kpy）与糖尿病足（DFU）

数据：GSE223964_allcell.h5ad（Theocharidis Nat Commun 2022 同篇；diabetic n=5 vs
      non-diabetic n=3 皮肤，27,998 细胞 × 18,287 基因，X=log1p，自带 Cell Type 注释）
定位：验证 KpyRN/枢纽基因在"糖尿病（非溃疡）皮肤"维度与细胞类型特异性；
      与 GSE165816（DFU vs 对照）互补。n=5 vs 3 → 结局对比仅探索级（如实标注）。

注意：无 raw counts（raw=False），X 为 log1p。伪 bulk 以"样本×细胞类型内 log1p
均值"实现（非 counts 求和），检验为探索性 Wilcoxon，稿件中如实声明。

输出：03_结果/单细胞GSE223964/
  module_scores_by_celltype.csv  模块打分（细胞类型 × 状态）
  hubgene_median_by_celltype.csv 枢纽基因中位表达
  pseudobulk_condition_stats.csv 样本×细胞类型伪 bulk 状态对比（探索性）
  fig_scores.png / fig_hubgenes.png
"""
import warnings

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")

BASE = _ROOT
OUT = f"{BASE}/03_结果/单细胞GSE223964"
import os
os.makedirs(OUT, exist_ok=True)

a = sc.read_h5ad(f"{BASE}/01_数据/原始数据/GSE223964/GSE223964_allcell.h5ad")
a.var_names_make_unique()
print(f"载入: {a.shape}")

K = pd.read_csv(f"{BASE}/01_数据/基因集/KpyRN_modules.csv")
K = K[K["scoring_include"] == "Yes"]
HUB = ["LDHA", "FAHD2B", "ADH1B", "DHX58", "ALDH3A2", "GBP1"]
present_hubs = [g for g in HUB if g in a.var_names]
print("枢纽基因在矩阵中:", present_hubs, "(缺:", set(HUB) - set(present_hubs), ")")

CT = "Cell Type"
cond = a.obs["condition"].astype(str).str.replace(" ", "_")
a.obs["cond"] = cond

# ---------- 1. 模块打分（与 GSE165816 一致：score_genes，ctrl=基因数×5，n_bins=25）----------
score_rows = []
for mod in sorted(K["module"].unique()):
    genes = [g for g in K.loc[K["module"] == mod, "gene_symbol"] if g in a.var_names]
    if len(genes) < 3:
        print(f"  !! {mod} 仅 {len(genes)} 基因在矩阵中，跳过")
        continue
    sc.tl.score_genes(a, gene_list=genes, score_name=f"s_{mod}",
                      ctrl_size=min(len(genes) * 5, 5000), n_bins=25,
                      random_state=42)
    sub = a.obs[[CT, "cond", f"s_{mod}"]]
    g = sub.groupby([CT, "cond"], observed=True)[f"s_{mod}"].agg(["mean", "median", "count"])
    for (ct, cd), r in g.iterrows():
        score_rows.append({"module": mod, "celltype": ct, "condition": cd,
                           "mean_score": r["mean"], "median_score": r["median"],
                           "n_cells": int(r["count"])})
    print(f"  {mod}: {len(genes)} 基因打分完成")
sc_df = pd.DataFrame(score_rows)
sc_df.to_csv(f"{OUT}/module_scores_by_celltype.csv", index=False)

# ---------- 2. 枢纽基因中位表达（细胞类型 × 状态）----------
raw_x = a[:, present_hubs].X
raw_x = raw_x.toarray() if hasattr(raw_x, "toarray") else np.asarray(raw_x)
hub_df = pd.DataFrame(raw_x, columns=present_hubs)
hub_df[CT] = a.obs[CT].values
hub_df["cond"] = a.obs["cond"].values
hg = (hub_df.groupby([CT, "cond"], observed=True)[present_hubs]
      .median().round(4).reset_index())
hg.to_csv(f"{OUT}/hubgene_median_by_celltype.csv", index=False)

# ---------- 3. 伪 bulk（样本 × 细胞类型 log1p 均值）+ 探索性 Wilcoxon ----------
obs = a.obs.copy()
pb_rows = []
genes_all = sorted(set(K["gene_symbol"]) & set(a.var_names) | set(present_hubs))
Xg = a[:, genes_all].X
Xg = Xg.toarray() if hasattr(Xg, "toarray") else np.asarray(Xg)
pb = pd.DataFrame(Xg, columns=genes_all)
pb["sample"] = a.obs["sample"].astype(str).values
pb["cond"] = a.obs["cond"].values
pb[CT] = a.obs[CT].values
for (ct, mod) in [(ct, m) for ct in a.obs[CT].unique() for m in sorted(K["module"].unique())]:
    genes_m = [g for g in K.loc[K["module"] == mod, "gene_symbol"] if g in pb.columns]
    if len(genes_m) < 3:
        continue
    sub = pb[pb[CT] == ct]
    if sub["sample"].nunique() < 4:      # 至少 4 个样本（两组）才比较
        continue
    # observed=True：sample 为分类变量，默认 observed=False 会展开 (sample × cond)
    # 笛卡尔积空组（8 样本 × 2 cond = 16 组，其中 8 组 0 细胞 → NaN 均值），
    # 导致 n_diabetic/n_nondiabetic 被虚报为 8/8、均值与 p 值全空（2026-10-25 修复）
    sm = (sub.groupby(["sample", "cond"], observed=True)[genes_m]
             .mean().mean(axis=1).reset_index())
    sm.columns = ["sample", "cond", "score"]
    d = sm.loc[sm["cond"] == "diabetic", "score"].values
    n = sm.loc[sm["cond"] == "non_diabetic", "score"].values
    if len(d) < 2 or len(n) < 2:
        continue
    # 独立两组样本（非配对）→ Mann-Whitney U；原用 wilcoxon（配对检验）方法不当
    try:
        p = mannwhitneyu(d, n, alternative="two-sided", method="exact").pvalue
    except ValueError:
        p = np.nan
    pb_rows.append({"celltype": ct, "module": mod, "n_diabetic": len(d),
                    "n_nondiabetic": len(n), "mean_diabetic": d.mean(),
                    "mean_nondiabetic": n.mean(),
                    "diff": d.mean() - n.mean(), "mwu_p": p})
pb_df = pd.DataFrame(pb_rows)
from statsmodels.stats.multitest import multipletests
pb_df["mwu_p_bf"] = multipletests(pb_df["mwu_p"], method="fdr_bh")[1]
pb_df["note"] = ("exploratory: n=5 vs 3 donors, pseudobulk=mean log1p per "
                 "sample x cell type, Mann-Whitney U (exact), 无 raw counts")
pb_df.to_csv(f"{OUT}/pseudobulk_condition_stats.csv", index=False)

# ---------- 4. 图 ----------
order = (a.obs[CT].value_counts().index.tolist())
mods = sorted(set(sc_df["module"]))
fig, axes = plt.subplots(1, len(mods), figsize=(3.1 * len(mods), 3.4),
                         sharey=False)
if len(mods) == 1:
    axes = [axes]
for ax, mod in zip(axes, mods):
    sub = sc_df[sc_df["module"] == mod]
    piv = (sub.pivot(index="celltype", columns="condition", values="mean_score")
           .reindex(order))
    ypos = np.arange(len(piv))
    ax.barh(ypos - 0.2, piv.get("diabetic", pd.Series(0, index=piv.index)),
            height=0.4, color="#C0392B", label="diabetic")
    ax.barh(ypos + 0.2, piv.get("non_diabetic", pd.Series(0, index=piv.index)),
            height=0.4, color="#2471A3", label="non-diabetic")
    ax.set_yticks(ypos); ax.set_yticklabels(piv.index, fontsize=7)
    ax.invert_yaxis()
    ax.set_title(mod.replace("_", " "), fontsize=8.5)
    ax.tick_params(labelsize=7)
axes[0].legend(fontsize=7)
fig.suptitle("GSE223964 KpyRN module scores by cell type (diabetic vs non-diabetic skin)",
             fontsize=9.5)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(f"{OUT}/fig_scores.png", dpi=200, bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(6.4, 3.6))
sub = hg[hg["cond"] == "diabetic"].set_index(CT).reindex(order)[present_hubs]
sub2 = hg[hg["cond"] == "non_diabetic"].set_index(CT).reindex(order)[present_hubs]
x = np.arange(len(present_hubs)); w = 0.38
for i, (ct) in enumerate(order):
    ax.plot(x, sub.loc[ct], "o-", color="#C0392B", alpha=0.65, lw=1)
    ax.plot(x, sub2.loc[ct], "s--", color="#2471A3", alpha=0.65, lw=1)
ax.plot([], [], "o-", color="#C0392B", label="diabetic (median)")
ax.plot([], [], "s--", color="#2471A3", label="non-diabetic (median)")
ax.set_xticks(x); ax.set_xticklabels(present_hubs, fontsize=8)
ax.set_ylabel("median log1p expression", fontsize=8)
ax.legend(fontsize=7); ax.tick_params(labelsize=7)
fig.suptitle("GSE223964 hub genes by cell type", fontsize=9.5)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(f"{OUT}/fig_hubgenes.png", dpi=200, bbox_inches="tight")
plt.close(fig)

print("== 伪 bulk 状态对比（探索性，FDR<0.10 标注）==")
sig = pb_df[pb_df["mwu_p_bf"] < 0.10]
print(sig.to_string(index=False) if len(sig) else "  无 FDR<0.10 的对比")
print("DONE — scrna_gse223964.py")
