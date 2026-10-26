# -*- coding: utf-8 -*-
"""
figures_gse230426.py — GSE230426 感染队列次级分析图件（任务 #11a 补充）

输入（均由 deg_gse230426.R 产出，本脚本不重新计算任何统计量，只做可视化）：
  03_结果/感染队列GSE230426/deg_main.csv         主模型全基因组 DEG
  03_结果/感染队列GSE230426/deg_paired.csv       推断配对敏感性模型
  03_结果/感染队列GSE230426/deg_after8.csv       仅 8wk 敏感性模型
  03_结果/感染队列GSE230426/camera_modules.csv   camera 模块检验（3 模型）
  03_结果/感染队列GSE230426/hub_genes.csv         6 个 ML 共识枢纽基因 DEG

图 1：火山图（主模型，全基因组）
图 2：KpyRN 128 基因 logFC 在主模型 vs 推断配对模型的一致性（Spearman）
图 3：camera 模块检验 -log10(FDR) 三模型并列（M4/M3 达显著）
图 4：6 个枢纽基因三模型 logFC 森林式点图（全部不显著，如实呈现）

输出：03_结果/感染队列GSE230426/figures/
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

BASE = _ROOT
RES = f"{BASE}/03_结果/感染队列GSE230426"
FIG = f"{RES}/figures"
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
                     "figure.dpi": 120, "savefig.dpi": 200, "font.family": "DejaVu Sans"})

deg_main = pd.read_csv(f"{RES}/deg_main.csv")
deg_pair = pd.read_csv(f"{RES}/deg_paired.csv")
cam = pd.read_csv(f"{RES}/camera_modules.csv")
hub = pd.read_csv(f"{RES}/hub_genes.csv")

K = pd.read_csv(f"{BASE}/01_数据/基因集/KpyRN_modules.csv")
K = K[K["scoring_include"] == "Yes"]
kpy = set(K["gene_symbol"])

# ---------- 图 1：火山图 ----------
fig, ax = plt.subplots(figsize=(4.4, 3.6))
lp = deg_main["logFC"].values
pv = -np.log10(deg_main["P.Value"].clip(lower=1e-300))
is_k = deg_main["gene"].isin(kpy).values
ax.scatter(lp[~is_k], pv[~is_k], s=5, c="#BFC8D2", alpha=0.55, linewidths=0,
           label=f"other genes (n={int((~is_k).sum())})")
ax.scatter(lp[is_k], pv[is_k], s=13, c="#2471A3", alpha=0.85, linewidths=0,
           label=f"KpyRN genes (n={int(is_k.sum())})")
ax.axhline(-np.log10(0.05), ls="--", lw=0.8, c="#C0392B")
ax.axvline(0, lw=0.6, c="#7F8C8D")
ax.set_xlabel("log$_2$FC (AFTER vs BEFORE)")
ax.set_ylabel("$-$log$_{10}$ nominal P")
ax.set_title("GSE230426 volcano (main model, n=18 vs 20)", fontsize=8.6)
ax.legend(fontsize=6.6, framealpha=0.92, loc="upper right")
ax.text(0.02, 0.03, "0 genes FDR<0.05", transform=ax.transAxes, fontsize=6.6,
        color="#C0392B")
fig.tight_layout()
fig.savefig(f"{FIG}/fig_gse230426_volcano.png", bbox_inches="tight")
plt.close(fig)
print("图 1 已生成：fig_gse230426_volcano.png")

# ---------- 图 2：主 vs 推断配对 一致性 ----------
m = deg_main[["gene", "logFC"]].rename(columns={"logFC": "logFC_main"})
p = deg_pair[["gene", "logFC"]].rename(columns={"logFC": "logFC_paired"})
j = m.merge(p, on="gene")
j["is_k"] = j["gene"].isin(kpy)
r_all, _ = spearmanr(j["logFC_main"], j["logFC_paired"])
r_k, _ = spearmanr(j.loc[j.is_k, "logFC_main"], j.loc[j.is_k, "logFC_paired"])
diragree = float((np.sign(j["logFC_main"]) == np.sign(j["logFC_paired"])).mean())
diragree_k = float((np.sign(j.loc[j.is_k, "logFC_main"]) ==
                    np.sign(j.loc[j.is_k, "logFC_paired"])).mean())

fig, ax = plt.subplots(figsize=(4.0, 4.0))
ax.scatter(j.loc[~j.is_k, "logFC_main"], j.loc[~j.is_k, "logFC_paired"],
           s=4, c="#C8D0DA", alpha=0.5, linewidths=0, label="other genes")
ax.scatter(j.loc[j.is_k, "logFC_main"], j.loc[j.is_k, "logFC_paired"],
           s=14, c="#1E8449", alpha=0.85, linewidths=0, label="KpyRN genes")
lim = np.nanmax(np.abs(j[["logFC_main", "logFC_paired"]].values)) * 1.05
ax.plot([-lim, lim], [-lim, lim], lw=0.7, ls="--", c="#7F8C8D")
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
ax.set_xlabel("log$_2$FC — main (unpaired voom)")
ax.set_ylabel("log$_2$FC — paired (inferred)")
ax.set_title(f"Model concordance\nall: Spearman r={r_all:.3f} | "
             f"KpyRN: r={r_k:.3f}", fontsize=8.6)
ax.legend(fontsize=6.8, loc="upper left", framealpha=0.92)
ax.text(0.03, 0.03, f"direction agreement — all genes: {diragree*100:.1f}%\n"
                    f"direction agreement — KpyRN 128: {diragree_k*100:.1f}%",
        transform=ax.transAxes, fontsize=6.6)
fig.tight_layout()
fig.savefig(f"{FIG}/fig_gse230426_model_concordance.png", bbox_inches="tight")
plt.close(fig)
print("图 2 已生成：fig_gse230426_model_concordance.png")

# ---------- 图 3：camera 模块检验 ----------
models = ["main_unpaired", "paired_inferred", "after8wk_only"]
mlab = {"main_unpaired": "main\n(unpaired)", "paired_inferred": "paired\n(inferred)",
        "after8wk_only": "8wk\nonly"}
mods = sorted(cam["module"].unique())
fig, ax = plt.subplots(figsize=(5.4, 3.2))
w = 0.26
x = np.arange(len(mods))
colors = {"main_unpaired": "#2471A3", "paired_inferred": "#1E8449", "after8wk_only": "#7D3C98"}
for i, mo in enumerate(models):
    sub = cam[cam["model"] == mo].set_index("module").reindex(mods)
    vals = -np.log10(sub["fdr"].clip(lower=1e-300)).values
    ax.bar(x + (i - 1) * w, vals, w, label=mlab[mo], color=colors[mo], alpha=0.9)
    for xi, v in zip(x + (i - 1) * w, vals):
        ax.text(xi, v + 0.06, f"{v:.2f}", ha="center", fontsize=5.4, rotation=90)
ax.axhline(-np.log10(0.05), ls="--", lw=0.8, c="#C0392B")
ax.text(len(mods) - 0.45, -np.log10(0.05) + 0.12, "FDR 0.05", fontsize=6.2, color="#C0392B")
ax.set_xticks(x)
ax.set_xticklabels([m.replace("_", "\n") for m in mods], fontsize=6.6)
ax.set_ylabel("$-$log$_{10}$ camera FDR")
ax.set_title("GSE230426 KpyRN module test (camera)\nAFTER (treatment) vs BEFORE "
             "(active infection), PBMC", fontsize=8.6)
ax.set_ylim(0, max(-np.log10(cam["fdr"].clip(lower=1e-300)).max() * 1.45, 2.6))
ax.legend(fontsize=6.4, loc="upper left", framealpha=0.95, ncol=3,
          columnspacing=0.9, handlelength=1.1)
fig.tight_layout()
fig.savefig(f"{FIG}/fig_gse230426_camera_modules.png", bbox_inches="tight")
plt.close(fig)
print("图 3 已生成：fig_gse230426_camera_modules.png")

# ---------- 图 4：枢纽基因三模型 logFC ----------
genes = sorted(hub["gene"].unique())
fig, ax = plt.subplots(figsize=(5.6, 3.2))
y = np.arange(len(genes))
for i, mo in enumerate(models):
    sub = hub[hub["model"] == mo].set_index("gene").reindex(genes)
    ax.scatter(sub["logFC"], y + (i - 1) * 0.18, s=26, color=colors[mo],
               label=mlab[mo], zorder=3)
ax.axvline(0, lw=0.8, c="#7F8C8D")
ax.set_yticks(y)
ax.set_yticklabels(genes, fontsize=7.4)
ax.set_xlabel("log$_2$FC (AFTER vs BEFORE)")
ax.set_title("GSE230426 ML-consensus hub genes in PBMC\nnone reaches FDR<0.05 in "
             "any model (negative result)", fontsize=8.6)
ax.legend(fontsize=6.4, loc="lower right", framealpha=0.92)
ax.set_ylim(-0.6, len(genes) - 0.2)
fig.tight_layout()
fig.savefig(f"{FIG}/fig_gse230426_hub_genes.png", bbox_inches="tight")
plt.close(fig)
print("图 4 已生成：fig_gse230426_hub_genes.png")

print("DONE — figures_gse230426.py")
