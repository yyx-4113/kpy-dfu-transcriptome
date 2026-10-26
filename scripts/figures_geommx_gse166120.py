# -*- coding: utf-8 -*-
"""
figures_geommx_gse166120.py — GSE166120 GeoMx DSP 图件（任务 #11c）

图 1：CTA 板对 KpyRN 的覆盖（按模块）
图 2：M4 IFN-I 模块评分在 ROI 位置上的分布（受试体级点 + ROI 级散点）
图 3：板内 3 个枢纽基因（LDHA / ADH1B / GBP1）的 ROI 位置分布
图 4：受试体级 vs ROI 级 p 值对照（展示 pseudoreplication 的影响）

所有统计口径以 region_stats.csv 为准（主口径=受试体级）。
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

BASE = _ROOT
RES = f"{BASE}/03_结果/空间GSE166120"
FIG = f"{RES}/figures"
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
                     "figure.dpi": 120, "savefig.dpi": 200, "font.family": "DejaVu Sans"})

cov = pd.read_csv(f"{RES}/target_coverage_modules.csv")
hubcov = pd.read_csv(f"{RES}/target_coverage_hubgenes.csv")
ms = pd.read_csv(f"{RES}/module_scores_by_roi.csv")
hub = pd.read_csv(f"{RES}/hubgene_expression.csv")
st = pd.read_csv(f"{RES}/region_stats.csv")

LOC_ORDER = ["epidermis", "non-ulcer", "ulcer edge", "ulcer"]
LOC_COL = {"epidermis": "#7D3C98", "non-ulcer": "#1E8449",
           "ulcer edge": "#D68910", "ulcer": "#C0392B"}

# ---------- 图 1：CTA 板覆盖 ----------
fig, ax = plt.subplots(figsize=(4.6, 3.0))
mods = cov.sort_values("pct", ascending=False)
y = np.arange(len(mods))
ax.barh(y, mods.n_in_board, color="#2471A3", alpha=0.9, label="genes present on CTA board")
ax.barh(y, mods.n_total - mods.n_in_board, left=mods.n_in_board, color="#D5D8DC",
        alpha=0.9, label="not on CTA board")
for i, (_, r) in enumerate(mods.iterrows()):
    ax.text(r.n_total + 1, i, f"{r.n_in_board}/{r.n_total} ({r.pct}%)",
            va="center", fontsize=6.8)
ax.set_yticks(y)
ax.set_yticklabels([m.replace("_", "\n") for m in mods.module], fontsize=6.8)
ax.set_xlabel("KpyRN scoring genes")
ax.set_xlim(0, 100)
ax.set_title("GeoMx CTA v2.0 board coverage of KpyRN\n(targeted panel, 1825 targets)",
             fontsize=8.6)
ax.legend(fontsize=6.4, loc="lower right", framealpha=0.92)
fig.tight_layout()
fig.savefig(f"{FIG}/fig_geommx_board_coverage.png", bbox_inches="tight")
plt.close(fig)
print("图 1 已生成：fig_geommx_board_coverage.png")

# ---------- 图 2：M4 分数按 ROI 位置 ----------
m4 = ms[ms.module == "M4_IFNI_readout"].copy()
fig, ax = plt.subplots(figsize=(5.4, 3.4))
rng = np.random.default_rng(42)
for i, lv in enumerate(LOC_ORDER):
    s = m4[m4.roi_location == lv]
    if len(s) == 0:
        continue
    x = np.full(len(s), i) + rng.uniform(-0.14, 0.14, len(s))
    ax.scatter(x, s.zscore, s=26, color=LOC_COL[lv], alpha=0.75, linewidths=0,
               label=f"{lv} (n={len(s)} ROI)")
    # 受试体级均值
    sm = s.groupby("subject", observed=True).zscore.mean()
    ax.scatter(np.full(len(sm), i), sm.values, s=110, marker="D",
               facecolors="none", edgecolors="#1B2631", linewidths=1.3, zorder=5,
               label=f"{lv} subject mean" if i == 0 else None)
ax.axhline(0, lw=0.7, ls="--", c="#7F8C8D")
ax.set_xticks(range(len(LOC_ORDER)))
ax.set_xticklabels(LOC_ORDER, fontsize=7.4)
ax.set_ylabel("M4 IFN-I readout z-score\n(board subset, 26/31 genes)")
ax.set_title("GSE166120 GeoMx: M4 IFN-I readout across DFU ROI regions\n"
             "ROI points (n=23) with subject-level means (n=4) overlaid", fontsize=8.6)
ax.legend(fontsize=6.2, loc="best", framealpha=0.92, ncol=2)
fig.tight_layout()
fig.savefig(f"{FIG}/fig_geommx_m4_by_region.png", bbox_inches="tight")
plt.close(fig)
print("图 2 已生成：fig_geommx_m4_by_region.png")

# ---------- 图 3：板内 3 个枢纽基因 ----------
genes_in = hubcov[hubcov.in_board].gene.tolist()
fig, axes = plt.subplots(1, len(genes_in), figsize=(3.0 * len(genes_in), 3.2),
                         sharey=False)
if len(genes_in) == 1:
    axes = [axes]
for ax, g in zip(axes, genes_in):
    s = hub[hub.gene == g]
    order = [lv for lv in LOC_ORDER if (s.roi_location == lv).any()]
    data = [s.loc[s.roi_location == lv, "rel_abund"].values * 1e3 for lv in order]
    bp = ax.boxplot(data, widths=0.55, patch_artist=True,
                    medianprops=dict(color="#1B2631", lw=1.1))
    for patch, lv in zip(bp["boxes"], order):
        patch.set_facecolor(LOC_COL.get(lv, "#AAB7B8"))
        patch.set_alpha(0.75)
    for j, d in enumerate(data):
        ax.scatter(np.full(len(d), j + 1) + rng.uniform(-0.12, 0.12, len(d)), d,
                   s=20, c="#1B2631", alpha=0.65, linewidths=0)
    ax.set_title(g, fontsize=8.6)
    ax.set_ylabel("relative abundance (x1e-3)" if g == genes_in[0] else "")
    ax.tick_params(labelsize=6.6)
    ax.set_xticklabels(order, fontsize=6.4, rotation=30, ha="right")
fig.suptitle("GSE166120 GeoMx: ML-consensus hub genes present on CTA board\n"
             "(3 of 6: FAHD2B / DHX58 / ALDH3A2 are not on the panel)",
             fontsize=8.8)
fig.tight_layout(rect=[0, 0, 1, 0.90])
fig.savefig(f"{FIG}/fig_geommx_hubgenes.png", bbox_inches="tight")
plt.close(fig)
print("图 3 已生成：fig_geommx_hubgenes.png")

# ---------- 图 4：主口径 vs 敏感性口径 p 值对照 ----------
fig, ax = plt.subplots(figsize=(5.2, 3.0))
x = np.arange(len(st))
ax.scatter(x, st.p_roi, s=42, c="#C0392B", alpha=0.85, label="ROI level (n=23, "
           "pseudoreplication-prone)", zorder=3)
ax.scatter(x, st.p_subject, s=42, c="#2471A3", alpha=0.9, marker="D",
           label="subject level (n=4, primary)", zorder=4)
ax.axhline(0.05, ls="--", lw=0.8, c="#7F8C8D")
ax.set_xticks(x)
ax.set_xticklabels([f"{r.module.split('_')[0]}\n{r.grouping[:8]}" for _, r in st.iterrows()],
                   fontsize=6.0)
ax.set_ylabel("nominal p")
ax.set_ylim(-0.03, 1.05)
ax.set_title("Pseudoreplication check: ROI-level vs subject-level p values\n"
             "no comparison is significant after FDR at either level", fontsize=8.6)
ax.legend(fontsize=6.4, loc="upper right", framealpha=0.92)
fig.tight_layout()
fig.savefig(f"{FIG}/fig_geommx_pseudoreplication_check.png", bbox_inches="tight")
plt.close(fig)
print("图 4 已生成：fig_geommx_pseudoreplication_check.png")

print("DONE — figures_geommx_gse166120.py")
