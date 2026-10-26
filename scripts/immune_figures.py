# -*- coding: utf-8 -*-
"""
immune_figures.py — KpyRN 模块评分 × 免疫浸润 + 免疫图表生成
项目：丙酮酰化（Kpy）与糖尿病足（DFU）

内容：
  1. KpyRN 五模块评分（M1 机器 / M2 丙酮酸供应 / M3 靶基因 / M4 IFN-I 读出；
     M5 组蛋白载体按既定规则不纳入评分）—— 每样本取模块内基因 log2 表达均值
  2. 模块评分 × 免疫细胞比例 Spearman 相关 + BH FDR
  3. 模块评分组间比较（GSE134431: Healer/Nonhealer/DFS；GSE199939: DFU/control）
  4. 图表：
     fig_immune_group_boxplots.png  — 组间差异显著的免疫指标箱线图（两队列）
     fig_hubgene_immune_heatmap.png — 枢纽基因 × 免疫指标 Spearman 相关热图（两队列）
     fig_method_concordance.png    — EPIC vs MCPcounter 一致性散点
     fig_kpyrn_module_immune.png   — KpyRN 模块评分 × 免疫指标相关热图
"""
from pathlib import Path

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))


import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(_ROOT)
PROC = BASE / "01_数据/处理后数据"
OUT = BASE / "03_结果/免疫浸润"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8, "figure.dpi": 300,
    "axes.spines.top": False, "axes.spines.right": False,
})

# 中国惯例：不使用红涨绿跌；分组箱线图用中性蓝橙
C_CTRL, C_DFU = "#4C72B0", "#DD8452"


def log(msg):
    print(msg, flush=True)
    with open(OUT / "immune_python.log", "a", encoding="utf-8") as fh:
        fh.write(str(msg) + "\n")


def spearman_fdr(x: np.ndarray, y: np.ndarray):
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan, np.nan
    rho, p = stats.spearmanr(x, y)
    return rho, p


# ---------------- 1. KpyRN 模块评分 ----------------
mods = pd.read_csv(BASE / "01_数据/基因集/KpyRN_modules.csv")
mods = mods[mods["scoring_include"].astype(str).str.strip().str.lower() == "yes"]
module_genes = {m: g["gene_symbol"].tolist() for m, g in mods.groupby("module")}
log(f"KpyRN 评分模块: { {k: len(v) for k, v in module_genes.items()} }")

cohort_info = {
    "GSE134431": {"contrasts": [("Healer", "Nonhealer", "healing"),
                                ("DFS", "DFU", "group2")]},
    "GSE199939": {"contrasts": [("DFU", "control", "group2")]},
}

module_scores, frac_all, ann_all = {}, {}, {}
for gse in cohort_info:
    ann = pd.read_csv(PROC / f"{gse}_annotation.csv").set_index("sample")
    expr = pd.read_csv(PROC / f"{gse}_log2_expr.csv", index_col=0)
    ms = pd.DataFrame(index=expr.columns)
    for m, genes in module_genes.items():
        present = [g for g in genes if g in expr.index]
        ms[m] = expr.loc[present].mean(axis=0).values
        log(f"[{gse}] {m}: {len(present)}/{len(genes)} 基因用于评分")
    ms.to_csv(OUT / f"{gse}_kpyrn_module_scores.csv")
    module_scores[gse] = ms
    frac_all[gse] = pd.read_csv(OUT / f"{gse}_estimate_scores.csv", index_col=0).join(
        pd.read_csv(OUT / f"{gse}_epic_fractions.csv", index_col=0).T.add_prefix("EPIC_")
    ).join(pd.read_csv(OUT / f"{gse}_mcpcounter_fractions.csv", index_col=0).add_prefix("MCP_"))
    ann_all[gse] = ann

# ---------------- 2-3. 模块评分 × 免疫相关 + 组间比较 ----------------
ms_corr_rows, ms_stat_rows = [], []
for gse in cohort_info:
    ms, frac, ann = module_scores[gse], frac_all[gse], ann_all[gse]
    samples = [s for s in ms.index if s in frac.index]
    for m in ms.columns:
        x = ms.loc[samples, m].values
        for met in frac.columns:
            y = frac.loc[samples, met].values
            rho, p = spearman_fdr(x, y)
            ms_corr_rows.append({"gse": gse, "module": m, "metric": met,
                                 "rho": rho, "p_raw": p})
    for g1, g2, gcol in cohort_info[gse]["contrasts"]:
        grp = ann[gcol].fillna("")
        s1 = [s for s in samples if grp.get(s) == g1]
        s2 = [s for s in samples if grp.get(s) == g2]
        if len(s1) < 3 or len(s2) < 3:
            continue
        for m in ms.columns:
            a, b = ms.loc[s1, m].values, ms.loc[s2, m].values
            try:
                p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
            except ValueError:
                p = np.nan
            ms_stat_rows.append({"gse": gse, "comparison": f"{g1}_vs_{g2}", "module": m,
                                 "mean_g1": a.mean(), "mean_g2": b.mean(), "p_raw": p})

ms_corr = pd.DataFrame(ms_corr_rows)
ok = ms_corr["p_raw"].notna()
ms_corr.loc[ok, "p_fdr"] = multipletests(ms_corr.loc[ok, "p_raw"], method="fdr_bh")[1]
ms_corr.to_csv(OUT / "kpyrn_module_immune_corr.csv", index=False)

ms_stat = pd.DataFrame(ms_stat_rows)
ok = ms_stat["p_raw"].notna()
ms_stat.loc[ok, "p_fdr"] = multipletests(ms_stat.loc[ok, "p_raw"], method="fdr_bh")[1]
ms_stat.to_csv(OUT / "kpyrn_module_group_stats.csv", index=False)

for gse in cohort_info:
    sub = ms_corr[(ms_corr["gse"] == gse) & (ms_corr["p_fdr"] < 0.05)]
    log(f"[{gse}] KpyRN 模块 × 免疫 FDR<0.05: {len(sub)} 对")
    for _, r in sub.sort_values("p_fdr").head(10).iterrows():
        log(f"    {r['module']:22s} ~ {r['metric']:22s} rho={r['rho']:+.3f} FDR={r['p_fdr']:.2e}")
    sub2 = ms_stat[(ms_stat["gse"] == gse) & (ms_stat["p_fdr"] < 0.05)]
    for _, r in sub2.iterrows():
        log(f"    组间[{r['comparison']}] {r['module']:22s} "
            f"{r['mean_g1']:.2f} vs {r['mean_g2']:.2f} FDR={r['p_fdr']:.3f}")

# ---------------- 4. 图表 ----------------
# 4a. 组间差异显著免疫指标箱线图
STATS_FILE = {"GSE134431": "GSE134431_DFS_vs_DFU_immune_stats.csv",
              "GSE199939": "GSE199939_DFU_vs_control_immune_stats.csv"}
fig, axes = plt.subplots(2, 6, figsize=(13, 5.2))
for col, (gse, g1, g2, gcol, color_pair) in enumerate([
        ("GSE134431", "DFS", "DFU", "group2", (C_CTRL, C_DFU)),
        ("GSE199939", "control", "DFU", "group2", (C_CTRL, C_DFU))]):
    ann = ann_all[gse]
    grp = ann[gcol].fillna("")
    frac = frac_all[gse]
    stat = pd.read_csv(OUT / STATS_FILE[gse])
    sig_metrics = stat[stat["p_fdr"] < 0.05]["metric"].tolist()[:6]
    for i, met in enumerate(sig_metrics):
        ax = axes[col, i]
        data = [frac.loc[grp[grp == g1].index, met].dropna().values,
                frac.loc[grp[grp == g2].index, met].dropna().values]
        bp = ax.boxplot(data, tick_labels=[f"{g1}\n(n={len(data[0])})", f"{g2}\n(n={len(data[1])})"],
                        patch_artist=True, widths=0.55, showfliers=False)
        for patch, c in zip(bp["boxes"], color_pair):
            patch.set_facecolor(c)
            patch.set_alpha(0.75)
        for j, d in enumerate(data):
            ax.scatter(np.full(len(d), j + 1) + np.random.default_rng(j).uniform(-.08, .08, len(d)),
                       d, s=6, color="black", alpha=0.5, zorder=3)
        prow = stat[stat["metric"] == met]["p_fdr"].iloc[0]
        ax.set_title(met.replace("EPIC_", "EPIC:").replace("MCP_", "MCP:")
                     + f"\nFDR={prow:.3f}", fontsize=7)
        ax.tick_params(labelsize=7)
        if i == 0:
            ax.set_ylabel(gse, fontsize=9)
fig.suptitle("Immune cell fractions: ulcer vs reference (Wilcoxon + BH FDR)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(FIG / "fig_immune_group_boxplots.png", bbox_inches="tight")
plt.close(fig)
log("fig_immune_group_boxplots.png saved")

# 4b. 枢纽基因 × 免疫相关热图
HUB = ["ADH1B", "ALDH3A2", "BPGM", "DHX58", "EP300", "FAHD2B", "GBP1",
       "HAT1", "LDHA", "PFKP", "PGAM1", "PKM2", "SIRT3", "STAT1", "TPI1"]
fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
for ax, gse in zip(axes, ["GSE134431", "GSE199939"]):
    corr = pd.read_csv(OUT / f"{gse}_hubgene_immune_corr.csv")
    piv = corr.pivot(index="gene", columns="metric", values="rho")
    piv = piv.reindex([g for g in HUB if g in piv.index])
    fdr = corr.pivot(index="gene", columns="metric", values="p_fdr").reindex(piv.index)[piv.columns]
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels([m.replace("EPIC_", "E:").replace("MCP_", "M:").replace("EST_", "S:")
                        for m in piv.columns], rotation=90, fontsize=6.5)
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels(piv.index, fontsize=7)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            if fdr.values[i, j] < 0.05:
                ax.text(j, i, "*", ha="center", va="center", fontsize=6, color="black")
    ax.set_title(f"{gse}: hub gene × immune fraction (Spearman; * FDR<0.05)", fontsize=8.5)
fig.colorbar(im, ax=axes, shrink=0.6, label="Spearman rho")
fig.savefig(FIG / "fig_hubgene_immune_heatmap.png", bbox_inches="tight")
plt.close(fig)
log("fig_hubgene_immune_heatmap.png saved")

# 4c. 方法一致性散点
conc = pd.read_csv(OUT / "method_concordance_epic_vs_mcp.csv")
show = conc[conc["epic"].isin(["EPIC_Bcells", "EPIC_CD8_Tcells", "EPIC_Monocytes",
                               "EPIC_Neutrophils", "EPIC_NKcells"])]
fig, axes = plt.subplots(2, 5, figsize=(13, 5.4))
for ax, (_, r) in zip(axes.T.reshape(-1), show.iterrows()):
    gse = r["gse"]
    frac = frac_all[gse]
    x = frac[r["epic"]].astype(float)
    y = frac[r["mcp"]].astype(float)
    ax.scatter(x, y, s=14, alpha=0.75, color="#4C72B0", edgecolor="none")
    ax.set_xlabel(r["epic"].replace("EPIC_", "EPIC "), fontsize=7)
    ax.set_ylabel(r["mcp"].replace("MCP_", "MCP "), fontsize=7)
    ax.set_title(f"{gse}\nrho={r['rho']:.2f}, p={r['p']:.1e}", fontsize=7.5)
    ax.tick_params(labelsize=6.5)
fig.suptitle("Method concordance: EPIC vs MCPcounter (per-sample fractions)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(FIG / "fig_method_concordance.png", bbox_inches="tight")
plt.close(fig)
log("fig_method_concordance.png saved")

# 4d. KpyRN 模块 × 免疫相关热图
fig, axes = plt.subplots(1, 2, figsize=(15, 3.2))
for ax, gse in zip(axes, ["GSE134431", "GSE199939"]):
    sub = ms_corr[ms_corr["gse"] == gse]
    piv = sub.pivot(index="module", columns="metric", values="rho")
    fdr = sub.pivot(index="module", columns="metric", values="p_fdr")
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels([m.replace("EPIC_", "E:").replace("MCP_", "M:").replace("EST_", "S:")
                        for m in piv.columns], rotation=90, fontsize=6.5)
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels(piv.index, fontsize=7)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            if fdr.values[i, j] < 0.05:
                ax.text(j, i, "*", ha="center", va="center", fontsize=6, color="black")
    ax.set_title(f"{gse}: KpyRN module score × immune fraction", fontsize=8.5)
fig.colorbar(im, ax=axes, shrink=0.6, label="Spearman rho")
fig.savefig(FIG / "fig_kpyrn_module_immune.png", bbox_inches="tight")
plt.close(fig)
log("fig_kpyrn_module_immune.png saved")

log("DONE — immune_figures.py")
