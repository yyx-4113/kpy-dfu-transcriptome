# -*- coding: utf-8 -*-
"""
immune_estimate_correlate.py — 免疫浸润第三方法（ESTIMATE）+ 组间比较 + 枢纽基因-免疫相关
项目：丙酮酰化（Kpy）与糖尿病足（DFU）

方法链（任务 #7，P6 阶段验收：≥2 种独立方法给出一致方向；枢纽基因与免疫细胞
Spearman 相关 + FDR 校正；相关性热图）：
  [R] EPIC 1.1.7（GfellerLab/EPIC，BRef，6 群 + otherCells）→ 03_结果/免疫浸润/<GSE>_epic_fractions.csv
  [R] MCPcounter（ebechet 官方算法直实现，10 群）→ 03_结果/免疫浸润/<GSE>_mcpcounter_fractions.csv
  [Python 本脚本] ESTIMATE（复刻 omnideconv/immunedeconv R/estimate.R 算法，签名取自其 extdata）
      + 组间差异（Wilcoxon + BH FDR）
      + 枢纽基因 × 免疫细胞 Spearman 相关 + BH FDR
      + KpyRN 模块评分 × 免疫细胞相关
      + EPIC vs MCPcounter 方法一致性

ESTIMATE 算法（忠实复刻 immunodeconv 版本）：
  1. 每样本内按表达秩次（ties.method=average），m = 10000 * rank / Ng
  2. 对每个签名：按 m 降序排列，RES = cumsum(TAG*w/sum_w) - cumsum((1-TAG)/Nm)，w=|m|^0.25
  3. ES = sum(RES)（注意：此为 immunodeconv 变体；原 Yoshihara 2013 用峰值 max/min(RES)）
  4. ESTIMATE score = stromal ES + immune ES
  注：肿瘤纯度公式仅适用于 Affymetrix 芯片原发肿瘤（Yoshihara 2013），DFU 非恶性组织不适用，
      故本分析不报告 TumorPurity。

输出：03_结果/免疫浸润/*.csv、03_结果/免疫浸润/figures/*.png
"""
import json

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import sys
from pathlib import Path

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
SIGDIR = BASE / "01_数据/基因集/免疫"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
})


def log(msg):
    print(msg, flush=True)
    with open(OUT / "immune_python.log", "a", encoding="utf-8") as fh:
        fh.write(str(msg) + "\n")


# ---------------- ESTIMATE（复刻 immunodeconv 算法） ----------------
def estimate_scores(expr: pd.DataFrame, sig_csv: Path) -> pd.DataFrame:
    """expr: genes x samples (log2 scale). 返回 sample x [StromalScore, ImmuneScore, ESTIMATEScore]"""
    sig = pd.read_csv(sig_csv)
    signatures = {name: sorted(g["gene"].tolist())
                  for name, g in sig.groupby("signature")}
    common = [g for g in expr.index if g in set(sum(signatures.values(), []))]
    log(f"[ESTIMATE] signature genes available: {len(common)} / "
        f"{len(set(sum(signatures.values(), [])))}")
    m = expr.loc[common].copy()
    Ng = len(m)
    # 每样本内秩次 -> 10000*rank/Ng
    m = m.rank(axis=0, method="average") * 10000.0 / Ng
    gene_names = list(m.index)
    sample_names = list(m.columns)
    Ns = len(sample_names)
    score = {}
    for gs_name, gene_set in signatures.items():
        overlap = [g for g in gene_set if g in gene_names]
        gene_set2 = [gene_names.index(g) for g in overlap]
        Nh, Nm = len(gene_set2), Ng - len(gene_set2)
        ES = np.zeros(Ns)
        for si, s in enumerate(sample_names):
            order = np.argsort(-m[s].values)          # 降序
            w = np.abs(m[s].values[order]) ** 0.25    # 权重
            tag = np.isin(order, gene_set2).astype(float)
            sum_w = w[tag == 1].sum()
            F0 = np.cumsum((1 - tag) / Nm)
            Fn = np.cumsum(tag * w / sum_w)
            RES = Fn - F0
            ES[si] = RES.sum()                        # immunodeconv 变体：sum(RES)
        score[gs_name] = ES
    out = pd.DataFrame(score, index=sample_names)
    out["ESTIMATEScore"] = out["StromalSignature"] + out["ImmuneSignature"]
    out.index.name = "sample"
    return out


# ---------------- 统计工具 ----------------
def wilcoxon_table(df: pd.DataFrame, group_col: pd.Series, g1: str, g2: str,
                   label: str) -> pd.DataFrame:
    """df: index=sample, columns=指标; group_col: sample->group"""
    rows = []
    s1 = df.index[group_col == g1]
    s2 = df.index[group_col == g2]
    for col in df.columns:
        a = df.loc[s1, col].dropna().values
        b = df.loc[s2, col].dropna().values
        if len(a) < 3 or len(b) < 3:
            p = np.nan
        else:
            try:
                p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
            except ValueError:
                p = np.nan
        rows.append({"comparison": label, "metric": col,
                     "mean_" + g1: a.mean(), "mean_" + g2: b.mean(),
                     "n_" + g1: len(a), "n_" + g2: len(b), "p_raw": p})
    res = pd.DataFrame(rows)
    ok = res["p_raw"].notna()
    res.loc[ok, "p_fdr"] = multipletests(res.loc[ok, "p_raw"], method="fdr_bh")[1]
    return res


def spearman_fdr(expr: pd.DataFrame, frac: pd.DataFrame,
                 genes: list) -> pd.DataFrame:
    """expr/frac: index=sample. 返回 gene x metric 的 rho 与 FDR（长表）"""
    samples = [s for s in expr.index if s in frac.index]
    rows = []
    for g in genes:
        if g not in expr.columns:
            continue
        x = expr.loc[samples, g].astype(float).values
        for met in frac.columns:
            y = frac.loc[samples, met].astype(float).values
            if np.std(y) == 0:
                continue
            rho, p = stats.spearmanr(x, y)
            rows.append({"gene": g, "metric": met, "rho": rho, "p_raw": p})
    res = pd.DataFrame(rows)
    if len(res):
        res["p_fdr"] = multipletests(res["p_raw"], method="fdr_bh")[1]
    return res


# ---------------- 主流程 ----------------
log("=" * 70)
log("immune_estimate_correlate.py — 启动")
log("=" * 70)

KPIRN_HUB = ["BPGM", "LDHA", "PFKP", "PGAM1", "PGM2", "TPI1"]      # WGCNA ∩ meta-DEG 糖酵解核心
ML_CONSENSUS = ["ADH1B", "ALDH3A2", "DHX58", "FAHD2B", "GBP1"]     # ML 共识（LDHA 已在核心）
M1_MACHINERY = ["HAT1", "EP300", "SIRT3", "STAT1", "PKM2"]          # Kpy 写入/擦除机器 + 底物/丙酮酸激酶
HUB_GENES = sorted(set(KPIRN_HUB + ML_CONSENSUS + M1_MACHINERY))
log(f"枢纽基因集合（{len(HUB_GENES)}）: {', '.join(HUB_GENES)}")

cohort_info = {
    "GSE134431": {"ann": PROC / "GSE134431_annotation.csv",
                  "group_col": "group2",
                  "contrasts": [("Healer", "Nonhealer", "healing"), ("DFS", "DFU", "group2")]},
    "GSE199939": {"ann": PROC / "GSE199939_annotation.csv",
                  "group_col": "group2",
                  "contrasts": [("DFU", "control", "group2")]},
}

all_est, all_frac = {}, {}
for gse, info in cohort_info.items():
    ann = pd.read_csv(info["ann"])
    ann = ann.set_index("sample")
    expr = pd.read_csv(PROC / f"{gse}_log2_expr.csv", index_col=0)
    epic = pd.read_csv(OUT / f"{gse}_epic_fractions.csv", index_col=0)
    mcp = pd.read_csv(OUT / f"{gse}_mcpcounter_fractions.csv", index_col=0)
    # EPIC CSV 由 R 端写出为 细胞类型×样本（行=细胞类型）；若索引不是样本名则转置
    if not set(epic.index) >= set(ann.index[:1]) or epic.shape[0] < epic.shape[1]:
        if set(epic.columns) >= set(ann.index):
            epic = epic.T
    log(f"[{gse}] EPIC fractions: {epic.shape}; MCPcounter: {mcp.shape}; expr: {expr.shape}")

    # ---- ESTIMATE ----
    est = estimate_scores(expr, SIGDIR / "estimate_signatures.csv")
    est.to_csv(OUT / f"{gse}_estimate_scores.csv")
    all_est[gse] = est
    log(f"[{gse}] ESTIMATE 完成: {est.shape[0]} samples; "
        f"median Immune={est['ImmuneSignature'].median():.1f}, "
        f"Stromal={est['StromalSignature'].median():.1f}")

    # ---- 合并免疫分数（EPIC + MCPcounter + ESTIMATE）----
    epic = epic.add_prefix("EPIC_")
    mcp = mcp.add_prefix("MCP_")
    est = est.add_prefix("EST_")
    frac = epic.join(mcp, how="outer").join(est, how="outer")
    all_frac[gse] = frac

    # ---- 组间比较 ----
    for g1, g2, gcol in info["contrasts"]:
        grp = ann[gcol].fillna("")
        n1, n2 = (grp == g1).sum(), (grp == g2).sum()
        if n1 < 3 or n2 < 3:
            log(f"[{gse}] {g1} vs {g2} ({gcol}): 样本不足 (n={n1}/{n2})，跳过")
            continue
        tab = wilcoxon_table(frac, grp, g1, g2, f"{gse}:{g1}_vs_{g2}")
        tab.to_csv(OUT / f"{gse}_{g1}_vs_{g2}_immune_stats.csv", index=False)
        sig = tab[tab["p_fdr"] < 0.05]
        log(f"[{gse}] {g1}(n={n1}) vs {g2}(n={n2}): FDR<0.05 的指标: "
            f"{', '.join(sig['metric'].tolist()) if len(sig) else '无'}")

    # ---- 枢纽基因 × 免疫相关 ----
    corr = spearman_fdr(expr.T, frac, HUB_GENES)
    corr.to_csv(OUT / f"{gse}_hubgene_immune_corr.csv", index=False)
    strong = corr[corr["p_fdr"] < 0.05].sort_values("p_fdr")
    log(f"[{gse}] hub gene × immune: FDR<0.05 相关 {len(strong)} 对（共 {len(corr)} 对）")
    for _, r in strong.head(12).iterrows():
        log(f"    {r['gene']:9s} ~ {r['metric']:22s} rho={r['rho']:+.3f} FDR={r['p_fdr']:.2e}")

# ---------------- 方法一致性（EPIC vs MCPcounter 跨样本相关） ----------------
log("---- 方法一致性：EPIC vs MCPcounter（全部样本合并）----")
conc_rows = []
for gse in cohort_info:
    frac = all_frac[gse]
    pairs = [("EPIC_Bcells", "MCP_B lineage"), ("EPIC_CD4_Tcells", "MCP_T cells"),
             ("EPIC_CD8_Tcells", "MCP_CD8 T cells"), ("EPIC_Monocytes", "MCP_Monocytic lineage"),
             ("EPIC_Neutrophils", "MCP_Neutrophils"), ("EPIC_NKcells", "MCP_NK cells")]
    for e_col, m_col in pairs:
        if e_col in frac.columns and m_col in frac.columns:
            sub = frac[[e_col, m_col]].dropna()
            if len(sub) > 5 and sub[m_col].std() > 0:
                rho, p = stats.spearmanr(sub[e_col], sub[m_col])
                conc_rows.append({"gse": gse, "epic": e_col, "mcp": m_col,
                                  "rho": rho, "p": p, "n": len(sub)})
conc = pd.DataFrame(conc_rows)
conc.to_csv(OUT / "method_concordance_epic_vs_mcp.csv", index=False)
log(conc.to_string(index=False))

log("DONE — immune_estimate_correlate.py")
