# -*- coding: utf-8 -*-
"""
scrna_composition_stats.py — GSE165816 单细胞 Stage 3：组成 + pseudobulk 统计
项目：丙酮酰化（Kpy）与糖尿病足（DFU）

输入：Stage 2 产物（cell_scores.csv + kpyrn_counts.npz + _genes.csv）

统计单位（诚信设计）：
  - 所有组间检验的单元 = 供体（样本），不是细胞（避免 pseudoreplication）。
  - 组成：每样本 × 谱系占比 → 组间 Wilcoxon（BH FDR）。
  - pseudobulk：每样本 × 谱系的 KpyRN 基因 counts 求和 → CPM → log2；
    仅纳入该谱系细胞数 ≥100 的样本×谱系（防小计数噪声）。

对比（预注册）：
  - 主对比（足皮肤）：DFU（healer+nonhealer）vs 非 DFU（非糖尿病 + 非溃疡糖尿病）
  - 次级（足皮肤）：4 组 Kruskal-Wallis；healer vs nonhealer（探索性）
  - 次级组织：前臂皮肤、PBMC 同样主对比（探索性）

输出：03_结果/单细胞/ 下组成表、pseudobulk 表、枢纽基因×谱系表、摘要 JSON
"""
import json

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

BASE = Path(_ROOT)
SCDIR = BASE / "01_数据/处理后数据/scRNA"
OUT = BASE / "03_结果/单细胞"
OUT.mkdir(parents=True, exist_ok=True)
MIN_CELLS_PB = 100  # pseudobulk 样本×谱系最小细胞数

t0 = time.time()


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    if ok.sum() > 0:
        q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return q


def wilson_effect(x, y):
    """rank-biserial correlation（Wilcoxon 效应量）"""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return np.nan
    u = mannwhitneyu(x, y, alternative="two-sided").statistic
    return 1 - 2 * u / (n1 * n2)


def group_compare(df, val_col, group_col, g1, g2, unit="gsm"):
    """donor 级两组比较：Wilcoxon + 效应量"""
    a = df.loc[df[group_col] == g1, [unit, val_col]].dropna()
    b = df.loc[df[group_col] == g2, [unit, val_col]].dropna()
    # 每供体去重（组成已是每供体一行；pseudobulk 同理）
    a = a.drop_duplicates(subset=[unit])
    b = b.drop_duplicates(subset=[unit])
    if len(a) < 3 or len(b) < 3:
        return {"n1": len(a), "n2": len(b), "p": np.nan, "q": np.nan,
                "effect_rb": np.nan, "median1": float(a[val_col].median()) if len(a) else np.nan,
                "median2": float(b[val_col].median()) if len(b) else np.nan}
    p = mannwhitneyu(a[val_col], b[val_col], alternative="two-sided").pvalue
    return {"n1": len(a), "n2": len(b), "p": float(p), "q": np.nan,
            "effect_rb": float(wilson_effect(a[val_col], b[val_col])),
            "median1": float(a[val_col].median()), "median2": float(b[val_col].median())}


# ---------- 读取 ----------
df = pd.read_csv(SCDIR / "scRNA_GSE165816_cell_scores.csv")
kp_counts = sp.load_npz(SCDIR / "scRNA_GSE165816_kpyrn_counts.npz").toarray().astype("float32")
kp_genes = pd.read_csv(SCDIR / "scRNA_GSE165816_kpyrn_genes.csv", header=None)[0].tolist()
log(f"cells: {len(df)}; kpyrn counts matrix: {kp_counts.shape}; genes: {len(kp_genes)}")
assert kp_counts.shape[0] == len(df), "row misalignment between scores CSV and counts npz"

kp = pd.read_csv(BASE / "01_数据/基因集/KpyRN_modules.csv")
scored = set(kp.loc[kp["scoring_include"] == "Yes", "gene_symbol"])  # noqa: E712
kp_scored_idx = [j for j, g in enumerate(kp_genes) if g in scored]
mod_of_gene = dict(zip(kp["gene_symbol"], kp["module"]))

df["group"] = np.where(df["disease"].isin(["DFU-healer", "DFU-nonhealer"]), "DFU", "nonDFU")
df["healing"] = df["disease"]

# ================= A. 组成分析（donor 级） =================
comp_rows = []
comp = (df.groupby(["gsm", "tissue", "disease", "group", "lineage_parent"])
          .size().rename("n_cells").reset_index())
tot = comp.groupby("gsm")["n_cells"].transform("sum")
comp["prop"] = comp["n_cells"] / tot

for tissue in ["Foot skin", "Forearm skin", "PBMCs"]:
    ct = comp[comp["tissue"] == tissue]
    n_by_group = ct.groupby("group")["gsm"].nunique().to_dict()
    log(f"[A] {tissue}: {n_by_group}")
    lineages_t = sorted(ct["lineage_parent"].unique())
    all_gsm = ct["gsm"].unique()
    gmap = ct.drop_duplicates("gsm").set_index("gsm")["group"]
    dmap = ct.drop_duplicates("gsm").set_index("gsm")["disease"]
    res = []
    for ln in lineages_t:
        # 每供体一行：该谱系占比（无该谱系 → 0），组别来自供体映射
        prop = (ct[ct["lineage_parent"] == ln].set_index("gsm")["prop"]
                  .reindex(all_gsm).fillna(0.0))
        d = pd.DataFrame({"prop": prop, "group": gmap})
        r = group_compare(d.reset_index(), "prop", "group", "DFU", "nonDFU")
        r.update({"tissue": tissue, "lineage": ln,
                  "mean_prop_DFU": float(d.loc[d["group"] == "DFU", "prop"].mean()),
                  "mean_prop_nonDFU": float(d.loc[d["group"] == "nonDFU", "prop"].mean())})
        res.append(r)
    res = pd.DataFrame(res)
    res["q"] = bh(res["p"].values)
    res = res.sort_values("q")
    res.to_csv(OUT / f"composition_{tissue.replace(' ', '_')}_DFU_vs_nonDFU.csv", index=False)
    sig = res[res["q"] < 0.05]
    log(f"[A] {tissue}: {len(sig)} lineages FDR<0.05: "
        f"{sig[['lineage','mean_prop_DFU','mean_prop_nonDFU','q']].to_dict('records') if len(sig) else 'none'}")

# 4 组 Kruskal-Wallis（足皮肤）
foot = comp[comp["tissue"] == "Foot skin"]
kw_rows = []
for ln in sorted(foot["lineage_parent"].unique()):
    sub = foot[foot["lineage_parent"] == ln]
    all_gsm = foot["gsm"].unique()
    pivot = (sub.pivot_table(index="gsm", columns="disease", values="prop")
                .reindex(all_gsm).fillna(0.0))
    groups = [pivot[c].dropna().values for c in pivot.columns]
    groups = [g for g in groups if len(g) >= 3]
    if len(groups) >= 3:
        p = kruskal(*groups).pvalue
        kw_rows.append({"lineage": ln, "p": float(p),
                        "n_groups": len(groups),
                        "median_prop_by_group": {c: float(pivot[c].median()) for c in pivot.columns}})
kw = pd.DataFrame(kw_rows)
kw["q"] = bh(kw["p"].values)
kw.sort_values("q").to_csv(OUT / "composition_Foot_skin_4group_kruskal.csv", index=False)
log(f"[A] Foot skin 4-group KW: {int((kw['q'] < 0.05).sum())} lineages FDR<0.05")

# healer vs nonhealer（足皮肤，探索性）
foot_gsm = foot["gsm"].unique()
hh = (foot.pivot_table(index="gsm", columns="lineage_parent", values="prop",
                       aggfunc="sum")
        .reindex(columns=sorted(foot["lineage_parent"].unique())).fillna(0.0))
hh["disease"] = hh.index.map(foot.drop_duplicates("gsm").set_index("gsm")["disease"])
hh_rows = []
for ln in [c for c in hh.columns if c != "disease"]:
    r = group_compare(hh.reset_index(), ln, "disease", "DFU-healer", "DFU-nonhealer")
    r["lineage"] = ln
    hh_rows.append(r)
hh_df = pd.DataFrame(hh_rows)[["lineage", "n1", "n2", "median1", "median2", "p", "effect_rb"]]
hh_df["q"] = bh(hh_df["p"].values)
hh_df.sort_values("q").to_csv(OUT / "composition_Foot_skin_healer_vs_nonhealer.csv", index=False)

# ================= B. Pseudobulk KpyRN（donor × lineage） =================
log("[B] pseudobulk aggregation...")
df["row"] = np.arange(len(df))
counts_sum = pd.DataFrame(kp_counts, columns=kp_genes)
counts_sum["gsm"] = df["gsm"].values
counts_sum["lineage_parent"] = df["lineage_parent"].values
counts_sum["n_counts"] = df["n_counts"].values
pb = counts_sum.groupby(["gsm", "lineage_parent"]).sum(numeric_only=True)
cell_n = counts_sum.groupby(["gsm", "lineage_parent"]).size().rename("n_cells")
tot_counts = counts_sum.groupby(["gsm", "lineage_parent"])["n_counts"].sum().rename("tot_counts")
pb = pb.join(cell_n).join(tot_counts).reset_index()
pb = pb[pb["n_cells"] >= MIN_CELLS_PB].copy()
meta_map = df.drop_duplicates("gsm").set_index("gsm")[["tissue", "disease", "group"]]
pb = pb.join(meta_map, on="gsm")
log(f"[B] pseudobulk rows (sample x lineage): {len(pb)}; "
    f"tissues: {pb['tissue'].value_counts().to_dict()}")

# CPM → log2（每行库大小归一）
for g in kp_genes:
    pb[g] = np.log2(pb[g] / pb["tot_counts"] * 1e6 + 1.0)

# 模块均值（仅评分基因）
for mod, gl in kp.groupby("module")["gene_symbol"].apply(list).items():
    use = [g for g in gl if g in pb.columns and g in scored]
    if use:
        pb[f"module__{mod}"] = pb[use].mean(axis=1)

pb.to_csv(OUT / "pseudobulk_donor_by_lineage.csv", index=False)

# 足皮肤 DFU vs nonDFU，逐谱系 × 基因（评分基因）+ 模块
foot_pb = pb[pb["tissue"] == "Foot skin"]
tests = []
test_genes = [g for g in kp_genes if g in scored] + [c for c in pb.columns if c.startswith("module__")]
for ln in sorted(foot_pb["lineage_parent"].unique()):
    sub = foot_pb[foot_pb["lineage_parent"] == ln]
    if sub["group"].nunique() < 2:
        continue
    for g in test_genes:
        r = group_compare(sub, g, "group", "DFU", "nonDFU")
        r.update({"lineage": ln, "feature": g})
        tests.append(r)
tdf = pd.DataFrame(tests)
# BH 逐谱系（谱系内基因相关，按谱系分层校正）
tdf["q"] = np.nan
for ln, idx in tdf.groupby("lineage").groups.items():
    tdf.loc[idx, "q"] = bh(tdf.loc[idx, "p"].values)
tdf = tdf.sort_values(["lineage", "q"])
tdf.to_csv(OUT / "pseudobulk_Foot_skin_DFU_vs_nonDFU_by_lineage.csv", index=False)
sig = tdf[(tdf["q"] < 0.05)]
log(f"[B] foot skin pseudobulk FDR<0.05 tests: {len(sig)}")
if len(sig):
    log(sig[["lineage", "feature", "median1", "median2", "q"]].to_string(index=False))

# ================= C. 枢纽基因 × 谱系（描述性，log1p CP10K） =================
HUB = ["LDHA", "PGAM1", "PGK1", "PFKP", "TPI1", "GAPDH", "ENO1", "ALDOA",
       "PGM2", "BPGM", "LDHB", "SLC16A1", "MPC2", "IFI44L", "IFI27", "STAT2",
       "DDX58", "DHX58", "ISG20", "GBP1", "HAT1", "EP300", "SIRT3", "STAT1",
       "FGFR4", "LGALS8", "RSPH4A", "ADH1B", "ALDH3A2", "FAHD2B"]
hub_cols = [g for g in HUB if g in df.columns]
hub = (df[df["tissue"] == "Foot skin"]
       .groupby("lineage_parent")[hub_cols].median().round(3))
hub.to_csv(OUT / "hubgene_median_by_lineage_footskin.csv")

# ================= 摘要 =================
n_by_grp = df[df["tissue"] == "Foot skin"].groupby(["group", "disease"])["gsm"].nunique()
n_by_grp = {f"{k[0]}|{k[1]}": int(v) for k, v in n_by_grp.items()}
summary = {
    "n_cells": int(len(df)),
    "foot_skin_donors_by_group": n_by_grp,
    "composition_files": ["composition_Foot_skin_DFU_vs_nonDFU.csv",
                          "composition_Forearm_skin_DFU_vs_nonDFU.csv",
                          "composition_PBMCs_DFU_vs_nonDFU.csv",
                          "composition_Foot_skin_4group_kruskal.csv",
                          "composition_Foot_skin_healer_vs_nonhealer.csv"],
    "pseudobulk_min_cells": MIN_CELLS_PB,
    "pseudobulk_rows": int(len(pb)),
    "pseudobulk_sig_tests": int(len(sig)),
    "runtime_sec": round(time.time() - t0, 1),
}
with open(OUT / "scrna_stage3_summary.json", "w", encoding="utf-8") as fh:
    json.dump(summary, fh, ensure_ascii=False, indent=2)
log(f"DONE stage3 ({time.time() - t0:.0f}s)")
