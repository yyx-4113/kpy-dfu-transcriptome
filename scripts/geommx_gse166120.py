# -*- coding: utf-8 -*-
"""
geommx_gse166120.py — GSE166120 GeoMx DSP 空间/区域验证（任务 #11c）

项目：丙酮酰化（Kpy）与糖尿病足（DFU）

数据核实（GEO series matrix + 原文 Theocharidis Nat Commun 2022）：
  - 技术：NanoString GeoMx **Digital Spatial Profiling (DSP)**，Cancer Transcriptome
    Atlas (CTA) **v2.0 靶向板，1825 target**（非计划书原写的 Visium —— 已纠正）
  - 组织：糖尿病足溃疡（DFU）皮肤，**Healers N=2 vs Non-Healers N=2**（12 周术后结局）
  - 样本：23 个 ROI（单个 GSM = 1 个 ROI），4 位受试者
    Subject 1/2 = Non-Healer（A02-A12, B01-B04），Subject 3/4 = Healer
  - ROI 位置分层：Ulcer / Ulcer edge / Non-ulcer / Epidermis

关键约束（决定本节能做什么、不能做什么，必须在稿件中如实声明）：
  1. **CTA 靶向板只覆盖 KpyRN 评分基因的 42.2%（54/128）**；M4 IFN-I 覆盖 83.9%
     （26/31）、M2 丙酮酸供应 27.8%（25/90）、M1 50%（2/4）、M3 33.3%（1/3）
  2. **6 个 ML 共识枢纽基因中仅 3 个在板内**：LDHA、ADH1B、GBP1
     （FAHD2B / DHX58 / ALDH3A2 确实不在 CTA 板 —— 已按 SystematicName 逐 probe 核实）
  3. **n=2 vs 2 受试者**（23 ROI 嵌套于 4 人）→ ROI 级检验存在 pseudoreplication，
     受试体级 n=4 功效极低；**本节为探索性空间定位证据，不作独立因果推断**

数据处理（DSP 官方流程的合理简化，如实记录）：
  - DCC 文件（GeoMx 原始格式，纯文本分层）<Code_Summary> 段给出 **每个 probe 的读段计数**
  - PKC 文件（Alpha_CTA_v2.0.pkc）给出 probe→target→gene 映射：
    **按 probe 级 SystematicName 映射**（不可用 DisplayName —— 存在 ADH1A target
    混编 ADH1A/B/C probe 的跨基因情况；也不可只取 SystematicName 首元素 —— 会索引错位漏基因）
  - 同基因多 probe 计数求和 → target 计数；DSQ2 上游 QC 不可得（无 .qcreport / .csv），
    故不做 GC/面积/核数归一化，改为**每个 ROI 内的相对丰度**（构成比）用于跨 ROI 比较
  - 模块评分 = 板内可用基因的 z-score 均值（**只代表板内子集**，不可与全 KpyRN 分数直接比较）

输出：03_结果/空间GSE166120/
  dcc_probe_counts.csv       23 ROI × 1825 target 计数矩阵
  roi_metadata.csv           ROI 元数据（受试者/结局/位置）
  target_coverage.csv        KpyRN/枢纽基因在 CTA 板的覆盖清单
  module_scores_by_roi.csv   板内子集模块评分（按 ROI）
  hubgene_expression.csv     板内 3 个枢纽基因的表达
  region_stats.csv           ROI 位置 / 结局 分层的探索性统计
  summary.json               关键数字汇总（供稿件引用）
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import re
import gzip
import json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, kruskal

BASE = _ROOT
RAW = f"{BASE}/01_数据/原始数据/GSE166120"
EXT = f"{RAW}/extracted"
OUT = f"{BASE}/03_结果/空间GSE166120"
os.makedirs(OUT, exist_ok=True)

HUB = ["LDHA", "FAHD2B", "ADH1B", "DHX58", "ALDH3A2", "GBP1"]

# ============ 1. 解析 series matrix 拿 ROI 元数据 ============
# 注意：行首标签是「字段名」，而列值随行变化；!Sample_title 行给的是孔位名(A02…)，
# GSM ID 在 !Sample_geo_accession 行 —— 必须分别取，不能混用。
# 关键：!Sample_characteristics_ch1 在文件中出现多次（每个 characteristic 一行，
# 行标签完全相同）→ 不能用 dict 按键存储（会互相覆盖丢字段），必须用 list 保留全部行。
rows_list = []
with gzip.open(f"{RAW}/GSE166120_series_matrix.txt.gz", "rt",
               encoding="utf-8", errors="ignore") as fh:
    for line in fh:
        if line.startswith("!Sample"):
            parts = line.rstrip("\n").split("\t")
            rows_list.append((parts[0], [v.strip('"') for v in parts[1:]]))
n_col = len(rows_list[0][1])
assert all(len(v) == n_col for _, v in rows_list), "元数据行列数不一致"

by_key = {}
for key, vals in rows_list:
    if key == "!Sample_characteristics_ch1":
        by_key.setdefault(key, []).append(vals)      # 保留全部 characteristics 行
    else:
        by_key[key] = vals

gsm = by_key["!Sample_geo_accession"]
titles = by_key["!Sample_title"]
fields = {}
for key, vals in by_key.items():
    if key in ("!Sample_title", "!Sample_geo_accession"):
        continue
    if key == "!Sample_characteristics_ch1":
        for row in vals:                              # 逐个 characteristic 行
            for i, v in enumerate(row):
                if ":" in v:
                    f, val = v.split(":", 1)
                    fields.setdefault(f.strip(), {})[i] = val.strip()
        continue
    for i, v in enumerate(vals):
        if ":" in v:
            f, val = v.split(":", 1)
            fields.setdefault(f.strip(), {})[i] = val.strip()

meta = pd.DataFrame({"gsm": gsm, "title": titles})
for f, d in fields.items():
    meta[f] = [d.get(i, "") for i in range(n_col)]
meta = meta.rename(columns={
    "subject id": "subject", "healing status": "healing",
    "roi location": "roi_location", "ulcer location": "ulcer_location",
    "tissue": "tissue", "Sex": "sex"})
# ROI 位置字段存在 'Non-ulcer' / 'Non-Ulcer' 两种大小写写法（A05 vs A07），
# 不统一会导致同一解剖分层被拆成两组 → 归一化
meta["roi_location"] = meta["roi_location"].str.strip().str.lower()
meta["healing"] = meta["healing"].str.strip()
meta["subject"] = meta["subject"].str.strip()
print("ROI 元数据:", meta.shape)
print(meta[["gsm", "title", "subject", "healing", "roi_location", "sex"]].to_string(index=False))
assert meta["gsm"].is_unique and set(meta["gsm"]) == set(gsm), "GSM 解析异常"

# ============ 2. 解析 PKC：probe 级 SystematicName 映射 ============
with gzip.open(f"{RAW}/GSE166120_Alpha_CTA_v2.0.pkc.txt.gz", "rt",
               encoding="utf-8", errors="ignore") as fh:
    pkc = json.load(fh)
probe_rows = []
for t in pkc["Targets"]:
    for pr in t["Probes"]:
        rts = pr.get("RTS_ID")
        sysl = pr.get("SystematicName") or []
        # 关键：逐 probe 用其自身的 SystematicName（可能多基因，取第一个非空）
        gene = next((s for s in sysl if s), None)
        probe_rows.append({"rts_id": rts, "display": t["DisplayName"],
                           "gene": gene, "n_sys": len(sysl)})
pm = pd.DataFrame(probe_rows)
print(f"\nPKC probe 数: {len(pm)} | 有 SystematicName 的: {pm.gene.notna().sum()}")
# 跨基因 target 检查（审计）
multi = pm.groupby("display")["gene"].nunique()
print("跨多基因的 target（DisplayName 混编）:", dict(multi[multi > 1]))
# 负对照探针（ERCC/NegProbe，49 个）必须排除在生物学分析之外
n_neg = int((pm.display == "NegProbe").sum())
pm_bio = pm[pm.display != "NegProbe"].copy()
print(f"排除 NegProbe 负对照探针: {n_neg} 个 → 剩余生物学探针 {len(pm_bio)}")

# ============ 3. 解析 23 个 DCC → 计数矩阵 ============
dcc_files = sorted(f for f in os.listdir(EXT) if f.endswith(".dcc.gz"))
print(f"\nDCC 文件数: {len(dcc_files)}")
probe_counts = {}
for fn in dcc_files:
    gsm_id = fn.split("_")[0]
    with gzip.open(f"{EXT}/{fn}", "rt", encoding="utf-8", errors="ignore") as fh:
        lines = fh.read().split("\n")
    i0 = lines.index("<Code_Summary>")
    i1 = lines.index("</Code_Summary>")
    d = {}
    for l in lines[i0 + 1:i1]:
        l = l.strip()
        if not l:
            continue
        parts = l.split(",")
        if len(parts) != 2:
            continue
        try:
            d[parts[0]] = int(parts[1])
        except ValueError:
            pass
    probe_counts[gsm_id] = d
    # 平台按 ROI 上报可用探针数不同（8668–8671，含 ERCC/NegProbe 差异），非解析错误
    assert 8660 <= len(d) <= 8680, f"{fn}: 解析到 {len(d)} 个 code（异常）"

all_codes = set().union(*[set(d.keys()) for d in probe_counts.values()])
print(f"probe code 并集: {len(all_codes)} | 映射表 probe: {pm.rts_id.nunique()}")
assert all_codes <= set(pm.rts_id), "存在未映射 code"

pc = pd.DataFrame(probe_counts).T  # rows = ROI(gsm), cols = rts_id
pc.index.name = "gsm"
print("probe 计数矩阵:", pc.shape)

# 仅保留生物学探针（剔除 NegProbe 负对照）
bio_codes = set(pm_bio.rts_id)
pc = pc[[c for c in pc.columns if c in bio_codes]]
print("剔除 NegProbe 后 probe 矩阵:", pc.shape)

# probe → gene 求和 → target 计数矩阵
pc2 = pc.copy()
pc2.columns = pm_bio.set_index("rts_id").loc[pc2.columns, "gene"].values
pc2 = pc2.T.groupby(level=0).sum().T  # 同基因多 probe 求和
pc2 = pc2.loc[:, pc2.columns.notna()]
print("target 计数矩阵:", pc2.shape)
pc2.to_csv(f"{OUT}/dcc_target_counts.csv")

meta = meta.set_index("gsm").loc[pc2.index].reset_index()
meta.to_csv(f"{OUT}/roi_metadata.csv", index=False)

# ============ 4. KpyRN 覆盖清单 ============
K = pd.read_csv(f"{BASE}/01_数据/基因集/KpyRN_modules.csv")
Ksc = K[K["scoring_include"] == "Yes"].copy()
board = set(pc2.columns)
cov_rows = []
for mod, grp in Ksc.groupby("module"):
    inb = [g for g in grp.gene_symbol if g in board]
    cov_rows.append({"module": mod, "n_total": len(grp), "n_in_board": len(inb),
                     "pct": round(100 * len(inb) / len(grp), 1),
                     "genes_in_board": ";".join(sorted(inb))})
cov = pd.DataFrame(cov_rows)
hub_cov = pd.DataFrame([{"gene": g, "in_board": g in board} for g in HUB])
print("\n== CTA 板覆盖 ==")
print(cov.to_string(index=False))
print(f"\n枢纽基因在板内: {hub_cov[hub_cov.in_board].gene.tolist()}")
print(f"枢纽基因不在板内: {hub_cov[~hub_cov.in_board].gene.tolist()}")
cov.to_csv(f"{OUT}/target_coverage_modules.csv", index=False)
hub_cov.to_csv(f"{OUT}/target_coverage_hubgenes.csv", index=False)

# ============ 5. 相对丰度 + 板内子集模块评分 ============
# 组成比归一（每个 ROI 内占比），用于跨 ROI 比较丰度
rel = pc2.div(pc2.sum(axis=1), axis=0)

mod_score_rows = []
for _, r in cov.iterrows():
    genes = [g for g in r.genes_in_board.split(";") if g]
    if len(genes) < 2:
        continue
    sub = rel[genes]
    # z-score 方向：逐基因跨 ROI 标准化（该基因在此 ROI 相对其他 ROI 的富集），
    # 再跨基因取均值 = 模块评分。
    # 注意：若做「逐 ROI 内标准化」再跨基因取均值，结果恒为 0（z 的行均值必为 0），
    # 这是一个真实的统计陷阱 —— 2026-10-25 修正。
    gene_mean = sub.mean(axis=0)                 # 每基因在 23 个 ROI 的均值
    gene_std = sub.std(axis=0, ddof=1).replace(0, np.nan)
    z = sub.sub(gene_mean, axis=1).div(gene_std, axis=1)
    mod_score_rows.append(pd.DataFrame({
        "gsm": pc2.index, "module": r.module, "n_genes_board": len(genes),
        "mean_rel_abund": sub.mean(axis=1).values,
        "zscore": z.mean(axis=1).values,
        "pos_frac": (sub > 0).mean(axis=1).values}))
ms = pd.concat(mod_score_rows, ignore_index=True)
ms = ms.merge(meta[["gsm", "subject", "healing", "roi_location"]], on="gsm")
ms.to_csv(f"{OUT}/module_scores_by_roi.csv", index=False)

# ============ 6. 板内枢纽基因表达 ============
hub_rows = []
for g in HUB:
    if g not in board:
        continue
    hub_rows.append(pd.DataFrame({
        "gsm": pc2.index, "gene": g,
        "raw_count": pc2[g].values,
        "rel_abund": rel[g].values,
        "log2_rel": np.log2(rel[g].values + 1e-6)}))
hub_df = pd.concat(hub_rows, ignore_index=True).merge(
    meta[["gsm", "subject", "healing", "roi_location"]], on="gsm")
hub_df.to_csv(f"{OUT}/hubgene_expression.csv", index=False)

# ============ 7. 探索性分层统计 ============
# 关键：ROI 嵌套于 4 位受试体（23 ROI / 4 人）→ ROI 级检验存在 pseudoreplication。
# 因此**主口径 = 受试体级**（先在受试体内对 ROI 取均值，n=4），ROI 级仅作敏感性分析。
stat_rows = []
for _, r in cov.iterrows():
    genes = [g for g in r.genes_in_board.split(";") if g]
    if len(genes) < 2:
        continue
    s = ms[ms.module == r.module].set_index("gsm")
    for grouping in ("roi_location", "healing"):
        levels = sorted(s[grouping].unique())
        if len(levels) < 2:
            continue
        # --- 受试体级（主口径）---
        sub_lvl = s.groupby(["subject", grouping], observed=True)["zscore"].mean().reset_index()
        gsub = [sub_lvl.loc[sub_lvl[grouping] == lv, "zscore"].values for lv in levels]
        gsub = [g for g in gsub if len(g) > 0]
        if len(gsub) >= 2 and min(len(g) for g in gsub) >= 2:
            try:
                if len(gsub) == 2:
                    p_sub = mannwhitneyu(gsub[0], gsub[1], alternative="two-sided",
                                        method="exact").pvalue
                    t_sub = "MWU exact"
                else:
                    p_sub = kruskal(*gsub).pvalue
                    t_sub = "Kruskal"
            except ValueError:
                p_sub, t_sub = np.nan, "NA"
        else:
            p_sub, t_sub = np.nan, "insufficient (n<2 per group)"
        # --- ROI 级（敏感性分析，含 pseudoreplication 风险）---
        groups = [s.loc[s[grouping] == lv, "zscore"].values for lv in levels]
        groups = [g for g in groups if len(g) > 0]
        if len(groups) >= 2:
            try:
                if len(groups) == 2:
                    p_roi = mannwhitneyu(groups[0], groups[1], alternative="two-sided",
                                         method="exact").pvalue
                    t_roi = "MWU exact"
                else:
                    p_roi = kruskal(*groups).pvalue
                    t_roi = "Kruskal"
            except ValueError:
                p_roi, t_roi = np.nan, "NA"
        else:
            p_roi, t_roi = np.nan, "NA"
        stat_rows.append({
            "module": r.module, "grouping": grouping, "levels": "|".join(levels),
            "subject_means": "|".join(f"{g.mean():.4f}" for g in gsub),
            "n_subjects": "|".join(str(len(g)) for g in gsub),
            "test_subject": t_sub, "p_subject": p_sub,
            "roi_means": "|".join(f"{g.mean():.4f}" for g in groups),
            "n_roi": "|".join(str(len(g)) for g in groups),
            "test_roi": t_roi, "p_roi": p_roi,
            "note": ("primary = subject-level (pseudoreplication-free); "
                     "ROI-level = sensitivity only (ROIs nested in 4 subjects)")})
st = pd.DataFrame(stat_rows)
if len(st):
    from statsmodels.stats.multitest import multipletests
    st["p_subject_bf"] = multipletests(st["p_subject"], method="fdr_bh")[1]
    st["p_roi_bf"] = multipletests(st["p_roi"], method="fdr_bh")[1]
st.to_csv(f"{OUT}/region_stats.csv", index=False)

# ============ 8. 汇总 ============
summary = {
    "technology": "GeoMx DSP, Alpha Cancer Transcriptome Atlas v2.0 (targeted, 1825 targets)",
    "n_roi": int(pc2.shape[0]),
    "n_subjects": int(meta.subject.nunique()),
    "healing_groups": {k: int(v) for k, v in
                       meta.groupby("healing").subject.nunique().items()},
    "roi_locations": {k: int(v) for k, v in meta.roi_location.value_counts().items()},
    "targets_parsed": int(pc2.shape[1]),
    "probes_per_roi": 8670,
    "kpyrn_scored_genes_total": int(len(Ksc)),
    "kpyrn_scored_genes_in_board": int(sum(cov.n_in_board)),
    "kpyrn_coverage_pct": round(100 * sum(cov.n_in_board) / len(Ksc), 1),
    "module_coverage": {r.module: {"n_in_board": int(r.n_in_board), "n_total": int(r.n_total),
                                   "pct": r.pct} for _, r in cov.iterrows()},
    "hub_genes_in_board": hub_cov[hub_cov.in_board].gene.tolist(),
    "hub_genes_not_in_board": hub_cov[~hub_cov.in_board].gene.tolist(),
    "n_significant_subject_level_fdr10": int((st.p_subject_bf < 0.10).sum()) if len(st) else 0,
    "n_significant_roi_level_fdr10_sensitivity": int((st.p_roi_bf < 0.10).sum()) if len(st) else 0,
    "note": ("exploratory only: n=2 vs 2 subjects; primary inference is SUBJECT-level "
             "(n=4, ROI nested within subject -> pseudoreplication-free); ROI-level tests "
             "reported as sensitivity only; module scores use ONLY board-covered genes "
             "and are NOT comparable to full-KpyRN scores"),
}
with open(f"{OUT}/summary.json", "w", encoding="utf-8") as fh:
    json.dump(summary, fh, ensure_ascii=False, indent=1, default=str)

print("\n== 区域/结局分层统计 ==")
if len(st):
    print(st[["module","grouping","levels","subject_means","n_subjects","test_subject","p_subject","p_subject_bf","roi_means","n_roi","p_roi","p_roi_bf"]]
          .round(4).to_string(index=False))
print("\n== summary ==")
print(json.dumps(summary, ensure_ascii=False, indent=1))
print("\nDONE — geommx_gse166120.py")
