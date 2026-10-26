# -*- coding: utf-8 -*-
"""
scrna_annotate_score.py — GSE165816 单细胞 Stage 2：谱系注释 + KpyRN 评分
项目：丙酮酰化（Kpy）与糖尿病足（DFU）

输入：Stage 1 落盘的逐样本 npz（cells × genes, CSR float32, 已 QC 过滤）
      + scRNA_GSE165816_cell_qc.csv（逐细胞 QC 元数据）
      + 01_数据/基因集/免疫/scRNA_lineage_markers.csv（19 谱系标记面板）
      + 01_数据/基因集/KpyRN_modules.csv（140 基因，128 纳入评分）

方法（如实记录，非参考基因组注释）：
  - 标记评分法（marker-scoring lineage assignment）：CP10K + log1p 标度后，
    谱系得分 = 该谱系全部标记基因的均值（缺失标记按 0 计，含 core/support 两层），
    逐细胞取 argmax；top1-top2 差值 < 0.25 判 "Ambiguous"（预注册阈值）。
  - 本方法为轻量启发式注释，非 SingleR/Azimuth 参考比对；
    以组织学 sanity check 验证（PBMC 血液谱系占比、皮肤实质谱系占比）。
  - KpyRN 模块评分：同法，模块内基因均值（M5 组蛋白按基因集定义不纳入评分）。
  - 枢纽基因面板：log1p(CP10K) 逐细胞表达值。
  - KpyRN 原始计数块另存 npz（供 Stage 3 pseudobulk；全 140 基因固定列序，
    样本缺失基因填 0，行序与 CSV 对齐）。

输出：
  01_数据/处理后数据/scRNA/scRNA_GSE165816_cell_scores.csv
  01_数据/处理后数据/scRNA/scRNA_GSE165816_kpyrn_counts.npz（+ _genes.csv）
  03_结果/QC/scRNA_GSE165816_stage2_summary.json
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

BASE = Path(_ROOT)
SCDIR = BASE / "01_数据/处理后数据/scRNA"
GSE165816 = SCDIR / "GSE165816"
QCDIR = BASE / "03_结果/QC"
MARKERS = BASE / "01_数据/基因集/免疫/scRNA_lineage_markers.csv"
KPYRN = BASE / "01_数据/基因集/KpyRN_modules.csv"

MARGIN_MIN = 0.25  # 预注册：top1-top2 < 0.25 → Ambiguous

t0 = time.time()


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------- 面板 ----------
mk = pd.read_csv(MARKERS)
lineages = sorted(mk["lineage"].unique())
lineage_genes = {ln: mk.loc[mk["lineage"] == ln, "gene_symbol"].tolist()
                 for ln in lineages}
# 亚谱系 → 母谱系（组成分析用）
LINEAGE_PARENT = {
    "Fibroblast_papillary": "Fibroblast", "Fibroblast_reticular": "Fibroblast",
    "Keratinocyte_basal": "Keratinocyte", "Keratinocyte_suprabasal": "Keratinocyte",
    "Endothelial_lymphatic": "Endothelial",
}
BLOOD_LINEAGES = {"T_CD4", "T_CD8", "NK", "B", "Plasma", "Mono_CD14",
                  "Mono_CD16", "Macrophage", "cDC", "pDC", "Neutrophil", "Mast"}
SKIN_PARENCHYMA = {"Keratinocyte", "Fibroblast", "Endothelial", "Pericyte_SMC",
                   "Schwann_neural", "Adipocyte", "Eccrine"}

kp = pd.read_csv(KPYRN)
kp_scored = kp.loc[kp["scoring_include"] == "Yes", "gene_symbol"].tolist()
kp_modules = (kp.loc[kp["scoring_include"] == "Yes"]  # noqa: E712
                 .groupby("module")["gene_symbol"].apply(list).to_dict())
kp_order = kp["gene_symbol"].tolist()  # 固定列序（全 140 基因）
log(f"lineages: {len(lineages)}; KpyRN scored genes: {len(kp_scored)}; "
    f"modules: {list(kp_modules)}")

HUB_PANEL = ["ADH1B", "ALDH3A2", "DHX58", "FAHD2B", "GBP1", "LDHA",
             "PGAM1", "PGK1", "PFKP", "TPI1", "GAPDH", "ENO1", "ALDOA",
             "PGM2", "BPGM", "LDHB", "SLC16A1", "PDHB", "MPC2",
             "IFI44L", "IFI27", "STAT2", "DDX58", "ISG20",
             "HAT1", "EP300", "SIRT3", "STAT1", "FGFR4", "LGALS8", "RSPH4A"]

# ---------- 基因全集（54 样本并集） ----------
gene_files = sorted(GSE165816.glob("*_genes.csv"))
universe = set()
for gf in gene_files:
    universe.update(pd.read_csv(gf, header=None)[0].tolist())
universe = sorted(universe)
log(f"gene universe (union of {len(gene_files)} samples): {len(universe)}")

mk_cov = {}
for ln in lineages:
    present = [g for g in lineage_genes[ln] if g in universe]
    mk_cov[ln] = {"n_markers": len(lineage_genes[ln]), "n_present": len(present),
                  "missing": [g for g in lineage_genes[ln] if g not in universe]}
    if len(present) < 2:
        log(f"WARN lineage {ln}: only {len(present)} markers present in dataset")
kp_present = [g for g in kp_order if g in universe]
log(f"KpyRN present in universe: {len(kp_present)}/{len(kp_order)}")

# ---------- 细胞元数据（QC 通过者） ----------
qc = pd.read_csv(QCDIR / "scRNA_GSE165816_cell_qc.csv")
qc_pass = qc[qc["pass_qc"]]
log(f"QC-pass cells: {len(qc_pass)}/{len(qc)}")

# ---------- 逐样本处理 ----------
recs = []
count_blocks = []
tissue_check = {}

for i, gf in enumerate(gene_files):
    gsm = gf.stem.replace("_genes", "")
    X = sp.load_npz(GSE165816 / f"{gsm}_counts.npz")  # cells × genes
    genes = pd.read_csv(gf, header=None)[0].tolist()
    gidx = {g: j for j, g in enumerate(genes)}
    n_cells, n_genes = X.shape
    n_counts = np.asarray(X.sum(axis=1)).ravel()
    # CP10K + log1p（稀疏安全）
    Xn = X.copy().astype("float32")
    scale = 1e4 / np.maximum(n_counts, 1.0)
    row_of_data = np.repeat(np.arange(n_cells), np.diff(Xn.indptr))
    Xn.data = np.log1p(Xn.data * scale[row_of_data])
    # 谱系得分（缺失标记按 0）
    scores = np.zeros((n_cells, len(lineages)), dtype="float32")
    for li, ln in enumerate(lineages):
        cols = [gidx[g] for g in lineage_genes[ln] if g in gidx]
        if cols:
            scores[:, li] = np.asarray(Xn[:, cols].mean(axis=1)).ravel()
    top1 = scores.max(axis=1)
    top2 = np.sort(scores, axis=1)[:, -2]
    assign = np.array([lineages[j] for j in scores.argmax(axis=1)], dtype=object)
    margin = top1 - top2
    assign = np.where(margin >= MARGIN_MIN, assign, "Ambiguous")
    parent = np.array([LINEAGE_PARENT.get(a, a) for a in assign], dtype=object)
    # KpyRN 模块评分（仅评分模块，缺失基因按 0）
    mod_scores = {}
    for mod, gl in kp_modules.items():
        cols = [gidx[g] for g in gl if g in gidx]
        mod_scores[mod] = (np.asarray(Xn[:, cols].mean(axis=1)).ravel() if cols
                           else np.zeros(n_cells, dtype="float32"))
    # 枢纽基因面板值
    panel_vals = {}
    for g in HUB_PANEL:
        panel_vals[g] = (np.asarray(Xn[:, gidx[g]].todense()).ravel()
                         if g in gidx else np.zeros(n_cells, dtype="float32"))
    # KpyRN 原始计数块（全 140 基因固定列序，缺失填 0）
    block = np.zeros((n_cells, len(kp_order)), dtype="float32")
    for j, g in enumerate(kp_order):
        if g in gidx:
            block[:, j] = np.asarray(X[:, gidx[g]].todense()).ravel()
    count_blocks.append(block)
    # 记录（列数组，避免逐细胞 Python 循环）
    sub = qc_pass[qc_pass["gsm"] == gsm].reset_index(drop=True)
    assert len(sub) == n_cells, f"{gsm}: qc rows {len(sub)} != matrix {n_cells}"
    rec = {"gsm": np.repeat(gsm, n_cells),
           "barcode": sub["barcode"].to_numpy(),
           "tissue": sub["tissue"].to_numpy(),
           "disease": sub["disease"].to_numpy(),
           "n_counts": n_counts.astype("float32"),
           "lineage": assign, "lineage_parent": parent,
           "margin": margin.astype("float32")}
    for li, ln in enumerate(lineages):
        rec[f"score_{ln}"] = scores[:, li]
    for mod, v in mod_scores.items():
        rec[f"kpyrn_{mod}"] = v.astype("float32")
    for g, v in panel_vals.items():
        rec[g] = v.astype("float32")
    recs.append(rec)
    # sanity check（组织学合理性）
    par = pd.Series(parent)
    tissue_check[gsm] = {"tissue": str(sub.at[0, "tissue"]),
                         "frac_blood": round(float(par.isin(BLOOD_LINEAGES).mean()), 4),
                         "frac_skin": round(float(par.isin(SKIN_PARENCHYMA).mean()), 4),
                         "frac_ambiguous": round(float((par == "Ambiguous").mean()), 4)}
    if (i + 1) % 10 == 0 or i == len(gene_files) - 1:
        log(f"scored {i + 1}/{len(gene_files)} samples ({time.time() - t0:.0f}s)")

# ---------- 输出 ----------
df = pd.concat([pd.DataFrame(r) for r in recs], ignore_index=True)
df.to_csv(SCDIR / "scRNA_GSE165816_cell_scores.csv", index=False)
counts_all = np.vstack(count_blocks) if count_blocks else np.zeros((0, 0))
sp.save_npz(SCDIR / "scRNA_GSE165816_kpyrn_counts.npz", sp.csr_matrix(counts_all))
pd.Series(kp_order).to_csv(SCDIR / "scRNA_GSE165816_kpyrn_genes.csv",
                           index=False, header=False)

tc = pd.DataFrame(tissue_check).T
tc["n_cells"] = tc.index.map(df["gsm"].value_counts())
summary = {
    "n_cells_scored": int(len(df)),
    "n_lineages": len(lineages),
    "lineage_marker_coverage": mk_cov,
    "kpyrn_present": f"{len(kp_present)}/{len(kp_order)}",
    "margin_threshold": MARGIN_MIN,
    "ambiguous_fraction": round(float((df["lineage"] == "Ambiguous").mean()), 4),
    "tissue_sanity": tc.reset_index().rename(columns={"index": "gsm"}).to_dict("records"),
    "runtime_sec": round(time.time() - t0, 1),
}
with open(QCDIR / "scRNA_GSE165816_stage2_summary.json", "w", encoding="utf-8") as fh:
    json.dump(summary, fh, ensure_ascii=False, indent=2)

log("--- tissue sanity (mean over samples) ---")
log(tc.groupby("tissue")[["frac_blood", "frac_skin", "frac_ambiguous"]].mean().round(3).to_string())
log("--- lineage composition (all cells) ---")
log(df["lineage_parent"].value_counts(normalize=True).head(12).round(4).to_string())
log(f"DONE: {len(df)} cells scored; ambiguous {summary['ambiguous_fraction']:.1%}; "
    f"{time.time() - t0:.0f}s")
