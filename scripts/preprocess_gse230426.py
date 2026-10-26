# -*- coding: utf-8 -*-
"""
preprocess_gse230426.py — GSE230426 感染队列预处理（任务 #11a）

数据：GSE230426_RAW.tar（38 个 tsv.gz，PBMC RNA-seq raw counts，GENCODE ENSG ID）
设计（Lavery et al. 2023, PMID 37434783；series matrix 2026-10-25 核实）：
  S. aureus 感染 DFU 患者 PBMC；Stage: BEFORE THERAPY(0wk, n=20) vs AFTER THERAPY
  (n=18: 15×8wk + 3×4wk)；愈合结局字段未在公开元数据提供 → 结局对比不构建（如实记录）。
  配对结构：文件名首数字为样本级 ID，无显式 donor 列；仅能按"连续 ID + 同年龄 +
  互补时点"推断配对（17 对 + 4 单例）→ 主分析采用非配对模型（保守），配对模型
  仅作敏感性分析（在 R 脚本中，标注推断性）。

输出：
  01_数据/处理后数据/GSE230426_counts.csv        （symbol × 样本 raw counts）
  01_数据/处理后数据/GSE230426_annotation.csv    （样本注释）
  03_结果/感染队列GSE230426/preprocess_report.json
"""
import gzip

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import json
import re
import tarfile
from pathlib import Path

import pandas as pd

BASE = Path(_ROOT)
RAW = BASE / "01_数据/原始数据/GSE230426/GSE230426_RAW.tar"
PROC = BASE / "01_数据/处理后数据"
OUT = BASE / "03_结果/感染队列GSE230426"
PROC.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

# ---------- 1. 样本注释：文件名 donor_id 年龄 时点 + series matrix 的 Stage ----------
smeta = {}
with gzip.open(BASE / "01_数据/元数据/GSE230426_series_matrix.txt.gz", "rt",
               encoding="utf-8", errors="ignore") as fh:
    titles, stages = None, None
    for line in fh:
        if line.startswith("!Sample_title"):
            titles = [s.strip('"') for s in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1") and "Stage:" in line:
            stages = [s.strip('"').split("Stage: ")[1] for s in
                      line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!series_id") and titles is not None and stages is not None:
            break
assert titles and stages and len(titles) == len(stages) == 38, \
    f"series matrix 解析异常: {len(titles)} titles, {len(stages)} stages"
for t, st in zip(titles, stages):
    smeta[t] = st

rows = []
for t in titles:
    m = re.match(r"^(\d+)_(\d+)_(\d+)wk$", t)
    assert m, f"文件名模式不匹配: {t}"
    sid, age, wk = m.group(1), int(m.group(2)), int(m.group(3))
    stage = "BEFORE" if "BEFORE" in smeta[t] else "AFTER"
    assert stage == ("BEFORE" if wk == 0 else "AFTER"), \
        f"时点与 Stage 不一致: {t} → {smeta[t]}"
    rows.append({"sample": t, "sample_id": sid, "age": age, "week": wk,
                 "stage": stage})
ann = pd.DataFrame(rows)
print(f"样本注释: {len(ann)}（BEFORE n={sum(ann.stage=='BEFORE')}, "
      f"AFTER n={sum(ann.stage=='AFTER')}）")
print(ann.groupby(["stage", "week"]).size().to_string())

# 配对推断：连续 sample_id（数值序）+ 同年龄 + 互补时点
ann["_id"] = ann["sample_id"].astype(int)
ann = ann.sort_values("_id").reset_index(drop=True)
donor = [None] * len(ann)
d = 0
i = 0
while i < len(ann):
    if (i + 1 < len(ann)
            and int(ann.loc[i, "sample_id"]) + 1 == int(ann.loc[i + 1, "sample_id"])
            and ann.loc[i, "age"] == ann.loc[i + 1, "age"]
            and ann.loc[i, "stage"] != ann.loc[i + 1, "stage"]):
        d += 1
        donor[i] = donor[i + 1] = f"D{d:02d}"
        i += 2
    else:
        donor[i] = "unpaired"
        i += 1
ann["donor_inferred"] = donor
n_pair = sum(1 for x in donor if x != "unpaired") // 2
print(f"配对推断（连续ID+同年龄+互补时点）：{n_pair} 对 + "
      f"{sum(1 for x in donor if x=='unpaired')} 单例（仅用于敏感性分析）")

# ---------- 2. counts 矩阵 ----------
ensg_map = json.load(open(BASE / "01_数据/元数据/ensg_symbol_map.json",
                          encoding="utf-8"))
counts = {}
with tarfile.open(RAW, "r") as tf:
    for mem in tf.getmembers():
        if not mem.name.endswith(".tsv.gz"):
            continue
        title = re.sub(r"^GSM\d+_", "", Path(mem.name).stem.replace(".tsv", ""))
        with tf.extractfile(mem) as fh:
            with gzip.open(fh, "rt", encoding="utf-8") as gz:
                df = pd.read_csv(gz, sep=",", index_col=0)
        s = df.iloc[:, 0]
        counts[title] = s
        assert len(s) == 58288, f"{title}: {len(s)} 基因（预期 58288）"

mat = pd.DataFrame(counts)
assert set(mat.columns) == set(ann["sample"]), "样本列与注释不对应"
mat = mat[ann["sample"].tolist()]  # 按注释顺序（ID 数值序）排列
print(f"counts 矩阵: {mat.shape[0]} ENSG × {mat.shape[1]} 样本")

# ENSG（去版本号）→ symbol；多对一取总和（与主分析口径一致则记录）
# 别名规范化：映射文件将 ENSG00000107201 记为别名 RIGI，官方 HGNC 符号为 DDX58
# （M4 IFN-I 读出成员 + ML 共识枢纽基因，必须回收）。其余 13 个"缺失"基因全部
# 为不参与评分的 M5 组蛋白载体模块（HGNC 改名/未映射），不影响模块分析，如实记录。
ALIAS_FIX = {"RIGI": "DDX58"}
mat.index = [g.split(".")[0] for g in mat.index]
sym = pd.Series([ALIAS_FIX.get(ensg_map.get(g, g), ensg_map.get(g, g))
                 for g in mat.index], index=mat.index)
dup = sym[sym.duplicated(keep=False)]
mat_sym = mat.groupby(sym).sum()
print(f"symbol 映射: {mat_sym.shape[0]} 唯一 symbol；"
      f"多对一基因 {dup.index.nunique()} 个（取 sum）")
n_unmapped = sum(1 for g in sym.index if g == ensg_map.get(g, g) and g.startswith("ENSG") and g not in ensg_map)
print(f"未映射 ENSG（保留原 ID）: {n_unmapped}")

mat_sym.to_csv(PROC / "GSE230426_counts.csv")
ann.set_index("sample").to_csv(PROC / "GSE230426_annotation.csv")

# ---------- 3. 预处理报告 ----------
KpyRN = pd.read_csv(BASE / "01_数据/基因集/KpyRN_modules.csv")
in_data = KpyRN[KpyRN["gene_symbol"].isin(mat_sym.index)]
report = {
    "n_samples": int(len(ann)),
    "stage_counts": {"BEFORE": int(sum(ann.stage == "BEFORE")),
                     "AFTER": int(sum(ann.stage == "AFTER"))},
    "after_week": {"8wk": int(sum((ann.stage == "AFTER") & (ann.week == 8))),
                   "4wk": int(sum((ann.stage == "AFTER") & (ann.week == 4)))},
    "paired_inferred": {"pairs": int(n_pair),
                        "unpaired": int(sum(1 for x in donor if x == "unpaired")),
                        "note": "连续ID+同年龄+互补时点推断，无显式donor列；主分析非配对"},
    "outcome_contrast": "NOT_POSSIBLE: 愈合结局字段未在公开元数据提供",
    "genes_ensg": int(mat.shape[0]),
    "genes_symbol": int(mat_sym.shape[0]),
    "kpyrn_genes_in_matrix": int(len(in_data)),
    "kpyrn_total": int(len(KpyRN)),
    "kpyrn_missing": sorted(set(KpyRN["gene_symbol"]) - set(mat_sym.index)),
}
(OUT / "preprocess_report.json").write_text(
    json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(report, ensure_ascii=False, indent=2))
print("DONE — preprocess_gse230426.py")
