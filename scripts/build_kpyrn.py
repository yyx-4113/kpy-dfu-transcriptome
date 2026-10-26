# -*- coding: utf-8 -*-
"""
build_kpyrn.py — 构建丙酮酰化调控网络（KpyRN）基因集
项目：丙酮酰化（Kpy）与糖尿病足（DFU）
溯源原则：每个基因必须能回溯到 (a) Kpy 原始文献或 (b) KEGG/Reactome 数据库条目；不编造。
输出：01_数据/基因集/KpyRN_modules.csv / .json / README.md
"""
import csv

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import json
import os
import re

BASE = _ROOT
META = f"{BASE}/01_数据/元数据"
OUT = f"{BASE}/01_数据/基因集"
os.makedirs(OUT, exist_ok=True)

CELL26 = "Cell 2026 (Zuo et al., doi:10.1016/j.cell.2026.01.023)"
NATM26 = "Nat Metab 2026 (Song et al., doi:10.1038/s42255-026-01556-2)"

# ---------- 1. 解析 KEGG 通路基因（带符号） ----------
def parse_kegg_genes(path: str) -> dict:
    """返回 {symbol: {'kegg_id':..., 'desc':..., 'ec':...}}"""
    out = {}
    with open(path, encoding="utf-8") as f:
        txt = f.read()
    gene_section = re.search(r"\nGENE\s+(.*?)(?:\n[A-Z]|\Z)", txt, re.S)
    if not gene_section:
        return out
    for line in gene_section.group(1).splitlines():
        m = re.match(r"\s+(\d+)\s+([A-Z0-9\-]+);\s+(.*)", line)
        if m:
            gid, sym, desc = m.group(1), m.group(2), m.group(3)
            ec = re.search(r"\[EC:(.*?)\]", desc)
            out[sym] = {
                "kegg_id": gid,
                "desc": re.sub(r"\s*\[KO:.*?\]", "", desc).strip(),
                "ec": ec.group(1) if ec else "",
            }
    return out


kegg_gly = parse_kegg_genes(f"{META}/kegg_hsa00010_entry.txt")   # Glycolysis / Gluconeogenesis
kegg_pyr = parse_kegg_genes(f"{META}/kegg_hsa00620_entry.txt")   # Pyruvate metabolism
kegg_jak = parse_kegg_genes(f"{META}/kegg_hsa04630_entry.txt")   # JAK-STAT signaling pathway
print(f"KEGG hsa00010: {len(kegg_gly)} genes; hsa00620: {len(kegg_pyr)}; hsa04630: {len(kegg_jak)}")

rows = []  # dict rows


def add(sym, module, submodule, role, source, evidence, scoring="Yes"):
    rows.append({
        "gene_symbol": sym, "module": module, "submodule": submodule,
        "role": role, "source": source, "evidence": evidence,
        "scoring_include": scoring,
    })


# ---------- 2. M1：Kpy 写入/擦除机器 + 已鉴定底物 ----------
M1 = [
    ("HAT1", "writer", NATM26, "组蛋白乙酰转移酶 HAT1 被证实催化丙酮酰辅酶A生成 Kpy（体外+细胞内验证）"),
    ("EP300", "writer", NATM26, "p300 被证实催化 Kpy 生成；p300/CBP 家族为经典酰基转移酶"),
    ("SIRT3", "eraser", NATM26, "NAD+依赖去酰化酶 SIRT3 负责去除 Kpy（可逆调控环路）"),
    ("STAT1", "substrate", CELL26, "STAT1-K201 为实验鉴定的丙酮酰化位点（+70.0468 Da）；PKM2 依赖；阻断 STAT1-STAT2 二聚化"),
]
for sym, role, src, ev in M1:
    add(sym, "M1_Kpy_machinery", "writer/eraser/substrate", role, src, ev)

# ---------- 3. M2：丙酮酸供应与分流（KEGG 溯源 + 文献补充） ----------
# 3a. KEGG Glycolysis/Gluconeogenesis
for sym, info in sorted(kegg_gly.items()):
    ev = f"KEGG hsa00010 Glycolysis/Gluconeogenesis (kegg:{info['kegg_id']}; {info['desc'][:80]})"
    role = "glycolysis enzyme"
    if sym in ("HK1", "HK2", "GCK"):
        role = "glycolysis enzyme (upstream of pyruvate; Cell 2026 shRNA screen)"
        ev += f"; {CELL26} shRNA screen: 敲低解除 IFN 抑制"
    elif sym in ("PFKM", "PFKL", "PFKP"):
        role = "glycolysis enzyme (upstream of pyruvate; Cell 2026 shRNA screen)"
        ev += f"; {CELL26} shRNA screen (PFK1)"
    elif sym in ("ENO1", "ENO2"):
        role = "glycolysis enzyme (upstream of pyruvate; Cell 2026 shRNA screen)"
        ev += f"; {CELL26} shRNA screen"
    elif sym == "PKM":
        role = "pyruvate-producing enzyme (PKM2 isoform)"
        ev += f"; {CELL26}: PKM2 催化产生丙酮酸为 Kpy 直接供体；敲低 PKM2 解除 IFN 抑制"
    elif sym in ("LDHA", "LDHB"):
        role = "pyruvate consumption (lactate branch)"
        ev += f"; {CELL26}: 敲低 LDHA 不解除 IFN 抑制 → LDH 位于 Kpy 节点下游"
    add(sym, "M2_Pyruvate_supply", "KEGG:hsa00010", role, "KEGG hsa00010", ev)

# 3b. KEGG Pyruvate metabolism
for sym, info in sorted(kegg_pyr.items()):
    if sym in {r["gene_symbol"] for r in rows}:
        # 交叉注释：同时属于糖酵解通路
        for r in rows:
            if r["gene_symbol"] == sym and r["module"] == "M2_Pyruvate_supply":
                r["submodule"] += " + KEGG:hsa00620"
                r["evidence"] += f"; 亦属 KEGG hsa00620 Pyruvate metabolism"
        continue
    add(sym, "M2_Pyruvate_supply", "KEGG:hsa00620", "pyruvate metabolism enzyme",
        "KEGG hsa00620", f"KEGG hsa00620 Pyruvate metabolism (kegg:{info['kegg_id']}; {info['desc'][:80]})")

# 3c. 文献补充：线粒体丙酮酸载体（MPC）复合体
for sym, ev in [
    ("MPC1", "线粒体丙酮酸载体亚基1（Bricker et al. 2012 / Herzig et al. 2012, Science）；Nat Cell Biol 2025 (Hsu et al.) 提示 ALDH4A1 为 MPC 复合体活性组分"),
    ("MPC2", "线粒体丙酮酸载体亚基2；同上"),
    ("SLC16A1", "质膜单羧酸转运体 MCT1，介导丙酮酸/乳酸跨膜转运；Nat Metab 2026 提及丙酮酸分配（routing）决定 Kpy 蓄积"),
    ("ALDH4A1", "Nat Cell Biol 2025 (Hsu et al.)：ALDH4A1 为 MPC 复合体活性组分，维持线粒体丙酮酸输入"),
]:
    if sym not in {r["gene_symbol"] for r in rows}:
        add(sym, "M2_Pyruvate_supply", "MPC/transport (literature)", "pyruvate transport", "Literature (Nat Cell Biol 2025; Nat Metab 2026)", ev)

# ---------- 4. M3：Kpy 转录输出靶基因 ----------
for sym in ["FGFR4", "LGALS8", "RSPH4A"]:
    add(sym, "M3_Kpy_targets", "CUT&Tag+RNA-seq targets", "Kpy target gene", NATM26,
        "启动子区 Kpy 富集且与转录激活正相关（CUT&Tag + RNA-seq 联合分析）")

# ---------- 5. M4：IFN-I 信号读出（STAT1-STAT2-IRF9-ISG 轴） ----------
M4_CORE = [
    ("IFNAR1", "IFN-I receptor subunit 1"),
    ("IFNAR2", "IFN-I receptor subunit 2"),
    ("JAK1", "receptor-associated kinase (IFN-I)"),
    ("TYK2", "receptor-associated kinase (IFN-I)"),
    ("STAT2", "STAT1 dimerization partner; K201 修饰阻断 STAT1-STAT2 结合"),
    ("IRF9", "ISGF3 component"),
    ("IRF1", "ISG transcription factor"),
]
M4_ISG = [
    "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM3", "MX1", "MX2",
    "OAS1", "OAS2", "OAS3", "OASL", "ISG15", "ISG20",
    "IFI6", "IFI16", "IFI27", "IFI44", "IFI44L", "RSAD2",
    "DDX58", "DHX58", "EIF2AK2", "CXCL10", "GBP1",
]
for sym, role in M4_CORE:
    jak_note = f"; KEGG hsa04630 JAK-STAT member" if sym in kegg_jak else ""
    add(sym, "M4_IFNI_readout", "IFN-I signaling", role, CELL26 + " + KEGG hsa04630",
        f"{CELL26}: IFN-I 信号通路核心组分，STAT1-K201 丙酮酰化的功能读出轴{jak_note}")
for sym in M4_ISG:
    jak_note = f"; KEGG hsa04630 JAK-STAT member" if sym in kegg_jak else ""
    ev = f"{CELL26}: 高糖条件下受抑制的干扰素刺激基因（ISG）类；Kpy-STAT1 修饰的功能读出"
    if sym in ("IFIT1", "ISG15", "RSAD2"):
        ev += "（Cell 2026 明确点名的 ISG）"
    add(sym, "M4_IFNI_readout", "ISG", "interferon-stimulated gene", CELL26 + " + canonical ISG list",
        ev + jak_note)

# ---------- 6. M5：组蛋白底物载体（不纳入评分） ----------
HISTONES = ["HIST1H1C", "HIST1H2AC", "HIST1H2AM", "HIST1H2BO", "HIST1H3B", "HIST1H4C",
            "HIST1H4E", "HIST2H2AC", "HIST2H2BE", "HIST2H3C", "HIST2H3D", "H2AFZ"]
for sym in HISTONES:
    add(sym, "M5_Histone_carriers", "histone substrate carrier", "substrate carrier", NATM26,
        "Nat Metab 2026 在哺乳动物细胞鉴定 10 个组蛋白 Kpy 位点（组蛋白 Kpy 富集于启动子）；"
        "本模块为组蛋白家族代表性基因。注意：组蛋白 mRNA 受细胞周期调控，不纳入表达评分（scoring_include=No），"
        "具体位点-组蛋白对应关系以原文为准。", scoring="No")

# ---------- 7. 输出 ----------
# STAT1 跨模块标注
for r in rows:
    if r["gene_symbol"] == "STAT1":
        r["evidence"] += "; 同时为 M4 IFN-I 读出轴核心（跨模块基因，主模块记 M1）"

with open(f"{OUT}/KpyRN_modules.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
with open(f"{OUT}/KpyRN_modules.json", "w", encoding="utf-8") as f:
    json.dump(rows, f, ensure_ascii=False, indent=2)

# 汇总
from collections import Counter
c = Counter(r["module"] for r in rows)
total = len(rows)
uniq = len({r["gene_symbol"] for r in rows})
print(f"\nTotal entries: {total} (unique genes: {uniq})")
for m, n in sorted(c.items()):
    scoring_yes = sum(1 for r in rows if r["module"] == m and r["scoring_include"] == "Yes")
    print(f"  {m:<24} n={n:<4} (scoring: {scoring_yes})")

# 与 M2 中 KEGG 基因一致性核对
m2 = {r["gene_symbol"] for r in rows if r["module"] == "M2_Pyruvate_supply"}
missing_from_kegg = {"MPC1", "MPC2", "SLC16A1", "ALDH4A1"} - m2
print(f"\nM2 total: {len(m2)}; literature-added present: {sorted({'MPC1','MPC2','SLC16A1','ALDH4A1'} & m2)}")
print(f"Output -> {OUT}/KpyRN_modules.csv")
