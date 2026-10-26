# -*- coding: utf-8 -*-
"""
docking_summary.py — 分子对接结果汇总与验收（任务 #10 配套）
输入：03_结果/药物预测/docking/docking_results.csv、positive_control_rmsd.csv
输出：03_结果/药物预测/docking/docking_summary.md  （稿件用结果表 + 验收判定）
纪律：结合能必须以阳性对照为基准（delta_vs_pos）；措辞仅"提示需实验验证的候选"
"""
import re

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

from pathlib import Path

import pandas as pd

BASE = Path(_ROOT)
DOCK = BASE / "03_结果/药物预测/docking"


def sani(s):
    return re.sub(r"[^A-Za-z0-9._-]", "", str(s))


def short_name(tag):
    if tag.startswith("POS_"):
        return f"{tag.split('_')[-1]} (co-crystal)"
    if tag.startswith("NEG_"):
        return "sucrose (negative)"
    if tag.startswith("LIG_"):
        return tag[4:]
    return tag


def main():
    df = pd.read_csv(DOCK / "docking_results.csv")
    rdf = pd.read_csv(DOCK / "positive_control_rmsd.csv")

    # ---- 验收 1：阳性对照重对接 RMSD ----
    rmsd_ok = rdf.dropna(subset=["rmsd"])
    n_pass = int((rmsd_ok["rmsd"] < 2.0).sum())
    n_total = len(rmsd_ok)
    print(f"[验收] 阳性对照重对接 RMSD <2 Å：{n_pass}/{n_total} 通过")

    # ---- 各体系：找最有利候选（相对阳性对照 delta 最负 = 结合更强）----
    out = ["# 分子对接结果汇总（AutoDock Vina 1.2.7）", ""]
    out.append("## 验收判定")
    out.append(f"- 阳性对照重对接 RMSD <2 Å（协议可靠）：**{n_pass}/{n_total} 通过**")
    out.append("- 结合能基准：每个体系以共晶配体（阳性对照）亲和力为基准，报告 Δ(kcal/mol)")
    out.append("- 措辞纪律：对接结果仅提示**需实验验证的候选**，不作为结合证据的唯一依据")
    out.append("")

    systems = list(dict.fromkeys(zip(df["target"], df["site"])))
    top_rows = []
    for tgt, site in systems:
        sub = df[(df["target"] == tgt) & (df["site"] == site)].copy()
        pos = sub[sub["role"] == "positive_control"]
        pos_aff = float(pos["affinity_kcal_mol"].iloc[0]) if len(pos) else float("nan")
        cands = sub[sub["role"] == "candidate"].dropna(
            subset=["affinity_kcal_mol"]).sort_values("affinity_kcal_mol")
        neg = sub[sub["role"] == "negative_control"]
        neg_aff = float(neg["affinity_kcal_mol"].iloc[0]) if len(neg) else float("nan")
        out.append(f"## {tgt} — {site}  (PDB {sub['pdb'].iloc[0]})")
        out.append(f"- 阳性对照（共晶配体）亲和力：**{pos_aff:.2f}** kcal/mol")
        out.append(f"- 负面对照（蔗糖）亲和力：{neg_aff:.2f} kcal/mol")
        out.append("")
        out.append("| 候选 | 亲和力 (kcal/mol) | Δ vs 共晶 | 备注 |")
        out.append("|---|---|---|---|")
        for _, r in cands.iterrows():
            dlt = r["delta_vs_pos"]
            dlt_s = f"{dlt:+.2f}" if pd.notna(dlt) else "—"
            out.append(f"| {short_name(r['ligand'])} | {r['affinity_kcal_mol']:.2f} | "
                       f"{dlt_s} | 提示需实验验证 |")
        out.append("")
        if len(cands):
            best = cands.iloc[0]
            best_dlt = best["delta_vs_pos"]
            flag = "更强（Δ<0）" if pd.notna(best_dlt) and best_dlt < 0 else "弱于共晶"
            top_rows.append((tgt, site, short_name(best["ligand"]),
                             best["affinity_kcal_mol"], best_dlt, flag))
    out.append("## 各体系最优候选（按亲和力）")
    out.append("| 靶点 | 位点 | 最优候选 | 亲和力 | Δ vs 共晶 | 相对共晶 |")
    out.append("|---|---|---|---|---|---|")
    for tgt, site, name, aff, dlt, flag in top_rows:
        dlt_s = f"{dlt:+.2f}" if pd.notna(dlt) else "—"
        out.append(f"| {tgt} | {site.split()[0]} | {name} | {aff:.2f} | {dlt_s} | {flag} |")
    out.append("")
    out.append("> 说明：Δ = 候选亲和力 − 共晶配体亲和力；负值表示对接打分更强于共晶配体。"
               "对接打分仅反映几何与经验势能，须实验（如 SPR/ITC/酶活/细胞实验）验证。")
    (DOCK / "docking_summary.md").write_text("\n".join(out), encoding="utf-8")
    print("汇总已写入 docking_summary.md")
    print(f"[概览] 共 {len(systems)} 个体系，{len(df)} 次对接")


if __name__ == "__main__":
    main()
