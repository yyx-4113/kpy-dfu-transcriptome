# -*- coding: utf-8 -*-
"""
统计每个对接配体的重原子数（AutoDock atom type 口径）与可旋转键数（BRANCH），
并计算配体效率 LE = Vina 亲和力 / 重原子数。

目的：稿件 §3.9 需要说明「Δ = 候选 − 共晶对照」未做分子尺寸归一，
以及按配体效率口径重新判读时结论是否反转。本脚本产出可追溯的产物表 S9c。

口径说明（写进产物 meta，避免与"元素列"口径混淆）：
  - 重原子 = AutoDock atom type 属于 {C/A, N/NA, O/OA, S/SA, P, F, Cl/CL, Br/BR, I, Mg, Zn, Mn, Ca, Fe} 的原子，
    即排除极性氢 H/HD 与虚拟原子。
  - 可旋转键 = PDBQT 中 BRANCH 记录数。
  - LE = affinity_kcal_per_mol / n_heavy_atoms（越负越好）。
"""
import csv
import glob
import json
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(BASE, "01_数据", "分子对接", "work")
LIG = os.path.join(WORK, "ligands")
OUT = os.path.join(BASE, "05_稿件", "附表")

HEAVY_TYPES = {
    "C", "A", "N", "NA", "NS", "O", "OA", "OS", "S", "SA", "P", "F",
    "Cl", "CL", "Br", "BR", "I", "Mg", "MG", "Zn", "ZN", "Mn", "MN",
    "Ca", "CA", "Fe", "FE",
}


def parse_pdbqt(path):
    """返回 (重原子数, BRANCH 数)。PDBQT 列不固定，按空白切分取末字段为 atom type。"""
    n_heavy = 0
    n_branch = 0
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("BRANCH"):
                n_branch += 1
            elif line.startswith(("ATOM", "HETATM")):
                tok = line.split()
                if len(tok) >= 2:
                    atype = tok[-1]
                    if atype in HEAVY_TYPES:
                        n_heavy += 1
    return n_heavy, n_branch


def main():
    aff = {}   # ligand_key -> affinity
    aff_path = os.path.join(BASE, "05_稿件", "附表", "S9a_docking_results.csv")
    with open(aff_path, encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            key = (r.get("ligand") or r.get("ligand_id") or "").strip()
            try:
                a = float(r["affinity_kcal_per_mol"])
            except (KeyError, TypeError, ValueError):
                for k, v in r.items():
                    if "aff" in k.lower():
                        try:
                            a = float(v)
                            break
                        except (TypeError, ValueError):
                            pass
                else:
                    continue
            aff.setdefault(key, a)
    # 若 S9a 无 affinity 列，退回产物目录
    if not aff:
        alt = os.path.join(BASE, "03_结果", "药物预测", "docking", "docking_results.csv")
        with open(alt, encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh):
                key = (r.get("ligand") or r.get("ligand_id") or "").strip()
                for k, v in r.items():
                    if "aff" in k.lower():
                        try:
                            aff.setdefault(key, float(v))
                            break
                        except (TypeError, ValueError):
                            pass

    rows = []
    for f in sorted(glob.glob(os.path.join(LIG, "*.pdbqt"))):
        name = os.path.basename(f)[:-6]
        n_h, n_b = parse_pdbqt(f)
        rows.append({
            "ligand_file": os.path.basename(f),
            "ligand_key": name,
            "n_heavy_atoms": n_h,
            "n_rotatable_bonds": n_b,
            "rigid": "yes" if n_b == 0 else "no",
        })

    # 关联亲和力：按 ligand_key 的尾部 token 匹配
    for r in rows:
        key = r["ligand_key"]
        a = aff.get(key)
        if a is None:
            tail = key.split("_")[-1]
            for k, v in aff.items():
                if k == tail or k.endswith("_" + tail) or tail in k:
                    a = v
                    break
        r["affinity_kcal_per_mol"] = "" if a is None else round(a, 2)
        r["ligand_efficiency"] = "" if (a is None or not r["n_heavy_atoms"]) \
            else round(a / r["n_heavy_atoms"], 3)

    dst = os.path.join(OUT, "S9c_ligand_size_metrics.csv")
    with open(dst, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ── 自检 ──
    pos = [r for r in rows if r["ligand_key"].startswith("POS_")]
    neg = [r for r in rows if r["ligand_key"].startswith("NEG_")]
    assert pos, "未找到共晶阳性对照配体"
    assert all(r["rigid"] == "yes" for r in pos), \
        f"共晶对照应为刚性，实际存在柔性：{[r['ligand_file'] for r in pos if r['rigid'] != 'yes']}"
    flex_cand = [r for r in rows if not r["ligand_key"].startswith(("POS_", "NEG_"))
                 and r["rigid"] == "no"]
    print(f"产物：{dst}（{len(rows)} 行）")
    print(f"自检1 共晶对照 {len(pos)} 个，全部刚性（BRANCH=0）："
          f"{'通过' if all(r['rigid'] == 'yes' for r in pos) else '失败'}")
    print(f"自检2 柔性候选 {len(flex_cand)} 个，阴性对照 {len(neg)} 个")
    print(f"自检3 重原子数均 > 0："
          f"{'通过' if all(r['n_heavy_atoms'] > 0 for r in rows) else '失败'}")
    for r in pos:
        print(f"  {r['ligand_key']:34s} heavy={r['n_heavy_atoms']:3d} "
              f"aff={r['affinity_kcal_per_mol']:>7} LE={r['ligand_efficiency']}")
    matched = sum(1 for r in rows if r["affinity_kcal_per_mol"] != "")
    print(f"自检4 亲和力匹配 {matched}/{len(rows)}"
          f"{'（未匹配者不参与稿件引用）' if matched < len(rows) else ''}")


if __name__ == "__main__":
    main()
