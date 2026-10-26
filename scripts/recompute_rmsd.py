# -*- coding: utf-8 -*-
"""
recompute_rmsd.py — 定点重算阳性对照重对接 RMSD（任务 #10 验证缺口修复）

背景：docking_pipeline.py 的 parse_pdbqt_atoms() 未区分 Vina 输出的多个 MODEL 块，
把 9 个 pose 的全部原子累加，导致 1OX/TF1(2 MODEL)、PV2/TF4(9 MODEL) 原子数翻倍、
与参考（单块）不匹配而无法计算 RMSD。FBP/ACO/AR6/ZHK 因恰好只有 1 个 MODEL 而
侥幸通过。

修复：仅取 pose 的 MODEL 1（最优 pose，Vina 已按能量排序）计算 Kabsch 重原子 RMSD，
与参考晶体坐标（单块）按索引配对。复用 docking_pipeline.py 已产出的 ref/pose 文件，
不重跑对接。

输出：覆盖 03_结果/药物预测/docking/positive_control_rmsd.csv
"""
import re

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path(_ROOT)
WORK = BASE / "01_数据/分子对接/work"
OUT = BASE / "03_结果/药物预测/docking"


def sani(s):
    """与 docking_pipeline.py 完全一致的标签消毒。"""
    return re.sub(r"[^A-Za-z0-9._-]", "", s)


# (target, 体系 site 全名, pdb, 共晶配体 resn) —— 文件名由 sani 派生，避免手写不匹配
SYSTEMS = [
    ("PKM2", "C-pocket (allosteric)", "4JPG", "1OX"),
    ("PKM2", "A-pocket (FBP site)",   "4JPG", "FBP"),
    ("EP300", "acetyl-CoA substrate site", "4PZS", "ACO"),
    ("SIRT3", "ADP-ribose/NAD+ C-pocket", "4BN4", "AR6"),
    ("LDHA", "active site", "4JNK", "ZHK"),
    ("PDK2", "ATP pocket", "4MPC", "PV2"),
    ("PDK2", "pyruvate site (R-domain)", "2BU8", "TF4"),
    ("PDK2", "R-domain N-terminal site", "2BU5", "TF1"),
]


def parse_pdbqt_atoms(path, first_model_only=False):
    atoms = []
    model_seen = 0
    for line in open(path, encoding="utf-8", errors="ignore"):
        if line.startswith("MODEL"):
            model_seen += 1
            if first_model_only and model_seen > 1:
                break
            continue
        if line.startswith("ENDMDL"):
            if first_model_only and model_seen >= 1:
                break
            continue
        if line.startswith(("ATOM", "HETATM")):
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                el = line[77:79].strip()
                nm = line[12:16].strip()
                atoms.append((x, y, z, el, nm))
            except ValueError:
                continue
    return atoms


def kabsch_rmsd(P, Q):
    P = np.asarray(P); Q = np.asarray(Q)
    assert P.shape == Q.shape, f"shape mismatch {P.shape} vs {Q.shape}"
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return float(np.sqrt(((Pc - Qc @ R) ** 2).sum() / len(P)))


rows = []
for tgt, site, pdb, resn in SYSTEMS:
    site_key = sani(site.split()[0])
    pos_tag = f"POS_{tgt}_{site_key}_{resn}"
    ref = WORK / "ligands" / f"{pos_tag}.pdbqt"
    pose = WORK / "poses" / f"{pdb}_{site_key}_{pos_tag}.pdbqt"
    if not (ref.exists() and pose.exists()):
        rows.append({"target": tgt, "site": site, "pos_lig": resn,
                     "n_atoms": None, "rmsd": None, "note": "file missing"})
        continue
    A = [a for a in parse_pdbqt_atoms(ref) if a[3] != "H"]                 # 参考：单块
    B = [a for a in parse_pdbqt_atoms(pose, first_model_only=True) if a[3] != "H"]  # pose：仅 MODEL 1
    nA, nB = len(A), len(B)
    if nA != nB or nA < 5:
        rows.append({"target": tgt, "site": site, "pos_lig": resn,
                     "n_atoms": min(nA, nB), "rmsd": None,
                     "note": f"atom count mismatch (ref={nA}, pose={nB})"})
        continue
    name_ok = all(A[i][4] == B[i][4] for i in range(nA))
    rmsd = kabsch_rmsd([a[:3] for a in A], [a[:3] for a in B])
    rows.append({"target": tgt, "site": site, "pos_lig": resn, "n_atoms": nA,
                 "rmsd": round(rmsd, 2),
                 "note": ("heavy-atom RMSD, model 1, index-paired"
                          + ("" if name_ok else " (atom-name mismatch, index-paired)"))})
    flag = "OK" if rmsd < 2.0 else "!! >2A"
    print(f"  {tgt:6s} {site[:24]:24s} {resn:4s}: RMSD={rmsd:.2f} Å (n={nA}) {flag}")

rdf = pd.DataFrame(rows)
rdf.to_csv(OUT / "positive_control_rmsd.csv", index=False)
print(f"\n重写 positive_control_rmsd.csv：{len(rdf)} 体系")
n_ok = (rdf["rmsd"] < 2.0).sum()
print(f"  RMSD<2 Å 体系数：{n_ok}/{len(rdf)}")
