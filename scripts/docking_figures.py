# -*- coding: utf-8 -*-
"""
docking_figures.py — 分子对接结果图表生成（任务 #10 配套）
输入：03_结果/药物预测/docking/docking_results.csv、positive_control_rmsd.csv
输出：figures/fig_docking_rmsd.png      阳性对照重对接 RMSD 验证图
      figures/fig_docking_affinity.png  各体系结合能图（以阳性对照为基准）
      figures/fig_docking_poses_*.png   前位候选姿态叠加图（三视图投影）
"""
import re

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(_ROOT)
DOCK = BASE / "03_结果/药物预测/docking"
WORK = BASE / "01_数据/分子对接/work"
FIG = DOCK / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def sani(s):
    """文件名消毒：去掉 / \\ : * ? " < > | + 等非法字符（与 docking_pipeline.py 一致）。"""
    return re.sub(r"[^A-Za-z0-9._-]", "", s)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
})

ROLE_COLOR = {"positive_control": "#C0392B", "negative_control": "#7F8C8D",
              "candidate": "#2471A3"}
ROLE_LABEL = {"positive_control": "Positive control (co-crystallized)",
              "negative_control": "Negative control (sucrose)",
              "candidate": "Candidate"}


def short_name(tag):
    """LIG_CHEMBL5590652 -> CHEMBL5590652; POS_PKM2_C-pocket_1OX -> 1OX (pos)"""
    t = tag
    if t.startswith("POS_"):
        parts = t.split("_")
        return f"{parts[-1]} (pos)"
    if t.startswith("NEG_"):
        return "sucrose (neg)"
    if t.startswith("LIG_"):
        return t[4:].replace("_", " ")
    return t


# ============ 图 1：阳性对照 RMSD ============
def fig_rmsd():
    rdf = pd.read_csv(DOCK / "positive_control_rmsd.csv")
    rdf = rdf.dropna(subset=["rmsd"])
    if rdf.empty:
        print("!! 无有效 RMSD 记录，跳过图 1")
        return
    rdf["label"] = rdf["target"] + " / " + rdf["site"].str.split(" ").str[0] \
        + " (" + rdf["pos_lig"] + ")"
    rdf = rdf.sort_values("rmsd")
    fig, ax = plt.subplots(figsize=(7.2, 0.42 * len(rdf) + 1.6))
    colors = ["#1E8449" if r < 2 else "#C0392B" for r in rdf["rmsd"]]
    ax.barh(rdf["label"], rdf["rmsd"], color=colors, height=0.62)
    ax.axvline(2.0, color="#C0392B", ls="--", lw=1)
    ax.text(2.02, len(rdf) - 0.4, "2 A protocol threshold", color="#C0392B",
            fontsize=8, va="top")
    for y, (r, n) in enumerate(zip(rdf["rmsd"], rdf["n_atoms"])):
        ax.text(r + 0.03, y, f"{r:.2f} (n={n})", va="center", fontsize=8)
    ax.set_xlabel("Heavy-atom RMSD of re-docked co-crystallized ligand (A)")
    ax.set_xlim(0, max(rdf["rmsd"].max() * 1.35, 2.3))
    ax.set_title("Positive-control re-docking validation (AutoDock Vina)",
                 fontsize=10, pad=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig_docking_rmsd.png", bbox_inches="tight")
    plt.close(fig)
    print("图 1 已生成：fig_docking_rmsd.png")


# ============ 图 2：各体系结合能 ============
def fig_affinity():
    df = pd.read_csv(DOCK / "docking_results.csv")
    systems = list(dict.fromkeys(zip(df["target"], df["site"])))
    n = len(systems)
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol,
                             figsize=(4.1 * ncol, 2.5 * nrow), squeeze=False)
    for idx, (tgt, site) in enumerate(systems):
        ax = axes[idx // ncol][idx % ncol]
        sub = df[(df["target"] == tgt) & (df["site"] == site)].copy()
        sub = sub.sort_values("affinity_kcal_mol")
        pos_aff = sub.loc[sub["role"] == "positive_control",
                          "affinity_kcal_mol"]
        pos_aff = pos_aff.iloc[0] if len(pos_aff) else np.nan
        labels = [short_name(t) for t in sub["ligand"]]
        colors = [ROLE_COLOR[r] for r in sub["role"]]
        y = np.arange(len(sub))
        ax.barh(y, sub["affinity_kcal_mol"], color=colors, height=0.66)
        if not np.isnan(pos_aff):
            ax.axvline(pos_aff, color="#C0392B", ls=":", lw=1)
        for yi, (aff, dlt) in enumerate(zip(sub["affinity_kcal_mol"],
                                            sub["delta_vs_pos"])):
            txt = f"{aff:.1f}" if np.isnan(dlt) else f"{aff:.1f} ({dlt:+.1f})"
            ax.text(aff + 0.08, yi, txt, va="center", fontsize=6.5)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=6.5)
        ax.invert_yaxis()
        pdb = sub["pdb"].iloc[0]
        ax.set_title(f"{tgt} - {site}\n(PDB {pdb})", fontsize=8)
        ax.set_xlabel("Affinity (kcal/mol)", fontsize=7.5)
        ax.tick_params(axis="x", labelsize=7)
        # 留出右侧标签空间
        lo, hi = ax.get_xlim()
        ax.set_xlim(lo, hi + (hi - lo) * 0.28)
    # 图例
    handles = [plt.Rectangle((0, 0), 1, 1, color=ROLE_COLOR[k]) for k in
               ("positive_control", "negative_control", "candidate")]
    fig.legend(handles, [ROLE_LABEL[k] for k in
                         ("positive_control", "negative_control", "candidate")],
               loc="lower center", ncol=3, frameon=False, fontsize=8,
               bbox_to_anchor=(0.5, -0.005))
    for i in range(n, nrow * ncol):
        axes[i // ncol][i % ncol].axis("off")
    fig.suptitle("Molecular docking affinities by system "
                 "(values in parentheses: delta vs positive control)",
                 fontsize=10.5)
    fig.tight_layout(rect=[0, 0.035, 1, 0.97])
    fig.savefig(FIG / "fig_docking_affinity.png", bbox_inches="tight")
    plt.close(fig)
    print("图 2 已生成：fig_docking_affinity.png")


# ============ 图 3：姿态叠加（三视图投影） ============
def parse_pose_atoms(path):
    """读取 PDBQT 的重原子坐标。

    Vina 输出含 MODEL/ENDMDL（取 MODEL 1）；晶体配体参考文件为干净的
    ROOT…ENDROOT 刚性 PDBQT（无 MODEL 标记），则读取全部 ATOM/HETATM。
    """
    atoms = []
    in_model = False
    model_seen = False
    for line in open(path, encoding="utf-8", errors="ignore"):
        if line.startswith("MODEL"):
            model_seen = True
            in_model = (line.split()[1] == "1")
            continue
        if line.startswith("ENDMDL"):
            if in_model:
                break
            in_model = False
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        el = line[77:79].strip()
        if el in ("H", "HS"):
            continue
        # 有 MODEL 标记的文件只取 MODEL 1；否则全取
        if model_seen and not in_model:
            continue
        atoms.append((float(line[30:38]), float(line[38:46]),
                      float(line[46:54])))
    return np.array(atoms)


def kabsch_align(P, Q):
    """把 Q 叠加到 P（返回变换后的 Q）。"""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return Qc @ R + P.mean(0)


def fig_poses(top_n=2):
    df = pd.read_csv(DOCK / "docking_results.csv")
    systems = list(dict.fromkeys(zip(df["target"], df["site"])))
    for tgt, site in systems:
        sub = df[(df["target"] == tgt) & (df["site"] == site)].copy()
        pos = sub[sub["role"] == "positive_control"]
        cands = sub[sub["role"] == "candidate"].dropna(
            subset=["affinity_kcal_mol"]).sort_values("affinity_kcal_mol")
        if pos.empty or cands.empty:
            continue
        pos_tag = pos["ligand"].iloc[0]
        pdb = sub["pdb"].iloc[0]
        site_key = sani(site.split()[0])
        ref = WORK / "ligands" / f"{pos_tag}.pdbqt"
        if not ref.exists():
            continue
        A = parse_pose_atoms(ref)
        if len(A) < 5:
            continue
        picks = cands.head(top_n)
        fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.4))
        views = [("XY", 0, 1), ("XZ", 0, 2), ("YZ", 1, 2)]
        for ax, (vname, i, j) in zip(axes, views):
            ax.scatter(A[:, i], A[:, j], s=26, c="#C0392B", marker="o",
                       edgecolors="none", alpha=0.85,
                       label=f"{pos_tag.replace('POS_', '')} crystal")
            for k, (_, row) in enumerate(picks.iterrows()):
                pose = WORK / "poses" / f"{pdb}_{site_key}_{row['ligand']}.pdbqt"
                if not pose.exists():
                    continue
                B = parse_pose_atoms(pose)
                if len(B) < 5:
                    continue
                # 共晶配体（A）与候选 pose（B）均处于同一受体坐标系（Vina 输出
                # 即为受体坐标帧），故直接按原始坐标叠加即可显示候选在口袋中的
                # 相对位置；不做跨分子 Kabsch 对齐（不同分子原子数不同，对齐非法
                # 且无意义）。仅作 2D 投影，非结构叠加证据。
                ax.scatter(B[:, i], B[:, j], s=20, c=["#2471A3", "#1E8449",
                                                      "#7D3C98"][k % 3],
                           marker="^", edgecolors="none", alpha=0.8,
                           label=f"{short_name(row['ligand'])} "
                                 f"({row['affinity_kcal_mol']:.1f})")
            ax.set_title(f"{vname} projection", fontsize=8.5)
            ax.set_xlabel(f"{'xyz'[i]} (A)", fontsize=7.5)
            ax.set_ylabel(f"{'xyz'[j]} (A)", fontsize=7.5)
            ax.tick_params(labelsize=7)
            ax.set_aspect("equal")
            ax.legend(fontsize=6.2, loc="best", framealpha=0.9)
        fig.suptitle(f"{tgt} / {site} (PDB {pdb}) - top candidate poses vs "
                     f"co-crystallized ligand (receptor-frame 2D projection)",
                     fontsize=9.5)
        fig.tight_layout(rect=[0, 0, 1, 0.93])
        safe = re.sub(r"[^A-Za-z0-9]+", "_", f"{tgt}_{site}")
        fig.savefig(FIG / f"fig_docking_poses_{safe}.png", bbox_inches="tight")
        plt.close(fig)
        print(f"图 3 已生成：fig_docking_poses_{safe}.png")


if __name__ == "__main__":
    fig_rmsd()
    fig_affinity()
    fig_poses(top_n=2)
    print("DONE — docking_figures.py")
