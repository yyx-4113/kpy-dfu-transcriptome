# -*- coding: utf-8 -*-
"""
assemble_submission_figures.py — 组装 BMC 投稿要求的"单文件复合图"（Fig 5/6/8/9/10）

背景：BMC/Springer Nature 投稿指南要求多面板图（parts a/b/c/d）必须以单一复合
文件提交（"Multi-panel figures ... should be submitted as a single composite file
that contains all parts of the figure"）。本稿 Fig 1–4、7 已有复合文件（04_图表/）；
Fig 5、6、8、10 此前仅以零散组件 PNG 存在于 03_结果/ 各子目录，Fig 9 为单面板图。

本脚本按稿件 Figure legends（manuscript_en.md，Fig 5/6/8/10 图注）的面板映射组装：
  fig5_ml_performance.png       A=重复5折CV AUC  B=GSE199939外验ROC / C=共识hub基因表达箱线
  fig6_immune_infiltration.png  A=免疫比例箱线  B=hub×免疫Spearman热图 / C=EPIC vs MCP一致性散点
  fig8_negative_findings.png    A=GSE223964模块得分 + MWU/置换P值矩阵
                                B=GSE166120 CTA覆盖 + hub基因覆盖 / C=ROI vs subject P值
                                D=GSE230426 camera模块 + hub基因阴性结果
  fig9_mr_forest.png            单面板MR森林图（重存300dpi元数据）
  fig10_docking.png             A=正控RMSD  B=亲和力总览 / C=四个最优系统构象2×2

纪律（与项目审计规范一致）：
- 组件→面板映射按稿件图注声明为常量，与 manuscript_en.md 图注一一对应；
- Fig 8A 的 P 值矩阵由 S6c/S6d 产物 CSV **动态渲染**，不硬编码任何数字；
- Fig 10C 的四个系统由 S9a 按"各系统最佳候选亲和力"**动态排序**取前四；
- 拼接为 1:1 像素白底填充，不做重采样、不改动任何数据内容；
- 全部输出 PNG 带 300 dpi 元数据（BMC 指南 ~300 dpi），逐文件校验 <10 MB。
"""
import os
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "03_结果")
TAB = os.path.join(BASE, "05_稿件", "附表")
FIG = os.path.join(BASE, "04_图表")

INK = (27, 38, 49)  # #1B2631，与既有主图标题色一致

_FONT_PATH = os.path.join(matplotlib.get_data_path(), "fonts", "ttf",
                          "DejaVuSans-Bold.ttf")
if not os.path.exists(_FONT_PATH):
    import matplotlib.font_manager as fm
    _FONT_PATH = fm.findfont(fm.FontProperties(family="DejaVu Sans", weight="bold"))
_FONT = ImageFont.truetype(_FONT_PATH, 52)


# ---------------------------------------------------------------- 基础拼接工具
def load(path):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return Image.open(path).convert("RGB")


def label_band(width, text):
    band = Image.new("RGB", (width, 80), "white")
    d = ImageDraw.Draw(band)
    d.text((16, 8), text, font=_FONT, fill=INK)
    return band


def with_label(img, text):
    return vstack([label_band(img.width, text), img])


def hstack(imgs, bg="white"):
    h = max(i.height for i in imgs)
    w = sum(i.width for i in imgs)
    out = Image.new("RGB", (w, h), bg)
    x = 0
    for i in imgs:
        out.paste(i, (x, 0))
        x += i.width
    return out


def vstack(imgs, bg="white"):
    w = max(i.width for i in imgs)
    h = sum(i.height for i in imgs)
    out = Image.new("RGB", (w, h), bg)
    y = 0
    for i in imgs:
        out.paste(i, ((w - i.width) // 2, y))
        y += i.height
    return out


def save300(img, path):
    img.save(path, dpi=(300, 300))
    kb = os.path.getsize(path) / 1024
    print(f"  -> {os.path.basename(path)}  {img.width}x{img.height}px  {kb:.0f} KB")
    assert kb < 10240, f"{path} 超过 BMC 单文件 10 MB 上限"


# ---------------------------------------------------- Fig 8A：P 值矩阵（动态）
def render_pmatrix(out_png, target_width):
    """由 S6c（Mann–Whitney U）/ S6d（置换检验）动态渲染供体级 P 值矩阵。"""
    c = pd.read_csv(os.path.join(TAB, "S6c_scrna223964_pseudobulk_stats.csv"))
    d = pd.read_csv(os.path.join(TAB, "S6d_scrna223964_permutation_check.csv"))
    mods = ["M1_Kpy_machinery", "M2_Pyruvate_supply",
            "M3_Kpy_targets", "M4_IFNI_readout"]
    mlab = ["M1\nKpy machinery", "M2\nPyruvate supply",
            "M3\nKpy targets", "M4\nIFN-I readout"]
    cts = list(dict.fromkeys(c["celltype"]))
    cmap = {(r["celltype"], r["module"]): r for _, r in c.iterrows()}
    dmap = {(r["celltype"], r["module"]): r for _, r in d.iterrows()}
    missing = [k for k in ((ct, m) for ct in cts for m in mods)
               if k not in cmap or k not in dmap]
    assert not missing, f"S6c/S6d 覆盖缺口: {missing}"

    q_mwu_min = float(c["mwu_p_bf"].min())
    q_perm_min = float(d["perm_p_bf"].min())
    assert q_mwu_min >= 0.05 and q_perm_min >= 0.05, (
        "存在 Bonferroni 后显著对比，与稿件 §3.6 阴性结论不符，需人工复核")

    nrow = len(cts)
    figw = 11.0
    figh = 0.42 * (nrow + 3) + 1.9
    fig, ax = plt.subplots(figsize=(figw, figh))
    ax.set_xlim(0, 9.0)
    ax.set_ylim(0, nrow + 3.4)
    ax.invert_yaxis()
    ax.axis("off")

    x_lab, x0, colw = 0.02, 2.05, 0.87

    def cell_xy(j):
        return x0 + j * colw

    # 表头
    for m_i, m in enumerate(mlab):
        ax.text(cell_xy(2 * m_i) + colw * 0.9, 0.55, m, ha="center", va="center",
                fontsize=8.2, fontweight="bold", color="#1B2631")
    for m_i in range(4):
        ax.text(cell_xy(2 * m_i) + colw * 0.5, 1.28, "MWU p", ha="center",
                va="center", fontsize=7.0, color="#5F5E5A")
        ax.text(cell_xy(2 * m_i + 1) + colw * 0.5, 1.28, "perm p", ha="center",
                va="center", fontsize=7.0, color="#5F5E5A")

    for r_i, ct in enumerate(cts):
        y = 2.05 + r_i
        n_d = int(cmap[(ct, mods[0])]["n_diabetic"])
        n_n = int(cmap[(ct, mods[0])]["n_nondiabetic"])
        ax.text(x_lab, y + 0.35, f"{ct}  ({n_d} v {n_n})", ha="left", va="center",
                fontsize=7.6, color="#1B2631")
        for m_i, m in enumerate(mods):
            p1 = float(cmap[(ct, m)]["mwu_p"])
            p2 = float(dmap[(ct, m)]["perm_p"])
            for j, p in enumerate((p1, p2)):
                x = cell_xy(2 * m_i + j)
                if p < 0.05:
                    ax.add_patch(plt.Rectangle((x, y), colw, 0.7,
                                               fc="#FADBD8", ec="none", zorder=1))
                ax.text(x + colw * 0.5, y + 0.35, f"{p:.3f}", ha="center",
                        va="center", fontsize=7.0,
                        color="#A32D2D" if p < 0.05 else "#2C2C2A", zorder=2)

    y_foot = 2.05 + nrow + 0.55
    fig.text(0.06, 1 - (y_foot + 0.75) / (nrow + 3.4),
             "Donor-level tests of KpyRN module scores, GSE223964 (diabetic vs non-diabetic skin). "
             "MWU: exact Mann–Whitney U on mean pseudobulk (path A); "
             "perm: 10,000-label permutation on median pseudobulk (sensitivity analysis).",
             fontsize=7.2, color="#5F5E5A", va="top")
    fig.text(0.06, 1 - (y_foot + 1.45) / (nrow + 3.4),
             f"Red cells: nominal p < 0.05. After Bonferroni adjustment no contrast is "
             f"significant: smallest MWU q = {q_mwu_min:.2f}, smallest permutation q = {q_perm_min:.2f}.",
             fontsize=7.2, color="#5F5E5A", va="top")

    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    fig.savefig(out_png, dpi=int(target_width / figw))
    plt.close(fig)


# ----------------------------------------------------------------- 复合图组装
def fig5():
    print("[Fig 5] ML performance / external validation (A|B / C)")
    A = with_label(load(os.path.join(FIG, "ML_performance_overview.png")), "A")
    B = with_label(load(os.path.join(FIG, "ML_ROC_external_GSE199939.png")), "B")
    C = with_label(load(os.path.join(FIG, "ML_consensus_gene_boxplots.png")), "C")
    save300(vstack([hstack([A, B]), C]),
            os.path.join(FIG, "fig5_ml_performance.png"))


def fig6():
    print("[Fig 6] Immune infiltration (A|B / C)")
    imp = os.path.join(RES, "免疫浸润", "figures")
    A = with_label(load(os.path.join(imp, "fig_immune_group_boxplots.png")), "A")
    B = with_label(load(os.path.join(imp, "fig_hubgene_immune_heatmap.png")), "B")
    C = with_label(load(os.path.join(imp, "fig_method_concordance.png")), "C")
    save300(vstack([hstack([A, B]), C]),
            os.path.join(FIG, "fig6_immune_infiltration.png"))


def fig8():
    print("[Fig 8] Negative and boundary findings (A|B / C|D)")
    scores = load(os.path.join(RES, "单细胞GSE223964", "fig_scores.png"))
    pm = os.path.join(FIG, "_tmp_fig8a_pmatrix.png")
    render_pmatrix(pm, target_width=scores.width)
    A = with_label(vstack([scores, load(pm)]), "A")
    os.remove(pm)
    B = with_label(hstack([
        load(os.path.join(RES, "空间GSE166120", "figures",
                          "fig_geommx_board_coverage.png")),
        load(os.path.join(RES, "空间GSE166120", "figures",
                          "fig_geommx_hubgenes.png"))]), "B")
    C = with_label(load(os.path.join(RES, "空间GSE166120", "figures",
                                     "fig_geommx_pseudoreplication_check.png")), "C")
    D = with_label(vstack([
        load(os.path.join(RES, "感染队列GSE230426", "figures",
                          "fig_gse230426_camera_modules.png")),
        load(os.path.join(RES, "感染队列GSE230426", "figures",
                          "fig_gse230426_hub_genes.png"))]), "D")
    save300(vstack([hstack([A, B]), hstack([C, D])]),
            os.path.join(FIG, "fig8_negative_findings.png"))


def fig9():
    print("[Fig 9] MR forest (single panel, re-saved with 300 dpi metadata)")
    src = os.path.join(RES, "MR", "figures", "fig_mr_forest.png")
    save300(load(src), os.path.join(FIG, "fig9_mr_forest.png"))


POSE_FILE = {
    ("PKM2", "4JPG", "C-pocket (allosteric)"):
        "fig_docking_poses_PKM2_C_pocket_allosteric_.png",
    ("PKM2", "4JPG", "A-pocket (FBP site)"):
        "fig_docking_poses_PKM2_A_pocket_FBP_site_.png",
    ("EP300", "4PZS", "acetyl-CoA substrate site"):
        "fig_docking_poses_EP300_acetyl_CoA_substrate_site.png",
    ("SIRT3", "4BN4", "ADP-ribose/NAD+ C-pocket"):
        "fig_docking_poses_SIRT3_ADP_ribose_NAD_C_pocket.png",
    ("LDHA", "4JNK", "active site"):
        "fig_docking_poses_LDHA_active_site.png",
    ("PDK2", "4MPC", "ATP pocket"):
        "fig_docking_poses_PDK2_ATP_pocket.png",
    ("PDK2", "2BU8", "pyruvate site (R-domain)"):
        "fig_docking_poses_PDK2_pyruvate_site_R_domain_.png",
    ("PDK2", "2BU5", "R-domain N-terminal site"):
        "fig_docking_poses_PDK2_R_domain_N_terminal_site.png",
}


def fig10():
    print("[Fig 10] Docking (A|B / C: 4 best systems by best candidate affinity)")
    s9 = pd.read_csv(os.path.join(TAB, "S9a_docking_results.csv"))
    cand = s9[s9["role"] == "candidate"]
    best = (cand.groupby(["target", "pdb", "site"], as_index=False)
                ["affinity_kcal_mol"].min()
                .sort_values("affinity_kcal_mol").head(4))
    print("  最优四系统（按最佳候选亲和力，动态计算）:")
    poses = []
    for r in best.itertuples():
        key = (r.target, r.pdb, r.site)
        fn = POSE_FILE[key]
        print(f"    {r.target} {r.pdb} {r.site}: {r.affinity_kcal_mol:.2f} kcal/mol -> {fn}")
        poses.append(load(os.path.join(RES, "药物预测", "docking",
                                       "figures", fn)))
    A = with_label(load(os.path.join(RES, "药物预测", "docking", "figures",
                                     "fig_docking_rmsd.png")), "A")
    B = with_label(load(os.path.join(RES, "药物预测", "docking", "figures",
                                     "fig_docking_affinity.png")), "B")
    C = with_label(vstack([hstack([poses[0], poses[1]]),
                           hstack([poses[2], poses[3]])]), "C")
    save300(vstack([hstack([A, B]), C]),
            os.path.join(FIG, "fig10_docking.png"))


def main():
    fig5()
    fig6()
    fig8()
    fig9()
    fig10()
    expected = ["fig1_workflow.png", "fig2_three_layer_intersection.png",
                "fig3_wgcna_module_trait.png", "fig4_kpyrn_meta_forest.png",
                "fig5_ml_performance.png", "fig6_immune_infiltration.png",
                "fig7_scrna_gse165816.png", "fig8_negative_findings.png",
                "fig9_mr_forest.png", "fig10_docking.png"]
    missing = [f for f in expected if not os.path.exists(os.path.join(FIG, f))]
    assert not missing, f"主图仍缺失: {missing}"
    print(f"DONE — 10/10 主图复合文件齐备于 04_图表/")


if __name__ == "__main__":
    main()
