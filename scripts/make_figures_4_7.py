# -*- coding: utf-8 -*-
"""
补齐稿件 Figure 4 与 Figure 7（正文引用了这两张图但仓库内无对应图件）。

Figure 4: META1 中显著的 KpyRN 基因的合并 log2FC 森林图（按模块着色，
          叠加两个队列各自的 logFC 以显示方向一致性）。
Figure 7: GSE165816 单细胞四面板（组成 / 糖酵解-IFN-I 热图 / 角质细胞 / 
          供体级推断示意）。

数据源全部为项目产物文件，脚本内含"数值来源自检"。
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "03_结果")
TAB = os.path.join(BASE, "05_稿件", "附表")
FIG = os.path.join(BASE, "04_图表")
os.makedirs(FIG, exist_ok=True)

MOD_COLOR = {
    "M1_Kpy_machinery": "#7F77DD",
    "M2_Pyruvate_supply": "#D85A30",
    "M3_Kpy_targets": "#1D9E75",
    "M4_IFNI_readout": "#378ADD",
    "M5_Histone_carriers": "#888780",
}


def fig4():
    meta = {}
    with open(os.path.join(TAB, "S3_meta_analysis_results.csv"), encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            if r["meta_analysis"] == "META1_DFU_vs_control":
                meta[r["gene"]] = r
    mod = {}
    with open(os.path.join(TAB, "S1_KpyRN_gene_set.csv"), encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            mod[r["gene_symbol"]] = r["module"]

    sig = [(g, r) for g, r in meta.items()
           if r["padj"] and float(r["padj"]) < 0.05 and g in mod]
    sig.sort(key=lambda t: float(t[1]["pooled_logFC"]))
    genes = [g for g, _ in sig]
    est = [float(r["pooled_logFC"]) for _, r in sig]
    lo = [float(r["pooled_logFC"]) - 1.96 * float(r["pooled_se"]) for _, r in sig]
    hi = [float(r["pooled_logFC"]) + 1.96 * float(r["pooled_se"]) for _, r in sig]
    c1 = [float(r["logFC_GSE134431"]) for _, r in sig]
    c2 = [float(r["logFC_GSE199939"]) for _, r in sig]

    # 自检：26 个显著基因、26/26 方向一致
    assert len(genes) == 26, f"显著 KpyRN 基因数应为 26，实际 {len(genes)}"
    concord = sum(1 for a, b in zip(c1, c2) if (a > 0) == (b > 0))
    assert concord == 26, f"方向一致应为 26/26，实际 {concord}/26"

    y = np.arange(len(genes))
    fig, ax = plt.subplots(figsize=(7.2, 8.4))
    for i, g in enumerate(genes):
        col = MOD_COLOR.get(mod[g], "#888780")
        ax.plot([lo[i], hi[i]], [i, i], color=col, lw=1.4, alpha=0.55)
        ax.plot(est[i], i, "o", color=col, ms=6.5, zorder=3)
        ax.plot(c1[i], i, marker="s", color="none", mec=col, mew=1.0, ms=4.5,
                alpha=0.85, zorder=3)
        ax.plot(c2[i], i, marker="^", color="none", mec=col, mew=1.0, ms=4.5,
                alpha=0.85, zorder=3)
    ax.axvline(0, color="#5F5E5A", lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(genes, fontsize=9)
    ax.set_xlabel("Pooled log2 fold change (DFU versus non-ulcer skin), "
                  "random effects", fontsize=9.5)
    ax.set_title("Figure 4. Meta-analytic effect sizes of the scored KpyRN genes\n"
                 "26 genes, FDR < 0.05; 95% CI from pooled SE",
                 fontsize=10.5, loc="left")
    handles = [plt.Line2D([], [], marker="o", ls="", color=c,
                          label=m.split("_", 1)[1].replace("_", " "))
               for m, c in MOD_COLOR.items() if m in {mod[g] for g in genes}]
    handles += [plt.Line2D([], [], marker="s", ls="", color="#444441", mfc="none",
                           label="GSE134431"),
                plt.Line2D([], [], marker="^", ls="", color="#444441", mfc="none",
                           label="GSE199939")]
    ax.legend(handles=handles, fontsize=8, loc="lower right", frameon=False)
    ax.tick_params(axis="x", labelsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    out = os.path.join(FIG, "fig4_kpyrn_meta_forest.png")
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[Fig 4] {out} | 26 genes, concordance {concord}/26 (自检通过)")


def fig7():
    comp = []
    with open(os.path.join(RES, "单细胞", "composition_Foot_skin_DFU_vs_nonDFU.csv"),
              encoding="utf-8-sig") as fh:
        comp = list(csv.DictReader(fh))
    pb = []
    with open(os.path.join(RES, "单细胞",
                           "pseudobulk_Foot_skin_DFU_vs_nonDFU_by_lineage.csv"),
              encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            pb.append(r)
    # median1 = DFU (n1), median2 = non-DFU (n2)  —— 与产物列语义一致
    comp_sorted = sorted(comp, key=lambda r: -float(r["median1"]))
    labels = [r["lineage"] for r in comp_sorted]
    m1 = [float(r["median1"]) for r in comp_sorted]
    m2 = [float(r["median2"]) for r in comp_sorted]
    qs = [float(r["q"]) for r in comp_sorted]
    sig_idx = [i for i, q in enumerate(qs) if q < 0.05]
    assert len(sig_idx) == 4, f"显著谱系应为 4 个，实际 {len(sig_idx)}"

    glyco = ["PGAM1", "TPI1", "LDHA", "ENO1", "PKM", "GAPDH", "HK2", "PGM2"]
    ifn = ["IFI44L", "IFI27", "DHX58", "DDX58", "STAT2", "ISG20", "GBP1"]
    feats = glyco + ifn
    lineages = ["Keratinocytes", "Fibroblasts", "Pericyte_SMC", "Mono_CD14",
                "T_CD8", "Endothelial"]
    lookup = {}
    for r in pb:
        lookup[(r["lineage"], r["feature"])] = (float(r["median1"]),
                                                float(r["median2"]))
    mat = []
    for f in feats:
        row = []
        for lg in lineages:
            v = lookup.get((lg, f))
            row.append((v[0] - v[1]) if v else np.nan)
        mat.append(row)
    mat = np.array(mat, dtype=float)

    kert = []
    for g in glyco:
        v = lookup.get(("Keratinocytes", g))
        if v:
            kert.append((g, v[0], v[1]))

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.0))

    ax = axes[0, 0]
    yy = np.arange(len(labels))
    ax.barh(yy, m1, color="#D85A30", height=0.42, label="DFU")
    ax.barh(yy + 0.42, m2, color="#378ADD", height=0.42, label="non-DFU")
    for i in sig_idx:
        ax.text(max(m1[i], m2[i]) * 1.05, i + 0.2, "*", fontsize=13,
                color="#A32D2D", va="center")
    ax.set_yticks(yy + 0.21)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xscale("symlog", linthresh=1e-4)
    ax.set_xlabel("Median proportion (donor level, symlog)", fontsize=9)
    ax.set_title("A. Foot-skin composition, DFU vs non-DFU\n"
                 "* FDR < 0.05 (donor-level Wilcoxon)", fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[0, 1]
    vmax = np.nanmax(np.abs(mat))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(lineages)))
    ax.set_xticklabels(lineages, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(feats)))
    ax.set_yticklabels(feats, fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if not np.isnan(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:+.2f}", ha="center", va="center",
                        fontsize=6.5, color="#2C2C2A")
    ax.axhline(len(glyco) - 0.5, color="#2C2C2A", lw=0.9)
    ax.set_title("B. Pseudobulk log2 CPM difference (DFU − non-DFU)\n"
                 "glycolysis (upper) and IFN-I read-out (lower)",
                 fontsize=10, loc="left")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)

    ax = axes[1, 0]
    kn = [k[0] for k in kert]
    kv1 = [k[1] for k in kert]
    kv2 = [k[2] for k in kert]
    x = np.arange(len(kn))
    ax.bar(x - 0.19, kv1, width=0.38, color="#D85A30", label="DFU")
    ax.bar(x + 0.19, kv2, width=0.38, color="#378ADD", label="non-DFU")
    ax.set_xticks(x)
    ax.set_xticklabels(kn, rotation=35, ha="right", fontsize=8.5)
    ax.set_ylabel("log2 CPM (pseudobulk, median)", fontsize=9)
    ax.set_title("C. Keratinocyte glycolytic enzymes", fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1, 1]
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    box = dict(boxstyle="round,pad=0.35", fc="#F1EFE8", ec="#888780", lw=0.8)
    ax.text(0.5, 0.86, "D. Inferential unit", fontsize=10, ha="center",
            weight="bold")
    ax.text(0.5, 0.62,
            "cells (n = 167,837)\n   nested within\n"
            "donors (14 DFU / 19 non-DFU)\n        ↓\n"
            "all tests at donor level",
            fontsize=9.5, ha="center", va="center", bbox=box)
    ax.text(0.5, 0.16,
            "Cell-level statistics are reported as descriptive only\n"
            "(pseudoreplication); no inference is drawn from them.",
            fontsize=8.5, ha="center", va="center", color="#5F5E5A")

    fig.suptitle("Figure 7. Single-cell localisation of the KpyRN signal "
                 "in GSE165816", fontsize=11.5, x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    out = os.path.join(FIG, "fig7_scrna_gse165816.png")
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[Fig 7] {out} | 显著谱系 {len(sig_idx)} 个、热图 {mat.shape} （自检通过）")


if __name__ == "__main__":
    fig4()
    fig7()
