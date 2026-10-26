# -*- coding: utf-8 -*-
"""
build_supplementary_tables.py — 生成稿件附表 S1–S11（任务 #12b）

全部附表**从产物文件直接导出**，不做任何重新计算或美化性省略；
每张表附 source 说明，列名保持产物原样以便审计复核。
S11（三层交集独立性检验）由 02_脚本/test_intersection_independence.py 计算，
本脚本仅转置导出。

输出：05_稿件/附表/S1..S11*.csv
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import glob
import json

import pandas as pd

BASE = _ROOT
OUT = f"{BASE}/05_稿件/附表"
os.makedirs(OUT, exist_ok=True)

written = []


def save(df, fname, note):
    path = f"{OUT}/{fname}"
    df.to_csv(path, index=False)
    written.append((fname, df.shape, note))
    print(f"  {fname:34s} {df.shape[0]:5d} 行 × {df.shape[1]:2d} 列  | {note}")


print("生成附表：")

# ── S1: KpyRN 基因集溯源表 ────────────────────────────────────────
kp = pd.read_csv(f"{BASE}/01_数据/基因集/KpyRN_modules.csv")
kp.insert(0, "gene_index", range(1, len(kp) + 1))
save(kp, "S1_KpyRN_gene_set.csv",
     "KpyRN 五模块基因集（140 条目，128 计分），含来源与证据说明")

# ── S2: 样本级元数据 ────────────────────────────────────────────
s2 = []
for f in sorted(glob.glob(f"{BASE}/01_数据/元数据/样本表/*_samples.csv")):
    gse = os.path.basename(f).replace("_samples.csv", "")
    d = pd.read_csv(f)
    d.insert(0, "GSE", gse)
    s2.append(d)
s2df = pd.concat(s2, ignore_index=True)
save(s2df, "S2_sample_metadata_all_GSE.csv",
     "全部 26 个 GEO 系列的样本级元数据（含未用于计算者）")

# ── S3: meta 分析结果 ────────────────────────────────────────────
m1 = pd.read_csv(f"{BASE}/03_结果/meta/META1_DFU_vs_control.csv")
m2 = pd.read_csv(f"{BASE}/03_结果/meta/META2_Healer_vs_Nonhealer.csv")
m1.insert(0, "meta_analysis", "META1_DFU_vs_control")
m2.insert(0, "meta_analysis", "META2_Healer_vs_Nonhealer")
m3 = pd.concat([m1, m2], ignore_index=True)
# 方向一致性按各自 meta 分析的两个构成队列分别判定：
#   META1 = GSE134431 vs GSE199939；META2 = GSE134431 vs GSE143735
# （旧逻辑误用 META2 的 GSE143735 去判 META1，导致 META1 行全为 NaN→False）
def _concord(df, c1, c2):
    return ((df[c1] * df[c2]) > 0).astype(bool)
cond = pd.Series(False, index=m3.index)
mask1 = m3.meta_analysis == "META1_DFU_vs_control"
mask2 = m3.meta_analysis == "META2_Healer_vs_Nonhealer"
cond[mask1] = _concord(m3.loc[mask1], "logFC_GSE134431", "logFC_GSE199939")
cond[mask2] = _concord(m3.loc[mask2], "logFC_GSE134431", "logFC_GSE143735")
m3["direction_concordant_2cohorts"] = cond
save(m3, "S3_meta_analysis_results.csv",
     "META1/META2 随机效应 meta 全表（pooled logFC/SE/z/p/FDR/I² + 队列特异）")

# ── S4: 机器学习性能 ────────────────────────────────────────────
s4 = []
for f, tag in [("TaskA_diagnosis_performance.csv", "TaskA_diagnosis"),
               ("TaskB_healing_performance.csv", "TaskB_healing")]:
    d = pd.read_csv(f"{BASE}/03_结果/ML/{f}")
    d.insert(0, "task", tag)
    s4.append(d)
s4df = pd.concat(s4, ignore_index=True)
save(s4df, "S4_ML_performance.csv",
     "Task A/B 各算法 CV AUC（10 种子）与外部队列 AUC + DeLong CI")

sel = pd.read_csv(f"{BASE}/03_结果/ML/TaskA_diagnosis_gene_selection.csv")
save(sel, "S4b_ML_gene_selection.csv",
     "26 基因候选池的算法选择频次与共识判定")

# ── S5: 免疫浸润 ────────────────────────────────────────────────
s5 = []
for f in sorted(glob.glob(f"{BASE}/03_结果/免疫浸润/*_immune_stats.csv")):
    s5.append(pd.read_csv(f))
save(pd.concat(s5, ignore_index=True), "S5a_immune_group_stats.csv",
     "三算法解卷积的组间统计（全部对比，含 FDR）")
s5b = []
for f in sorted(glob.glob(f"{BASE}/03_结果/免疫浸润/*_hubgene_immune_corr.csv")):
    s5b.append(pd.read_csv(f))
save(pd.concat(s5b, ignore_index=True), "S5b_hubgene_immune_correlation.csv",
     "枢纽基因 × 免疫指标 Spearman 相关（含 FDR）")
conc = pd.read_csv(f"{BASE}/03_结果/免疫浸润/method_concordance_epic_vs_mcp.csv")
save(conc, "S5c_method_concordance.csv", "EPIC 与 MCP-counter 跨样本一致性")

# ── S6: 单细胞与空间 ────────────────────────────────────────────
s6 = []
for f in sorted(glob.glob(f"{BASE}/03_结果/单细胞/composition_*.csv")):
    s6.append(pd.read_csv(f))
save(pd.concat(s6, ignore_index=True), "S6a_scrna_composition.csv",
     "GSE165816 受试体级细胞组成统计（足皮肤/前臂/PBMC/愈合结局/4 组）")
pb = pd.read_csv(f"{BASE}/03_结果/单细胞/pseudobulk_Foot_skin_DFU_vs_nonDFU_by_lineage.csv")
save(pb, "S6b_scrna_pseudobulk_footskin.csv",
     "GSE165816 足皮肤 pseudobulk 基因×谱系统计（median1/median2/q/effect_rb）")
p223 = pd.read_csv(f"{BASE}/03_结果/单细胞GSE223964/pseudobulk_condition_stats.csv")
save(p223, "S6c_scrna223964_pseudobulk_stats.csv",
     "GSE223964 供体级 细胞类型×模块 统计（Mann-Whitney U exact）")
perm = pd.read_csv(f"{BASE}/03_结果/单细胞GSE223964/xcheck_pseudobulk_perm.csv")
save(perm, "S6d_scrna223964_permutation_check.csv",
     "GSE223964 独立路径：供体级中位 pseudobulk + 10000 次标签置换")
roi = pd.read_csv(f"{BASE}/03_结果/空间GSE166120/roi_metadata.csv")
save(roi, "S6e_geommx_roi_metadata.csv", "GSE166120 23 个 ROI 的元数据（受试体/位置/结局）")
ms = pd.read_csv(f"{BASE}/03_结果/空间GSE166120/module_scores_by_roi.csv")
save(ms, "S6f_geommx_module_scores.csv", "GSE166120 ROI 级模块评分")
rs = pd.read_csv(f"{BASE}/03_结果/空间GSE166120/region_stats.csv")
save(rs, "S6g_geommx_region_stats.csv",
     "GSE166120 区域统计（受试体级为主口径 + ROI 级敏感性）")

# ── S7: 感染队列 ────────────────────────────────────────────────
s7 = []
for f in ["deg_main.csv", "deg_paired.csv", "deg_after8.csv"]:
    s7.append(pd.read_csv(f"{BASE}/03_结果/感染队列GSE230426/{f}"))
save(pd.concat(s7, ignore_index=True), "S7a_gse230426_deg_all_models.csv",
     "GSE230426 三模型全基因组 DEG（各 57,059 基因，全阴性）")
cam = pd.read_csv(f"{BASE}/03_结果/感染队列GSE230426/camera_modules.csv")
save(cam, "S7b_gse230426_camera_modules.csv", "GSE230426 camera 模块检验（三模型 × 四模块）")
hg = pd.read_csv(f"{BASE}/03_结果/感染队列GSE230426/hub_genes.csv")
save(hg, "S7c_gse230426_hub_genes.csv", "GSE230426 六个枢纽基因三模型结果（全不显著）")

# ── S8: MR ──────────────────────────────────────────────────────
mr = pd.read_csv(f"{BASE}/03_结果/MR/mr_results.csv")
save(mr, "S8a_MR_results.csv", "MR 全部 34 对 gene×outcome 结果（含 BH-FDR）")
inst = pd.read_csv(f"{BASE}/03_结果/MR/mr_instruments.csv")
save(inst, "S8b_MR_instruments.csv", "cis-eQTL 工具变量与 F 统计")
harm = pd.read_csv(f"{BASE}/03_结果/MR/mr_harmonization_audit.csv")
save(harm, "S8c_MR_harmonization_audit.csv", "等位基因协调审计（含全部排除记录）")

# ── S9: 分子对接 ────────────────────────────────────────────────
dk = pd.read_csv(f"{BASE}/03_结果/药物预测/docking/docking_results.csv")
save(dk, "S9a_docking_results.csv",
     "全部 55 次对接（target/site/PDB/ligand/role/affinity/Δ/box）")
rmsd = pd.read_csv(f"{BASE}/03_结果/药物预测/docking/positive_control_rmsd.csv")
save(rmsd, "S9b_positive_control_rmsd.csv",
     "阳性对照重对接 RMSD 验证（8/8 = 0.00 Å，验收阈值 < 2 Å）")

# ── S10: RT-qPCR 设计与验证 ─────────────────────────────────────
qp = pd.read_csv(f"{BASE}/03_结果/实验验证设计/qpcr_primers.csv")
save(qp, "S10a_qpcr_primers.csv", "12 对引物（6 目标 + 6 内参候选）")
pv = pd.read_csv(f"{BASE}/03_结果/实验验证设计/primer_verification.csv")
save(pv, "S10b_primer_verification.csv", "引物二级结构/自互补/GC/Tm 独立复算")
ss = pd.read_csv(f"{BASE}/03_结果/实验验证设计/sample_size_calculation.csv")
save(ss, "S10c_sample_size_calculation.csv",
     "基于真实效应量的样本量计算（配对/独立/保守情景）")
# S10d: BLAST 特异性（优先 v3 = core_nt 最终版；回退 v2 仅供历史比对）
# 修复记录：旧版写 `f"{OUT}/../primer_blast_specificity_v2.csv"`，OUT 为
# `05_稿件/附表`，`..` 解析到 `05_稿件/`，而实际文件在 `03_结果/实验验证设计/`
# → 路径不存在 → S10d 被**静默跳过**（附表索引中无 S10d）。现改为绝对路径，
# 且文件缺失时**显式报错**，不再静默跳过。
_bf_cands = [
    f"{BASE}/03_结果/实验验证设计/primer_blast_specificity_v3.csv",
    f"{BASE}/03_结果/实验验证设计/primer_blast_specificity_v2.csv",
]
_bf = next((p for p in _bf_cands if os.path.exists(p)), None)
if _bf is None:
    raise FileNotFoundError(
        "S10d 需要 BLAST 特异性结果文件，但以下候选均不存在：\n  "
        + "\n  ".join(_bf_cands)
        + "\n请先运行 02_脚本/verify_primers_blast_v3.py 生成。")
b = pd.read_csv(_bf)
# v3 与 v2 的列名不同，按版本选择对应列（缺失列自动剔除）
if "primer_blast_specificity_v3" in _bf:
    _cols_pref = [
        "gene", "role",
        "forward_target_human_hits", "reverse_target_human_hits",
        "forward_target_max_identity", "reverse_target_max_identity",
        "forward_target_top_acc", "reverse_target_top_acc",
        "forward_n_human_refseq_tx", "reverse_n_human_refseq_tx",
        "forward_n_offtarget_refseq_tx", "reverse_n_offtarget_refseq_tx",
        "forward_offtarget_tx_genes", "reverse_offtarget_tx_genes",
        "forward_n_offtarget_genomic", "reverse_n_offtarget_genomic",
        "forward_paralog_genes", "reverse_paralog_genes",
        "forward_pseudogene_genes", "reverse_pseudogene_genes",
        "forward_human_hit_genes", "reverse_human_hit_genes",
        "offtarget_refseq_tx_total", "offtarget_genomic_total",
        "paralog_total", "specific",
    ]
    _note = ("BLAST 单靶特异性 v3（core_nt 数据库，XML 完整描述；"
             "真脱靶=人源 RefSeq 转录本命中，同源旁系同源/假基因/基因组层面命中分层单列）")
else:
    _cols_pref = [
        "gene", "role",
        "forward_human_hit_genes", "reverse_human_hit_genes",
        "forward_random_offtarget_genes", "reverse_random_offtarget_genes",
        "forward_paralog_genes", "reverse_paralog_genes",
        "forward_pseudogene_genes", "reverse_pseudogene_genes",
        "forward_target_max_identity", "reverse_target_max_identity",
        "random_offtarget_total", "paralog_total", "specific",
    ]
    _note = ("BLAST 单靶特异性 v2（refseq_rna 数据库；真脱靶/旁系同源/假基因三层判定）")
cols = [c for c in _cols_pref if c in b.columns]
save(b[cols], "S10d_primer_blast_specificity.csv", _note)

# ── S11: 三层交集独立性检验（全为糖酵解酶是否为独立证据）─────────────
# 稿件 §3.2 与 §4.1 引用了本表。产物由 02_脚本/test_intersection_independence.py
# 生成；此处仅做**转置导出**，不重算，保证"表内数字 = 脚本输出"。
_iij = f"{BASE}/03_结果/WGCNA/intersection_independence.json"
if os.path.exists(_iij):
    with open(_iij, encoding="utf-8") as _fh:
        _d = json.load(_fh)
    _rows = []
    for _name, _v in (_d.get("T2_hypergeometric") or {}).items():
        _rows.append({
            "subset": _name,
            "total_genes": _v.get("total"),
            "genes_in_wgcna_background": _v.get("in_universe"),
            "observed_in_DFU_pool": _v.get("observed_in_pool"),
            "expected": _v.get("expected"),
            "obs_over_exp": _v.get("ratio_obs_exp"),
            "direction": _v.get("direction", "depleted" if _name == "nonM2_combined" else ""),
            "p_value": _v.get("p_two_sided_lesser_tail", _v.get("p_depleted")),
        })
    _s11 = pd.DataFrame(_rows)
    _note11 = ("三层交集独立性检验（10000 次置换 + 超几何精确检验，seed 20261025）。"
               "抽样数一律取「该子集 ∩ WGCNA 5000 高变基因背景」的大小；"
               "pool = META1-padj<0.05 ∩ DFU 显著相关模块"
               "（magenta/yellow/red/blue）共 "
               f"{_d['layers']['pool_S_size']} 基因。")
    save(_s11, "S11_intersection_independence.csv", _note11)
else:
    print(f"  [SKIP] S11: {_iij} 不存在 —— 请先运行 "
          f"02_脚本/test_intersection_independence.py")

# ── 附表索引 ────────────────────────────────────────────────────
idx = pd.DataFrame([
    {"table": w[0], "n_rows": w[1][0], "n_cols": w[1][1], "content": w[2]}
    for w in written])
idx.to_csv(f"{OUT}/附表索引.csv", index=False)

with open(f"{OUT}/附表说明.md", "w", encoding="utf-8") as fh:
    fh.write("# 附表索引与溯源说明\n\n")
    fh.write("本目录所有附表均由 `02_脚本/build_supplementary_tables.py` "
             "**从产物文件直接导出**，未做任何重新计算或美化性省略；"
             "列名保持产物原样以便审计复核。\n\n")
    fh.write("| 附表 | 行×列 | 内容 |\n|---|---|---|\n")
    for w in written:
        fh.write(f"| {w[0]} | {w[1][0]}×{w[1][1]} | {w[2]} |\n")
    fh.write("\n所有附表对应的主图与正文引用见 `05_稿件/manuscript_en.md` "
             "的 Figure legends 与 Supplementary material 章节。\n")

print(f"\n共 {len(written)} 张附表 -> 05_稿件/附表/")
print("DONE — build_supplementary_tables.py")
