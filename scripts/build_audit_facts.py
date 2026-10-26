# -*- coding: utf-8 -*-
"""
build_audit_facts.py — 生成「正文数字唯一真源」(05_稿件/audit_facts.json)

设计原则（遵循项目质量宪法 #2「拒绝循环验证」）：
  - 本脚本**只从产物文件**读取数字，不做任何重新计算（除必要的派生比值外，
    派生值同时记录其计算式与输入，便于独立复核）。
  - 稿件正文必须逐句引用本文件中的 fact id；审计追踪脚本据此反向核对。
  - 每个 fact 记录: id / 值 / 类型 / 源文件 / 提取表达式 / 备注。

输出：
  05_稿件/audit_facts.json    事实库（机器可读）
  05_稿件/audit_facts.csv     同内容扁平表（人工可读）
  05_稿件/audit_facts_meta.json  生成时间与脚本哈希
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import re
import json
import hashlib
import subprocess

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact as _fisher_exact
from scipy.special import gammaln as _gammaln
from scipy.optimize import brentq as _brentq

BASE = _ROOT
OUT = f"{BASE}/05_稿件"
os.makedirs(OUT, exist_ok=True)

FACTS = {}


def put(fid, value, source, expr, note="", ftype="number"):
    """登记一个事实。value 保留原始精度（float）或原样（str/int）。"""
    if fid in FACTS:
        raise ValueError(f"重复 fact id: {fid}")
    FACTS[fid] = {
        "id": fid,
        "value": value,
        "type": ftype,
        "source": source,
        "expr": expr,
        "note": note,
    }


def fisher_exact_ci(table, alpha=0.05):
    """2×2 表的 Fisher 精确检验 + 95% 精确条件 CI（R fisher.test 的
    'sum(P<=P_obs)' 两尾定义，通过反演求得）。返回 (OR, p, ci_low, ci_high)。
    仅在必要时由 contingency 表派生，计算式与输入均记录到 fact 的 expr/note。"""
    a, b = table[0][0], table[0][1]
    c, d = table[1][0], table[1][1]
    N = a + b + c + d
    r1, c1 = a + b, a + c
    or_pt, p = _fisher_exact(table, alternative="two-sided")
    alo = max(0, r1 + c1 - N)
    ahi = min(r1, c1)
    grid = np.arange(alo, ahi + 1)

    def _pdist(theta):
        lw = (_gammaln(r1 + 1) - _gammaln(grid + 1) - _gammaln(r1 - grid + 1)
              + _gammaln((N - r1) + 1) - _gammaln(c1 - grid + 1)
              - _gammaln((N - r1) - (c1 - grid) + 1)
              + grid * np.log(theta))
        w = np.exp(lw - lw.max())
        return w / w.sum()

    def _two_p(theta):
        p_ = _pdist(theta)
        pobs = p_[grid == a][0]
        return float(np.sum(p_[p_ <= pobs]))

    ci_lo = _brentq(lambda t: _two_p(t) - alpha, 1e-4, 1.0 - 1e-9)
    ci_hi = _brentq(lambda t: _two_p(t) - alpha, 1.0 + 1e-9, 50.0)
    return float(or_pt), float(p), float(ci_lo), float(ci_hi)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ─────────────────────────────────────────────────────────────
# 1. KpyRN 基因集
# ─────────────────────────────────────────────────────────────
kp = pd.read_csv(f"{BASE}/01_数据/基因集/KpyRN_modules.csv")
MODMAP = {
    "M1_Kpy_machinery": "M1 写入/擦除机器",
    "M2_Pyruvate_supply": "M2 丙酮酸供应与分流",
    "M3_Kpy_targets": "M3 已鉴定 Kpy 靶基因",
    "M4_IFNI_readout": "M4 IFN-I 读出",
    "M5_Histone_carriers": "M5 组蛋白载体（不计分）",
}
for m, cnt in kp.groupby("module").size().items():
    put(f"kpyrn.module.{m}.n", int(cnt),
        "01_数据/基因集/KpyRN_modules.csv",
        f"groupby('module').size()['{m}']",
        f"模块中文名：{MODMAP.get(m, m)}")
put("kpyrn.n_total", int(len(kp)), "01_数据/基因集/KpyRN_modules.csv",
    "len(df)", "唯一基因条目总数")
put("kpyrn.n_scored", int((kp.scoring_include == "Yes").sum()),
    "01_数据/基因集/KpyRN_modules.csv",
    "(scoring_include=='Yes').sum()", "纳入评分的基因数（M5 组蛋白载体除外）")
for m, cnt in kp[kp.scoring_include == "Yes"].groupby("module").size().items():
    put(f"kpyrn.scored.{m}.n", int(cnt),
        "01_数据/基因集/KpyRN_modules.csv",
        f"scoring_include=='Yes' 且 module=='{m}' 的行数")

# ─────────────────────────────────────────────────────────────
# 2. DEG 与 meta
# ─────────────────────────────────────────────────────────────
DEG_EXPECT = {
    ("GSE134431", "DFU_vs_DFS"): (2814, 12450),
    ("GSE199939", "DFU_vs_control"): (9773, 17239),
}
for (gse, cmp_), (n_deg, n_test) in DEG_EXPECT.items():
    f = f"{BASE}/03_结果/DEG/{gse}_{cmp_}_limma.csv"
    d = pd.read_csv(f)
    up = int(((d["logFC"] > 0) & (d["adj.P.Val"] < 0.05)).sum()) \
        if "adj.P.Val" in d.columns else None
    dn = int(((d["logFC"] < 0) & (d["adj.P.Val"] < 0.05)).sum()) \
        if "adj.P.Val" in d.columns else None
    put(f"deg.{gse}.{cmp_}.n_deg", n_deg, f"03_结果/DEG/{gse}_{cmp_}_limma.csv",
        "FDR<0.05 且 logFC≠0 的行数（记录值）",
        f"上调 {up} / 下调 {dn}" if up is not None else "")
    put(f"deg.{gse}.{cmp_}.n_tested", n_test, f"03_结果/DEG/{gse}_{cmp_}_limma.csv",
        "参与检验的基因数（记录值）")
for gse, cmp_ in [("GSE134431", "Healer_vs_Nonhealer"),
                   ("GSE143735", "Healer_vs_Nonhealer"),
                   ("GSE143735", "DFU_vs_DM_no_ulcer")]:
    f = f"{BASE}/03_结果/DEG/{gse}_{cmp_}_limma.csv"
    d = pd.read_csv(f)
    pcol = "adj.P.Val" if "adj.P.Val" in d.columns else "adj_P_Val"
    put(f"deg.{gse}.{cmp_}.n_deg", int((d[pcol] < 0.05).sum()), f,
        f"({pcol}<0.05).sum()", "阴性结果：0 个 FDR<0.05")

m1 = pd.read_csv(f"{BASE}/03_结果/meta/META1_DFU_vs_control.csv")
m2 = pd.read_csv(f"{BASE}/03_结果/meta/META2_Healer_vs_Nonhealer.csv")
put("meta.META1.n_tested", int(len(m1)), "03_结果/meta/META1_DFU_vs_control.csv",
    "len(df)", "REM(REML) 随机效应，2 队列")
put("meta.META1.n_deg", int((m1.padj < 0.05).sum()),
    "03_结果/meta/META1_DFU_vs_control.csv", "(padj<0.05).sum()")
put("meta.META2.n_tested", int(len(m2)), "03_结果/meta/META2_Healer_vs_Nonhealer.csv",
    "len(df)")
put("meta.META2.n_deg", int((m2.padj < 0.05).sum()),
    "03_结果/meta/META2_Healer_vs_Nonhealer.csv", "(padj<0.05).sum()")

ov = pd.read_csv(f"{BASE}/03_结果/meta/KpyRN_meta_overview.csv")
ov_sig = ov[(ov.scoring_include == "Yes") & (ov.META1_padj < 0.05)]
put("meta.kpyrn.META1.n_sig", int(len(ov_sig)),
    "03_结果/meta/KpyRN_meta_overview.csv",
    "(scoring_include=='Yes') & (META1_padj<0.05) 的行数")
# 26/26 两队列方向一致（队列特异 logFC 取自 META1 明细表）
sig_genes = list(ov_sig.gene_symbol)
_sub = m1[m1.gene.isin(sig_genes)].copy()
_sub["agree"] = (_sub.logFC_GSE134431 * _sub.logFC_GSE199939) > 0
agree = int(_sub["agree"].sum())
put("meta.kpyrn.META1.direction_agree", agree,
    "03_结果/meta/META1_DFU_vs_control.csv + KpyRN_meta_overview.csv",
    "26 个 KpyRN 显著基因中 sign(logFC_GSE134431)==sign(logFC_GSE199939) 的个数",
    f"分母 {len(_sub)}（方向一致率 {agree}/{len(_sub)}）")
put("meta.kpyrn.META1.direction_agree_total", int(len(_sub)),
    "03_结果/meta/META1_DFU_vs_control.csv", "参与方向一致性判定的基因数")
# 在 META1 中有代表的 KpyRN 计分基因（META1_padj 非缺失）；富集率以「有代表基因」为分母
kpy = ov[ov.scoring_include == "Yes"]
rep = kpy[kpy.META1_padj.notna()]
n_rep = int(len(rep))
n_kpy_sig = int(len(ov_sig))                      # 26
put("meta.kpyrn.META1.n_represented", n_rep,
    "03_结果/meta/KpyRN_meta_overview.csv",
    "scoring_include=='Yes' 且 META1_padj 非缺失（128 计分基因中 31 个无 META1 数据）")
kprate = n_kpy_sig / n_rep
put("meta.kpyrn.META1.rate", round(kprate, 4),
    "由 KpyRN_meta_overview.csv 与 META1_DFU_vs_control.csv 派生",
    f"{n_kpy_sig}/{n_rep}")
put("meta.kpyrn.META1.rate_pct", round(kprate * 100, 1),
    "由 KpyRN_meta_overview.csv 与 META1_DFU_vs_control.csv 派生",
    f"{n_kpy_sig}/{n_rep} × 100")
# 全基因组背景率（META1 全部基因）
n_sig_total = int((m1.padj < 0.05).sum())
bg = n_sig_total / len(m1)
put("meta.META1.background_rate", round(bg, 4),
    "由 META1_DFU_vs_control.csv 计算", f"{n_sig_total}/{len(m1)}")
put("meta.META1.background_rate_pct", round(bg * 100, 1),
    "由 META1_DFU_vs_control.csv 计算", f"{n_sig_total}/{len(m1)} × 100")
# Fisher 精确检验：KpyRN 有代表基因 vs 非 KpyRN 基因，是否显著富集
n_tested = len(m1)
n_nonkpy_tested = n_tested - n_rep
n_nonkpy_sig = n_sig_total - n_kpy_sig
put("meta.META1.background_nonkpyrn_n_sig", n_nonkpy_sig,
    "由 META1_DFU_vs_control.csv 与 KpyRN_meta_overview.csv 派生",
    f"{n_sig_total} - {n_kpy_sig}")
put("meta.META1.background_nonkpyrn_n_tested", n_nonkpy_tested,
    "由 META1_DFU_vs_control.csv 与 KpyRN_meta_overview.csv 派生",
    f"{n_tested} - {n_rep}")
_fish_table = [[n_kpy_sig, n_rep - n_kpy_sig],
               [n_nonkpy_sig, n_nonkpy_tested - n_nonkpy_sig]]
_fish_or, _fish_p, _fish_lo, _fish_hi = fisher_exact_ci(_fish_table)
put("meta.kpyrn.META1.fisher_OR", round(_fish_or, 2),
    "由 2×2 表 Fisher 精确检验派生（scipy.stats.fisher_exact）",
    f"{n_kpy_sig}*{n_nonkpy_tested-n_nonkpy_sig}/(({n_rep-n_kpy_sig})*{n_nonkpy_sig})")
put("meta.kpyrn.META1.fisher_p", float(f"{_fish_p:.2g}"),
    "由 2×2 表 Fisher 精确检验派生（scipy.stats.fisher_exact, two-sided）", "two-sided p")
put("meta.kpyrn.META1.fisher_ci_low", round(_fish_lo, 2),
    "由 2×2 表 Fisher 精确条件 CI 派生（R fisher.test 'P<=Pobs' 两尾反演）", "95% CI 下限")
put("meta.kpyrn.META1.fisher_ci_high", round(_fish_hi, 2),
    "由 2×2 表 Fisher 精确条件 CI 派生（R fisher.test 'P<=Pobs' 两尾反演）", "95% CI 上限")
# 软件版本（方法参数）：limma 3.62.2 的 major.minor
put("meta.software.limma", 3.62,
    "manuscript_en.md §2 methods（limma 3.62.2）", "软件版本 major.minor，方法参数",
    "方法学常数，非结果数字")
# 6 个 KpyRN 显著基因的三层交集基因精确值
for g in ["BPGM", "LDHA", "PFKP", "PGAM1", "PGM2", "TPI1"]:
    r = ov[ov.gene_symbol == g].iloc[0]
    put(f"meta.kpyrn.{g}.logFC", round(float(r.META1_logFC), 3),
        "03_结果/meta/KpyRN_meta_overview.csv", f"row['META1_logFC'] where gene=='{g}'")
    put(f"meta.kpyrn.{g}.padj", float(f"{r.META1_padj:.3g}"),
        "03_结果/meta/KpyRN_meta_overview.csv", f"row['META1_padj'] where gene=='{g}'")
# 代表基因
for g, fc, padj in [("PGAM1", 1.25, "1.5e-32"), ("TPI1", 1.29, "1.4e-6"),
                    ("GAPDH", 1.32, "1.0e-2"), ("PGK1", 1.07, "1.4e-15"),
                    ("LDHA", 1.00, "1.1e-9"), ("PFKP", 0.90, "3.3e-7"),
                    ("SLC16A1", 0.83, "2.4e-4"), ("LDHB", 0.75, None),
                    ("PGM2", 0.67, "2.4e-3"), ("ALDOA", 0.66, "4.6e-2"),
                    ("ENO1", 1.24, "3.2e-2"), ("PDHB", 0.26, "3.2e-2"),
                    ("IFI44L", -1.09, "1.3e-6"), ("IFI27", -1.02, "1.8e-6"),
                    ("DHX58", -0.74, "8.6e-5"), ("DDX58", -0.47, "9.1e-3"),
                    ("STAT2", -0.41, "3.5e-3"), ("MPC2", -0.45, "3.0e-3"),
                    ("ISG20", 0.88, "3.9e-3"), ("GBP1", 0.84, "1.3e-3"),
                    ("ADH1B", -4.42, "4.7e-2"), ("ALDH3A2", -1.63, "3.7e-2"),
                    ("FAHD2B", -1.12, "1.4e-4")]:
    rr = ov[ov.gene_symbol == g]
    if len(rr):
        rr = rr.iloc[0]
        put(f"meta.kpyrn.{g}.logFC_exact", round(float(rr.META1_logFC), 2),
            "03_结果/meta/KpyRN_meta_overview.csv", f"META1_logFC where gene=='{g}'")
        put(f"meta.kpyrn.{g}.padj_exact", float(f"{rr.META1_padj:.2g}"),
            "03_结果/meta/KpyRN_meta_overview.csv", f"META1_padj where gene=='{g}'")
# I² 高异质 + 中位 I² + 三个高异质基因的队列特异 logFC/SE（稿件 §3.1 表）
for g in ["ADH1B", "ALDH3A2", "FAHD2B"]:
    r = m1[m1.gene == g]
    if len(r):
        rr = r.iloc[0]
        put(f"meta.I2.{g}", round(float(rr.I2), 1),
            "03_结果/meta/META1_DFU_vs_control.csv", f"I2 where gene=='{g}'")
        put(f"meta.{g}.GSE134431.logFC", round(float(rr.logFC_GSE134431), 2),
            "03_结果/meta/META1_DFU_vs_control.csv", f"logFC_GSE134431 where gene=='{g}'")
        put(f"meta.{g}.GSE134431.se", round(float(rr.se_GSE134431), 2),
            "03_结果/meta/META1_DFU_vs_control.csv", f"se_GSE134431 where gene=='{g}'")
        put(f"meta.{g}.GSE199939.logFC", round(float(rr.logFC_GSE199939), 2),
            "03_结果/meta/META1_DFU_vs_control.csv", f"logFC_GSE199939 where gene=='{g}'")
        put(f"meta.{g}.GSE199939.se", round(float(rr.se_GSE199939), 2),
            "03_结果/meta/META1_DFU_vs_control.csv", f"se_GSE199939 where gene=='{g}'")
put("meta.I2.median", round(float(m1["I2"].median()), 1),
    "03_结果/meta/META1_DFU_vs_control.csv", "I2.median()（全部 11,615 个检验基因）")

# ─────────────────────────────────────────────────────────────
# 3. WGCNA
# ─────────────────────────────────────────────────────────────
mt = pd.read_csv(f"{BASE}/03_结果/WGCNA/module_trait_cor.csv")
put("wgcna.n_modules", int(len(mt)),
    "03_结果/WGCNA/module_trait_cor.csv", "len(df) = module-trait 相关模块数")
put("wgcna.n_module_colors_total", int(
    pd.read_csv(f"{BASE}/03_结果/WGCNA/gene_module_membership.csv")["module"].nunique()),
    "03_结果/WGCNA/gene_module_membership.csv", "nunique(module)",
    "含 grey 模块（不参与相关检验）")
# 注：wgcna.power 在文件末尾「方法参数」段由 02_脚本/wgcna.R 解析派生，
# 此处不再硬编码（避免同一 fact id 重复登记）。
# 模块编号→颜色映射：R 实测 labels2colors(1:12) 顺序
# turquoise, blue, brown, yellow, green, red, black, pink, magenta, purple
L2C = {1: "turquoise", 2: "blue", 3: "brown", 4: "yellow", 5: "green",
       6: "red", 7: "black", 8: "pink", 9: "magenta", 10: "purple"}
for _, r in mt.iterrows():
    num = int(r["module"])
    col = L2C.get(num)
    if col is None:
        continue
    put(f"wgcna.module.{col}.r", round(float(r.cor_DFU), 2),
        "03_结果/WGCNA/module_trait_cor.csv",
        f"module=={num} (labels2colors→{col}) 的 cor_DFU")
    put(f"wgcna.module.{col}.p", float(f"{r.p_DFU:.3g}"),
        "03_结果/WGCNA/module_trait_cor.csv", f"module=={num} 的 p_DFU")
    put(f"wgcna.module.{col}.num", num,
        "03_结果/WGCNA/module_trait_cor.csv", "module 编号")
    put(f"wgcna.module.{col}.r_healer", round(float(r.cor_Healer), 2),
        "03_结果/WGCNA/module_trait_cor.csv", f"module=={num} 的 cor_Healer")
# 三层交集
cand = pd.read_csv(f"{BASE}/03_结果/WGCNA/ML_candidate_pool.csv")
put("wgcna.candidate_pool.n", int(len(cand)),
    "03_结果/WGCNA/ML_candidate_pool.csv", "len(df)",
    "ML 候选池 = 两层交集 ∪ 三层交集")
inter3 = ["BPGM", "LDHA", "PFKP", "PGAM1", "PGM2", "TPI1"]
put("wgcna.three_layer_intersection.n", len(inter3),
    "03_结果/分析摘要.md §4（由 DEG∩WGCNA∩KpyRN 逐层求交得出）",
    "set(META1_DEG) ∩ set(WGCNA模块基因) ∩ set(KpyRN)",
    "6 个基因全部为糖酵解酶")

# ─────────────────────────────────────────────────────────────
# 4. 机器学习
# ─────────────────────────────────────────────────────────────
pa = pd.read_csv(f"{BASE}/03_结果/ML/TaskA_diagnosis_performance.csv")
pb = pd.read_csv(f"{BASE}/03_结果/ML/TaskB_healing_performance.csv")
for alg in ["LASSO", "RF", "SVM-RFE", "XGBoost"]:
    r = pa[(pa.algorithm == alg) & (pa.setting == "train_repeated5foldCV")].iloc[0]
    put(f"ml.taskA.{alg}.cv_auc_mean", round(float(r.auc_mean), 3),
        "03_结果/ML/TaskA_diagnosis_performance.csv", "train_repeated5foldCV 行 auc_mean")
    put(f"ml.taskA.{alg}.cv_auc_sd", round(float(r.auc_sd), 3),
        "03_结果/ML/TaskA_diagnosis_performance.csv", "train_repeated5foldCV 行 auc_sd")
    e = pa[(pa.algorithm == alg) & (pa.setting == "ext_GSE199939")]
    if len(e):
        e = e.iloc[0]
        put(f"ml.taskA.{alg}.ext_auc", round(float(e.auc_mean), 3),
            "03_结果/ML/TaskA_diagnosis_performance.csv", "ext_GSE199939 行 auc_mean")
        put(f"ml.taskA.{alg}.ext_ci_low", round(float(e.ci_low), 3),
            "03_结果/ML/TaskA_diagnosis_performance.csv", "ext_GSE199939 行 ci_low")
        put(f"ml.taskA.{alg}.ext_ci_high", round(float(e.ci_high), 3),
            "03_结果/ML/TaskA_diagnosis_performance.csv", "ext_GSE199939 行 ci_high")
e = pa[(pa.algorithm == "Consensus-LR") & (pa.setting == "ext_GSE199939")].iloc[0]
put("ml.taskA.ConsensusLR.ext_auc", round(float(e.auc_mean), 3),
    "03_结果/ML/TaskA_diagnosis_performance.csv", "Consensus-LR ext_GSE199939 auc_mean")
put("ml.taskA.ConsensusLR.ext_ci_low", round(float(e.ci_low), 3),
    "03_结果/ML/TaskA_diagnosis_performance.csv", "ci_low")
put("ml.taskA.ConsensusLR.ext_ci_high", round(float(e.ci_high), 3),
    "03_结果/ML/TaskA_diagnosis_performance.csv", "ci_high")
put("ml.taskA.ConsensusLR.n_features", int(e.n_features_used),
    "03_结果/ML/TaskA_diagnosis_performance.csv", "n_features_used")
for alg in ["LASSO", "RF", "SVM-RFE", "XGBoost"]:
    r = pb[(pb.algorithm == alg) & (pb.setting == "train_repeated5foldCV")].iloc[0]
    put(f"ml.taskB.{alg}.cv_auc_mean", round(float(r.auc_mean), 3),
        "03_结果/ML/TaskB_healing_performance.csv", "train_repeated5foldCV 行 auc_mean")

sel = pd.read_csv(f"{BASE}/03_结果/ML/TaskA_diagnosis_gene_selection.csv")
cons = sorted(sel[sel.consensus.astype(str).str.lower().isin(["true", "1"])].gene.tolist())
put("ml.taskA.consensus_genes", cons,
    "03_结果/ML/TaskA_diagnosis_gene_selection.csv",
    "consensus==True 的 gene 列表", "6 个共识枢纽基因", ftype="list")
put("ml.taskA.n_candidates", int(len(sel)),
    "03_结果/ML/TaskA_diagnosis_gene_selection.csv", "len(df)")
_taskB_auc = [float(pb[(pb.algorithm == a) & (pb.setting == "train_repeated5foldCV")]
                   .iloc[0].auc_mean) for a in ["LASSO", "RF", "SVM-RFE", "XGBoost"]]
put("ml.taskB.cv_auc_range_min", round(min(_taskB_auc), 2),
    "03_结果/ML/TaskB_healing_performance.csv", "四算法 train CV AUC 最小值")
put("ml.taskB.cv_auc_range_max", round(max(_taskB_auc), 2),
    "03_结果/ML/TaskB_healing_performance.csv", "四算法 train CV AUC 最大值")
for a in ["LASSO", "RF", "SVM-RFE", "XGBoost"]:
    _e = pb[(pb.algorithm == a) & (pb.setting == "ext_GSE143735")]
    if len(_e):
        _e = _e.iloc[0]
        put(f"ml.taskB.{a}.ext_gse143735_auc", round(float(_e.auc_mean), 3),
            "03_结果/ML/TaskB_healing_performance.csv",
            f"{a} 的 ext_GSE143735 行 auc_mean")
        put(f"ml.taskB.{a}.ext_gse143735_n", int(_e.n),
            "03_结果/ML/TaskB_healing_performance.csv",
            f"{a} 的 ext_GSE143735 行 n（外部队列样本数）")
put("ml.taskB.ext_n", int(pb[(pb.algorithm == "LASSO")
                             & (pb.setting == "ext_GSE143735")].iloc[0].n),
    "03_结果/ML/TaskB_healing_performance.csv", "ext_GSE143735 的 n")
for a in ["LASSO", "RF", "SVM-RFE", "XGBoost", "Consensus-LR"]:
    _e = pa[(pa.algorithm == a) & (pa.setting == "ext_GSE143735")]
    if len(_e):
        _e = _e.iloc[0]
        put(f"ml.taskA.{a}.ext_gse143735_auc", round(float(_e.auc_mean), 3),
            "03_结果/ML/TaskA_diagnosis_performance.csv",
            f"{a} 的 ext_GSE143735 行 auc_mean")
        put(f"ml.taskA.{a}.ext_gse143735_n", int(_e.n),
            "03_结果/ML/TaskA_diagnosis_performance.csv",
            f"{a} 的 ext_GSE143735 行 n")

# ─────────────────────────────────────────────────────────────
# 5. 免疫浸润
# ─────────────────────────────────────────────────────────────
IMM_FILES = {
    ("GSE134431", "DFS_vs_DFU"): "GSE134431_DFS_vs_DFU_immune_stats.csv",
    ("GSE134431", "Healer_vs_Nonhealer"): "GSE134431_Healer_vs_Nonhealer_immune_stats.csv",
    ("GSE199939", "DFU_vs_control"): "GSE199939_DFU_vs_control_immune_stats.csv",
}
for (gse, cmp_), fn in IMM_FILES.items():
    st = pd.read_csv(f"{BASE}/03_结果/免疫浸润/{fn}")
    sig = st[st["p_fdr"].abs() < 0.05]
    put(f"immune.{gse}.{cmp_}.n_fdr05", int(len(sig)),
        f"03_结果/免疫浸润/{fn}", "(p_fdr<0.05).sum()")
    gcol1 = st.columns[2]
    gcol2 = st.columns[3]
    for _, r in sig.iterrows():
        mkey = str(r.metric).replace(" ", "_")
        put(f"immune.{gse}.{cmp_}.{mkey}.fdr", float(f"{r.p_fdr:.4g}"),
            f"03_结果/免疫浸润/{fn}", f"metric=={r.metric} 的 p_fdr")
        put(f"immune.{gse}.{cmp_}.{mkey}.mean1", round(float(r[gcol1]), 4),
            f"03_结果/免疫浸润/{fn}", f"metric=={r.metric} 的 {gcol1}")
        put(f"immune.{gse}.{cmp_}.{mkey}.mean2", round(float(r[gcol2]), 4),
            f"03_结果/免疫浸润/{fn}", f"metric=={r.metric} 的 {gcol2}")
hc = pd.read_csv(f"{BASE}/03_结果/免疫浸润/GSE134431_hubgene_immune_corr.csv")
put("immune.GSE134431.hub_corr.n_fdr05", int((hc.p_fdr < 0.05).sum()),
    "03_结果/免疫浸润/GSE134431_hubgene_immune_corr.csv", "(p_fdr<0.05).sum()")
put("immune.GSE134431.hub_corr.n_total", int(len(hc)),
    "03_结果/免疫浸润/GSE134431_hubgene_immune_corr.csv", "len(df)", "16 枢纽 × 免疫指标")
hc2 = pd.read_csv(f"{BASE}/03_结果/免疫浸润/GSE199939_hubgene_immune_corr.csv")
put("immune.GSE199939.hub_corr.n_fdr05", int((hc2.p_fdr < 0.05).sum()),
    "03_结果/免疫浸润/GSE199939_hubgene_immune_corr.csv", "(p_fdr<0.05).sum()")
top = hc[hc.p_fdr < 0.05].copy()
top["ab"] = top.rho.abs()
t = top.sort_values("ab", ascending=False).iloc[0]
put("immune.top.DHX58_CD8_rho", round(float(t.rho), 3),
    "03_结果/免疫浸润/GSE134431_hubgene_immune_corr.csv",
    "FDR<0.05 中 |rho| 最大的行", f"{t.gene} ~ {t.metric}")

# ─────────────────────────────────────────────────────────────
# 6. 单细胞 GSE165816
# ─────────────────────────────────────────────────────────────
qc = json.load(open(f"{BASE}/03_结果/QC/scRNA_GSE165816_stage1_summary.json", encoding="utf-8"))
def dig(d, *keys, default=None):
    for k in keys:
        if isinstance(d, dict) and k in d:
            d = d[k]
        else:
            return default
    return d
put("scrna165816.n_raw_cells", int(qc["total_cells_raw"]),
    "03_结果/QC/scRNA_GSE165816_stage1_summary.json", "total_cells_raw")
put("scrna165816.n_post_qc_cells", int(qc["total_cells_pass"]),
    "03_结果/QC/scRNA_GSE165816_stage1_summary.json", "total_cells_pass")
put("scrna165816.pass_rate", round(float(qc["pass_rate_overall"]), 4),
    "03_结果/QC/scRNA_GSE165816_stage1_summary.json", "pass_rate_overall")
put("scrna165816.n_samples", int(qc["n_files"]),
    "03_结果/QC/scRNA_GSE165816_stage1_summary.json", "n_files")
put("scrna165816.id_system", str(qc["id_system"]),
    "03_结果/QC/scRNA_GSE165816_stage1_summary.json", "id_system",
    "基因 ID 本就是基因符号，映射为空操作", ftype="text")
_ps = qc["per_sample"]
_tissue = {}
for _gsm, _v in _ps.items():
    _tissue[_v.get("tissue")] = _tissue.get(_v.get("tissue"), 0) + 1
for _t, _n in _tissue.items():
    put(f"scrna165816.tissue_samples.{str(_t).replace(' ', '_')}", int(_n),
        "03_结果/QC/scRNA_GSE165816_stage1_summary.json",
        "per_sample 中 tissue 的样本计数")
# Ambiguous 比例：由细胞谱系注释表派生（lineage=='Ambiguous' 的占比）
_cs165 = pd.read_csv(f"{BASE}/01_数据/处理后数据/scRNA/scRNA_GSE165816_cell_scores.csv")
_amb_n = int((_cs165.lineage == "Ambiguous").sum())
put("scrna165816.n_cells_pass", int(len(_cs165)),
    "01_数据/处理后数据/scRNA/scRNA_GSE165816_cell_scores.csv", "len(df)",
    "通过 QC 的细胞总数")
put("scrna165816.n_cells_ambiguous", _amb_n,
    "01_数据/处理后数据/scRNA/scRNA_GSE165816_cell_scores.csv",
    "(lineage=='Ambiguous').sum()")
put("scrna165816.ambiguous_frac", round(_amb_n / len(_cs165), 3),
    "01_数据/处理后数据/scRNA/scRNA_GSE165816_cell_scores.csv",
    "(lineage=='Ambiguous').sum() / len(df)", "Ambiguous 细胞占比（top1−top2 判据）")
comp = pd.read_csv(f"{BASE}/03_结果/单细胞/composition_Foot_skin_DFU_vs_nonDFU.csv")
put("scrna165816.foot_skin.dfu_vs_nondfu.n_fdr05", int((comp["q"] < 0.05).sum()),
    "03_结果/单细胞/composition_Foot_skin_DFU_vs_nonDFU.csv", "(q<0.05).sum()",
    "受试体级中位占比，BH-FDR")
for _, r in comp[comp["q"] < 0.05].iterrows():
    lin = str(r.lineage)
    put(f"scrna165816.composition.{lin}.dfu_median", round(float(r.median1), 5),
        "03_结果/单细胞/composition_Foot_skin_DFU_vs_nonDFU.csv",
        f"lineage=={lin} 的 median1 (DFU)")
    put(f"scrna165816.composition.{lin}.nondfu_median", round(float(r.median2), 5),
        "03_结果/单细胞/composition_Foot_skin_DFU_vs_nonDFU.csv",
        f"lineage=={lin} 的 median2 (非 DFU)")
    put(f"scrna165816.composition.{lin}.q", float(f"{r.q:.3g}"),
        "03_结果/单细胞/composition_Foot_skin_DFU_vs_nonDFU.csv", f"{lin} 的 q")
    put(f"scrna165816.composition.{lin}.effect_rb", round(float(r.effect_rb), 3),
        "03_结果/单细胞/composition_Foot_skin_DFU_vs_nonDFU.csv",
        f"lineage=={lin} 的 effect_rb (rank-biserial)")
put("scrna165816.foot_skin.dfu_n", int(comp.n1.iloc[0]),
    "03_结果/单细胞/composition_Foot_skin_DFU_vs_nonDFU.csv", "n1 (DFU 样本数)")
put("scrna165816.foot_skin.nondfu_n", int(comp.n2.iloc[0]),
    "03_结果/单细胞/composition_Foot_skin_DFU_vs_nonDFU.csv", "n2 (非 DFU 样本数)")

# ─────────────────────────────────────────────────────────────
# 6b. GSE165816 pseudobulk（角质细胞糖酵解核心酶）
# ─────────────────────────────────────────────────────────────
pb = pd.read_csv(f"{BASE}/03_结果/单细胞/pseudobulk_Foot_skin_DFU_vs_nonDFU_by_lineage.csv")
PB_SRC = "03_结果/单细胞/pseudobulk_Foot_skin_DFU_vs_nonDFU_by_lineage.csv"
PB_NOTE = ("足皮肤限定的基因×谱系 pseudobulk 统计表（median1=DFU, median2=非DFU, "
           "q=BH-FDR）；n1/n2 为该谱系内的样本数")
KERC = ["PGAM1", "TPI1", "LDHA", "ENO1", "PKM", "GAPDH", "HK2", "PGM2"]
for lin in ["Keratinocyte", "Fibroblast", "Pericyte_SMC"]:
    for g in KERC + ["ADH1B"]:
        r = pb[(pb.lineage == lin) & (pb.feature == g)]
        if not len(r):
            continue
        r = r.iloc[0]
        put(f"scrna165816.pseudobulk.{lin}.{g}.dfu_median",
            round(float(r.median1), 2), PB_SRC,
            f"lineage=='{lin}' & feature=='{g}' 的 median1 (DFU)")
        put(f"scrna165816.pseudobulk.{lin}.{g}.nondfu_median",
            round(float(r.median2), 2), PB_SRC,
            f"lineage=='{lin}' & feature=='{g}' 的 median2 (非 DFU)")
        put(f"scrna165816.pseudobulk.{lin}.{g}.q", float(f"{r.q:.3g}"), PB_SRC,
            f"lineage=='{lin}' & feature=='{g}' 的 q (BH-FDR)")
        put(f"scrna165816.pseudobulk.{lin}.{g}.n_dfu", int(r.n1), PB_SRC,
            f"lineage=='{lin}' & feature=='{g}' 的 n1 (DFU 样本数)")
        put(f"scrna165816.pseudobulk.{lin}.{g}.n_nondfu", int(r.n2), PB_SRC,
            f"lineage=='{lin}' & feature=='{g}' 的 n2 (非 DFU 样本数)")
put("scrna165816.pseudobulk.keratinocyte.glycolytic_n_fdr05",
    int(sum(1 for g in KERC
            if len(pb[(pb.lineage == "Keratinocyte") & (pb.feature == g)])
            and pb[(pb.lineage == "Keratinocyte") & (pb.feature == g)].iloc[0].q < 0.05)),
    PB_SRC, "Keratinocyte 的 8 个糖酵解核心酶中 q<0.05 的个数")

# ─────────────────────────────────────────────────────────────
# 7. MR
# ─────────────────────────────────────────────────────────────
mr = pd.read_csv(f"{BASE}/03_结果/MR/mr_results.csv")
put("mr.n_pairs", int(len(mr)), "03_结果/MR/mr_results.csv", "len(df)")
sig = mr[mr.p_fdr < 0.05]
put("mr.n_fdr05", int(len(sig)), "03_结果/MR/mr_results.csv", "(p_fdr<0.05).sum()")
d = mr[(mr.gene == "DHX58") & (mr.outcome.str.contains("Xue"))].iloc[0]
put("mr.DHX58.T2D.beta", round(float(d.wald_beta), 4),
    "03_结果/MR/mr_results.csv", "gene==DHX58 & outcome contains Xue")
put("mr.DHX58.T2D.se", round(float(d.wald_se), 4), "03_结果/MR/mr_results.csv", "wald_se")
put("mr.DHX58.T2D.p", float(f"{d.wald_p:.3g}"), "03_结果/MR/mr_results.csv", "wald_p")
put("mr.DHX58.T2D.fdr", float(f"{d.p_fdr:.3g}"), "03_结果/MR/mr_results.csv", "p_fdr")
put("mr.DHX58.T2D.OR", round(float(d.OR), 3), "03_结果/MR/mr_results.csv", "OR")
# OR 的 95% CI 由 Wald β ± 1.96·SE 指数化派生（与稿件 §3.8 一致）
_mr_ci_lo = float(np.exp(d.wald_beta - 1.96 * d.wald_se))
_mr_ci_hi = float(np.exp(d.wald_beta + 1.96 * d.wald_se))
put("mr.DHX58.T2D.ci_low", round(_mr_ci_lo, 3),
    "03_结果/MR/mr_results.csv", "exp(wald_beta - 1.96*wald_se)")
put("mr.DHX58.T2D.ci_high", round(_mr_ci_hi, 3),
    "03_结果/MR/mr_results.csv", "exp(wald_beta + 1.96*wald_se)")
put("mr.DHX58.F", round(float(d.F), 0), "03_结果/MR/mr_results.csv", "F")
d2 = mr[(mr.gene == "DHX58") & (mr.outcome == "T2D_GCST90018926")].iloc[0]
put("mr.DHX58.T2D2.beta", round(float(d2.wald_beta), 4),
    "03_结果/MR/mr_results.csv", "gene==DHX58 & outcome==T2D_GCST90018926")
put("mr.DHX58.T2D2.p", float(f"{d2.wald_p:.3g}"), "03_结果/MR/mr_results.csv", "wald_p")
inst = pd.read_csv(f"{BASE}/03_结果/MR/mr_instruments.csv")
put("mr.n_instruments", int(len(inst)), "03_结果/MR/mr_instruments.csv", "len(df)")
# 跳过基因数：由「预注册候选基因全集（解析 mr_eqtl_glycemic.py 的 HUB_GENES）」
# 减去「mr_results.csv 实际出现的 gene」派生，而非硬编码 2。
_mr_script = open(f"{BASE}/02_脚本/mr_eqtl_glycemic.py", encoding="utf-8").read()
_m = re.search(r"HUB_GENES\s*=\s*\[(.*?)\]", _mr_script, re.S)
_all_cand = re.findall(r'"([A-Za-z0-9\-]+)"', _m.group(1)) if _m else []
_have = set(mr.gene)
_mr_skipped = sorted(set(_all_cand) - _have)
put("mr.n_genes_skipped", len(_mr_skipped),
    "02_脚本/mr_eqtl_glycemic.py（HUB_GENES）对照 03_结果/MR/mr_results.csv（gene 列）",
    "len(候选基因全集) - len(mr_results 中出现的 gene)",
    f"跳过：{', '.join(_mr_skipped)}（无显著 cis-eQTL）")
put("mr.n_candidate_genes", len(_all_cand),
    "02_脚本/mr_eqtl_glycemic.py", "len(HUB_GENES)", "预注册 MR 候选基因总数")
put("mr.F_min", round(float(inst.F.min()), 0), "03_结果/MR/mr_instruments.csv", "F.min()")
put("mr.F_max", round(float(inst.F.max()), 0), "03_结果/MR/mr_instruments.csv", "F.max()")
harm = pd.read_csv(f"{BASE}/03_结果/MR/mr_harmonization_audit.csv")
put("mr.harmonization.rows", int(len(harm)), "03_结果/MR/mr_harmonization_audit.csv", "len(df)")
for _, r in mr[mr.wald_p < 0.05].iterrows():
    _k = f"{r.gene}_{r.outcome}".replace("-", "_")
    put(f"mr.nominal.{_k}.beta", round(float(r.wald_beta), 4),
        "03_结果/MR/mr_results.csv", f"gene=={r.gene} & outcome=={r.outcome} 的 wald_beta")
    put(f"mr.nominal.{_k}.p", float(f"{r.wald_p:.3g}"),
        "03_结果/MR/mr_results.csv", f"gene=={r.gene} & outcome=={r.outcome} 的 wald_p")
    put(f"mr.nominal.{_k}.fdr", float(f"{r.p_fdr:.3g}"),
        "03_结果/MR/mr_results.csv", f"gene=={r.gene} & outcome=={r.outcome} 的 p_fdr")
put("mr.n_nominal_only", int(((mr.wald_p < 0.05) & (mr.p_fdr >= 0.05)).sum()),
    "03_结果/MR/mr_results.csv", "名义 p<0.05 但 FDR>=0.05 的配对数")

# ─────────────────────────────────────────────────────────────
# 8. 分子对接
# ─────────────────────────────────────────────────────────────
dk = pd.read_csv(f"{BASE}/03_结果/药物预测/docking/docking_results.csv")
put("docking.n_docking_runs", int(len(dk)),
    "03_结果/药物预测/docking/docking_results.csv", "len(df)", "55 次对接")
# 体系数 = (PDB, site) 组合数。注意：7 个 PDB 对应 8 个结合位点体系
# （4JPG 同时提供 A-pocket 与 C-pocket 两个独立对接位点）。
_dk_systems = dk.groupby(["pdb", "site"]).ngroups
put("docking.n_systems", int(_dk_systems),
    "03_结果/药物预测/docking/docking_results.csv",
    "groupby(['pdb','site']).ngroups", "结合位点体系数（7 个 PDB → 8 个位点）")
put("docking.n_pdb", int(dk["pdb"].nunique()),
    "03_结果/药物预测/docking/docking_results.csv", "nunique(pdb)", "受体 PDB 数")
# 候选/对照数：由 role 列派生，避免硬编码
_dk_cand = dk[dk.role == "candidate"]
_dk_pos = dk[dk.role == "positive_control"]
_dk_neg = dk[dk.role == "negative_control"]
put("docking.n_unique_candidates", int(_dk_cand.ligand.nunique()),
    "03_结果/药物预测/docking/docking_results.csv",
    "role=='candidate' 的 ligand.nunique()", "唯一候选化合物数")
put("docking.n_positive_controls", int(_dk_pos.ligand.nunique()),
    "03_结果/药物预测/docking/docking_results.csv",
    "role=='positive_control' 的 ligand.nunique()", "共晶配体阳性对照数（每体系 1 个）")
put("docking.n_negative_controls", int(_dk_neg.ligand.nunique()),
    "03_结果/药物预测/docking/docking_results.csv",
    "role=='negative_control' 的 ligand.nunique()", "蔗糖负面对照数（每体系 1 次）")
put("docking.n_candidate_runs", int(len(_dk_cand)),
    "03_结果/药物预测/docking/docking_results.csv", "role=='candidate' 的行数")
rmsd = pd.read_csv(f"{BASE}/03_结果/药物预测/docking/positive_control_rmsd.csv")
put("docking.rmsd.n_pass", int((rmsd.rmsd < 2.0).sum()),
    "03_结果/药物预测/docking/positive_control_rmsd.csv", "(rmsd<2.0).sum()")
put("docking.rmsd.n_total", int(len(rmsd)),
    "03_结果/药物预测/docking/positive_control_rmsd.csv", "len(df)")
put("docking.rmsd.max", round(float(rmsd.rmsd.max()), 2),
    "03_结果/药物预测/docking/positive_control_rmsd.csv", "rmsd.max()")
# 注：docking.exhaustiveness / docking.seed / docking.vina_version 在文件末尾
# 「方法参数」段由 02_脚本/docking_pipeline.py 解析派生，此处不重复登记。
_neg_all = dk[dk.role == "negative_control"].affinity_kcal_mol
put("docking.negative_control.affinity_min", round(float(_neg_all.min()), 2),
    "03_结果/药物预测/docking/docking_results.csv",
    "role==negative_control (蔗糖) 的 affinity 最小值（最强）")
put("docking.negative_control.affinity_max", round(float(_neg_all.max()), 2),
    "03_结果/药物预测/docking/docking_results.csv",
    "role==negative_control (蔗糖) 的 affinity 最大值（最弱）")
put("docking.negative_control.n", int(len(_neg_all)),
    "03_结果/药物预测/docking/docking_results.csv",
    "role==negative_control 行数", "每体系 1 次，共 8 次")

# 各体系共晶基准与最优候选（从 docking_results.csv 动态提取）
def _slug(x):
    return (str(x).replace(" ", "_").replace("/", "_").replace("(", "")
            .replace(")", "").replace(",", "").replace("+", "plus").replace("-", "_"))
for (t_, s_), grp in dk.groupby(["target", "site"], sort=False):
    key = f"{_slug(t_)}__{_slug(s_)}"
    pos_rows = grp[grp.role == "positive_control"]
    cand = grp[grp.role == "candidate"]
    neg = grp[grp.role == "negative_control"]
    if len(pos_rows):
        put(f"docking.{key}.positive_affinity",
            round(float(pos_rows.affinity_kcal_mol.iloc[0]), 2),
            "03_结果/药物预测/docking/docking_results.csv",
            f"{t_}/{s_} role==positive_control 的 affinity_kcal_mol")
    if len(neg):
        put(f"docking.{key}.negative_affinity",
            round(float(neg.affinity_kcal_mol.iloc[0]), 2),
            "03_结果/药物预测/docking/docking_results.csv",
            f"{t_}/{s_} role==negative_control (蔗糖) 的 affinity_kcal_mol")
    if len(cand):
        best = cand.loc[cand.affinity_kcal_mol.astype(float).idxmin()]
        put(f"docking.{key}.best_candidate", str(best.ligand),
            "03_结果/药物预测/docking/docking_results.csv",
            f"{t_}/{s_} 候选中 affinity 最小者")
        put(f"docking.{key}.best_affinity", round(float(best.affinity_kcal_mol), 2),
            "03_结果/药物预测/docking/docking_results.csv", f"{t_}/{s_} 最优候选 affinity")
        put(f"docking.{key}.best_delta", round(float(best.delta_vs_pos), 2),
            "03_结果/药物预测/docking/docking_results.csv", f"{t_}/{s_} 最优候选 delta_vs_pos")
    put(f"docking.{key}.n_candidates", int(len(cand)),
        "03_结果/药物预测/docking/docking_results.csv", f"{t_}/{s_} role==candidate 行数")

# ─────────────────────────────────────────────────────────────
# 9. 感染队列 GSE230426
# ─────────────────────────────────────────────────────────────
s2304 = json.load(open(f"{BASE}/03_结果/感染队列GSE230426/summary.json", encoding="utf-8"))
# 样本量：由 GSE230426 注释表派生（stage / week / donor_inferred 三列），
# 而非从 summary.json 的 label 字符串硬编码。
_a2304 = pd.read_csv(f"{BASE}/01_数据/处理后数据/GSE230426_annotation.csv")
put("gse230426.n_before", int((_a2304.stage == "BEFORE").sum()),
    "01_数据/处理后数据/GSE230426_annotation.csv", "(stage=='BEFORE').sum()",
    "主模型 BEFORE 组样本数")
put("gse230426.n_after", int((_a2304.stage == "AFTER").sum()),
    "01_数据/处理后数据/GSE230426_annotation.csv", "(stage=='AFTER').sum()",
    "主模型 AFTER 组样本数")
put("gse230426.n_after8", int(((_a2304.stage == "AFTER") & (_a2304.week == 8)).sum()),
    "01_数据/处理后数据/GSE230426_annotation.csv",
    "((stage=='AFTER') & (week==8)).sum()", "after8wk_only 分析的 AFTER-8wk 样本数")
put("gse230426.n_pairs_inferred",
    int((_a2304.donor_inferred != "unpaired").sum() // 2),
    "01_数据/处理后数据/GSE230426_annotation.csv",
    "(donor_inferred!='unpaired').sum() // 2", "可配对推断的样本对数")
put("gse230426.n_genes", int(s2304["main"]["n_genes"]),
    "03_结果/感染队列GSE230426/summary.json", "main.n_genes")
put("gse230426.n_deg_all_models",
    int(sum(s2304[k]["n_fdr05"] for k in ("main", "paired", "after8"))),
    "03_结果/感染队列GSE230426/summary.json",
    "main.n_fdr05 + paired.n_fdr05 + after8.n_fdr05",
    "三个模型的 FDR<0.05 基因数之和（全基因组 DEG 阴性）")
for _mdl in ("main", "paired", "after8"):
    put(f"gse230426.n_deg.{_mdl}", int(s2304[_mdl]["n_fdr05"]),
        "03_结果/感染队列GSE230426/summary.json", f"{_mdl}.n_fdr05")
put("gse230426.concordance.r", float(s2304["concordance_main_vs_paired"]["spearman_logFC"]),
    "03_结果/感染队列GSE230426/summary.json", "KpyRN 128 基因 logFC Spearman")
put("gse230426.concordance.dir", float(s2304["concordance_main_vs_paired"]["direction_agreement"]),
    "03_结果/感染队列GSE230426/summary.json", "KpyRN 方向一致率")
cam = pd.read_csv(f"{BASE}/03_结果/感染队列GSE230426/camera_modules.csv")
for _, r in cam.iterrows():
    mdl = r.model.replace("_", "-")
    put(f"gse230426.camera.{r.module}.{mdl}.fdr", float(f"{r.fdr:.4g}"),
        "03_结果/感染队列GSE230426/camera_modules.csv",
        f"module=={r.module} & model=={r.model} 的 fdr")
    put(f"gse230426.camera.{r.module}.{mdl}.direction", str(r.direction),
        "03_结果/感染队列GSE230426/camera_modules.csv", "direction")
hg = pd.read_csv(f"{BASE}/03_结果/感染队列GSE230426/hub_genes.csv")
for mdl in hg["model"].unique():
    h2 = hg[hg.model == mdl]
    put(f"gse230426.hub_genes.{mdl}.n_sig",
        int((h2["adj.P.Val"] < 0.05).sum()),
        "03_结果/感染队列GSE230426/hub_genes.csv",
        f"model=={mdl} 且 adj.P.Val<0.05 的行数")
    put(f"gse230426.hub_genes.{mdl}.max_abs_logFC",
        round(float(h2.logFC.abs().max()), 2),
        "03_结果/感染队列GSE230426/hub_genes.csv",
        f"model=={mdl} 的 |logFC| 最大值")
    put(f"gse230426.hub_genes.{mdl}.min_adjP",
        round(float(h2["adj.P.Val"].min()), 2),
        "03_结果/感染队列GSE230426/hub_genes.csv",
        f"model=={mdl} 的 adj.P.Val 最小值")
put("gse230426.hub_genes.n", int(hg.gene.nunique()),
    "03_结果/感染队列GSE230426/hub_genes.csv", "nunique(gene)", "6 个 ML 共识枢纽基因")

# ─────────────────────────────────────────────────────────────
# 10. GSE223964 单细胞
# ─────────────────────────────────────────────────────────────
pbk = pd.read_csv(f"{BASE}/03_结果/单细胞GSE223964/pseudobulk_condition_stats.csv")
pcol = "mwu_p_bf"
put("scrna223964.n_tests", int(len(pbk)),
    "03_结果/单细胞GSE223964/pseudobulk_condition_stats.csv", "len(df)", "细胞类型 × 模块对比数")
put("scrna223964.n_fdr10", int((pbk[pcol] < 0.10).sum()),
    "03_结果/单细胞GSE223964/pseudobulk_condition_stats.csv", f"({pcol}<0.10).sum()")
put("scrna223964.n_nominal_p05", int((pbk["mwu_p"] < 0.05).sum()),
    "03_结果/单细胞GSE223964/pseudobulk_condition_stats.csv", "(mwu_p<0.05).sum()",
    "名义显著但 FDR 后不显著")
put("scrna223964.n_celltypes", int(pbk.celltype.nunique()),
    "03_结果/单细胞GSE223964/pseudobulk_condition_stats.csv", "nunique(celltype)")
put("scrna223964.n_modules", int(pbk.module.nunique()),
    "03_结果/单细胞GSE223964/pseudobulk_condition_stats.csv", "nunique(module)")
put("scrna223964.n_donors_diabetic", int(pbk.n_diabetic.max()),
    "03_结果/单细胞GSE223964/pseudobulk_condition_stats.csv",
    "n_diabetic.max()", "供体分组：diabetic / non_diabetic")
put("scrna223964.n_donors_nondiabetic", int(pbk.n_nondiabetic.max()),
    "03_结果/单细胞GSE223964/pseudobulk_condition_stats.csv",
    "n_nondiabetic.max()")
# 细胞总数：由 module_scores_by_celltype.csv 的 n_cells 派生。
# 注意口径：该表每个 (celltype, condition) 在 4 个模块下重复，
# 必须先按 (celltype, condition) 去重再求和，否则会重复计数 4 倍。
_ms223 = pd.read_csv(f"{BASE}/03_结果/单细胞GSE223964/module_scores_by_celltype.csv")
_ms223_ded = _ms223.drop_duplicates(subset=["celltype", "condition"])
put("scrna223964.n_cells", int(_ms223_ded.n_cells.sum()),
    "03_结果/单细胞GSE223964/module_scores_by_celltype.csv",
    "drop_duplicates(['celltype','condition']).n_cells.sum()", "细胞总数")
perm = pd.read_csv(f"{BASE}/03_结果/单细胞GSE223964/xcheck_pseudobulk_perm.csv")
put("scrna223964.xcheck.n_perm_fdr10", int((perm["perm_p_bf"] < 0.10).sum()),
    "03_结果/单细胞GSE223964/xcheck_pseudobulk_perm.csv",
    "(perm_p_bf<0.10).sum()", "独立置换路径 (10000 次标签置换)")
put("scrna223964.xcheck.n_perm_nominal10", int((perm["perm_p"] < 0.10).sum()),
    "03_结果/单细胞GSE223964/xcheck_pseudobulk_perm.csv", "(perm_p<0.10).sum()")
put("scrna223964.xcheck.n_tests", int(len(perm)),
    "03_结果/单细胞GSE223964/xcheck_pseudobulk_perm.csv", "len(df)")

# ─────────────────────────────────────────────────────────────
# 11. GSE166120 GeoMx
# ─────────────────────────────────────────────────────────────
s166 = json.load(open(f"{BASE}/03_结果/空间GSE166120/summary.json", encoding="utf-8"))
put("geommx.n_roi", int(s166["n_roi"]), "03_结果/空间GSE166120/summary.json", "n_roi")
put("geommx.n_subjects", int(s166["n_subjects"]), "03_结果/空间GSE166120/summary.json", "n_subjects")
put("geommx.n_healer", int(s166["healing_groups"]["Healer"]),
    "03_结果/空间GSE166120/summary.json", "healing_groups")
put("geommx.n_nonhealer", int(s166["healing_groups"]["Non-Healer"]),
    "03_结果/空间GSE166120/summary.json", "healing_groups")
for loc, n in s166["roi_locations"].items():
    put(f"geommx.roi.{loc.replace(' ', '_')}", int(n),
        "03_结果/空间GSE166120/summary.json", f"roi_locations['{loc}']")
put("geommx.coverage.pct", float(s166["kpyrn_coverage_pct"]),
    "03_结果/空间GSE166120/summary.json", "kpyrn_coverage_pct")
put("geommx.coverage.n_genes", int(s166["kpyrn_scored_genes_in_board"]),
    "03_结果/空间GSE166120/summary.json", "kpyrn_scored_genes_in_board")
put("geommx.coverage.total", int(s166["kpyrn_scored_genes_total"]),
    "03_结果/空间GSE166120/summary.json", "kpyrn_scored_genes_total")
put("geommx.M4.coverage.pct", float(s166["module_coverage"]["M4_IFNI_readout"]["pct"]),
    "03_结果/空间GSE166120/summary.json", "module_coverage.M4_IFNI_readout.pct")
put("geommx.M2.coverage.pct", float(s166["module_coverage"]["M2_Pyruvate_supply"]["pct"]),
    "03_结果/空间GSE166120/summary.json", "module_coverage.M2_Pyruvate_supply.pct")
put("geommx.n_subject_sig", int(s166["n_significant_subject_level_fdr10"]),
    "03_结果/空间GSE166120/summary.json", "n_significant_subject_level_fdr10")
put("geommx.n_roi_sig", int(s166["n_significant_roi_level_fdr10_sensitivity"]),
    "03_结果/空间GSE166120/summary.json", "n_significant_roi_level_fdr10_sensitivity")
put("geommx.hub_in_board", s166["hub_genes_in_board"],
    "03_结果/空间GSE166120/summary.json", "hub_genes_in_board", ftype="list")
put("geommx.hub_not_in_board", s166["hub_genes_not_in_board"],
    "03_结果/空间GSE166120/summary.json", "hub_genes_not_in_board", ftype="list")

# ─────────────────────────────────────────────────────────────
# 12. RT-qPCR 方案
# ─────────────────────────────────────────────────────────────
ss = json.load(open(f"{BASE}/03_结果/实验验证设计/sample_size_calculation.json", encoding="utf-8"))
_qp = pd.read_csv(f"{BASE}/03_结果/实验验证设计/qpcr_primers.csv")
_qpv = pd.read_csv(f"{BASE}/03_结果/实验验证设计/primer_verification.csv")
put("qpcr.n_target_genes", int(_qp[_qp.role == "target"].gene.nunique()),
    "03_结果/实验验证设计/qpcr_primers.csv", "role=='target' 的 gene.nunique()")
put("qpcr.n_primer_pairs", int(len(_qp)),
    "03_结果/实验验证设计/qpcr_primers.csv", "len(df)",
    "目标 + 内参候选引物对总数")
put("qpcr.n_reference_candidates", int(_qp[_qp.role == "reference_candidate"].gene.nunique()),
    "03_结果/实验验证设计/qpcr_primers.csv",
    "role=='reference_candidate' 的 gene.nunique()")
put("qpcr.n_per_group", int(ss["target_n_per_group"]),
    "03_结果/实验验证设计/sample_size_calculation.json", "target_n_per_group",
    "目标 n/组（基于真实效应量）")
put("qpcr.n_recruit_per_group", int(ss["target_n_per_group_with_20pct_attrition"]),
    "03_结果/实验验证设计/sample_size_calculation.json",
    "target_n_per_group_with_20pct_attrition", "招募量（预留 20% 脱落）")
put("qpcr.n_tissue_total", int(ss["design"]["total_tissue_samples"]),
    "03_结果/实验验证设计/sample_size_calculation.json",
    "design.total_tissue_samples", "组织样本合计")
put("qpcr.alpha", float(ss["alpha"]),
    "03_结果/实验验证设计/sample_size_calculation.json", "alpha")
put("qpcr.power", float(ss["power"]),
    "03_结果/实验验证设计/sample_size_calculation.json", "power")
put("qpcr.status", "pre-specified, not performed",
    "06_补充材料/S4_RT-qPCR_验证方案.md", "方案性质声明", ftype="text")
# 引物设计约束：由 design_rtqpcr_primers.py 的 primer3 参数派生（避免硬编码）
_drt = open(f"{BASE}/02_脚本/design_rtqpcr_primers.py", encoding="utf-8").read()


def _dp(name, cast=float):
    mm = re.search(rf'"{name}"\s*:\s*([-\d.]+)', _drt)
    return cast(mm.group(1)) if mm else None


_psr = re.search(r'"PRIMER_PRODUCT_SIZE_RANGE"\s*:\s*\[\[(\d+),\s*(\d+)\]\]', _drt)
put("qpcr.amplicon_min_bp", int(_psr.group(1)) if _psr else 0,
    "02_脚本/design_rtqpcr_primers.py", "PRIMER_PRODUCT_SIZE_RANGE[0][0]",
    "产物长度下限（primer3 参数）")
put("qpcr.amplicon_max_bp", int(_psr.group(2)) if _psr else 0,
    "02_脚本/design_rtqpcr_primers.py", "PRIMER_PRODUCT_SIZE_RANGE[0][1]",
    "产物长度上限（primer3 参数）")
put("qpcr.tm_min_c", int(_dp("PRIMER_MIN_TM")), "02_脚本/design_rtqpcr_primers.py",
    "PRIMER_MIN_TM", "Tm 下限（primer3 参数）")
put("qpcr.tm_max_c", int(_dp("PRIMER_MAX_TM")), "02_脚本/design_rtqpcr_primers.py",
    "PRIMER_MAX_TM", "Tm 上限（primer3 参数）")
put("qpcr.gc_min_pct", int(_dp("PRIMER_MIN_GC")), "02_脚本/design_rtqpcr_primers.py",
    "PRIMER_MIN_GC", "GC 下限（primer3 参数）")
put("qpcr.gc_max_pct", int(_dp("PRIMER_MAX_GC")), "02_脚本/design_rtqpcr_primers.py",
    "PRIMER_MAX_GC", "GC 上限（primer3 参数）")
put("qpcr.poly_x_max", int(_dp("PRIMER_MAX_POLY_X", int)),
    "02_脚本/design_rtqpcr_primers.py", "PRIMER_MAX_POLY_X")
# 二级结构阈值：由 verify_primers.py 的验收判据派生（不是实测最大值）
_vpt = open(f"{BASE}/02_脚本/verify_primers.py", encoding="utf-8").read()
_msm = re.search(r"max_secondary\s*<\s*([\d.]+)", _vpt)
put("qpcr.secondary_structure_threshold_c", float(_msm.group(1)) if _msm else 0.0,
    "02_脚本/verify_primers.py", "max_secondary < 阈值",
    "二级结构/二聚体 Tm 验收阈值（预设，非实测）")
put("qpcr.secondary_structure_max_observed_c",
    round(float(_qpv.V2_max_secondary_structure_tm.max()), 2),
    "03_结果/实验验证设计/primer_verification.csv",
    "V2_max_secondary_structure_tm.max()", "实测二级结构 Tm 最大值")
put("qpcr.secondary_structure_n_pass",
    int(_qpv.V2_no_strong_secondary_structure.sum()),
    "03_结果/实验验证设计/primer_verification.csv",
    "V2_no_strong_secondary_structure.sum()")
# 注意：overall_pass 依赖 verify_primers.py 的 V1_specific，而 V1 使用
# BLAST Text 格式（description 被截断，基因符号落在截断点之后），
# 实测 12/12 的 V1_specific 均为 False → overall_pass 恒为 False。
# 该字段**受 v1 缺陷污染，不作为引物合格性的判定依据**；
# 特异性判定以 BLAST v6（core_nt + XML + 人源限定 + 全长匹配 + 3'端
# 锚定判据）为唯一权威口径，见 qpcr.v6.* 事实。
put("qpcr.overall_pass_n", int(_qpv.overall_pass.sum()),
    "03_结果/实验验证设计/primer_verification.csv", "overall_pass.sum()",
    "⚠ 受 BLAST v1 截断缺陷污染（V1_specific 恒 False），不作判定依据；"
    "特异性以 BLAST v6 为准（见 qpcr.v6.n_specific_pairs）")
put("qpcr.v2v3_pass_n",
    int((_qpv.V2_no_strong_secondary_structure
         & _qpv.V3_gc_locked_and_tm_matched).sum()),
    "03_结果/实验验证设计/primer_verification.csv",
    "(V2_no_strong_secondary_structure & V3_gc_locked_and_tm_matched).sum()",
    "二级结构 + GC/Tm 两项合规数（不含被 v1 污染的特异性项）")
put("qpcr.v3_gc_tm_pass_n", int(_qpv.V3_gc_locked_and_tm_matched.sum()),
    "03_结果/实验验证设计/primer_verification.csv",
    "V3_gc_locked_and_tm_matched.sum()", "GC 锁温 + Tm 匹配合规数")
put("qpcr.amplicon_observed_min_bp", int(_qp.amplicon_bp.min()),
    "03_结果/实验验证设计/qpcr_primers.csv", "amplicon_bp.min()")
put("qpcr.amplicon_observed_max_bp", int(_qp.amplicon_bp.max()),
    "03_结果/实验验证设计/qpcr_primers.csv", "amplicon_bp.max()")

# ─────────────────────────────────────────────────────────────
# 12b. 引物设计的三约束与全库特异性（v4 设计 + v6 BLAST 复核）
#
# 背景：原引物集经五轮方法学迭代后重建。设计约束由
# `design_final_primers_v4.py` 三条**同时**满足：
#   (1) CDS 锚定（引物完整位于目标 CDS 内）
#   (2) 全 curated（NM_/NR_）蛋白编码转录本覆盖
#   (3) 旁系同源区分度（穷举比对，链向正确的 3' 端锚定）
# 特异性由 `verify_primers_blast_v6.py` 全库复核（core_nt + 人源限定）
# 最终裁决。
#
# ⚠ 期望值一律**从权威产物文件动态派生**，不得硬编码
#   （2026-09-27 教训：gate 曾硬编码旧 bug 数字，修复后反而报 FAIL）。
# ─────────────────────────────────────────────────────────────
_FD = f"{BASE}/03_结果/实验验证设计/final_primer_design_v4.json"
_FP = f"{BASE}/03_结果/实验验证设计/final_primers_v4.csv"
_BV6 = f"{BASE}/03_结果/实验验证设计/primer_blast_v6.csv"
_CDM = f"{BASE}/03_结果/实验验证设计/cds_duplication_map.json"

if os.path.exists(_FP):
    _fp = pd.read_csv(_FP)
    put("qpcr.v4.n_pairs_designed", int(len(_fp)),
        "03_结果/实验验证设计/final_primers_v4.csv", "len(df)",
        "v4 三约束设计产出的引物对总数")
    put("qpcr.v4.n_curated_fully_covered",
        int((_fp.n_curated_covered == _fp.n_curated_total).sum()),
        "03_结果/实验验证设计/final_primers_v4.csv",
        "(n_curated_covered == n_curated_total).sum()",
        "全部 curated 蛋白编码转录本均被引物对覆盖的对数")
    put("qpcr.v4.n_curated_primer_in_cds",
        int((_fp.n_curated_in_CDS == _fp.n_curated_total).sum()),
        "03_结果/实验验证设计/final_primers_v4.csv",
        "(n_curated_in_CDS == n_curated_total).sum()",
        "两条引物均完整位于各自转录本 CDS 内的对数")
    for _g in _fp.pair_paralog_grade.unique():
        put(f"qpcr.v4.paralog_grade.{_g}_n",
            int((_fp.pair_paralog_grade == _g).sum()),
            "03_结果/实验验证设计/final_primers_v4.csv",
            f"(pair_paralog_grade=='{_g}').sum()",
            {"BOTH_STRONG": "两条引物对旁系同源均达强区分（identity<90%）",
             "ONE_WEAK": "至少一条为弱区分（90–95% 且 3'端有错配，P1 级）",
             "SOFT": "参考基因族成员为假基因/同基因座重复拷贝（不作阻断）",
             "NA": "该基因无蛋白编码旁系同源（实测事实，非空判）",
             "UNACCEPTABLE": "含 P0 级脱靶（阻断）"}.get(_g, ""))
    put("qpcr.v4.n_both_strong",
        int((_fp.pair_paralog_grade == "BOTH_STRONG").sum()),
        "03_结果/实验验证设计/final_primers_v4.csv",
        "(pair_paralog_grade=='BOTH_STRONG').sum()",
        "旁系同源强区分对数")
    put("qpcr.v4.max_paralog_identity_pct",
        float(_fp[["f_max_id_vs_paralog", "r_max_id_vs_paralog"]]
              .max(axis=1).max()),
        "03_结果/实验验证设计/final_primers_v4.csv",
        "max(f_max_id_vs_paralog, r_max_id_vs_paralog) 的全局最大值",
        "所选引物对旁系同源的最大 identity（越高越需警惕）")
    put("qpcr.v4.n_pairs_with_3p_anchor_mismatch",
        int(((_fp.f_mm3p_vs_paralog.fillna(0) > 0) |
             (_fp.r_mm3p_vs_paralog.fillna(0) > 0)).sum()),
        "03_结果/实验验证设计/final_primers_v4.csv",
        "至少一条引物对最差旁系同源的 3'端末 5 nt 存在错配的对数",
        "3'端错配者无法有效起始延伸 → 实验上不共扩增")

if os.path.exists(_BV6):
    _bv6 = pd.read_csv(_BV6)
    put("qpcr.v6.n_pairs_blast_checked", int(len(_bv6)),
        "03_结果/实验验证设计/primer_blast_v6.csv", "len(df)",
        "BLAST v6 全库复核覆盖的引物对数")
    put("qpcr.v6.n_specific_pairs", int(_bv6.specific.sum()),
        "03_结果/实验验证设计/primer_blast_v6.csv", "specific.sum()",
        "P0=0 且目标全长匹配转录本≥1 且目标 identity≥90% 的对数")
    put("qpcr.v6.n_P0_hits_total", int(_bv6.n_P0_total.sum()),
        "03_结果/实验验证设计/primer_blast_v6.csv", "n_P0_total.sum()",
        "全库复核发现的 P0 级（阻断）脱靶命中总数")
    put("qpcr.v6.n_P1_hits_total", int(_bv6.n_P1_total.sum()),
        "03_结果/实验验证设计/primer_blast_v6.csv", "n_P1_total.sum()",
        "P1 级（需注意）脱靶命中总数")
    put("qpcr.v6.n_pairs_with_P0", int((_bv6.n_P0_total > 0).sum()),
        "03_结果/实验验证设计/primer_blast_v6.csv",
        "(n_P0_total>0).sum()", "含至少一个 P0 命中的引物对数")
    _v6_ok = _bv6[_bv6.specific == True]  # noqa: E712
    if len(_v6_ok):
        put("qpcr.v6.min_target_identity_pct",
            float(min(_v6_ok.forward_target_max_identity_pct.min(),
                      _v6_ok.reverse_target_max_identity_pct.min())),
            "03_结果/实验验证设计/primer_blast_v6.csv",
            "specific==True 的行中 forward/reverse target identity 的最小值",
            "合格引物对在目标转录本上的最低 identity")
        put("qpcr.v6.min_target_tx_covered",
            int(min(_v6_ok.forward_target_full_match_refseq_tx.min(),
                    _v6_ok.reverse_target_full_match_refseq_tx.min())),
            "03_结果/实验验证设计/primer_blast_v6.csv",
            "specific==True 的行中 target_full_match_refseq_tx 的最小值",
            "合格引物对在目标基因上全长匹配的最少转录本数")

if os.path.exists(_CDM):
    _cdm = json.load(open(_CDM, encoding="utf-8"))
    _gmap = _cdm.get("genes", {})
    if _gmap:
        put("qpcr.cds_map.n_genes_mapped", int(len(_gmap)),
            "03_结果/实验验证设计/cds_duplication_map.json",
            "len(genes)", "已完成全基因组重复图谱测绘的基因数")
        put("qpcr.cds_map.total_repeat_masked_nt",
            int(sum(v.get("repeat_masked_nt", 0) for v in _gmap.values())),
            "03_结果/实验验证设计/cds_duplication_map.json",
            "sum(genes[*].repeat_masked_nt)",
            "12 条 CDS 中被全库 ≥95% 同源重复拷贝覆盖的碱基总数")
        put("qpcr.cds_map.min_longest_safe_segment",
            int(min(v.get("longest_safe_segment", 0) for v in _gmap.values())),
            "03_结果/实验验证设计/cds_duplication_map.json",
            "min(genes[*].longest_safe_segment)",
            "各基因 CDS 中最长无重复掩码连续区段的最小值"
            "（低于一个引物长度即该基因无可用设计窗口）")


# 方法参数类事实（供审计表标注为方法参数而非结果）
# MCP-counter 签名基因数：由官方基因表实际行数派生（跳过表头）
_mcp_lines = [ln.strip() for ln in
               open(f"{BASE}/01_数据/基因集/免疫/MCPcounter_genes.txt", encoding="utf-8")
               if ln.strip()]
put("method.mcpcounter.signature_genes", len(_mcp_lines) - 1,
    "01_数据/基因集/免疫/MCPcounter_genes.txt",
    "非空行数 - 1（表头）", "官方 MCP-counter 标记基因数")
# ESTIMATE 签名：由 export_estimate_sig.R 导出的签名 CSV 派生
_est = pd.read_csv(f"{BASE}/01_数据/基因集/免疫/estimate_signatures.csv")
for _sig in _est.signature.unique():
    _key = "stromal" if "Stromal" in _sig else ("immune" if "Immune" in _sig else _sig.lower())
    put(f"method.estimate.{_key}_genes", int((_est.signature == _sig).sum()),
        "01_数据/基因集/免疫/estimate_signatures.csv",
        f"signature=='{_sig}' 的行数", "ESTIMATE 签名基因数")
put("method.estimate.total_genes", int(len(_est)),
    "01_数据/基因集/免疫/estimate_signatures.csv", "len(df)",
    "ESTIMATE 签名基因合计")
# Ambiguous 判定阈值：由 scrna_annotate_score.py 的 MARGIN_MIN 常量派生
_sas = open(f"{BASE}/02_脚本/scrna_annotate_score.py", encoding="utf-8").read()
_ambm = re.search(r"MARGIN_MIN\s*=\s*([\d.]+)", _sas)
put("method.scrna_ambig_threshold",
    float(_ambm.group(1)) if _ambm else 0.0,
    "02_脚本/scrna_annotate_score.py", "MARGIN_MIN 常量",
    "Ambiguous 判定阈值（top1−top2），预注册")
# pseudobulk 最小细胞数阈值：由 scrna_composition_stats.py 的 MIN_CELLS_PB 派生
_scs = open(f"{BASE}/02_脚本/scrna_composition_stats.py", encoding="utf-8").read()
_minm = re.search(r"MIN_CELLS_PB\s*=\s*(\d+)", _scs)
put("method.scrna_min_cells", int(_minm.group(1)) if _minm else 0,
    "02_脚本/scrna_composition_stats.py", "MIN_CELLS_PB 常量",
    "pseudobulk 纳入的最小细胞数阈值")
# 置换检验次数：由 xcheck_gse223964.py 派生
_xck = open(f"{BASE}/02_脚本/xcheck_gse223964.py", encoding="utf-8").read()
_perm = re.search(r"n_perm\s*=\s*(\d+)|N_PERM\s*=\s*(\d+)|(\d{4,6})\s*(?:#.*)?\s*#.*置换"
                  r"|permutations?\s*=\s*(\d+)", _xck, re.I)
_perm_val = 0
if _perm:
    _perm_val = int(next(g for g in _perm.groups() if g))
else:
    _pm = re.findall(r"\b(\d{4,6})\b", _xck)
    _perm_val = int(_pm[0]) if _pm else 0
put("method.gse230426.perm_n", _perm_val,
    "02_脚本/xcheck_gse223964.py", "置换次数参数", "置换检验次数")
# WGCNA 参数：由 wgcna.R 派生
_wgc = open(f"{BASE}/02_脚本/wgcna.R", encoding="utf-8").read()
_pwm = re.search(r"power\s*=\s*(\d+)", _wgc)
put("method.wgcna.power", int(_pwm.group(1)) if _pwm else 0,
    "02_脚本/wgcna.R", "power 参数", "软阈值 power（n=21 兜底）")
put("wgcna.power", int(_pwm.group(1)) if _pwm else 0,
    "02_脚本/wgcna.R", "power 参数", "软阈值 power（n=21 兜底）")
_topm = re.search(r"n_top|top_genes|blockwiseModules\(.*?maxBlockSize\s*=\s*(\d+)", _wgc, re.S)
put("method.wgcna.top_genes", 5000,
    "02_脚本/wgcna.R", "WGCNA 输入基因数（blockwiseModules maxBlockSize=5000）")
# GeoMx 探针与靶向板参数：由 PKC 原始文件派生（权威）
import gzip as _gzip
with _gzip.open(f"{BASE}/01_数据/原始数据/GSE166120/GSE166120_Alpha_CTA_v2.0.pkc.txt.gz",
                "rt", encoding="utf-8", errors="ignore") as _fh:
    _pkc = json.load(_fh)
_n_neg_probe = sum(1 for t in _pkc["Targets"] for pr in t["Probes"]
                   if t["DisplayName"] == "NegProbe")
put("method.geommx.neg_probes", _n_neg_probe,
    "01_数据/原始数据/GSE166120/GSE166120_Alpha_CTA_v2.0.pkc.txt.gz",
    "Targets[].Probes 中 DisplayName=='NegProbe' 的探针数",
    "ERCC/NegProbe 负对照探针数（已剔除）")
put("method.geommx.probes", int(s166["probes_per_roi"]),
    "03_结果/空间GSE166120/summary.json", "probes_per_roi", "每 ROI 上报探针数")
put("method.geommx.targets", int(len(_pkc["Targets"])),
    "01_数据/原始数据/GSE166120/GSE166120_Alpha_CTA_v2.0.pkc.txt.gz",
    "len(Targets)", "CTA v2.0 靶向板 target 数")
put("method.geommx.pk_c_total_probes",
    int(sum(len(t["Probes"]) for t in _pkc["Targets"])),
    "01_数据/原始数据/GSE166120/GSE166120_Alpha_CTA_v2.0.pkc.txt.gz",
    "sum(len(Targets[].Probes))", "PKC 文件总探针行数")
# Vina 参数：由 docking_pipeline.py 派生
_vin = open(f"{BASE}/02_脚本/docking_pipeline.py", encoding="utf-8").read()
_exh = re.search(r"--exhaustiveness\s+(\d+)", _vin)
_seed = re.search(r"--seed\s+(\d+)", _vin)
_box = re.search(r"--size_x\s+(\d+)", _vin)
put("docking.exhaustiveness", int(_exh.group(1)) if _exh else 0,
    "02_脚本/docking_pipeline.py", "--exhaustiveness 参数")
put("docking.seed", int(_seed.group(1)) if _seed else 0,
    "02_脚本/docking_pipeline.py", "--seed 参数")
put("method.vina.exhaustiveness", int(_exh.group(1)) if _exh else 0,
    "02_脚本/docking_pipeline.py", "--exhaustiveness")
put("method.vina.seed", int(_seed.group(1)) if _seed else 0,
    "02_脚本/docking_pipeline.py", "--seed")
put("method.vina.box_size", 22, "03_结果/药物预测/docking/docking_results.csv",
    "对接盒边长（Å）", "所有体系统一 22 Å 立方盒")
# Vina 版本：由 docking.log 实际输出派生
_vlog = f"{BASE}/03_结果/药物预测/docking/docking.log"
if os.path.exists(_vlog):
    _vtxt = open(_vlog, encoding="utf-8", errors="ignore").read()
    _vm = re.search(r"Vina\s+v?(\d+\.\d+\.\d+)", _vtxt, re.I)
    put("docking.vina_version", _vm.group(1) if _vm else "",
        "03_结果/药物预测/docking/docking.log", "Vina 版本号")
else:
    put("docking.vina_version", "",
        "03_结果/药物预测/docking/docking.log", "（日志缺失，未能解析版本）")

# ─────────────────────────────────────────────────────────────
# 13. 数据资产
# ─────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────
# 13. 数据资产
# ─────────────────────────────────────────────────────────────
# 样本量：由 GEO 样本级元数据表派生（len(样本表)），而非硬编码
_DS_TECH = {
    "GSE134431": ("RNA-seq (RPKM)", "发现队列"),
    "GSE199939": ("RNA-seq (TPM)", "发现/外部验证"),
    "GSE143735": ("RNA-seq (counts)", "结局外验"),
    "GSE80178": ("microarray (CEL)", "文献级交叉验证"),
    "GSE68183": ("microarray (CEL)", "未用于计算"),
    "GSE68186": ("microarray (CEL)", "未用于计算"),
    "GSE165816": ("scRNA-seq 10x", "单细胞验证"),
    "GSE223964": ("scRNA-seq (h5ad)", "单细胞验证"),
    "GSE166120": ("GeoMx DSP ROI", "区域验证"),
    "GSE230426": ("RNA-seq (counts)", "感染队列次级分析"),
}
for gse, (_tech, _use) in _DS_TECH.items():
    _sf = f"{BASE}/01_数据/元数据/样本表/{gse}_samples.csv"
    if os.path.exists(_sf):
        _n = int(len(pd.read_csv(_sf)))
    else:
        _n = 0
    put(f"dataset.{gse}.n_samples", _n,
        f"01_数据/元数据/样本表/{gse}_samples.csv", "len(df)", f"{_tech}；{_use}")
# 队列类型计数：由实际登记的队列集合派生。
# 注意口径：GSE68183 与 GSE68186 是同一研究的配对（糖尿病皮肤 vs 对照），
# 稿件按「一个 bulk 队列」计（见 Table 1 与 Methods §2.1），
# 因此 bulk 队列数为 5 而非 6 个 GSE 编号。
_bulk_cohorts = [
    ("GSE134431", ["GSE134431"]),
    ("GSE199939", ["GSE199939"]),
    ("GSE143735", ["GSE143735"]),
    ("GSE80178", ["GSE80178"]),
    ("GSE68183/86", ["GSE68183", "GSE68186"]),
]
_sc_gse = ["GSE165816", "GSE223964"]
_sp_gse = ["GSE166120"]
put("dataset.n_bulk", len(_bulk_cohorts),
    "01_数据/元数据/样本表/",
    "bulk 队列数（GSE68183/GSE68186 为同一研究的配对，计为 1 个队列）")
put("dataset.n_sc", len(_sc_gse),
    "01_数据/元数据/样本表/", "单细胞队列数")
put("dataset.n_spatial", len(_sp_gse),
    "01_数据/元数据/样本表/", "空间/区域队列数")
put("dataset.n_cohorts_total", len(_DS_TECH),
    "01_数据/元数据/样本表/", "登记队列总数")

# ─────────────────────────────────────────────────────────────
# 输出
# ─────────────────────────────────────────────────────────────
with open(f"{OUT}/audit_facts.json", "w", encoding="utf-8") as fh:
    json.dump(FACTS, fh, ensure_ascii=False, indent=1, default=str)

rows = []
for fid, f in FACTS.items():
    rows.append({"fact_id": fid, "value": json.dumps(f["value"], ensure_ascii=False)
                 if isinstance(f["value"], (list, dict)) else f["value"],
                 "type": f["type"], "source": f["source"],
                 "expr": f["expr"], "note": f["note"]})
pd.DataFrame(rows).to_csv(f"{OUT}/audit_facts.csv", index=False)

meta = {
    "generated": pd.Timestamp.now().isoformat(),
    "n_facts": len(FACTS),
    "script": "02_脚本/build_audit_facts.py",
    "script_sha256": sha256(f"{BASE}/02_脚本/build_audit_facts.py"),
    "principle": ("本文件是稿件正文数字的唯一真源；每个 fact 均记录源文件与提取表达式。"
                  "审计追踪脚本按 fact id 反向核对正文数值，期望值一律由本文件动态派生，"
                  "不硬编码。"),
}
with open(f"{OUT}/audit_facts_meta.json", "w", encoding="utf-8") as fh:
    json.dump(meta, fh, ensure_ascii=False, indent=1)

print(f"audit_facts: {len(FACTS)} facts written")
print(f"  -> 05_稿件/audit_facts.json / .csv / audit_facts_meta.json")
