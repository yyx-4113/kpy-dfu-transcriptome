# -*- coding: utf-8 -*-
"""
compute_qpcr_sample_size.py — RT-qPCR 验证的样本量计算（任务 #11d）

原则：样本量**不得拍脑袋**。本脚本的效应量直接取自本项目已完成的 meta 分析产物
`03_结果/meta/META1_DFU_vs_control.csv` 的 pooled_logFC 与 pooled_se，
Cohen's d = pooled_logFC / pooled_se；并同时给出「保守情景」（效应量 ×0.5），
以覆盖 bulk RNA-seq 效应量在 RT-qPCR 平台上的不确定性。

设计：
  - 配对（DFU 溃疡组织 vs 同一患者配对非溃疡皮肤）→ 配对 t 检验
  - 独立（DFU vs 健康对照皮肤）→ 两样本 t 检验
α=0.05 双侧，power=0.80

输出：03_结果/实验验证设计/sample_size_calculation.csv / .json
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import json

import numpy as np
import pandas as pd
from scipy import stats

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
os.makedirs(OUT, exist_ok=True)

HUB = ["LDHA", "FAHD2B", "ADH1B", "DHX58", "ALDH3A2", "GBP1"]
ALPHA, POWER = 0.05, 0.80


def n_paired_t(d, alpha=ALPHA, power=POWER):
    za = stats.norm.ppf(1 - alpha / 2)
    zb = stats.norm.ppf(power)
    return int(np.ceil(((za + zb) / d) ** 2 + (za + zb) ** 2 / 2))


def n_indep_t(d, alpha=ALPHA, power=POWER):
    za = stats.norm.ppf(1 - alpha / 2)
    zb = stats.norm.ppf(power)
    return int(np.ceil(2 * ((za + zb) / d) ** 2))


def main():
    m1 = pd.read_csv(f"{BASE}/03_结果/meta/META1_DFU_vs_control.csv")
    h = m1[m1.gene.isin(HUB)].copy()
    h["cohen_d"] = h.pooled_logFC / h.pooled_se
    rows = []
    for _, r in h.iterrows():
        d = abs(r.cohen_d)
        d_half = d * 0.5
        rows.append({
            "gene": r.gene,
            "META1_pooled_logFC": round(r.pooled_logFC, 4),
            "META1_pooled_se": round(r.pooled_se, 4),
            "META1_z": round(r.z, 4),
            "META1_padj": r.padj,
            "META1_I2_pct": round(r.I2, 2),
            "n_cohorts_pos": r.n_cohorts_pos,
            "n_cohorts_neg": r.n_cohorts_neg,
            "cohen_d": round(r.cohen_d, 3),
            "n_paired_primary": n_paired_t(d),
            "n_per_group_independent": n_indep_t(d),
            "n_paired_conservative_50pct_effect": n_paired_t(d_half),
            "n_per_group_independent_conservative_50pct_effect": n_indep_t(d_half),
            "heterogeneity_flag": ("HIGH (I2>80) — pooled effect unstable, "
                                   "replication failure would NOT falsify the "
                                   "core claim" if r.I2 > 80 else
                                   "LOW/MODERATE — reliable verification target"),
        })
    df = pd.DataFrame(rows).sort_values("cohen_d", key=lambda s: s.abs(),
                                         ascending=False)
    df.to_csv(f"{OUT}/sample_size_calculation.csv", index=False)

    # 目标样本量（计划书要求 n>=8/组）
    need_paired = int(df.n_paired_conservative_50pct_effect.max())
    need_indep = int(df.n_per_group_independent_conservative_50pct_effect.max())
    target = max(8, need_paired, need_indep)
    recruit = int(np.ceil(target * 1.2))       # 预留 20% 脱落/技术失败

    summary = {
        "date": "2026-10-25",
        "effect_size_source": "03_结果/meta/META1_DFU_vs_control.csv "
                              "(pooled_logFC / pooled_se)",
        "alpha": ALPHA, "power": POWER, "test": "two-sided",
        "n_paired_required_worst_gene_conservative": need_paired,
        "n_per_group_independent_required_worst_gene_conservative": need_indep,
        "target_n_per_group": target,
        "target_n_per_group_with_20pct_attrition": recruit,
        "design": {
            "DFU_ulcer_tissue": recruit,
            "paired_non_ulcer_skin_same_patient": recruit,
            "healthy_control_skin": max(8, int(np.ceil(8 * 1.2))),
            "paired_peripheral_blood": recruit,
            "total_tissue_samples": recruit * 2 + max(8, int(np.ceil(8 * 1.2))),
        },
        "caveats": [
            "effect sizes are RNA-seq pooled (log2) and not directly comparable "
            "to RT-qPCR absolute expression; only direction is pre-specified",
            "three genes have I2>87% (ADH1B, ALDH3A2, FAHD2B) -> pooled effects "
            "unstable; failure to replicate them does NOT falsify the core claim",
            "no assumption about Kpy modification itself: no Kpy antibody or MS "
            "is available, so only upstream pyruvate supply and downstream "
            "IFN-I readout can be validated",
        ],
        "per_gene": rows,
    }
    with open(f"{OUT}/sample_size_calculation.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1, default=str)

    print("== 样本量计算（效应量来自 META1_DFU_vs_control.csv）==")
    print(df[["gene", "META1_pooled_logFC", "cohen_d", "META1_I2_pct",
              "n_paired_primary", "n_paired_conservative_50pct_effect",
              "n_per_group_independent_conservative_50pct_effect"]]
          .to_string(index=False))
    print(f"\n最保守基因所需: 配对 n={need_paired} | 独立 n/组={need_indep}")
    print(f"目标 n/组 = {target}（满足计划书 n>=8）；招募 {recruit}（含 20% 脱落）")
    print(f"组织样本合计: DFU {recruit} + 配对非溃疡 {recruit} + "
          f"健康对照 {max(8, int(np.ceil(8*1.2)))} = "
          f"{summary['design']['total_tissue_samples']}")
    print("DONE — compute_qpcr_sample_size.py")


if __name__ == "__main__":
    main()
