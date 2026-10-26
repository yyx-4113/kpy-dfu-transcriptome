# -*- coding: utf-8 -*-
"""
test_intersection_independence.py — 三层交集"全为糖酵解酶"的独立性检验

═══ 为什么必须做这个检验 ═══
稿件 §3.2 称"三层交集 = 6 个基因，全部为糖酵解酶"，并以此作为结构证据。
但 KpyRN 的 M2 模块本身就是 KEGG hsa00010（糖酵解/糖异生）+ hsa00620（丙酮酸代谢）
共 90 个基因——**交集必然落在糖酵解酶上**。若不做独立性检验，"all glycolytic enzymes"
就是把柄：审稿人一句"这是先验注定的"即可拆掉该节。

═══ 三项检验 ═══
T1  剔除 M2 后重跑三层交集（剩 M1+M3+M4 = 38 基因），看是否仍有非糖酵解基因汇聚
T2  超几何检验：给定交集大小 k，k 个全部落在糖酵解子集内的概率
    （零假设：KpyRN 成员在 DFU 相关池中随机分布）
T3  结构置换检验（10000 次，固定种子）：从 WGCNA 5000 基因背景随机抽 128 个作伪 KpyRN，
    与同一池求交，得交集大小零分布 → 判断真实交集 6 是否超出随机
    * T3 同时报告：零分布中"交集全为糖酵解"的频率，用于判断"全糖酵解"是否为额外发现

判据（预注册于本文件）：
  - T3 交集大小 p<0.05  → KpyRN 与 DFU 相关层存在非随机关联（正面结论成立）
  - T3 "全糖酵解"频率 ≥0.05 → "all glycolytic" 不构成独立证据，§3.2/§4.1 措辞必须降级
  - T1 交集为空 → 必须明确写出"三层结构证据完全由 M2 模块驱动"

输出：03_结果/WGCNA/intersection_independence.json
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import csv
import json
import random
from math import comb

BASE = _ROOT
OUT = f"{BASE}/03_结果"
WGCNA = f"{OUT}/WGCNA"

SEED = 20261025
N_PERM = 10000
PADJ_CUT = 0.05

# 模块编号→颜色：R WGCNA labels2colors(1:12) 顺序（与 build_audit_facts.py 的 L2C 同源）
# module_trait_cor.csv 的 color 列存的是**编号数字串**而非颜色名，直接用会让池为空。
L2C = {1: "turquoise", 2: "blue", 3: "brown", 4: "yellow", 5: "green",
       6: "red", 7: "black", 8: "pink", 9: "magenta", 10: "purple"}


def read_csv(p, enc="utf-8-sig"):
    return list(csv.DictReader(open(p, encoding=enc)))


def main():
    # ── 层 1：META1 显著 ──
    meta = read_csv(f"{OUT}/meta/META1_DFU_vs_control.csv")
    sig = {r["gene"] for r in meta if float(r["padj"]) < PADJ_CUT}
    print(f"[L1] META1 padj<0.05: {len(sig)} / {len(meta)}")

    # ── 层 2：DFU 显著相关模块 ──
    trait = read_csv(f"{WGCNA}/module_trait_cor.csv")
    sig_nums = [int(r["module"]) for r in trait if float(r["p_DFU"]) < PADJ_CUT]
    sig_colors = {L2C[n] for n in sig_nums if n in L2C}
    print(f"[L2] DFU 显著相关模块: 编号 {sorted(sig_nums)} → 颜色 {sorted(sig_colors)}")
    for r in trait:
        n = int(r["module"])
        if n in sig_nums:
            print(f"     module {n} = {L2C.get(n)}: cor_DFU={float(r['cor_DFU']):+.3f}, "
                  f"p={float(r['p_DFU']):.2e}")

    mem = read_csv(f"{WGCNA}/gene_module_membership.csv")
    mod2gene = {r["gene"]: r["module"] for r in mem}
    universe = set(mod2gene)                      # WGCNA 5000 基因背景
    in_sig_modules = {g for g, c in mod2gene.items() if c in sig_colors}
    print(f"[L2] 背景 {len(universe)} 基因；落入 DFU 显著模块: {len(in_sig_modules)}")

    # ── 层 3：KpyRN ──
    kpy_rows = read_csv(f"{BASE}/01_数据/基因集/KpyRN_modules.csv")
    kpy = {r["gene_symbol"] for r in kpy_rows if r["scoring_include"] == "Yes"}
    kpy_module = {r["gene_symbol"]: r["module"] for r in kpy_rows}
    kpy_sub = {r["gene_symbol"]: r["submodule"] for r in kpy_rows}
    # 糖酵解标签：客观取自 KEGG hsa00010 子模块标注（非人工判定）
    glycolysis = {r["gene_symbol"] for r in kpy_rows if "hsa00010" in r["submodule"]}
    print(f"[L3] KpyRN 评分基因 {len(kpy)}；其中 hsa00010(糖酵解) 标注: {len(glycolysis)}")

    # ── 交集池 S = L1 ∩ L2 ──
    S = sig & in_sig_modules
    print(f"[S ] 池 S = META1显著 ∩ DFU显著模块 = {len(S)}")

    three = sorted(kpy & S)
    print(f"\n[三层交集] k = {len(three)}: {three}")
    print(f"  其中糖酵解: {sorted(set(three) & glycolysis)}  "
          f"({len(set(three) & glycolysis)}/{len(three)})")

    # ══════ T1：剔除 M2 后重跑 ══════
    kpy_noM2 = {g for g in kpy if kpy_module.get(g) != "M2_Pyruvate_supply"}
    three_noM2 = sorted(kpy_noM2 & S)
    print(f"\n[T1] 剔除 M2 后 KpyRN = {len(kpy_noM2)} 基因；三层交集 = "
          f"{len(three_noM2)}: {three_noM2}")

    # ══════ T2：超几何精确检验（按"背景内子集"逐个检验，避免抽样基数错配）══════
    # ★ 关键约束：WGCNA 只用了 5000 个高变基因，KpyRN 128 个中仅 56 个在该背景内。
    #   若直接以 128 为抽样数做置换，会把"不在背景内"误算成"未富集"。
    #   故每个子集一律以「该子集 ∩ universe」为抽样数 n，以 S 为成功态做超几何检验。
    N_uni, K_pool = len(universe), len(S & universe)

    def hyper_ge(n, k):
        """P(X >= k)，X ~ Hypergeom(N=背景, K=池, n=抽样数)"""
        if n == 0:
            return float("nan")
        lo = min(n, K_pool)
        return sum(comb(K_pool, x) * comb(N_uni - K_pool, n - x)
                   for x in range(k, lo + 1)) / comb(N_uni, n)

    def hyper_le(n, k):
        """P(X <= k)"""
        if n == 0:
            return float("nan")
        return sum(comb(K_pool, x) * comb(N_uni - K_pool, n - x)
                   for x in range(0, k + 1)) / comb(N_uni, n)

    mod_groups = {
        "M1_Kpy_machinery": {g for g in kpy if kpy_module.get(g) == "M1_Kpy_machinery"},
        "M2_Pyruvate_supply": {g for g in kpy if kpy_module.get(g) == "M2_Pyruvate_supply"},
        "M3_Kpy_targets": {g for g in kpy if kpy_module.get(g) == "M3_Kpy_targets"},
        "M4_IFNI_readout": {g for g in kpy if kpy_module.get(g) == "M4_IFNI_readout"},
    }
    print(f"\n[T2] 背景 universe={N_uni}，池 S={K_pool}")
    print(f"{'subset':22s} {'total':>5s} {'inUni':>5s} {'obs_inS':>7s} "
          f"{'exp':>6s} {'ratio':>6s} {'p':>10s}")
    t2 = {}
    subsets = {"KpyRN_all": kpy, "glycolysis_hsa00010": glycolysis, **mod_groups}
    for name, gs in subsets.items():
        n = len(gs & universe)
        obs = len(gs & S)
        exp = n * K_pool / N_uni
        p = hyper_ge(n, obs) if obs > exp else hyper_le(n, obs)
        t2[name] = {"total": len(gs), "in_universe": n, "observed_in_pool": obs,
                    "expected": round(exp, 2),
                    "ratio_obs_exp": round(obs / exp, 3) if exp else None,
                    "p_two_sided_lesser_tail": p,
                    "direction": "enriched" if obs > exp else "depleted"}
        print(f"{name:22s} {len(gs):5d} {n:5d} {obs:7d} {exp:6.2f} "
              f"{obs/exp if exp else 0:6.2f} {p:10.4f}")

    # 非 M2 合并子集（T1 的直接检验）
    nonM2 = {g for g in kpy if kpy_module.get(g) != "M2_Pyruvate_supply"}
    n_nm2, obs_nm2 = len(nonM2 & universe), len(nonM2 & S)
    p_nm2 = hyper_le(n_nm2, obs_nm2)
    t2["nonM2_combined"] = {"total": len(nonM2), "in_universe": n_nm2,
                            "observed_in_pool": obs_nm2,
                            "expected": round(n_nm2 * K_pool / N_uni, 2),
                            "p_depleted": p_nm2}
    print(f"{'nonM2_combined':22s} {len(nonM2):5d} {n_nm2:5d} {obs_nm2:7d} "
          f"{n_nm2*K_pool/N_uni:6.2f} {0:6.2f} {p_nm2:10.4f}")

    # ══════ T3：结构置换（抽样数 = 各子集在背景内的大小，交叉验证 T2）══════
    rng = random.Random(SEED)
    uni = sorted(universe)
    perm_note = {}
    for name, gs in [("KpyRN_all", kpy), ("glycolysis_hsa00010", glycolysis),
                     ("nonM2_combined", nonM2)]:
        n = len(gs & universe)
        obs = len(gs & S)
        cnt = 0
        for _ in range(N_PERM):
            if len(set(rng.sample(uni, n)) & S) >= obs:
                cnt += 1
        perm_note[name] = {"n": n, "obs": obs,
                           "p_perm_ge": (cnt + 1) / (N_PERM + 1)}
        print(f"[T3] {name:22s} n={n:3d} obs={obs:2d} "
              f"p_perm(>=obs)={(cnt+1)/(N_PERM+1):.4f}")

    # ══════ 结论判定 ══════
    verdict = {
        "T1_three_layer_driven_by_M2_only": len(three_noM2) == 0,
        "T2_KpyRN_whole_enriched": t2["KpyRN_all"]["direction"] == "enriched"
        and t2["KpyRN_all"]["p_two_sided_lesser_tail"] < 0.05,
        "T2_glycolysis_enriched":
            t2["glycolysis_hsa00010"]["direction"] == "enriched"
            and t2["glycolysis_hsa00010"]["p_two_sided_lesser_tail"] < 0.05,
        "T2_M4_depleted": t2["M4_IFNI_readout"]["direction"] == "depleted"
        and t2["M4_IFNI_readout"]["p_two_sided_lesser_tail"] < 0.05,
        "all_glycolytic_is_independent_evidence":
            t2["glycolysis_hsa00010"]["direction"] == "enriched"
            and t2["glycolysis_hsa00010"]["p_two_sided_lesser_tail"] < 0.05,
    }
    print("\n[判定]")
    for kk, vv in verdict.items():
        print(f"  {kk} = {vv}")

    out = {
        "script": "02_脚本/test_intersection_independence.py",
        "seed": SEED, "n_perm": N_PERM, "padj_cut": PADJ_CUT,
        "design_note": ("WGCNA 仅使用 5000 个高变基因，故每个子集的抽样数一律取"
                        "「该子集 ∩ universe」的大小；若误用 128 会把背景外基因"
                        "算作未富集，导致假阴性/假阳性。"),
        "layers": {
            "L1_META1_sig": len(sig),
            "L2_sig_module_colors": sorted(sig_colors),
            "L2_sig_module_numbers": sorted(sig_nums),
            "L2_genes_in_sig_modules": len(in_sig_modules),
            "L3_kpyrn_scored": len(kpy),
            "L3_kpyrn_in_universe": len(kpy & universe),
            "L3_glycolysis_hsa00010": len(glycolysis),
            "L3_glycolysis_in_universe": len(glycolysis & universe),
            "pool_S_size": K_pool,
            "wgcna_universe": N_uni,
            "pool_fraction_of_universe": round(K_pool / N_uni, 4),
        },
        "three_layer": {"genes": three, "k": len(three),
                        "n_glycolysis": len(set(three) & glycolysis),
                        "module_of_each": {g: kpy_module.get(g) for g in three},
                        "submodule_of_each": {g: kpy_sub.get(g) for g in three}},
        "T1_exclude_M2": {"kpyrn_size_without_M2": len(kpy_noM2),
                          "kpyrn_without_M2_in_universe": len(kpy_noM2 & universe),
                          "three_layer_genes": three_noM2, "k": len(three_noM2)},
        "T2_hypergeometric": t2,
        "T3_permutation": perm_note,
        "verdict": verdict,
    }
    os.makedirs(WGCNA, exist_ok=True)
    with open(f"{WGCNA}/intersection_independence.json", "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(f"\n-> {WGCNA}/intersection_independence.json")
    print("DONE — test_intersection_independence.py")


if __name__ == "__main__":
    main()
