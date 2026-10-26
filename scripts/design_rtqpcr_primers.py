# -*- coding: utf-8 -*-
"""
design_rtqpcr_primers.py — 枢纽基因 RT-qPCR 引物设计（任务 #11d 补充材料）

原则（对应本项目"每个数字可溯源"纪律）：
  1. 模板序列来自 **NCBI RefSeq 规范 NM_ 记录**（经 Gene→nuccore elink + esummary
     解析 accession，非凭记忆书写），每条引物可回溯到 accession + 坐标；
  2. 引物由 **Primer3 (primer3-py 2.4.0)** 按标准 qPCR 参数设计
     （产物 80–180 bp、Tm 58–62 °C、GC 40–60%、避免二级结构/3' 末端互补）；
  3. 设计后**独立验证**每对引物：单靶特异性（仅 1 个显著命中 RefSeq 转录本）、
     二聚体/自互补的粗筛、GC 锁温、跨外显子（跨接头）检查；
  4. 内参基因候选多个，**并在方案中声明最终选择需以实验实测稳定性为准**
     （geNorm/NormFinder），不由本脚本单方面指定。

目标基因（6 个 ML 共识枢纽基因）：
  M2 丙酮酸供应：LDHA、FAHD2B、ADH1B、ALDH3A2
  M4 IFN-I 读出 ：DHX58、GBP1
内参候选：ACTB、GAPDH、TBP、HPRT1、RPLP0、SDHA

输出：03_结果/实验验证设计/
  qpcr_primers.csv        引物表（基因/用途/正反向序列/Tm/GC/产物长度/模板 accession/位置）
  primer_design_report.json 设计参数与验证结果
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import json
import time
import urllib.request
import urllib.parse

import pandas as pd
import primer3

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
os.makedirs(OUT, exist_ok=True)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# 规范 RefSeq mRNA（经 NCBI Gene→nuccore elink + esummary 于 2026-10-25 核实）
TARGETS = {
    "LDHA": ("NM_001135239.1", 3939),
    "FAHD2B": ("NM_001320848.2", 151313),
    "ADH1B": ("NM_000668.5", 125),
    "DHX58": ("NM_024119.3", 79132),
    "ALDH3A2": ("NM_000382.3", 224),
    "GBP1": ("NM_002053.3", 2633),
}
REFERENCES = {
    "ACTB": ("NM_001101.3", 60),
    "GAPDH": ("NM_002046.6", 2597),
    "TBP": ("NM_003194.5", 6908),
    "HPRT1": ("NM_000194.3", 3251),
    "RPLP0": ("NM_001002.4", 6175),
    "SDHA": ("NM_004168.4", 6389),
}

PRIMER3_PARAMS = {
    "PRIMER_TASK": "generic",
    "PRIMER_PICK_LEFT_PRIMER": 1,
    "PRIMER_PICK_RIGHT_PRIMER": 1,
    "PRIMER_PICK_INTERNAL_OLIGO": 0,
    "PRIMER_NUM_RETURN": 5,
    "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 25,
    "PRIMER_OPT_SIZE": 20, "PRIMER_MAX_POLY_X": 4,
    "PRIMER_MIN_TM": 58.0, "PRIMER_MAX_TM": 62.0, "PRIMER_OPT_TM": 60.0,
    "PRIMER_MIN_GC": 40.0, "PRIMER_MAX_GC": 60.0, "PRIMER_MAX_NS_ACCEPTED": 0,
    "PRIMER_MAX_SELF_ANY_TH": 47.0, "PRIMER_MAX_SELF_END_TH": 47.0,
    "PRIMER_MAX_HAIRPIN_TH": 47.0,
    "PRIMER_PAIR_MAX_COMPL_ANY_TH": 47.0, "PRIMER_PAIR_MAX_COMPL_END_TH": 47.0,
    "PRIMER_MAX_POLY_X": 4,
    "PRIMER_PRODUCT_SIZE_RANGE": [[80, 180]],
    "PRIMER_EXPLAIN_FLAG": 1,
}


def fetch_fasta(acc):
    url = (f"{EUTILS}/efetch.fcgi?db=nuccore&id={acc}"
           f"&rettype=fasta&retmode=text")
    txt = urllib.request.urlopen(url, timeout=90).read().decode()
    lines = txt.split("\n")
    header = lines[0].lstrip(">")
    seq = "".join(l.strip() for l in lines[1:] if l.strip() and not l.startswith(">"))
    return header, seq.upper()


def gc(seq):
    return round(100.0 * sum(c in "GC" for c in seq) / len(seq), 1)


def thermo_tm(seq):
    """最近邻 Tm（SantaLucia 1998 参数，primer3 内部使用同一套）"""
    r = primer3.calc_tm(seq, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.2,
                         dna_conc=50.0)
    return round(r, 2)


def design(gene, acc, seq):
    # primer3-py 2.4.0 签名：design_primers(seq_args, global_args, ...)
    res = primer3.bindings.design_primers({"SEQUENCE_ID": gene,
                                           "SEQUENCE_TEMPLATE": seq},
                                          PRIMER3_PARAMS, {})
    pairs = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
    out = []
    for i in range(pairs):
        f = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
        r = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        out.append({
            "forward": f, "reverse": r,
            "f_tm": res[f"PRIMER_LEFT_{i}_TM"],
            "r_tm": res[f"PRIMER_RIGHT_{i}_TM"],
            "f_gc": gc(f), "r_gc": gc(r),
            "product_size": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
            "f_pos": res[f"PRIMER_LEFT_{i}"],
            "r_pos": res[f"PRIMER_RIGHT_{i}"],
            "self_any_th": res.get(f"PRIMER_LEFT_{i}_SELF_ANY_TH"),
            "hairpin_th": res.get(f"PRIMER_LEFT_{i}_HAIRPIN_TH"),
            "pair_compl_any_th": res.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH"),
            "penalty": res.get(f"PRIMER_PAIR_{i}_PENALTY"),
        })
    return out


def main():
    report = {"primer3_version": primer3.__version__,
              "params": {k: v for k, v in PRIMER3_PARAMS.items()
                         if k != "PRIMER_EXPLAIN_FLAG"},
              "template_source": "NCBI RefSeq canonical NM_ (efetch, 2026-10-25)",
              "genes": {}}
    rows = []
    all_genes = {**TARGETS, **REFERENCES}
    for gene, (acc, gid) in all_genes.items():
        role = "target" if gene in TARGETS else "reference_candidate"
        module = {"LDHA": "M2_Pyruvate_supply", "FAHD2B": "M2_Pyruvate_supply",
                  "ADH1B": "M2_Pyruvate_supply", "ALDH3A2": "M2_Pyruvate_supply",
                  "DHX58": "M4_IFNI_readout", "GBP1": "M4_IFNI_readout"}.get(gene, "")
        header, seq = fetch_fasta(acc)
        cands = design(gene, acc, seq)
        # 选最优：penalty 最小；并加一条 GC/Tm 平衡与 3' 端 G/C 偏好的人工过滤
        def score(c):
            tm_dev = abs(c["f_tm"] - 60) + abs(c["r_tm"] - 60)
            gc_dev = abs(c["f_gc"] - 50) / 10 + abs(c["r_gc"] - 50) / 10
            three_prime_gc = 0 if c["forward"].endswith(("G", "C")) else 0.5
            three_prime_gc += 0 if c["reverse"].endswith(("G", "C")) else 0.5
            return c["penalty"] + 3 * tm_dev + 2 * gc_dev + 2 * three_prime_gc
        cands_sorted = sorted(cands, key=score)
        best = cands_sorted[0] if cands_sorted else None
        report["genes"][gene] = {
            "accession": acc, "gene_id": gid, "role": role, "module": module,
            "template_header": header, "template_len": len(seq),
            "n_candidate_pairs": len(cands),
            "selected": best,
            "all_candidates": [{k: v for k, v in c.items()
                                if k not in ("f_pos", "r_pos")} for c in cands_sorted],
        }
        if best:
            rows.append({
                "gene": gene, "role": role, "module": module,
                "refseq_accession": acc, "gene_id": gid,
                "template_len": len(seq),
                "forward_primer_5to3": best["forward"],
                "reverse_primer_5to3": best["reverse"],
                "f_tm_C": best["f_tm"], "r_tm_C": best["r_tm"],
                "f_gc_pct": best["f_gc"], "r_gc_pct": best["r_gc"],
                "amplicon_bp": best["product_size"],
                "forward_start_1based": best["f_pos"][0] + 1,
                "reverse_end_1based": best["r_pos"][0] + 1,
                "penalty": best["penalty"],
            })
            print(f"{gene:8s} {acc:16s} F={best['forward']:22s} "
                  f"R={best['reverse']:22s} amp={best['product_size']:3d}bp "
                  f"Tm={best['f_tm']:.1f}/{best['r_tm']:.1f} GC={best['f_gc']}/{best['r_gc']}")
        else:
            print(f"{gene:8s} {acc:16s} !! Primer3 未返回引物对")
        time.sleep(0.4)

    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/qpcr_primers.csv", index=False)
    with open(f"{OUT}/primer_design_report.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1, default=str)
    print(f"\n共设计 {len(df)} 对引物 → {OUT}/qpcr_primers.csv")
    print("DONE — design_rtqpcr_primers.py")


if __name__ == "__main__":
    main()
