# -*- coding: utf-8 -*-
"""
verify_primers.py — 引物独立验证（v2，2026-10-25 修订）

═══ v1 → v2 的修订（均为实测发现的缺陷）═══

**V1 单靶特异性**：v1 用 NCBI E-utilities `esearch` 精确序列检索
（`"seq"[Sequence] AND srcdb_refseq[PROP]`）取命中基因集。该路径存在
**命中列表截断**（retmax=50）且检索语义为"精确序列出现在转录本中"，
无法区分全长匹配与部分匹配、无法评估 3' 端锚定。
→ v2 改为**直接读取 `primer_blast_v5.csv`**（BLAST v5 全库复核结果：
人源物种限定 + XML 完整描述 + 全长匹配 + 3' 端锚定分层），
V1_specific 定义为：`P0 == 0` 且 `target_full_match_refseq_tx ≥ 1` 且
`target_max_identity ≥ 90%`。若 v5 结果文件不存在则**显式报错**，
不再静默回退到有缺陷的 esearch 路径。

**V4 跨外显子**：v1 用 RefSeq **mRNA** 记录的 `exon` 特征判定，而该坐标是
**mRNA 自身坐标**（已拼接），相邻 exon 连续 → 判定在数学上不成立。
→ v2 改为基于 **RefSeqGene（NG_）基因组坐标**的判定：
引物在基因组上的落点若横跨该转录本的外显子边界，则 gDNA 亦可扩增
（cDNA 不受影响）→ 作为**信息性指标**记录并在方案中声明；
不作为否决项（RT-qPCR 跨接头引物本身是允许且常用的策略）。

**V2 二聚体**：不变（primer3.bindings 独立复算路径，与设计阶段无关）。
**V3 GC 锁温**：不变（MIQE 风格：GC 差 ≤10%、Tm 差 ≤2 °C）。

输出：03_结果/实验验证设计/primer_verification.csv + .json
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import json

import pandas as pd
import primer3

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
PROV = f"{OUT}/provenance"
SEC_THRESHOLD = 47.0

# 基因 → RefSeqGene（基因组坐标，用于 V4 跨接头判定）
NG_BY_GENE = {
    "LDHA": "NG_008185.1",
    "ALDH3A2": "NG_007095.2",
    "SDHA": "NG_012339.1",
}


def gc(seq):
    return round(100.0 * sum(c in "GC" for c in seq) / len(seq), 1)


def load_v5():
    p = f"{OUT}/primer_blast_v5.csv"
    if not os.path.exists(p):
        raise FileNotFoundError(
            f"缺少 BLAST v5 复核结果：{p}\n"
            "V1 单靶特异性必须由 BLAST v5 提供（v1 的 esearch 路径存在"
            "命中列表截断且无法评估全长匹配与 3' 端锚定）。\n"
            "请先运行：python 02_脚本/verify_primers_blast_v5.py")
    return pd.read_csv(p).set_index("gene")


def load_genome_anchor():
    p = f"{PROV}/refseqgene_records.json"
    if not os.path.exists(p):
        return {}
    return json.load(open(p, encoding="utf-8"))


def locate_genome(gseq, primer):
    rc_map = str.maketrans("ACGT", "TGCA")
    hits = []
    for s, strand in ((primer, "+"), (primer.translate(rc_map)[::-1], "-")):
        st = 0
        while True:
            i = gseq.find(s, st)
            if i < 0:
                break
            hits.append((i, i + len(s), strand))
            st = i + 1
    return hits


def main():
    df = pd.read_csv(f"{OUT}/qpcr_primers.csv")
    v5 = load_v5()
    genome = load_genome_anchor()

    report = {"verification_date": "2026-10-25",
              "primer3_version": primer3.__version__,
              "revision": "v2（2026-10-25）：V1 改由 BLAST v5 供数；"
                          "V4 改用 RefSeqGene 基因组坐标",
              "method": {
                  "V1": "BLAST v5 全库复核（core_nt + txid9606 人源限定 + XML）："
                        "specific = (P0 == 0) 且 目标基因存在全长匹配的人源 "
                        "RefSeq 转录本 且 目标位点 identity ≥ 90%",
                  "V2": "primer3.bindings calc_hairpin / calc_homodimer / "
                        "calc_heterodimer 独立复算（与设计阶段参数无关的路径）",
                  "V3": "MIQE 风格 GC 锁温（差 ≤10%）与 Tm 匹配（差 ≤2 °C）",
                  "V4": "RefSeqGene（NG_）基因组坐标下的外显子边界判定："
                        "引物落点若横跨该转录本外显子边界则 gDNA 可扩增"
                        "（信息性指标，不作否决项）",
              },
              "results": []}
    rows = []
    for _, r in df.iterrows():
        gene = r["gene"]
        f, rev = r["forward_primer_5to3"], r["reverse_primer_5to3"]
        # ---- V1 单靶特异性（来自 BLAST v5）----
        if gene not in v5.index:
            raise KeyError(f"BLAST v5 结果缺少基因 {gene}")
        v = v5.loc[gene]
        p0 = int(v.get("forward_n_P0", 0)) + int(v.get("reverse_n_P0", 0))
        p1 = int(v.get("forward_n_P1", 0)) + int(v.get("reverse_n_P1", 0))
        ft = int(v.get("forward_target_full_match_refseq_tx", 0))
        rt = int(v.get("reverse_target_full_match_refseq_tx", 0))
        fi = float(v.get("forward_target_max_identity_pct", 0))
        ri = float(v.get("reverse_target_max_identity_pct", 0))
        specific = bool(p0 == 0 and ft >= 1 and rt >= 1 and fi >= 90 and ri >= 90)
        # ---- V2 二聚体独立复算 ----
        hp_f = primer3.bindings.calc_hairpin(f)
        hd_f = primer3.bindings.calc_homodimer(f)
        hp_r = primer3.bindings.calc_hairpin(rev)
        hd_r = primer3.bindings.calc_homodimer(rev)
        het = primer3.bindings.calc_heterodimer(f, rev)
        max_secondary = max(hp_f.tm, hp_r.tm, hd_f.tm, hd_r.tm, het.tm)
        dimer_ok = max_secondary < SEC_THRESHOLD
        # ---- V3 GC 锁温 / Tm 匹配 ----
        gc_diff = abs(r["f_gc_pct"] - r["r_gc_pct"])
        tm_diff = abs(r["f_tm_C"] - r["r_tm_C"])
        gc_tm_ok = (gc_diff <= 10.0) and (tm_diff <= 2.0)
        # ---- V4 跨外显子（基因组坐标；信息性）----
        v4_note = "no RefSeqGene record available"
        v4_span = ""
        ng = genome.get(gene)
        if ng:
            gseq = ng.get("genome_seq", "")
            joins = ng.get("mrna_joins", [])
            if gseq and joins:
                segs = joins[0]["segments"]
                spans = []
                for primer, lab in ((f, "F"), (rev, "R")):
                    hit_spans = []
                    for s0, s1, strand in locate_genome(gseq, primer):
                        crossed = False
                        for k in range(len(segs) - 1):
                            a_end = segs[k][1]
                            b_start = segs[k + 1][0]
                            if s0 < b_start - 1 and s1 > a_end:
                                crossed = True
                        hit_spans.append((strand, crossed))
                    spans.append((lab, hit_spans))
                any_span = any(any(c for _, c in hs) for _, hs in spans)
                v4_span = any_span
                v4_note = (f"genome hits: "
                           + "; ".join(f"{lab}=" + ",".join(
                               f"{st}{'*' if cr else ''}"
                               for st, cr in hs) for lab, hs in spans))
            else:
                v4_note = "RefSeqGene record fetched but no genome sequence/join parsed"
        rec = {
            "gene": gene, "role": r["role"], "module": r["module"],
            "refseq_accession": r["refseq_accession"],
            "forward_primer_5to3": f, "reverse_primer_5to3": rev,
            "amplicon_bp": int(r["amplicon_bp"]),
            "V1_source": "BLAST v5 (core_nt, human-restricted, full-length + 3' anchor)",
            "V1_n_P0_offtarget": p0, "V1_n_P1_offtarget": p1,
            "V1_target_full_match_refseq_tx": ft + rt,
            "V1_target_max_identity_pct": max(fi, ri),
            "V1_P0_genes": (str(v.get("forward_P0_genes", "")) + ";"
                            + str(v.get("reverse_P0_genes", ""))).strip(";"),
            "V1_specific": specific,
            "V2_hairpin_f_tm": round(hp_f.tm, 2),
            "V2_hairpin_r_tm": round(hp_r.tm, 2),
            "V2_homodimer_f_tm": round(hd_f.tm, 2),
            "V2_homodimer_r_tm": round(hd_r.tm, 2),
            "V2_heterodimer_tm": round(het.tm, 2),
            "V2_max_secondary_structure_tm": round(max_secondary, 2),
            "V2_no_strong_secondary_structure": bool(dimer_ok),
            "V3_gc_diff_pct": gc_diff, "V3_tm_diff_C": round(tm_diff, 2),
            "V3_gc_locked_and_tm_matched": bool(gc_tm_ok),
            "V4_spans_exon_junction": v4_span,
            "V4_note": v4_note,
        }
        rows.append(rec)
        report["results"].append(rec)
        print(f"{gene:9s} V1_specific={specific} P0={p0} P1={p1} "
              f"tgtTx={ft + rt} V2_maxTm={max_secondary:.1f} "
              f"V3_gcDiff={gc_diff:.0f} tmDiff={tm_diff:.1f}")

    out = pd.DataFrame(rows)
    out["overall_pass"] = (out["V1_specific"]
                           & out["V2_no_strong_secondary_structure"]
                           & out["V3_gc_locked_and_tm_matched"])
    out.to_csv(f"{OUT}/primer_verification.csv", index=False)
    with open(f"{OUT}/primer_verification.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1, default=str)
    n = len(out)
    print(f"\n验证完成: {n} 对引物 | V1 单靶 {int(out.V1_specific.sum())}/{n} "
          f"| V2 无强二级结构 {int(out.V2_no_strong_secondary_structure.sum())}/{n} "
          f"| V3 GC/Tm 合规 {int(out.V3_gc_locked_and_tm_matched.sum())}/{n} "
          f"| 总通过 {int(out.overall_pass.sum())}/{n}")
    print("DONE — verify_primers.py (v2)")


if __name__ == "__main__":
    main()
