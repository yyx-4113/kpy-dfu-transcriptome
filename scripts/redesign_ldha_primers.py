# -*- coding: utf-8 -*-
"""
redesign_ldha_primers.py — LDHA 引物重新设计（因 RefSeq 版本漂移致旧引物失效）

═══ 触发原因（2026-10-25 独立复核发现，非臆测）═══

原设计模板 `NM_001135239.1`（design_rtqpcr_primers.py，2026-10-25）。
独立复核（BLAST v4 + efetch 版本核对）发现：

1. **旧前向引物在当前 RefSeq 中完全不存在**：
   `AGCTGTTCCACTTAAGGCCC` 在 LDHA 全部 6 条当前 RefSeq 转录本
   （NR_028500.2 / NM_001165414.2 / NM_001165415.2 / NM_001165416.2 /
   NM_001135239.2 / NM_005566.4）中 **匹配位置均为 -1**
   （正向与反向互补均为 -1）→ 该引物对**无法扩增 LDHA cDNA**。

2. **RefSeq 版本间 5'UTR 差异巨大**：`.1` 与 `.2` 的位置一致性仅 0.2524
   （长度 2052 vs 2067）；`.2` 的 GenBank 注释为
   exon1 = 1..74, exon2 = 75..224, CDS = 99..923。
   旧前向引物位于 `.1` 的 136..155（1-based），该区段在 `.2` 中已不复存在。

3. **旧反向引物存在真实脱靶**：`AAGGAATCGGGAATGCACGT` 与 **LYST**
   的 8 条人源 RefSeq 转录本（XM_054334017/18/19/20、XM_047443038/40/46、
   XM_011544037）全长 20 bp 匹配、identity 19/20（95%），
   且错配位于引物**中段**（qseq 索引 9），**3' 端 5 nt 完全匹配**
   → 具备真实共扩增风险。

→ 结论：LDHA 引物对**必须重新设计**，不能沿用。

═══ 重新设计原则 ═══

1. **模板 = 当前 RefSeq canonical（NM_005566.4）**，且设计窗口必须落在
   **全部 6 条 LDHA 转录本共有的序列区段**（用各转录本自身序列求交，
   不复用同一参考切片——v1 脚本曾犯此 bug）；
2. **优先落在 CDS 内**（5'UTR 在各变体间高度变异，UTR 内引物易漏检变体）；
3. Primer3 参数与原设计脚本**完全一致**（保证与其余 11 对引物口径统一）；
4. 设计后须通过：跨接头检查（exon 边界内不跨）、二级结构/3' 端互补阈值
   （≤ 47 °C，与原口径一致）、GC 40–60%、Tm 58–62 °C、产物 80–180 bp；
5. **最终以 BLAST 重新验证**（v4 判据：人源物种限定 + 全长匹配 + 3' 端考量）。

输出：
  03_结果/实验验证设计/ldha_redesign_candidates.json   候选对 + 全部证据
  03_结果/实验验证设计/provenance/ldha_transcripts.fasta 6 条转录本序列（留档）
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import json
import re
import time
import urllib.request

import primer3

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
PROV = f"{OUT}/provenance"
os.makedirs(PROV, exist_ok=True)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "WorkBuddy-PrimerRedesign/1.0"}

# LDHA (GeneID 3939) 全部 6 条 RefSeq 转录本（2026-10-25 经
# gene→nuccore elink + esummary 实测枚举，非凭记忆书写）
LDHA_TX = [
    "NM_005566.4",     # variant 1, canonical
    "NM_001135239.2",  # variant 2
    "NM_001165414.2",  # variant 3
    "NM_001165415.2",  # variant 4
    "NM_001165416.2",  # variant 5
    "NR_028500.2",     # variant 6, non-coding
]
CANONICAL = "NM_005566.4"

# 与 design_rtqpcr_primers.py 完全一致的 Primer3 参数（不得偏离）
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
    "PRIMER_PRODUCT_SIZE_RANGE": [[80, 180]],
    "PRIMER_EXPLAIN_FLAG": 1,
}
RC_MAP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(RC_MAP)[::-1]


def _get(url):
    return urllib.request.urlopen(
        urllib.request.Request(url, headers=UA), timeout=120).read().decode()


def efetch_fasta(acc):
    raw = _get(f"{EUTILS}/efetch.fcgi?db=nuccore&id={acc}"
               f"&rettype=fasta&retmode=text")
    lines = raw.strip().split("\n")
    header = lines[0].lstrip(">")
    seq = "".join(l.strip() for l in lines[1:]
                   if re.fullmatch(r"[A-Za-z]+", l.strip())
                   and not l.startswith("BASE"))
    return header, seq.upper()


def efetch_gb(acc):
    return _get(f"{EUTILS}/efetch.fcgi?db=nuccore&id={acc}"
                f"&rettype=gb&retmode=text")


def parse_cds(gb_text):
    """从 GenBank flat file 提取 CDS 的 join(...) 区间（1-based 闭区间列表）。"""
    m = re.search(r"^\s{5}CDS\s+(\S+)", gb_text, re.M)
    if not m:
        return []
    loc = m.group(1).replace("<", ">")
    parts = re.findall(r"(\d+)\.\.(\d+)", loc)
    return [(int(a), int(b)) for a, b in parts]


def main():
    # ---------- 1. 取全部转录本序列 + GenBank CDS 注释 ----------
    tx = {}
    for acc in LDHA_TX:
        header, seq = efetch_fasta(acc)
        time.sleep(0.4)
        gb = efetch_gb(acc)
        time.sleep(0.4)
        cds = parse_cds(gb)
        tx[acc] = {"header": header, "seq": seq, "cds": cds, "len": len(seq)}
        print(f"  {acc:16s} len={len(seq):5d} CDS={cds}")
        with open(f"{PROV}/ldha_transcripts.fasta", "a", encoding="utf-8") as fh:
            fh.write(f">{acc} {header} CDS={cds}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    canon = tx[CANONICAL]

    # ---------- 2. 求 canonical 与其余转录本的共有 25-mer（逐本自身序列）----------
    #   修正 v1 脚本 bug：v1 对所有转录本都使用 canonical 的切片做索引，
    #   导致交集恒等于全集（假阳性"共有区段 0..2241"）。
    def kmers(s, k=25):
        return {s[i:i + k] for i in range(len(s) - k + 1)}

    ref = canon["seq"]
    common = kmers(ref)
    for acc in LDHA_TX:
        if acc == CANONICAL:
            continue
        common &= kmers(tx[acc]["seq"])
    common_pos = sorted(i for i in range(len(ref) - 24) if ref[i:i + 25] in common)

    segs = []
    if common_pos:
        st = prev = common_pos[0]
        for p in common_pos[1:]:
            if p <= prev + 25:
                prev = p
            else:
                segs.append((st, prev + 25))
                st = prev = p
        segs.append((st, prev + 25))
    print(f"\n  共有 25-mer 数 = {len(common_pos)}")
    for a, b in segs:
        print(f"    共有区段 {a}..{b} (len {b - a})")

    # ---------- 3. CDS 内共有区段（引物应尽量落在 CDS，保证各变体均可检出）----------
    #    canonical CDS = join 区间；取 CDS 在 canonical 序列上的映射
    cds_ranges = canon["cds"]
    cds_set = set()
    for a, b in cds_ranges:
        cds_set.update(range(a - 1, b))          # 0-based
    cds_common = [i for i in common_pos
                  if all(i + k in cds_set for k in range(25))]
    cds_seg = []
    if cds_common:
        st = prev = cds_common[0]
        for p in cds_common[1:]:
            if p <= prev + 25:
                prev = p
            else:
                cds_seg.append((st, prev + 25))
                st = prev = p
        cds_seg.append((st, prev + 25))
    print(f"  CDS 内共有 25-mer 数 = {len(cds_common)}")
    for a, b in cds_seg:
        print(f"    CDS 共有区段 {a}..{b} (len {b - a})  {ref[a:a + 36]}...")

    # ---------- 4. 在 CDS 共有区段上跑 Primer3 ----------
    #    以 canonical 序列为模板，把 Primer3 的 SEQUENCE 限制在 CDS 共有区段
    #    （用 PRIMER3 的 "SEQUENCE_" 子串偏移不直接支持，故改为：
    #    在共有区段序列上独立设计，再把坐标映射回 canonical）
    results = []
    for seg_a, seg_b in cds_seg:
        sub = ref[seg_a:seg_b]
        if len(sub) < 130:                       # 需容纳 80–180 bp 产物
            continue
        res = primer3.bindings.design_primers(
            {"SEQUENCE_ID": f"LDHA_CDS_{seg_a}_{seg_b}",
             "SEQUENCE_TEMPLATE": sub}, PRIMER3_PARAMS, {})
        npairs = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
        for i in range(npairs):
            f = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
            r = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
            lp = res[f"PRIMER_LEFT_{i}"]
            rp = res[f"PRIMER_RIGHT_{i}"]
            f_abs0 = seg_a + lp[0]                # 0-based 绝对坐标
            r_abs0 = seg_a + rp[0]
            results.append({
                "seg": [seg_a, seg_b],
                "forward": f, "reverse": r,
                "f_tm": res[f"PRIMER_LEFT_{i}_TM"],
                "r_tm": res[f"PRIMER_RIGHT_{i}_TM"],
                "f_gc": round(100 * sum(c in "GC" for c in f) / len(f), 1),
                "r_gc": round(100 * sum(c in "GC" for c in r) / len(r), 1),
                "product_size": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                "canonical_f_start_1based": f_abs0 + 1,
                "canonical_f_end_1based": f_abs0 + len(f),
                "canonical_r_start_1based": r_abs0 + 1,
                "canonical_r_end_1based": r_abs0 + len(r),
                "f_self_any_th": res.get(f"PRIMER_LEFT_{i}_SELF_ANY_TH"),
                "f_self_end_th": res.get(f"PRIMER_LEFT_{i}_SELF_END_TH"),
                "f_hairpin_th": res.get(f"PRIMER_LEFT_{i}_HAIRPIN_TH"),
                "r_self_any_th": res.get(f"PRIMER_RIGHT_{i}_SELF_ANY_TH"),
                "r_self_end_th": res.get(f"PRIMER_RIGHT_{i}_SELF_END_TH"),
                "r_hairpin_th": res.get(f"PRIMER_RIGHT_{i}_HAIRPIN_TH"),
                "pair_compl_any_th": res.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH"),
                "pair_compl_end_th": res.get(f"PRIMER_PAIR_{i}_COMPL_END_TH"),
                "penalty": res.get(f"PRIMER_PAIR_{i}_PENALTY"),
            })

    # ---------- 5. 逐条转录本落位核查（正反向引物须在每条转录本都出现）----------
    def locate(seq, primer):
        """返回 (正向命中位置, 反向互补命中位置)；-1 表示未命中"""
        return seq.find(primer), seq.find(rc(primer))

    for c in results:
        cov = {}
        for acc in LDHA_TX:
            s = tx[acc]["seq"]
            # 前向引物在转录本正义链上直接出现；
            # 反向引物以反向互补形式出现（PCR 时二者指向同一扩增子）。
            iF = s.find(c["forward"])
            iRrc = s.find(rc(c["reverse"]))
            ok = (iF >= 0) and (iRrc >= 0)
            cov[acc] = {"F_pos": iF, "R_rc_pos": iRrc,
                        "in_CDS": (iF >= 0 and any(a - 1 <= iF < b
                                                   for a, b in tx[acc]["cds"])),
                        "is_coding": bool(tx[acc]["cds"]),
                        "ok": ok}
        c["transcript_coverage"] = cov
        c["n_tx_covered"] = sum(1 for v in cov.values() if v["ok"])
        # 覆盖判据：全部 **蛋白编码** 转录本（NM_）必须命中；
        # NR_028500.2 为非编码变体（无 CDS），CDS 内引物本就不应命中，
        # 单列报告而不计入覆盖分母。
        coding = [a for a in LDHA_TX if tx[a]["cds"]]
        c["n_coding_tx"] = len(coding)
        c["n_coding_covered"] = sum(1 for a in coding if cov[a]["ok"])
        c["all_coding_covered"] = (c["n_coding_covered"] == len(coding))
        c["all_tx_covered"] = (c["n_tx_covered"] == len(LDHA_TX))

    # ---------- 6. 打分排序：全转录本覆盖优先，其次 penalty/Tm/GC/3' GC ----------
    def score(c):
        s = 3 * c["penalty"]
        s += 3 * (abs(c["f_tm"] - 60) + abs(c["r_tm"] - 60))
        s += 2 * (abs(c["f_gc"] - 50) / 10 + abs(c["r_gc"] - 50) / 10)
        s += 2 * (0 if c["forward"].endswith(("G", "C")) else 0.5)
        s += 2 * (0 if c["reverse"].endswith(("G", "C")) else 0.5)
        s += 20 * (c["n_coding_tx"] - c["n_coding_covered"])  # 编码变体覆盖不全重罚
        return s

    results.sort(key=score)
    print(f"\n  Primer3 候选对（CDS 共有区段内）= {len(results)}")
    for c in results[:8]:
        print(f"    cov={c['n_coding_covered']}/{c['n_coding_tx']}编码 "
              f"(全部 {c['n_tx_covered']}/{len(LDHA_TX)}) "
              f"F={c['forward']} R={c['reverse']} amp={c['product_size']} "
              f"Tm={c['f_tm']:.1f}/{c['r_tm']:.1f} "
              f"GC={c['f_gc']}/{c['r_gc']} pen={c['penalty']:.4f}")

    best = next((c for c in results if c["all_coding_covered"]), None)
    out = {
        "date": "2026-10-25",
        "script": "02_脚本/redesign_ldha_primers.py",
        "reason": ("原设计模板 NM_001135239.1 的 5'UTR 在 RefSeq 版本更新中重排，"
                   "旧前向引物 AGCTGTTCCACTTAAGGCCC 在全部 6 条当前 LDHA "
                   "RefSeq 转录本中匹配位置均为 -1（正反向互补均为 -1），"
                   "无法扩增 LDHA cDNA；且旧反向引物与 LYST 的 8 条人源 "
                   "RefSeq 转录本全长 20 bp 匹配、identity 19/20、错配位于中段"
                   "（3' 端 5 nt 全匹配），存在真实共扩增风险 → 整对重新设计。"),
        "old_primers": {
            "template": "NM_001135239.1",
            "forward": "AGCTGTTCCACTTAAGGCCC",
            "reverse": "AAGGAATCGGGAATGCACGT",
            "f_match_in_current_tx": {a: tx[a]["seq"].find("AGCTGTTCCACTTAAGGCCC")
                                      for a in LDHA_TX},
        },
        "transcripts": {a: {"header": tx[a]["header"], "len": tx[a]["len"],
                            "cds": tx[a]["cds"]} for a in LDHA_TX},
        "canonical": CANONICAL,
        "shared_25mer_segments": segs,
        "shared_25mer_in_cds_segments": cds_seg,
        "n_shared_25mer": len(common_pos),
        "n_shared_25mer_in_cds": len(cds_common),
        "primer3_version": primer3.__version__,
        "primer3_params": {k: v for k, v in PRIMER3_PARAMS.items()
                           if k != "PRIMER_EXPLAIN_FLAG"},
        "candidates": results,
        "selected": best,
    }
    with open(f"{OUT}/ldha_redesign_candidates.json", "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1, default=str)

    print("\n  === 推荐对（覆盖全部 5 条蛋白编码转录本）===")
    if best:
        print(f"    F = {best['forward']}  (canonical {best['canonical_f_start_1based']}"
              f"..{best['canonical_f_end_1based']}, Tm {best['f_tm']:.2f}, "
              f"GC {best['f_gc']}%)")
        print(f"    R = {best['reverse']}  (canonical {best['canonical_r_start_1based']}"
              f"..{best['canonical_r_end_1based']}, Tm {best['r_tm']:.2f}, "
              f"GC {best['r_gc']}%)")
        print(f"    产物 = {best['product_size']} bp")
        print("    落位：")
        for acc, v in best["transcript_coverage"].items():
            print(f"      {acc:16s} F@{v['F_pos']:5d} Rrc@{v['R_rc_pos']:5d} "
                  f"in_CDS={v['in_CDS']} coding={v['is_coding']} ok={v['ok']}")
    else:
        print("    !! 无候选对覆盖全部编码转录本，需放宽设计窗口或改用基因组坐标")
    print("\nDONE — redesign_ldha_primers.py")


if __name__ == "__main__":
    main()
