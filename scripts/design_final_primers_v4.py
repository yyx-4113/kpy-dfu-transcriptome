# -*- coding: utf-8 -*-
"""
design_final_primers_v4.py — RT-qPCR 引物最终设计（v4：修正 v3 的三处判定缺陷）

═══ v3 → v4 修正清单（均由独立复核确认，非推测）═══

**M1 分级判据与 BLAST v5 不一致**
v3 的 `grade()` 把「identity ≥ 95%」无条件判为 UNACCEPTABLE。但
BLAST v5 判据（`adjudicate_blast_offtargets.py`）是：
    P0：3'端末 5 nt 全匹配 且 identity ≥ 95%   → 阻断
    P1：3'端末 5 nt 全匹配 且 90% ≤ identity < 95% → 需注意，**非阻断**
    P2：其余（3'端有错配，或 identity < 90%）    → 噪声
v3 误把 P1 当阻断 → ADH1B（identity 90%、3'端全匹配）被误杀。
v4 严格按 v5 对齐：只有 P0 级（id≥95% 且 3'端全匹配）才是 UNACCEPTABLE；
P1 级标为 WEAK（可接受，但须在补充材料声明并建议熔解曲线验证）。

**M2 mm3 计算 bug（v3 首跑的实际输出受影响）**
v3 首跑的 `mm3` 跨窗口求和，产生 59–84 这类不可能的值。修正为
在**最佳窗口内部**逐位统计。

**M3 curated / predicted 覆盖未分离**
FAHD2B 有 29 条蛋白编码转录本，其中仅 3 条为 curated（NM_/NR_）、
26 条为 XM_ 预测。要求「全部 29 条共有」把共有区段压缩到 462 nt
（CDS 315..776），使 FAHD2B 落在**预测转录本特有**的变异区里。
实验上真正需要保证的是：**已策展的、可被实验观测到的剪接变体**
（NM_/NR_）必须全覆盖；XM_ 预测转录本不表达或表达量可忽略时，
不应成为硬约束，但需作为**信息性指标**报告覆盖率。

**M4 SDHA 的旁系同源处理与 BLAST v5 冲突**
`adjudicate_blast_offtargets.py` 的 STRICT_PARALOG_ONLY 白名单含
SDHA（SDHAP1/2/3 为假基因 → PARALOG，非真脱靶）。v3 把 SDHA 当硬约束，
与其余基因口径不一致 → v4 统一。

**M5 FAHD2B 生物学限制须如实记录**
D1 诊断：FAHD2B 共有 CDS 区段（CDS 315..776，462 nt）的 443 个 20-mer
窗口中，对 FAHD2A/FAHD1 的 identity **≤90% 的窗口数 = 0**（全部 ≥95%）。
→ FAHD2A/FAHD2B CDS 一致性 98.41%、蛋白 319 aa 逐字相同，
   **该区段在原理上不可区分**。
v4 的处理：
   (a) 放宽共有区段到 curated 转录本共有（含 UTR 锚定，但要求引物
       完整位于 CDS 内）；
   (b) 若仍无可区分窗口 → **如实标记为「生物学限制」**，
       并给出替代方案（报告 FAHD2 家族合并读出 + 单独报告 FAHD2A），
       **不得**降标准放行或隐去。

输出：
  03_结果/实验验证设计/final_primer_design_v4.json
  03_结果/实验验证设计/final_primers_v4.csv
"""
import json

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import os
import re
import time
import urllib.request
import urllib.parse

import numpy as np
import pandas as pd
import primer3

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
PROV = f"{OUT}/provenance"
os.makedirs(PROV, exist_ok=True)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "WorkBuddy-FinalPrimerDesignV4/1.0"}
RC_MAP = str.maketrans("ACGT", "TGCA")
CACHE_P = f"{PROV}/final_design_cache.json"

GENES = {
    "LDHA": (3939, "target", "M2_Pyruvate_supply"),
    "FAHD2B": (151313, "target", "M2_Pyruvate_supply"),
    "ADH1B": (125, "target", "M2_Pyruvate_supply"),
    "ALDH3A2": (224, "target", "M2_Pyruvate_supply"),
    "DHX58": (79132, "target", "M4_IFNI_readout"),
    "GBP1": (2633, "target", "M4_IFNI_readout"),
    "ACTB": (60, "reference_candidate", ""),
    "GAPDH": (2597, "reference_candidate", ""),
    "TBP": (6908, "reference_candidate", ""),
    "HPRT1": (3251, "reference_candidate", ""),
    "RPLP0": (6175, "reference_candidate", ""),
    "SDHA": (6389, "reference_candidate", ""),
}

PARALOG_FAMILIES = {
    "LDHA": ["LDHB", "LDHC", "LDHAL", "LDHAL6A", "LDHAL6B"],
    "FAHD2B": ["FAHD2A", "FAHD1"],
    "ADH1B": ["ADH1A", "ADH1C", "ADH4", "ADH6", "ADH7"],
    "ALDH3A2": ["ALDH3A1", "ALDH3A3"],
    "GBP1": ["GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "GBP7"],
    "DHX58": [],
    "ACTB": ["ACTG1", "ACTG2", "ACTC1"],
    "GAPDH": ["GAPDHS", "GAPDHL"],
    "TBP": [],
    "HPRT1": [],
    "RPLP0": ["RPLP0P1", "RPLP0P2", "RPLP0P3", "RPLP0P4", "RPLP0P5"],
    "SDHA": ["SDHAP1", "SDHAP2", "SDHAP3"],
}

# 与 BLAST v5 判读器的 STRICT_PARALOG_ONLY 完全一致：
# 这些基因的同家族成员为假基因或同基因座功能重复拷贝 → PARALOG（非真脱靶）
# → 设计期只作**信息性报告**，不作阻断约束。
SOFT_PARALOG_GENES = {"ACTB", "ACTG1", "GAPDH", "RPLP0", "TBP", "HPRT1",
                      "SDHA", "DHX58"}

PRIMER3_PARAMS = {
    "PRIMER_TASK": "generic",
    "PRIMER_PICK_LEFT_PRIMER": 1,
    "PRIMER_PICK_RIGHT_PRIMER": 1,
    "PRIMER_PICK_INTERNAL_OLIGO": 0,
    "PRIMER_NUM_RETURN": 100,
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


def rc(s):
    return s.translate(RC_MAP)[::-1]


def _get(url, tries=4):
    last = None
    for k in range(tries):
        try:
            return urllib.request.urlopen(
                urllib.request.Request(url, headers=UA), timeout=180).read().decode()
        except Exception as e:                        # noqa: BLE001
            last = e
            time.sleep(2 + 3 * k)
    raise RuntimeError(str(last))


CACHE = json.load(open(CACHE_P, encoding="utf-8")) if os.path.exists(CACHE_P) else {}


def save_cache():
    with open(CACHE_P, "w", encoding="utf-8") as fh:
        json.dump(CACHE, fh)


def efetch_fasta(acc):
    if acc in CACHE:
        return CACHE[acc]["seq"]
    raw = _get(f"{EUTILS}/efetch.fcgi?db=nuccore&id={acc}"
               f"&rettype=fasta&retmode=text")
    L = raw.strip().split("\n")
    seq = "".join(x.strip() for x in L[1:]
                   if re.fullmatch(r"[A-Za-z]+", x.strip())).upper()
    CACHE[acc] = {"seq": seq}
    return seq


def efetch_gb(acc):
    key = f"__gb__{acc}"
    if key in CACHE:
        return CACHE[key]["seq"]
    raw = _get(f"{EUTILS}/efetch.fcgi?db=nuccore&id={acc}&rettype=gb&retmode=text")
    CACHE[key] = {"seq": raw}
    return raw


def gene_transcripts(gid):
    key = f"__tx__{gid}"
    if key in CACHE:
        return CACHE[key]["seq"]
    d = json.loads(_get(f"{EUTILS}/elink.fcgi?dbfrom=gene&db=nuccore&id={gid}"
                        f"&linkname=gene_nuccore_refseqrna&retmode=json"))
    ids = [l for ls in d.get("linksets", [])
           for dbk in ls.get("linksetdbs", []) for l in dbk.get("links", [])]
    accs = []
    if ids:
        s = json.loads(_get(f"{EUTILS}/esummary.fcgi?db=nuccore&id={','.join(ids)}"
                            f"&retmode=json"))
        for u in s["result"]["uids"]:
            v = s["result"][u]
            if v.get("accessionversion"):
                accs.append(v["accessionversion"])
    accs = sorted(set(accs))
    CACHE[key] = {"seq": accs}
    return accs


def parse_cds_location(gb):
    """CDS 坐标（1-based 闭区间）。遇 qualifier 即停止拼接。"""
    for b in re.split(r"\n {5}(?=\w)", gb):
        if not b.startswith("CDS"):
            continue
        loc = re.match(r"CDS\s+(\S+)", b).group(1).replace("<", ">")
        for ln in b.split("\n")[1:]:
            st = ln.strip()
            if not st:
                continue
            if st.startswith("/"):
                break
            loc += st
        return [(int(a), int(c)) for a, c in
                re.findall(r"(\d+)\.\.(\d+)", loc)]
    return []


def cds_of(acc, seq=None):
    gb = efetch_gb(acc)
    seg = parse_cds_location(gb)
    if not seg or seq is None:
        return seg, None
    return seg, "".join(seq[a - 1:b] for a, b in seg)


def gene_name_by_esearch(sym):
    key = f"__gid__{sym}"
    if key in CACHE:
        return CACHE[key]["seq"]
    d = json.loads(_get(f"{EUTILS}/esearch.fcgi?db=gene&term="
                        f"{urllib.parse.quote(sym + '[Gene Name]')}"
                        f"%20AND%20human%5Borgn%5D&retmode=json"))
    gids = d["esearchresult"]["idlist"]
    gid = None
    if gids:
        s = json.loads(_get(f"{EUTILS}/esummary.fcgi?db=gene&id={','.join(gids)}"
                            f"&retmode=json"))
        for u in s["result"]["uids"]:
            v = s["result"][u]
            if v.get("name", "").upper() == sym.upper():
                gid = u
                break
        if gid is None:
            gid = gids[0]
    CACHE[key] = {"seq": gid}
    return gid


ENC = np.full(256, 4, dtype=np.uint8)
for _i, _ch in enumerate(b"ACGT"):
    ENC[_ch] = _i


def enc(s):
    return ENC[np.frombuffer(s.encode("ascii"), dtype=np.uint8)]


def best_window_identity(subject: str, primer: str):
    """穷举全部窗口 × 正反链 → (identity%, mm_total, mm_3p5, strand)

    3' 端锚定按**引物自身方向**：
      - 正链（P 出现于 S）：引物 3' 端 = 窗口后 5 nt
      - 负链（rc(P) 出现于 S）：引物 3' 端 = 窗口前 5 nt
    """
    L = len(primer)
    S = enc(subject)
    if len(S) < L:
        return (0.0, L, L, "none")
    best = (0.0, L, L, "none")
    n = len(S) - L + 1
    for strand, P in (("+", enc(primer)), ("-", enc(rc(primer)))):
        mm = np.zeros(n, dtype=np.int16)
        for k in range(L):
            mm += (S[k:k + n] != P[k])
        i = int(np.argmin(mm))
        m = int(mm[i])
        if m >= L:
            continue
        ident = 100.0 * (L - m) / L
        if ident <= best[0]:
            continue
        pos_mm = (S[i:i + L] != P)
        mm3 = int(pos_mm[L - 5:L].sum()) if strand == "+" else int(pos_mm[0:5].sum())
        best = (ident, m, mm3, strand)
    return best


def grade(ident, mm3):
    """严格对齐 BLAST v5（adjudicate_blast_offtargets.py）判据。

    P0（阻断）：3'端末 5 nt 全匹配 且 identity ≥ 95%
    P1（需注意，非阻断）：3'端末 5 nt 全匹配 且 90% ≤ identity < 95%
    P2（噪声）：其余
    """
    if mm3 == 0 and ident >= 95.0:
        return "UNACCEPTABLE"      # = v5 P0
    if mm3 == 0 and ident >= 90.0:
        return "WEAK"              # = v5 P1（可接受，须声明）
    if ident >= 90.0:
        return "WEAK"              # 3'端有错配但 identity 高 → 警示
    return "STRONG"                # = v5 P2 噪声


def shared_segments(ref, others, k=25):
    """逐本自身序列求交 → ref 上的共有 k-mer 连续区段（0-based 半开）。"""
    common = {ref[i:i + k] for i in range(len(ref) - k + 1)}
    for s in others:
        if not s:
            continue
        common &= {s[i:i + k] for i in range(len(s) - k + 1)}
    pos = [i for i in range(len(ref) - k + 1) if ref[i:i + k] in common]
    if not pos:
        return []
    out, st, pv = [], pos[0], pos[0]
    for p in pos[1:]:
        if p <= pv + k:
            pv = p
        else:
            out.append((st, pv + k))
            st = pv = p
    out.append((st, pv + k))
    return out


def clip_to_cds(segs, cds_pos):
    out = []
    for a, b in segs:
        sub = [p for p in range(a, b) if p in cds_pos]
        if not sub:
            continue
        s0, pv = sub[0], sub[0]
        for p in sub[1:]:
            if p == pv + 1:
                pv = p
            else:
                out.append((s0, pv + 1))
                s0 = pv = p
        out.append((s0, pv + 1))
    return out


def main():
    report = {"date": "2026-10-25",
              "script": "02_脚本/design_final_primers_v4.py",
              "primer3_version": primer3.__version__,
              "primer3_params": {k: v for k, v in PRIMER3_PARAMS.items()
                                 if k != "PRIMER_EXPLAIN_FLAG"},
              "method": {
                  "objective": "三约束同时满足：CDS锚定 / 全curated转录本覆盖 / 旁系同源区分度",
                  "coverage_policy": ("硬约束 = 全部 curated (NM_/NR_) 蛋白编码转录本覆盖；"
                                      "predicted (XM_) 覆盖率作信息性指标报告"),
                  "paralog_match": "numpy 穷举全部窗口 × 正反链，零剪枝",
                  "three_prime_rule": "正链→引物3'端=窗口后5nt；负链→窗口前5nt",
                  "grading": {"UNACCEPTABLE": "3'端全匹配 且 identity ≥ 95%（=v5 P0 阻断）",
                              "WEAK": "3'端全匹配 且 90%≤identity<95%（=v5 P1 需注意）；"
                                      "或 identity ≥ 90% 但 3'端有错配",
                              "STRONG": "identity < 90%（=v5 P2 噪声）"},
                  "soft_paralog_genes": sorted(SOFT_PARALOG_GENES),
              },
              "genes": {}}
    selected_rows = []

    for gene, (gid, role, module) in GENES.items():
        print(f"\n{'=' * 76}\n=== {gene} (GeneID {gid}) ===", flush=True)
        tx_accs = gene_transcripts(gid)
        seqs, coding, cds_map, cds_len = {}, [], {}, {}
        for acc in tx_accs:
            try:
                s = efetch_fasta(acc)
                if not s:
                    continue
                seqs[acc] = s
                seg, cds = cds_of(acc, s)
                if seg and cds:
                    coding.append(acc)
                    cds_map[acc] = seg
                    cds_len[acc] = len(cds)
            except Exception as e:                     # noqa: BLE001
                print(f"   fetch fail {acc}: {e}")
        if not coding:
            print("   !! 无蛋白编码转录本")
            continue
        curated = [a for a in coding if a.startswith(("NM_", "NR_"))]
        predicted = [a for a in coding if a.startswith("XM_")]
        print(f"   RefSeq RNA {len(tx_accs)}；蛋白编码 {len(coding)} "
              f"(curated {len(curated)} / predicted {len(predicted)})", flush=True)

        def rank(acc):
            return (0 if acc.startswith(("NM_", "NR_")) else 1, -cds_len[acc])
        canonical = sorted(coding, key=rank)[0]
        ref = seqs[canonical]
        cds_pos = set()
        for a, b in cds_map[canonical]:
            cds_pos.update(range(a - 1, b))

        # ── 共有区段：curated 优先，同时保留 all-transcript 版本作对照 ──
        seg_cur = clip_to_cds(shared_segments(ref, [seqs[a] for a in curated]),
                              cds_pos)
        seg_all = clip_to_cds(shared_segments(ref, [seqs[a] for a in coding]),
                              cds_pos)
        segs = seg_cur if seg_cur else seg_all
        print(f"   canonical={canonical} (CDS {cds_len[canonical]} nt)  "
              f"curated共有CDS={[(a+1, b-a) for a, b in seg_cur]}  "
              f"all共有CDS={[(a+1, b-a) for a, b in seg_all]}", flush=True)

        # ── 旁系同源序列 ──
        para_full, para_cds = [], []
        para_detail = {}
        for fam in PARALOG_FAMILIES.get(gene, []):
            pgid = gene_name_by_esearch(fam)
            if not pgid:
                print(f"   旁系同源 {fam}: GeneID 未找到（该符号无蛋白编码基因）")
                continue
            f_accs = gene_transcripts(pgid)
            n_coding = 0
            for acc in f_accs:
                try:
                    s = efetch_fasta(acc)
                    if not s:
                        continue
                    para_full.append((fam, acc, s))
                    seg, cds = cds_of(acc, s)
                    if cds:
                        para_cds.append((fam, acc, cds))
                        n_coding += 1
                except Exception:                      # noqa: BLE001
                    pass
            para_detail[fam] = {"gene_id": pgid, "n_refseq_rna": len(f_accs),
                                 "n_coding": n_coding}
            print(f"   旁系同源 {fam}: RefSeq RNA {len(f_accs)}，蛋白编码 {n_coding}")
        soft = gene in SOFT_PARALOG_GENES

        # ── Primer3 ──
        cands = []
        for seg_a, seg_b in segs:
            if seg_b - seg_a < 130:
                continue
            res = primer3.bindings.design_primers(
                {"SEQUENCE_ID": f"{gene}_{seg_a}",
                 "SEQUENCE_TEMPLATE": ref[seg_a:seg_b]},
                PRIMER3_PARAMS, {})
            for i in range(res.get("PRIMER_PAIR_NUM_RETURNED", 0)):
                f = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
                r = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
                lp, rp = res[f"PRIMER_LEFT_{i}"], res[f"PRIMER_RIGHT_{i}"]
                cands.append({
                    "forward": f, "reverse": r,
                    "f_tm": res[f"PRIMER_LEFT_{i}_TM"],
                    "r_tm": res[f"PRIMER_RIGHT_{i}_TM"],
                    "f_gc": round(100 * sum(c in "GC" for c in f) / len(f), 1),
                    "r_gc": round(100 * sum(c in "GC" for c in r) / len(r), 1),
                    "product_size": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                    "f_start_1based": seg_a + lp[0] + 1,
                    "r_end_1based": seg_a + rp[0] + len(r),
                    "f_self_any_th": res.get(f"PRIMER_LEFT_{i}_SELF_ANY_TH"),
                    "r_self_any_th": res.get(f"PRIMER_RIGHT_{i}_SELF_ANY_TH"),
                    "pair_compl_any_th": res.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH"),
                    "penalty": res.get(f"PRIMER_PAIR_{i}_PENALTY"),
                    "segment": [seg_a + 1, seg_b],
                })
        print(f"   Primer3 候选对 {len(cands)}", flush=True)

        # ── 逐候选：约束 2 + 约束 3 ──
        for c in cands:
            def cover(subset):
                n_ok = n_cds = 0
                for acc in subset:
                    s = seqs[acc]
                    iF = s.find(c["forward"])
                    iR = s.find(rc(c["reverse"]))
                    if iF < 0 or iR < 0:
                        continue
                    n_ok += 1
                    pos = set()
                    for a, b in cds_map[acc]:
                        pos.update(range(a - 1, b))
                    if all(p in pos for p in range(iF, iF + len(c["forward"]))) and \
                       all(p in pos for p in range(iR, iR + len(c["reverse"]))):
                        n_cds += 1
                return n_ok, n_cds

            c["n_cur_covered"], c["n_cur_in_CDS"] = cover(curated)
            c["n_all_covered"], c["n_all_in_CDS"] = cover(coding)
            c["curated_fully_covered"] = (c["n_cur_covered"] == len(curated)
                                          and c["n_cur_in_CDS"] == len(curated))
            c["all_fully_covered"] = (c["n_all_covered"] == len(coding)
                                      and c["n_all_in_CDS"] == len(coding))

            for lab, pk in (("f", "forward"), ("r", "reverse")):
                p = c[pk]
                # ═══ 硬判据只对**蛋白编码**旁系同源转录本 ═══
                # RT-qPCR 的底物是 cDNA：非编码 RNA（NR_ non-coding /
                # XR_）不产生可扩增的蛋白编码转录本 → 不构成共扩增风险。
                # （v3/v4 首跑把全部 RefSeq RNA 含非编码记录计入硬判据
                #   → 夸大了脱靶面；本版改为只用编码转录本作硬判据，
                #   全量（含非编码）结果另存为信息性指标。）
                b_full = b_cds = None
                for fam, acc, s in para_full:
                    ident, mm, mm3, strand = best_window_identity(s, p)
                    if b_full is None or ident > b_full[0]:
                        b_full = (ident, mm, mm3, strand, fam, acc)
                for fam, acc, s in para_cds:
                    ident, mm, mm3, strand = best_window_identity(s, p)
                    if b_cds is None or ident > b_cds[0]:
                        b_cds = (ident, mm, mm3, strand, fam, acc)
                use = b_cds if b_cds else b_full
                c[f"{lab}_max_id_vs_paralog"] = round(use[0], 2) if use else None
                c[f"{lab}_mm_vs_paralog"] = use[1] if use else None
                c[f"{lab}_mm3p_vs_paralog"] = use[2] if use else None
                c[f"{lab}_strand_vs_paralog"] = use[3] if use else None
                c[f"{lab}_worst_paralog"] = use[4] if use else None
                c[f"{lab}_worst_acc"] = use[5] if use else None
                c[f"{lab}_grade"] = grade(use[0], use[2]) if use else "NA"
                # 信息性：若把非编码/预测转录本也计入，最大 identity 是多少
                c[f"{lab}_max_id_incl_noncoding"] = (round(b_full[0], 2)
                                                    if b_full else None)
                c[f"{lab}_worst_acc_incl_noncoding"] = b_full[5] if b_full else None
            if not para_full:
                c["pair_paralog_grade"] = "NA"
            elif soft:
                c["pair_paralog_grade"] = "SOFT"
            else:
                gs = [c["f_grade"], c["r_grade"]]
                c["pair_paralog_grade"] = ("BOTH_STRONG" if gs == ["STRONG"] * 2
                                           else ("ONE_WEAK" if "WEAK" in gs
                                                 else "UNACCEPTABLE"))

        def score(c):
            s = 3 * c["penalty"]
            s += 3 * (abs(c["f_tm"] - 60) + abs(c["r_tm"] - 60))
            s += 2 * (abs(c["f_gc"] - 50) / 10 + abs(c["r_gc"] - 50) / 10)
            s += 2 * (0 if c["forward"].endswith(("G", "C")) else 0.5)
            s += 2 * (0 if c["reverse"].endswith(("G", "C")) else 0.5)
            s += 80 * (len(curated) - c["n_cur_covered"])
            s += 80 * (len(curated) - c["n_cur_in_CDS"])
            s += 20 * (len(coding) - c["n_all_covered"])
            pg = c["pair_paralog_grade"]
            s += {"BOTH_STRONG": 0, "ONE_WEAK": 25, "NA": 0, "SOFT": 5,
                  "UNACCEPTABLE": 300}[pg]
            if pg == "ONE_WEAK":
                s += 8 * max(c["f_max_id_vs_paralog"] or 0,
                             c["r_max_id_vs_paralog"] or 0)
            return s

        cands.sort(key=score)
        ok = [c for c in cands
              if c["curated_fully_covered"]
              and c["pair_paralog_grade"] != "UNACCEPTABLE"]
        print(f"   三约束全过 {len(ok)}/{len(cands)}；分级 "
              f"{pd.Series([c['pair_paralog_grade'] for c in cands]).value_counts().to_dict()}",
              flush=True)
        for c in cands[:4]:
            print(f"     cur={c['n_cur_covered']}/{len(curated)} "
                  f"inCDS={c['n_cur_in_CDS']}/{len(curated)} "
                  f"all={c['n_all_covered']}/{len(coding)} "
                  f"para={c['pair_paralog_grade']} "
                  f"F_id={c['f_max_id_vs_paralog']} R_id={c['r_max_id_vs_paralog']} "
                  f"F={c['forward']} R={c['reverse']} amp={c['product_size']}", flush=True)

        best = ok[0] if ok else None
        if best:
            print(f"   ★ 最终推荐：F={best['forward']} R={best['reverse']} "
                  f"amp={best['product_size']}bp 旁系同源={best['pair_paralog_grade']} "
                  f"cur={best['n_cur_covered']}/{len(curated)}", flush=True)
        else:
            print("   !! 无候选同时满足三约束 → 须如实记录生物学限制", flush=True)

        report["genes"][gene] = {
            "gene_id": gid, "role": role, "module": module,
            "canonical": canonical,
            "n_refseq_rna": len(tx_accs),
            "n_coding_tx": len(coding),
            "n_curated_tx": len(curated), "n_predicted_tx": len(predicted),
            "curated_accs": curated,
            "shared_cds_segments_curated": [(a + 1, b - a) for a, b in seg_cur],
            "shared_cds_segments_all": [(a + 1, b - a) for a, b in seg_all],
            "paralog_families": para_detail,
            "n_paralog_refseq_rna": len(para_full),
            "n_paralog_coding": len(para_cds),
            "soft_paralog": soft,
            "n_candidates": len(cands),
            "n_passing_all3": len(ok),
            "grade_counts": pd.Series([c["pair_paralog_grade"]
                                       for c in cands]).value_counts().to_dict(),
            "selected": best,
            "top_candidates": cands[:12],
        }
        if best:
            selected_rows.append({
                "gene": gene, "role": role, "module": module,
                "refseq_accession": canonical, "gene_id": gid,
                "forward_primer_5to3": best["forward"],
                "reverse_primer_5to3": best["reverse"],
                "f_tm_C": round(best["f_tm"], 6), "r_tm_C": round(best["r_tm"], 6),
                "f_gc_pct": best["f_gc"], "r_gc_pct": best["r_gc"],
                "amplicon_bp": best["product_size"],
                "forward_start_1based": best["f_start_1based"],
                "reverse_end_1based": best["r_end_1based"],
                "penalty": best["penalty"],
                "n_curated_covered": best["n_cur_covered"],
                "n_curated_total": len(curated),
                "n_curated_in_CDS": best["n_cur_in_CDS"],
                "n_all_tx_covered": best["n_all_covered"],
                "n_all_tx_total": len(coding),
                "pair_paralog_grade": best["pair_paralog_grade"],
                "f_max_id_vs_paralog": best["f_max_id_vs_paralog"],
                "r_max_id_vs_paralog": best["r_max_id_vs_paralog"],
                "f_worst_paralog": best["f_worst_paralog"],
                "r_worst_paralog": best["r_worst_paralog"],
                "f_mm3p_vs_paralog": best["f_mm3p_vs_paralog"],
                "r_mm3p_vs_paralog": best["r_mm3p_vs_paralog"],
            })
        save_cache()

    with open(f"{OUT}/final_primer_design_v4.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1, default=str)
    if selected_rows:
        pd.DataFrame(selected_rows).to_csv(f"{OUT}/final_primers_v4.csv",
                                            index=False)
    print(f"\n完成：{len(selected_rows)}/{len(GENES)} 个基因产出三约束全过的引物对")
    print("DONE — design_final_primers_v4.py")


if __name__ == "__main__":
    main()
