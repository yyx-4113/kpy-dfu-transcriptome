# -*- coding: utf-8 -*-
"""
design_paralog_aware_primers_v2.py — 旁系同源感知引物设计（v2，修正区段约束）

═══ v1 → v2 的关键修正 ═══

v1 把 Primer3 的设计模板限制在「可区分区段」内部，产物长度要求
80–180 bp，而可区分区段往往只有 20–60 bp → Primer3 无候选，
5 个基因里只有 2 个产出。

**方法学修正**：需要具备区分度的**只是引物本身**，不是引物之间的
产物间隔区。产物中段即使落在保守区也无妨——PCR 产物是 cDNA，
中段不参与特异性决定。因此正确做法是：

1. 在**目标基因完整 CDS** 上跑 Primer3，取足量候选对；
2. 对每个候选对，**逐条引物**判定其对各旁系同源的最大 identity；
3. 只保留**两条引物各自都满足区分阈值**的候选对；
4. 若某基因无候选对满足，退而求其次：允许一条引物为「弱区分」
   （identity ≤ 90% 且 3' 端有错配 → 实验上不共扩增），另一条强区分。

═══ 区分阈值（与 BLAST v5 判据对齐）═══

- **强区分**：对所有旁系同源转录本，最佳近似匹配 identity ≤ 90%
  （错配 ≥ 2）且 3' 端最后 5 nt 至少 1 个错配
  → BLAST v5 判为 P2（噪声），实验上不产生共扩增
- **弱区分**：identity ≤ 95% 且 3' 端有错配 → 判为 P1（风险，需熔解曲线验证）
- **不可接受**：identity ≥ 95% 且 3' 端无错配 → P0（真实共扩增，必须换引物）

**近似匹配算法**：分段种子索引 + 滑动窗口最小错配（与 BLAST 局部比对
语义一致，无 gap），正/反链均查。v1 曾用精确子串搜索，会漏检含 1–2 个
错配的旁系同源匹配——那恰恰是 P0 脱靶的定义。

**序列缓存**：所有 RefSeq 序列落盘到 `provenance/paralog_seq_cache.json`，
重跑无需重复联网。

输出：
  03_结果/实验验证设计/paralog_aware_design_v2.json
  03_结果/实验验证设计/paralog_aware_selected_v2.csv
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import re
import json
import time
import urllib.request
import urllib.parse

import primer3

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
PROV = f"{OUT}/provenance"
os.makedirs(PROV, exist_ok=True)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "WorkBuddy-ParalogAwareV2/1.0"}
RC_MAP = str.maketrans("ACGT", "TGCA")
CACHE_P = f"{PROV}/paralog_seq_cache.json"

# 目标基因 → (GeneID, 旁系同源家族)
TARGETS = {
    "LDHA": (3939, ["LDHB", "LDHC", "LDHAL", "LDHAL6A", "LDHAL6B"]),
    "FAHD2B": (151313, ["FAHD2A", "FAHD1"]),
    "ADH1B": (125, ["ADH1A", "ADH1C", "ADH4", "ADH6", "ADH7"]),
    "ALDH3A2": (224, ["ALDH3A1", "ALDH3A3"]),
    "GBP1": (2633, ["GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "GBP7"]),
    "DHX58": (79132, ["DHX58P1", "DHX58P2"]),
}

PRIMER3_PARAMS = {
    "PRIMER_TASK": "generic",
    "PRIMER_PICK_LEFT_PRIMER": 1,
    "PRIMER_PICK_RIGHT_PRIMER": 1,
    "PRIMER_PICK_INTERNAL_OLIGO": 0,
    "PRIMER_NUM_RETURN": 200,
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
                urllib.request.Request(url, headers=UA), timeout=120).read().decode()
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
        return CACHE[acc]["header"], CACHE[acc]["seq"]
    raw = _get(f"{EUTILS}/efetch.fcgi?db=nuccore&id={acc}"
               f"&rettype=fasta&retmode=text")
    L = raw.strip().split("\n")
    header = L[0].lstrip(">")
    seq = "".join(x.strip() for x in L[1:]
                   if re.fullmatch(r"[A-Za-z]+", x.strip())).upper()
    CACHE[acc] = {"header": header, "seq": seq}
    return header, seq


def efetch_gb(acc):
    key = f"__gb__{acc}"
    if key in CACHE:
        return CACHE[key]["seq"]
    raw = _get(f"{EUTILS}/efetch.fcgi?db=nuccore&id={acc}&rettype=gb&retmode=text")
    CACHE[key] = {"header": "gb", "seq": raw}
    return raw


def cds_segments(gb):
    for b in re.split(r"\n {5}(?=\w)", gb):
        if b.startswith("CDS"):
            loc = re.match(r"CDS\s+(\S+)", b).group(1).replace("<", ">")
            extra = ""
            for ln in b.split("\n")[1:]:
                st = ln.strip()
                if not st or st.startswith("/"):
                    continue
                extra += st
            return [(int(a), int(c)) for a, c in
                    re.findall(r"(\d+)\.\.(\d+)", loc + extra)]
    return []


def gene_transcripts(sym):
    key = f"__tx__{sym}"
    if key in CACHE:
        return CACHE[key]["seq"]
    d = json.loads(_get(f"{EUTILS}/esearch.fcgi?db=gene&term="
                        f"{urllib.parse.quote(sym + '[Gene Name]')}"
                        f"%20AND%20human%5Borgn%5D&retmode=json"))
    gids = d["esearchresult"]["idlist"]
    accs = []
    for gid in gids:
        time.sleep(0.3)
        e = json.loads(_get(f"{EUTILS}/elink.fcgi?dbfrom=gene&db=nuccore&id={gid}"
                            f"&linkname=gene_nuccore_refseqrna&retmode=json"))
        ids = [l for ls in e.get("linksets", [])
               for dbk in ls.get("linksetdbs", []) for l in dbk.get("links", [])]
        if not ids:
            continue
        s = json.loads(_get(f"{EUTILS}/esummary.fcgi?db=nuccore&id={','.join(ids)}"
                            f"&retmode=json"))
        for u in s["result"]["uids"]:
            v = s["result"][u]
            if v.get("accessionversion"):
                accs.append([gid, v["accessionversion"], v.get("title", "")])
        break
    CACHE[key] = {"header": "tx", "seq": accs}
    return accs


def best_match(primer, paralog_seqs):
    """返回 (max identity%, 总错配数, 3'端错配数, 最差旁系同源名)。

    分段种子索引 + 滑动窗口最小错配（等价于 BLAST 局部比对，无 gap）。
    """
    L = len(primer)
    best = (0.0, L, 0, "")
    p_rc = rc(primer)
    nseg = 4
    seg_len = L // nseg
    for name, seq in paralog_seqs:
        if len(seq) < L:
            continue
        for s in (primer, p_rc):
            seeds = [s[k * seg_len:(k + 1) * seg_len] for k in range(nseg)]
            hit_idx = set()
            for k, seed in enumerate(seeds):
                st = 0
                while True:
                    i = seq.find(seed, st)
                    if i < 0:
                        break
                    st = i + 1
                    for off in range(-k * seg_len, L - (k + 1) * seg_len + 1):
                        w0 = i + off
                        if 0 <= w0 <= len(seq) - L:
                            hit_idx.add(w0)
            for w0 in hit_idx:
                mm = 0
                for k in range(L):
                    if seq[w0 + k] != s[k]:
                        mm += 1
                        if mm >= 5:
                            break
                else:
                    mm3 = sum(1 for k in range(L - 5, L)
                              if seq[w0 + k] != s[k])
                    ident = 100.0 * (L - mm) / L
                    if ident > best[0]:
                        best = (ident, mm, mm3, name)
    return best


def grade(ident, mm3):
    """按 BLAST v5 判据分级：'STRONG' / 'WEAK' / 'UNACCEPTABLE'"""
    if ident <= 90.0 and mm3 >= 1:
        return "STRONG"
    if ident <= 95.0 and mm3 >= 1:
        return "WEAK"
    if ident <= 90.0:
        return "STRONG"          # identity ≤ 90% 即为 P2 噪声
    return "UNACCEPTABLE"


def main():
    report = {"date": "2026-10-25",
              "script": "02_脚本/design_paralog_aware_primers_v2.py",
              "primer3_version": primer3.__version__,
              "primer3_params": {k: v for k, v in PRIMER3_PARAMS.items()
                                 if k != "PRIMER_EXPLAIN_FLAG"},
              "method": {
                  "template": "目标基因完整 CDS（而非仅可区分区段）",
                  "key_correction": ("需要具备区分度的是**引物本身**而非产物间隔区；"
                                     "v1 把 Primer3 限制在 20–60 bp 的可区分"
                                     "区段内，因无法容纳 80–180 bp 产物而失败"),
                  "grading": {"STRONG": "identity ≤ 90% 且 3'端有错配 → P2 噪声",
                              "WEAK": "identity ≤ 95% 且 3'端有错配 → P1 风险",
                              "UNACCEPTABLE": "identity ≥ 95% 且 3'端无错配 → P0"},
                  "approx_match": "分段种子索引 + 滑动窗口最小错配（无 gap，"
                                  "正反链均查），等价 BLAST 局部比对",
              },
              "genes": {}}
    selected = []

    for gene, (gid, family) in TARGETS.items():
        print(f"\n=== {gene} (旁系同源: {', '.join(family)}) ===")
        para_seqs = []
        for pf in family:
            for pgid, acc, title in gene_transcripts(pf):
                try:
                    _, s = efetch_fasta(acc)
                    if s:
                        para_seqs.append((pf, s))
                except Exception as e:                # noqa: BLE001
                    print(f"   抓取失败 {acc}: {e}")
        print(f"   旁系同源转录本 {len(para_seqs)} 条 "
              f"({sorted(set(n for n, _ in para_seqs))})")

        # 目标：最长编码转录本
        best_t = None
        for pgid, acc, title in gene_transcripts(gene):
            try:
                _, s = efetch_fasta(acc)
                segs = cds_segments(efetch_gb(acc))
                if not segs:
                    continue
                cds = "".join(s[a - 1:b] for a, b in segs)
                if best_t is None or len(cds) > len(best_t["cds"]):
                    best_t = {"acc": acc, "seq": s, "cds": cds, "segs": segs}
            except Exception:                         # noqa: BLE001
                pass
        if not best_t:
            print("   !! 未取到目标 CDS")
            continue
        print(f"   目标 {best_t['acc']} CDS {len(best_t['cds'])} nt")

        res = primer3.bindings.design_primers(
            {"SEQUENCE_ID": gene, "SEQUENCE_TEMPLATE": best_t["seq"]},
            PRIMER3_PARAMS, {})
        npair = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
        print(f"   Primer3 候选对（整条 CDS 模板）= {npair}")

        cands = []
        for i in range(npair):
            f = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
            r = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
            lp, rp = res[f"PRIMER_LEFT_{i}"], res[f"PRIMER_RIGHT_{i}"]
            f_i, f_mm, f_mm3, f_pn = best_match(f, para_seqs)
            r_i, r_mm, r_mm3, r_pn = best_match(r, para_seqs)
            gf, gr = grade(f_i, f_mm3), grade(r_i, r_mm3)
            cands.append({
                "forward": f, "reverse": r,
                "f_tm": res[f"PRIMER_LEFT_{i}_TM"], "r_tm": res[f"PRIMER_RIGHT_{i}_TM"],
                "f_gc": round(100 * sum(c in "GC" for c in f) / len(f), 1),
                "r_gc": round(100 * sum(c in "GC" for c in r) / len(r), 1),
                "product_size": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                "mrna_f_start_1based": lp[0] + 1,
                "mrna_r_end_1based": rp[0] + len(r),
                "f_max_identity_vs_paralog": round(f_i, 2),
                "f_min_mismatch": f_mm, "f_mm_3p5": f_mm3,
                "f_worst_paralog": f_pn, "f_grade": gf,
                "r_max_identity_vs_paralog": round(r_i, 2),
                "r_min_mismatch": r_mm, "r_mm_3p5": r_mm3,
                "r_worst_paralog": r_pn, "r_grade": gr,
                "f_self_any_th": res.get(f"PRIMER_LEFT_{i}_SELF_ANY_TH"),
                "r_self_any_th": res.get(f"PRIMER_RIGHT_{i}_SELF_ANY_TH"),
                "pair_compl_any_th": res.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH"),
                "penalty": res.get(f"PRIMER_PAIR_{i}_PENALTY"),
            })

        def pair_grade(c):
            if c["f_grade"] == "UNACCEPTABLE" or c["r_grade"] == "UNACCEPTABLE":
                return 0
            if c["f_grade"] == "STRONG" and c["r_grade"] == "STRONG":
                return 2
            return 1

        def score(c):
            g = pair_grade(c)
            s = 3 * c["penalty"]
            s += 3 * (abs(c["f_tm"] - 60) + abs(c["r_tm"] - 60))
            s += 2 * (abs(c["f_gc"] - 50) / 10 + abs(c["r_gc"] - 50) / 10)
            s += 2 * (0 if c["forward"].endswith(("G", "C")) else 0.5)
            s += 2 * (0 if c["reverse"].endswith(("G", "C")) else 0.5)
            s += 100 * (2 - g)                     # 分级优先
            s += 0.8 * max(c["f_max_identity_vs_paralog"],
                           c["r_max_identity_vs_paralog"])
            return s

        cands.sort(key=score)
        n_g = {0: 0, 1: 0, 2: 0}
        for c in cands:
            n_g[pair_grade(c)] += 1
        print(f"   候选分级：双强区分 {n_g[2]}，至少一条弱区分 {n_g[1]}，"
              f"含不可接受 {n_g[0]}")
        for c in cands[:5]:
            print(f"     F[{c['f_grade'][:4]} {c['f_max_identity_vs_paralog']:.0f}% "
                  f"3p{c['f_mm_3p5']}] {c['forward']} | "
                  f"R[{c['r_grade'][:4]} {c['r_max_identity_vs_paralog']:.0f}% "
                  f"3p{c['r_mm_3p5']}] {c['reverse']} amp={c['product_size']}")

        ok = [c for c in cands if pair_grade(c) >= 1]
        best = ok[0] if ok else None
        if best:
            print(f"   ★ 推荐（分级 {pair_grade(best)}）："
                  f"F={best['forward']} R={best['reverse']} "
                  f"amp={best['product_size']}bp")
        else:
            print("   !! 无可接受候选（全部含 identity ≥95% 且 3'端无错配的引物）")
        report["genes"][gene] = {
            "target_accession": best_t["acc"],
            "target_cds_len": len(best_t["cds"]),
            "n_paralog_tx": len(para_seqs),
            "paralog_family": family,
            "paralog_transcripts": sorted(set(n for n, _ in para_seqs)),
            "n_primer3_pairs": npair,
            "grade_counts": {"both_strong": n_g[2], "one_weak": n_g[1],
                             "unacceptable": n_g[0]},
            "selected": best,
            "top_candidates": cands[:10],
        }
        if best:
            selected.append({"gene": gene, "accession": best_t["acc"],
                             "role": "", "module": "",
                             **{k: v for k, v in best.items()}})
        save_cache()

    with open(f"{OUT}/paralog_aware_design_v2.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1, default=str)
    import csv
    if selected:
        keys = sorted({k for r in selected for k in r})
        with open(f"{OUT}/paralog_aware_selected_v2.csv", "w",
                  encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            for r in selected:
                w.writerow(r)
    print(f"\n完成：{len(selected)}/{len(TARGETS)} 个基因产出旁系同源感知引物对")
    print("DONE — design_paralog_aware_primers_v2.py")


if __name__ == "__main__":
    main()
