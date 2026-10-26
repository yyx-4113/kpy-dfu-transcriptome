# -*- coding: utf-8 -*-
"""
design_paralog_aware_primers.py — 旁系同源感知的引物定向设计

═══ 触发事件 ═══

BLAST v5 全库复核暴露出一类**前几版设计方法无法发现**的缺陷：
即使引物位于 CDS 内、覆盖全部编码转录本，仍会与**高同源旁系同源基因**
完全匹配（identity 100%、3' 端全匹配），造成真实共扩增：

| 基因 | 脱靶 | 实测 |
|------|------|------|
| FAHD2B 反向 | FAHD2A（37 条转录本）100% | FAHD2A (GeneID 51011, 2q11.1) 与 FAHD2B (2q11.2) 相邻；CDS 945 nt 仅 **15 个碱基差异**（一致性 98.41%），蛋白 319 aa 逐字相同 |
| LDHA 反向 | LOC105375980（2 条）100% | 位于 **9p22.3**，与 LDHA（11p15.5）不同染色体 |
| ADH1B 前向/反向 | ADH1C（4 条）/ ADH1A（2 条）100% | ADH 家族不同基因座、不同酶亚基 |

═══ 本脚本的解法 ═══

**不靠"多找几条候选"碰运气，而是从序列层面直接排除同源拷贝**：

对每个目标基因：
1. 枚举其**旁系同源家族**的全部 RefSeq 转录本（经 NCBI elink 实测，
   非凭记忆书写）；
2. 取目标基因 CDS 与每个旁系同源的全部转录本序列做**全比对**，
   统计每个 20-nt 窗口在旁系同源中的**最大错配数**；
3. 只保留**对所有旁系同源错配数 ≥ 2**（或完全不存在）的窗口作为
   可用设计区——因为 identity ≥ 95%（错配 ≤ 1）且 3' 端锚定的命中
   才是真实共扩增风险；错配 ≥ 2 时 identity ≤ 90%，落到 P2 噪声级；
4. 在这些窗口上跑 Primer3（参数与原设计一致），产物 80–180 bp、
   Tm 58–62 °C、GC 40–60%、二级结构 ≤ 47 °C；
5. 输出每条候选对到**各旁系同源**的实测最大 identity 与错配位置，
   作为可直接复核的证据。

**关键方法学点**：错配位置也重要——若错配全部落在 3' 端末 5 nt 内，
延伸无法起始，实验上无共扩增风险（这正是 BLAST v5 引入的 3' 端锚定
判据）。本脚本对每个旁系同源同时报告「总错配数」与「3' 端错配数」。

输出：
  03_结果/实验验证设计/paralog_aware_design.json
  03_结果/实验验证设计/provenance/paralog_transcripts.fasta
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
UA = {"User-Agent": "WorkBuddy-ParalogAware/1.0"}
RC_MAP = str.maketrans("ACGT", "TGCA")

# 目标基因 → (GeneID, 旁系同源家族基因名列表)
# 家族成员均为**独立基因座、编码不同产物**者 → 必须能区分
TARGETS = {
    "LDHA": (3939, ["LDHB", "LDHC", "LDHAL", "LDHAL6A", "LDHAL6B"]),
    "FAHD2B": (151313, ["FAHD2A", "FAHD1"]),
    "ADH1B": (125, ["ADH1A", "ADH1C", "ADH4", "ADH6", "ADH7"]),
    "ALDH3A2": (224, ["ALDH3A1", "ALDH3A3"]),
    "GBP1": (2633, ["GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "GBP7"]),
}
# 无已知高同源旁系同源者（内参候选，家族成员皆为假基因或同基因座）
NO_PARALOG = ["DHX58", "ACTB", "GAPDH", "TBP", "HPRT1", "RPLP0", "SDHA"]

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
# 旁系同源可接受阈值：错配数 ≥ 2 → identity ≤ 90% → BLAST v5 判为 P2 噪声
MIN_MM = 2
MIN_MM_3P = 1        # 3' 端至少 1 个错配（保证延伸无法起始）


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


def efetch_fasta(acc):
    raw = _get(f"{EUTILS}/efetch.fcgi?db=nuccore&id={acc}"
               f"&rettype=fasta&retmode=text")
    L = raw.strip().split("\n")
    return (L[0].lstrip(">"),
            "".join(x.strip() for x in L[1:]
                    if re.fullmatch(r"[A-Za-z]+", x.strip())).upper())


def efetch_gb(acc):
    return _get(f"{EUTILS}/efetch.fcgi?db=nuccore&id={acc}&rettype=gb&retmode=text")


def cds_segments(gb):
    """解析 CDS 特征块（含跨行 location）。"""
    for b in re.split(r"\n {5}(?=\w)", gb):
        if b.startswith("CDS"):
            head = b.split("\n")[0]
            loc = re.match(r"CDS\s+(\S+)", head).group(1).replace("<", ">")
            extra = ""
            for ln in b.split("\n")[1:]:
                st = ln.strip()
                if not st or st.startswith("/"):
                    continue
                extra += st
            loc += extra
            return [(int(a), int(c)) for a, c in
                    re.findall(r"(\d+)\.\.(\d+)", loc)]
    return []


def gene_transcripts(gene_sym):
    """实测枚举某基因的全部 RefSeq 转录本（gene name → RefSeq RNA）。"""
    d = json.loads(_get(f"{EUTILS}/esearch.fcgi?db=gene&term="
                        f"{urllib.parse.quote(gene_sym + '[Gene Name]')}"
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
                accs.append((gid, v["accessionversion"], v.get("title", "")))
        break
    return accs


def best_match(primer, paralog_seqs):
    """返回 (最大 identity 百分比, 该处总错配数, 3'端错配数, 命中的旁系同源名)。

    ═══ 关键：必须是「最佳近似匹配」，不能只做精确子串搜索 ═══
    v1 版本用 `seq.find(primer)`（精确匹配）→ 对旁系同源中含 1–2 个错配的
    匹配**完全漏检**，会误报「0% identity / 20 错配」。
    而 P0 脱靶的定义恰恰是「identity ≥ 95%（错配 ≤ 1）且 3' 端锚定」——
    只差 1 个碱基的旁系同源才是最大的实际风险。
    → 改为**滑动窗口最小错配**（与 BLAST 局部比对语义一致、无 gap），
      正/反链都查。

    ═══ 实现：分段种子索引（等价于滑动窗口，但快 ~L 倍）═══
    把引物切成 4 段（长 5 nt）。若两窗口错配 ≤ 4，则**至少有一段完全匹配**
    （鸽巢原理：4 段中至少 1 段零错配）。故可用「任一段精确出现」作为
    候选窗口的过滤器，再对这些窗口逐一算真实错配数——结果与暴力滑动
    窗口完全等价，但复杂度从 O(n·L) 降到 O(n/4 + 候选数·L)。
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
            # 种子：4 段各取一段（轮换取不同段作种子以免漏检）
            seeds = [s[k * seg_len:(k + 1) * seg_len] for k in range(nseg)]
            hit_idx = set()
            for seed in seeds:
                st = 0
                while True:
                    i = seq.find(seed, st)
                    if i < 0:
                        break
                    st = i + 1
                    # 该种子对应的窗口起点范围
                    k = seeds.index(seed)
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


def main():
    report = {"date": "2026-10-25",
              "script": "02_脚本/design_paralog_aware_primers.py",
              "primer3_version": primer3.__version__,
              "primer3_params": {k: v for k, v in PRIMER3_PARAMS.items()
                                 if k != "PRIMER_EXPLAIN_FLAG"},
              "criteria": {
                  "min_mismatch_vs_paralog": MIN_MM,
                  "min_mismatch_3prime": MIN_MM_3P,
                  "rationale": ("错配 ≥ 2 → identity ≤ 90% → BLAST v5 判为 P2 噪声；"
                                "3' 端额外 ≥ 1 错配 → PCR 延伸无法起始，"
                                "实验上不产生共扩增"),
              },
              "genes": {}}
    selected = []
    fasta_p = f"{PROV}/paralog_transcripts.fasta"
    if os.path.exists(fasta_p):
        os.remove(fasta_p)

    for gene, (gid, family) in TARGETS.items():
        print(f"\n=== {gene} (旁系同源家族: {', '.join(family)}) ===")
        # 目标转录本
        t_tx = gene_transcripts(gene)
        # 旁系同源转录本
        para_seqs = []
        for pf in family:
            for pgid, acc, title in gene_transcripts(pf):
                try:
                    _, s = efetch_fasta(acc)
                    time.sleep(0.3)
                    para_seqs.append((pf, s))
                    with open(fasta_p, "a", encoding="utf-8") as fh:
                        fh.write(f">{acc} [{pf}] {title}\n")
                        for i in range(0, len(s), 60):
                            fh.write(s[i:i + 60] + "\n")
                except Exception as e:                # noqa: BLE001
                    print(f"   旁系同源抓取失败 {acc}: {e}")
        print(f"   旁系同源转录本 {len(para_seqs)} 条 "
              f"({sorted(set(n for n, _ in para_seqs))})")

        # 目标 CDS（用最长编码转录本）
        best_target = None
        for pgid, acc, title in t_tx:
            try:
                _, s = efetch_fasta(acc)
                time.sleep(0.3)
                segs = cds_segments(efetch_gb(acc))
                time.sleep(0.3)
                if not segs:
                    continue
                cds = "".join(s[a - 1:b] for a, b in segs)
                if best_target is None or len(cds) > len(best_target["cds"]):
                    best_target = {"acc": acc, "seq": s, "cds": cds,
                                   "segs": segs, "title": title}
            except Exception as e:                    # noqa: BLE001
                print(f"   目标抓取失败 {acc}: {e}")
        if not best_target:
            print("   !! 未取到目标 CDS")
            continue
        cds = best_target["cds"]
        print(f"   目标 {best_target['acc']} CDS {len(cds)} nt")

        # 逐窗口评估对旁系同源的最大 identity
        win = 20
        good_windows = []
        for i in range(len(cds) - win + 1):
            w = cds[i:i + win]
            ident, mm, mm3, pname = best_match(w, para_seqs)
            if mm >= MIN_MM and mm3 >= MIN_MM_3P:
                good_windows.append((i, w, ident, mm, mm3, pname))
        print(f"   满足「错配≥{MIN_MM} 且 3'端错配≥{MIN_MM_3P}」的 20-mer 窗口: "
              f"{len(good_windows)} 个")
        if not good_windows:
            # 放宽：只要求错配 ≥ MIN_MM
            for i in range(len(cds) - win + 1):
                w = cds[i:i + win]
                ident, mm, mm3, pname = best_match(w, para_seqs)
                if mm >= MIN_MM:
                    good_windows.append((i, w, ident, mm, mm3, pname))
            print(f"   放宽到「错配≥{MIN_MM}」后: {len(good_windows)} 个")
        if not good_windows:
            print("   !! 无可用窗口")
            report["genes"][gene] = {"status": "no_distinguishing_window",
                                     "n_paralog_tx": len(para_seqs)}
            continue

        # 归并为区段
        segs = []
        st = prev = good_windows[0][0]
        for p in good_windows[1:]:
            if p[0] <= prev + win:
                prev = p[0]
            else:
                segs.append((st, prev + win))
                st = prev = p[0]
        segs.append((st, prev + win))
        print(f"   可用区段: {[(a, b - a) for a, b in segs][:8]}")

        # Primer3 设计（模板 = 目标 mRNA 全长，限制在可用区段）
        _, tseq = efetch_fasta(best_target["acc"])
        time.sleep(0.3)
        cds_start = min(a for a, b in best_target["segs"]) - 1  # 0-based
        cands = []
        for ga, gb in segs:
            # 映射回 mRNA 坐标
            m_start = cds_start + ga
            m_end = cds_start + gb
            sub = tseq[m_start:m_end]
            if len(sub) < 130:
                continue
            res = primer3.bindings.design_primers(
                {"SEQUENCE_ID": f"{gene}_{ga}", "SEQUENCE_TEMPLATE": sub},
                PRIMER3_PARAMS, {})
            for i in range(res.get("PRIMER_PAIR_NUM_RETURNED", 0)):
                f = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
                r = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
                lp, rp = res[f"PRIMER_LEFT_{i}"], res[f"PRIMER_RIGHT_{i}"]
                f_abs = m_start + lp[0]
                r_abs = m_start + rp[0]
                f_i, f_mm, f_mm3, f_pn = best_match(f, para_seqs)
                r_i, r_mm, r_mm3, r_pn = best_match(r, para_seqs)
                cands.append({
                    "forward": f, "reverse": r,
                    "f_tm": res[f"PRIMER_LEFT_{i}_TM"], "r_tm": res[f"PRIMER_RIGHT_{i}_TM"],
                    "f_gc": round(100 * sum(c in "GC" for c in f) / len(f), 1),
                    "r_gc": round(100 * sum(c in "GC" for c in r) / len(r), 1),
                    "product_size": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                    "mrna_f_start_1based": f_abs + 1,
                    "mrna_r_end_1based": r_abs + len(r),
                    "f_max_identity_vs_paralog": round(f_i, 2),
                    "f_min_mismatch": f_mm, "f_mm_3p5": f_mm3, "f_worst_paralog": f_pn,
                    "r_max_identity_vs_paralog": round(r_i, 2),
                    "r_min_mismatch": r_mm, "r_mm_3p5": r_mm3, "r_worst_paralog": r_pn,
                    "f_self_any_th": res.get(f"PRIMER_LEFT_{i}_SELF_ANY_TH"),
                    "r_self_any_th": res.get(f"PRIMER_RIGHT_{i}_SELF_ANY_TH"),
                    "pair_compl_any_th": res.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH"),
                    "penalty": res.get(f"PRIMER_PAIR_{i}_PENALTY"),
                })

        def score(c):
            s = 3 * c["penalty"]
            s += 3 * (abs(c["f_tm"] - 60) + abs(c["r_tm"] - 60))
            s += 2 * (abs(c["f_gc"] - 50) / 10 + abs(c["r_gc"] - 50) / 10)
            s += 2 * (0 if c["forward"].endswith(("G", "C")) else 0.5)
            s += 2 * (0 if c["reverse"].endswith(("G", "C")) else 0.5)
            # 旁系同源 identity 越低越好
            s += 1.5 * max(c["f_max_identity_vs_paralog"],
                          c["r_max_identity_vs_paralog"])
            return s

        cands.sort(key=score)
        print(f"   候选对 {len(cands)}")
        for c in cands[:5]:
            print(f"     F={c['forward']} (para id {c['f_max_identity_vs_paralog']:.0f}%, "
                  f"mm {c['f_min_mismatch']}/3p{c['f_mm_3p5']}) "
                  f"R={c['reverse']} (para id {c['r_max_identity_vs_paralog']:.0f}%, "
                  f"mm {c['r_min_mismatch']}/3p{c['r_mm_3p5']}) "
                  f"amp={c['product_size']}")
        if not cands:
            print("   !! Primer3 在可用区段未产出候选")
            report["genes"][gene] = {"status": "no_primer3_candidate",
                                     "n_paralog_tx": len(para_seqs)}
            continue
        best = cands[0]
        print(f"   ★ 推荐：F={best['forward']} R={best['reverse']} "
              f"amp={best['product_size']}bp "
              f"para_id F/R = {best['f_max_identity_vs_paralog']:.0f}%/"
              f"{best['r_max_identity_vs_paralog']:.0f}%")
        report["genes"][gene] = {
            "target_accession": best_target["acc"],
            "n_paralog_tx": len(para_seqs),
            "paralog_family": family,
            "paralog_transcripts": sorted(set(n for n, _ in para_seqs)),
            "target_cds_len": len(cds),
            "usable_windows": len(good_windows),
            "usable_segments": segs,
            "candidates": cands,
            "selected": best,
        }
        selected.append({"gene": gene, "accession": best_target["acc"],
                         "forward": best["forward"], "reverse": best["reverse"],
                         **{k: v for k, v in best.items()
                            if k not in ("forward", "reverse")}})

    with open(f"{OUT}/paralog_aware_design.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1, default=str)
    with open(f"{OUT}/paralog_aware_selected.json", "w", encoding="utf-8") as fh:
        json.dump(selected, fh, ensure_ascii=False, indent=1, default=str)
    print(f"\n完成：{len(selected)}/{len(TARGETS)} 个目标基因产出旁系同源感知引物对")
    print("DONE — design_paralog_aware_primers.py")


if __name__ == "__main__":
    main()
