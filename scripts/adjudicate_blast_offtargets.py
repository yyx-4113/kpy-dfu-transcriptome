# -*- coding: utf-8 -*-
"""
adjudicate_blast_offtargets.py — BLAST v4 命中的严格脱靶判读（独立代码路径）

═══ 为什么需要独立判读器 ═══

`verify_primers_blast_v4.py` 的 `classify()` 把「全长匹配（align_len == 引物长度）
的人源非目标 RefSeq 转录本命中」全部计为真脱靶。该判据**过于粗糙**，
在两个方向上都会出错：

1. **假阳性（把噪声当脱靶）**：ALDH3A2 前向引物在人源库中有 3,111 条
   全长命中、反向 4,713 条——这是 Alu/LINE 重复序列散布拷贝的特征，
   而非 4,713 个"会扩增的脱靶"。这类命中绝大多数 identity 仅 18–19/20，
   且 3' 端存在错配。
2. **假阴性（漏掉真脱靶）**：未考虑 **3' 端错配位置**。PCR 延伸依赖
   3' 端锚定：3' 端末 5 nt 完全匹配 + 总 identity ≥ 95% 的非目标转录本
   是**真实的共扩增风险**；反之 3' 端有错配者几乎无法起始延伸。

═══ 本判读器判据（每条命中逐项记录，可复算）═══

对每条全长（align_len == 引物长度）人源 RefSeq 转录本命中：

| 指标 | 定义 |
|------|------|
| `identity_n` | 匹配碱基数（Hsp_identity，注意**不是**百分比） |
| `identity_pct` | identity_n / align_len × 100 |
| `mm_3p5` | 3' 端最后 5 nt 的错配数（PCR 延伸关键区） |
| `mm_total` | 全长错配数 |
| `is_pseudogene` | title 含 pseudogene |
| `is_paralog` | 基因属目标基因的旁系同源家族 |
| `is_target` | 基因 == 目标基因 |

分级：
- **P0 真脱靶（阻断级）**：非目标、非旁系同源、非假基因，且
  `mm_3p5 == 0` 且 `identity_pct ≥ 95`  → 必须更换引物
- **P1 风险脱靶（需谨慎）**：`mm_3p5 == 0` 且 `90 ≤ identity_pct < 95`
  → 报告并建议熔解曲线/测序验证
- **P2 噪声**：其余（3' 端有错配，或 identity < 90%）→ 不计脱靶

同时独立复核：
- 目标基因是否存在**全长匹配的人源 RefSeq 转录本**（不存在则引物不可用）
- 目标位点的最大 identity

输出：
  03_结果/实验验证设计/blast_offtarget_adjudication.csv  逐命中明细
  03_结果/实验验证设计/blast_offtarget_adjudication.json  逐引物汇总
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import re
import json
import csv
import glob

import pandas as pd

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
RAW = f"{OUT}/provenance/blast_xml_raw_v4"
RC_MAP = str.maketrans("ACGT", "TGCA")

PARALOG_FAMILIES = {
    "ACTB": {"ACTB", "ACTG1", "ACTG2", "ACTC1"},
    "ACTG1": {"ACTB", "ACTG1", "ACTG2", "ACTC1"},
    "GAPDH": {"GAPDH", "GAPDHS", "GAPDHL"},
    "RPLP0": {"RPLP0", "RPLP0P1", "RPLP0P2", "RPLP0P3", "RPLP0P4", "RPLP0P5"},
    "TBP": {"TBP"},
    "HPRT1": {"HPRT1"},
    "SDHA": {"SDHA", "SDHAP1", "SDHAP2", "SDHAP3"},
    "GBP1": {"GBP1", "GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "GBP7"},
    "DHX58": {"DHX58", "DHX58P1", "DHX58P2"},
    "LDHA": {"LDHA", "LDHAL", "LDHAL6A", "LDHAL6B", "LDHB", "LDHC"},
    "ALDH3A2": {"ALDH3A2", "ALDH3A1", "ALDH3A3"},
    "ADH1B": {"ADH1B", "ADH1A", "ADH1C", "ADH4", "ADH1B", "ADH6", "ADH7"},
    "FAHD2B": {"FAHD2B", "FAHD2A", "FAHD1"},
}
# ═══ 旁系同源家族处理策略（逐基因显式声明，不靠隐式推断）═══
#
# 判据语义：只有 **同一基因座的功能重复拷贝**（真正的旁系同源，实验上
# 通常可接受并需在文中说明）才归为 PARALOG；其余同家族但独立基因座的
# 成员一律按真脱靶处理，因为交叉扩增会污染读出。
#
#   ACTB   → ACTG1/ACTG2/ACTC1：肌动蛋白家族，各基因座独立编码，
#            但作为内参属经典高一致性区域 → 计入 PARALOG（需说明）
#   GAPDH  → GAPDHS/GAPDHL：同上
#   RPLP0  → RPLP0P*：核糖体蛋白假基因/假基因家族 → PARALOG
#   TBP/HPRT1/SDHA：单拷贝，家族内无其他成员
#   LDHA   → LDHAL/LDHAL6A/LDHAL6B/LDHB/LDHC：乳酸脱氢酶家族，
#            各成员为独立基因座且均表达 → **计入真脱靶**（LDHB 与 LDHA
#            共用代谢通量，交叉扩增会直接污染读出）
#   ALDH3A2 → ALDH3A1/ALDH3A3：各自独立编码（Schiff 碱基形成酶 vs
#            脂肪酸醛脱氢酶）→ **真脱靶**
#   ADH1B  → ADH1A/ADH1C/ADH6：不同酶亚基基因座 → **真脱靶**
#   GBP1   → GBP2..GBP7：各 GBP 亚型为独立基因座且均可被干扰素诱导
#            → **真脱靶**（IFN-I 读出的核心风险）
#   DHX58  → DHX58P1/P2：假基因 → PARALOG/PSEUDO
#   FAHD2B → FAHD2A/FAHD1：FAHD2A 与 FAHD2B 属 FAHD2 家族，
#            FAHD2A 为独立基因座且在组织中广泛表达 → **真脱靶**
STRICT_PARALOG_ONLY = {"ACTB", "ACTG1", "GAPDH", "RPLP0", "TBP", "HPRT1",
                        "SDHA", "DHX58"}


def rc(s):
    return s.translate(RC_MAP)[::-1]


def blk(s, t):
    return re.findall(rf"<{t}>(.*?)</{t}>", s, re.S)


def parse_xml(path):
    x = open(path, encoding="utf-8", errors="ignore").read()
    out = []
    for hb in re.findall(r"<Hit>(.*?)</Hit>", x, re.S):
        acc = (blk(hb, "Hit_accession") or [""])[0]
        title = (blk(hb, "Hit_def") or [""])[0].strip()
        for hp in re.findall(r"<Hsp>(.*?)</Hsp>", hb, re.S):
            g = lambda t: (blk(hp, t) or [""])[0]      # noqa: E731
            idn = int(g("Hsp_identity") or 0)
            al = int(g("Hsp_align-len") or 0)
            qs, hs = g("Hsp_qseq"), g("Hsp_hseq")
            if not qs or not hs:
                continue
            out.append({
                "acc": acc, "title": title,
                "idn": idn, "al": al,
                "identity_pct": round(idn * 100.0 / al, 2) if al else 0.0,
                "qseq": qs, "hseq": hs,
                "human": "Homo sapiens" in title,
                "is_refseq_tx": bool(re.match(r"(NM_|XM_|NR_|XR_)", acc)),
                "pseudo": "pseudogene" in title.lower(),
            })
    return out


def gene_of(title):
    m = re.search(r"\(([A-Za-z0-9\-]{1,20})\)\s*(?:,|\s|$)", title)
    if m and not m.group(1).isdigit():
        return m.group(1)
    m2 = re.search(r"\b([A-Z][A-Z0-9\-]{1,14})\s*\(", title)
    if m2:
        return m2.group(1)
    m3 = re.search(r"Homo sapiens\s+(?:[A-Za-z0-9\-]+\s+){0,3}?"
                   r"([A-Z][A-Z0-9\-]{1,14})\b", title)
    return m3.group(1) if m3 else ""


def classify_hit(hit, gene, primer_len, fam, strict_flag):
    """返回 (级别, 目标基因名)。级别 ∈ {TARGET, PARALOG, PSEUDO, P0, P1, P2}"""
    if hit["al"] != primer_len:
        return ("PARTIAL", "")
    g = gene_of(hit["title"])
    qs, hs = hit["qseq"], hit["hseq"]
    L = min(len(qs), len(hs))
    mm_pos = [i for i in range(L) if qs[i] != hs[i] and hs[i] not in "- "]
    mm_total = len(mm_pos)
    mm_3p5 = sum(1 for i in mm_pos if i >= L - 5)
    hit["mm_total"] = mm_total
    hit["mm_3p5"] = mm_3p5
    hit["gene"] = g
    if g and g == gene:
        return ("TARGET", g)
    if hit["pseudo"]:
        return ("PSEUDO", g)
    # 只有声明为"严格旁系同源"的基因，其家族成员才归 PARALOG；
    # 其余基因的同家族成员按真脱靶（P0/P1）处理。
    if strict_flag and g and g in fam:
        return ("PARALOG", g)
    if not hit["is_refseq_tx"]:
        return ("P2", g)
    if mm_3p5 == 0 and hit["identity_pct"] >= 95:
        return ("P0", g)
    if mm_3p5 == 0 and hit["identity_pct"] >= 90:
        return ("P1", g)
    return ("P2", g)


def main():
    primers = {r["gene"]: r for _, r in
               pd.read_csv(f"{OUT}/qpcr_primers.csv").iterrows()}
    seq2p = {}
    for g, r in primers.items():
        seq2p[r["forward_primer_5to3"]] = (g, "forward")
        seq2p[r["reverse_primer_5to3"]] = (g, "reverse")

    detail_rows, per_primer = [], {}
    for path in sorted(glob.glob(f"{RAW}/core_nt_*.xml")):
        bname = os.path.basename(path)
        if "independent" in bname:
            continue
        hits = parse_xml(path)
        if not hits:
            continue
        q = hits[0]["qseq"]
        if q not in seq2p:
            continue
        gene, lab = seq2p[q]
        r = primers[gene]
        plen = len(q)
        fam = PARALOG_FAMILIES.get(gene, {gene})
        strict = gene in STRICT_PARALOG_ONLY
        counts = {"TARGET": 0, "PARALOG": 0, "PSEUDO": 0, "P0": 0,
                  "P1": 0, "P2": 0, "PARTIAL": 0}
        p0_genes, p1_genes = {}, {}
        tgt_full_tx = 0
        tgt_max_id = 0.0
        for h in hits:
            if not h["human"]:
                continue
            lvl, g = classify_hit(h, gene, plen, fam, strict)
            h["level"] = lvl
            counts[lvl] = counts.get(lvl, 0) + 1
            if lvl == "TARGET":
                if h["is_refseq_tx"]:
                    tgt_full_tx += 1
                tgt_max_id = max(tgt_max_id, h["identity_pct"])
            if lvl == "P0":
                p0_genes.setdefault(g or "?", []).append(h["acc"])
            if lvl == "P1":
                p1_genes.setdefault(g or "?", []).append(h["acc"])
            if lvl in ("P0", "P1", "TARGET", "PARALOG", "PSEUDO"):
                detail_rows.append({
                    "gene": gene, "primer": lab, "query": q,
                    "accession": h["acc"], "hit_gene": g,
                    "identity_n": h["idn"], "align_len": h["al"],
                    "identity_pct": h["identity_pct"],
                    "mm_total": h.get("mm_total"), "mm_3p5": h.get("mm_3p5"),
                    "level": lvl, "title": h["title"][:150],
                })
        key = f"{gene}_{lab}"
        per_primer[key] = {
            "gene": gene, "primer": lab, "query": q,
            "n_hits_total": len(hits),
            "n_hits_human": sum(1 for h in hits if h["human"]),
            "n_full_length_human": counts["TARGET"] + counts["PARALOG"]
                                    + counts["PSEUDO"] + counts["P0"]
                                    + counts["P1"] + counts["P2"],
            "n_partial_human": counts["PARTIAL"],
            "counts": counts,
            "target_full_match_refseq_tx": tgt_full_tx,
            "target_max_identity_pct": tgt_max_id,
            "P0_genes": {k: len(v) for k, v in sorted(p0_genes.items())},
            "P1_genes": {k: len(v) for k, v in sorted(p1_genes.items())},
            "P0_accessions": p0_genes,
            "usable": (tgt_full_tx >= 1 and counts["P0"] == 0),
        }

    df = pd.DataFrame(detail_rows)
    df.to_csv(f"{OUT}/blast_offtarget_adjudication.csv", index=False)
    with open(f"{OUT}/blast_offtarget_adjudication.json", "w", encoding="utf-8") as fh:
        json.dump({"date": "2026-10-25",
                   "script": "02_脚本/adjudicate_blast_offtargets.py",
                   "criteria": {
                       "P0": "全长匹配 + 人源 RefSeq 转录本 + 非目标 + 非旁系同源"
                             " + 非假基因 + 3'端末5nt 全匹配 + identity ≥ 95%",
                       "P1": "同上但 90% ≤ identity < 95%",
                       "P2": "其余（3'端有错配或 identity < 90%）→ 判为噪声",
                       "PARTIAL": "align_len < 引物长度 → 无法完全退火，判为噪声",
                   },
                   "per_primer": per_primer},
                  fh, ensure_ascii=False, indent=1, default=str)

    print(f"{'primer':16s} {'full':>6s} {'part':>6s} {'tgt_tx':>7s} "
          f"{'P0':>4s} {'P1':>4s} {'para':>5s}  verdict")
    for k, v in per_primer.items():
        c = v["counts"]
        print(f"{k:16s} {v['n_full_length_human']:6d} "
              f"{v['n_partial_human']:6d} {v['target_full_match_refseq_tx']:7d} "
              f"{c['P0']:4d} {c['P1']:4d} {c['PARALOG']:5d}  "
              f"{'USABLE' if v['usable'] else ('NO_TARGET' if v['target_full_match_refseq_tx']==0 else 'OFFTARGET')}"
              f"  P0genes={v['P0_genes']} P1genes={v['P1_genes']}")
    print(f"\n明细行数 {len(df)} → {OUT}/blast_offtarget_adjudication.csv")
    print("DONE — adjudicate_blast_offtargets.py")


if __name__ == "__main__":
    main()
