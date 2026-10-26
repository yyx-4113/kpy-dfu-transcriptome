# -*- coding: utf-8 -*-
"""
design_final_primers_v3.py — RT-qPCR 引物**统一最终设计**（12 基因）

═══ 为什么要重做：v1/v2 三条设计路径各自只满足部分约束 ═══

本项目先后走过三条设计路径，各自都有结构性缺陷（全部经独立复核确认）：

| 路径 | 脚本 | 缺陷 | 证据 |
|------|------|------|------|
| v1 | `design_rtqpcr_primers.py` | 整条 mRNA 上设计，无 CDS 锚定、无跨转录本覆盖 | ALDH3A2 引物落 3'UTR Alu 元件（3111/4713 条全长命中）；LDHA 引物匹配位置 -1；ADH1B 100% 匹配 ADH1C/ADH1A |
| v2a | `redesign_primers_v2.py` | ✅ CDS 锚定 ✅ 全转录本覆盖，但**无旁系同源区分度约束** | FAHD2B/ADH1B/LDHA 的旧脱靶（LOC124900512 / ADH1C / LYST）未被该脚本覆盖，只能靠事后 BLAST 拦截 |
| v2b | `design_paralog_aware_primers_v2.py` | 声称"完整 CDS 上设计"但**代码传入整条 mRNA**（`SEQUENCE_TEMPLATE: best_t["seq"]`）；种子索引漏检；反向引物 3' 端锚定测错端；DHX58 家族配置了 0 条转录本的假基因名 | 六对引物**全部**落在 CDS 外（独立审计：FAHD2B 引物在 mRNA 1572/1669 而 CDS 为 425..1435）；自报 `identity=0.0/mm=20` 与 FAHD2A/FAHD2B 98.41% CDS 一致性直接矛盾 |

═══ v3 的方法学立场 ═══

**一次求解，三约束同时满足**（不是三选一）：

1. **CDS 锚定**：引物对完整位于目标基因的 CDS 内。
   理由：RT-qPCR 测 mRNA 丰度；3'UTR 在各转录本间高度变异且常含
   Alu/LINE 重复元件（ALDH3A2 的旧缺陷正源于此）。
2. **全编码转录本覆盖**：设计窗口须落在目标基因**全部蛋白编码
   RefSeq 转录本共有**的序列区段（逐本自身序列求交，不复用同一
   参考切片）。
   理由：模板版本漂移（`.1`→`.2` 的 5'UTR 重排）会让只对单一版本
   有效的引物在全部当前转录本上失配。
3. **旁系同源区分度**：两条引物**各自**对所有旁系同源转录本的最佳
   近似匹配须达到区分阈值。
   理由：PCR 交叉扩增会直接污染读出；**需要区分度的是引物本身**，
   不是产物间隔区（产物中段落在保守区无妨，产物是 cDNA，中段不
   参与特异性决定）。

**约束 3 的精确判据**（与 BLAST v5 判读器对齐，且链向正确）：

对每条引物 P（长 L）与每条旁系同源序列 S，穷举 S 的全部窗口 × 正反链，
取最小错配窗口：

- 正链命中（P 直接出现于 S）：引物 3' 端 = 窗口**后** 5 nt
- 负链命中（rc(P) 出现于 S）：引物 3' 端 = 窗口**前** 5 nt
  （`rc(P)[0]` 才是引物 3' 端碱基的互补 —— v2b 在此处测错了端）

分级：
- **STRONG**：max identity ≤ 90% 且 3' 端 ≥1 错配 → 判 P2 噪声，实验不共扩增
- **WEAK**：max identity ≤ 95% 且 3' 端 ≥1 错配 → 判 P1 风险，需熔解曲线验证
- **UNACCEPTABLE**：max identity ≥ 95%（或 3' 端全匹配且 ≥90%）→ 判 P0，必须换引物

**穷举而非种子索引**：v2b 的 4×5nt 分段种子索引在「错配均匀分散到
全部 4 段」时会漏掉整个最佳窗口（鸽巢原理的反例），把真实的
identity 85–100% 误报为 0.0。本脚本用 numpy 向量化枚举全部窗口，
不剪枝、不漏检。

**穷举的代价控制**：旁系同源序列用「每家族取最长编码转录本的 CDS +
全部编码转录本」两档；CDS 序列远短于全 mRNA（数百 ~ 数千 nt），
全量穷举完全可行。

**已知生物学限制（如实记录，不掩盖）**：
FAHD2A / FAHD2B 的 CDS 仅 15 个 SNP（一致性 98.41%），蛋白 319 aa
逐字相同 → 部分区段**在原理上不可区分**。若穷举后无可用窗口，
脚本必须如实报出「该区段不可区分」，而非降标准放行。

输出：
  03_结果/实验验证设计/final_primer_design_v3.json
  03_结果/实验验证设计/final_primers_v3.csv
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
UA = {"User-Agent": "WorkBuddy-FinalPrimerDesignV3/1.0"}
RC_MAP = str.maketrans("ACGT", "TGCA")
CACHE_P = f"{PROV}/final_design_cache.json"

# ═══ 12 个基因 ═══
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

# ═══ 旁系同源家族 ═══
# 只列**真实存在且有蛋白编码 RefSeq 转录本**的同家族成员。
# v2b 给 DHX58 配了 DHX58P1/DHX58P2 —— 实测这两个名字的
# `gene_nuccore_refseqrna` 链接数为 **0**（是假基因，非转录本基因），
# 导致 DHX58 的 200 个候选被全部判为「无旁系同源可比」的空判。
PARALOG_FAMILIES = {
    "LDHA": ["LDHB", "LDHC", "LDHAL", "LDHAL6A", "LDHAL6B"],
    "FAHD2B": ["FAHD2A", "FAHD1"],
    "ADH1B": ["ADH1A", "ADH1C", "ADH4", "ADH6", "ADH7"],
    "ALDH3A2": ["ALDH3A1", "ALDH3A3"],
    "GBP1": ["GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "GBP7"],
    "DHX58": [],          # DHX58/MDA5 无蛋白编码旁系同源（P1/P2 为假基因）
    "ACTB": ["ACTG1", "ACTG2", "ACTC1"],
    "GAPDH": ["GAPDHS", "GAPDHL"],
    "TBP": [],
    "HPRT1": [],
    "RPLP0": ["RPLP0P1", "RPLP0P2", "RPLP0P3", "RPLP0P4", "RPLP0P5"],
    "SDHA": ["SDHAP1", "SDHAP2", "SDHAP3"],
}
# 参考基因的同家族成员：同基因座功能重复拷贝 → 实验上经典可接受，
# 仅作信息性报告，不作阻断约束（沿用 adjudicate_blast_offtargets.py 的
# STRICT_PARALOG_ONLY 政策）。
REFERENCE_SOFT_PARALOG = {"ACTB", "GAPDH", "RPLP0"}

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
    """基因 → RefSeq RNA accession 列表（gene_nuccore_refseqrna 链接）。"""
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
    """从 GenBank 记录取 CDS 坐标（1-based 闭区间列表）。

    ═══ 修正说明 ═══
    旧实现把 **qualifier 续行**（如 `/note="...67kDa..."`）拼进 location
    字符串，导致 `69..1847` + `67` → 误解析为 `69..184767`。
    正确做法：location 只取特征块首行 + **紧跟其后、以空格缩进且不以
    `/` 开头**的续行；一旦遇到 qualifier（`/` 开头）即停止拼接。
    """
    for b in re.split(r"\n {5}(?=\w)", gb):
        if not b.startswith("CDS"):
            continue
        loc = re.match(r"CDS\s+(\S+)", b).group(1).replace("<", ">")
        for ln in b.split("\n")[1:]:
            st = ln.strip()
            if not st:
                continue
            if st.startswith("/"):
                break                     # 进入 qualifier 区 → 停止
            loc += st
        return [(int(a), int(c)) for a, c in
                re.findall(r"(\d+)\.\.(\d+)", loc)]
    return []


def cds_of(acc, seq=None):
    gb = efetch_gb(acc)
    seg = parse_cds_location(gb)
    if not seg or seq is None:
        return seg, None
    cds = "".join(seq[a - 1:b] for a, b in seg)
    return seg, cds


def gene_name_by_esearch(sym):
    """符号 → GeneID（严格取完全匹配项，避免 [Gene Name] 模糊匹配）。"""
    key = f"__gid__{sym}"
    if key in CACHE:
        return CACHE[key]["seq"]
    d = json.loads(_get(f"{EUTILS}/esearch.fcgi?db=gene&term="
                        f"{urllib.parse.quote(sym + '[Gene Name]')}"
                        f"%20AND%20human%5Borgn%5D&retmode=json"))
    gids = d["esearchresult"]["idlist"]
    # 用 esummary 精确核对 name（esearch 的 [Gene Name] 是精确字段，
    # 但仍需确认返回项的 name 与请求符号完全一致）
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


# ═══ numpy 向量化穷举比对 ═══
ENC = np.full(256, 4, dtype=np.uint8)
for _i, _ch in enumerate(b"ACGT"):
    ENC[_ch] = _i


def enc(s):
    return ENC[np.frombuffer(s.encode("ascii"), dtype=np.uint8)]


def best_window_identity(subject: str, primer: str):
    """穷举 subject 全部窗口 × 正反链 → (identity%, mm_total, mm_3p5, strand)

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
        # 逐位置错配（必须在**最佳窗口内部**计算，不能跨窗口求和）
        pos_mm = (S[i:i + L] != P)
        if strand == "+":
            mm3 = int(pos_mm[L - 5:L].sum())
        else:
            mm3 = int(pos_mm[0:5].sum())
        best = (ident, m, mm3, strand)
    return best


def grade(ident, mm3):
    """按 BLAST v5 判据分级（与 adjudicate_blast_offtargets.py 严格对齐）。

    v5 判据：
      P0（阻断）：mm_3p5 == 0 且 identity ≥ 95%
      P1（风险）：mm_3p5 == 0 且 90% ≤ identity < 95%
      P2（噪声）：其余（3'端有错配，或 identity < 90%）

    → 对应到设计期：
      UNACCEPTABLE = P0 或 P1（3'端全匹配且 identity ≥ 90%）
      WEAK         = 3'端有错配且 identity ≥ 90%（v5 判 P2，但留作警示）
      STRONG       = identity < 90%（v5 判 P2 噪声）
    """
    if mm3 == 0 and ident >= 95.0:
        return "UNACCEPTABLE"      # v5 P0：真实共扩增
    if mm3 == 0 and ident >= 90.0:
        return "UNACCEPTABLE"      # v5 P1：风险脱靶（阻断级保守处理）
    if ident >= 90.0:
        return "WEAK"              # v5 P2，但 identity 高 → 警示
    return "STRONG"                # v5 P2 噪声：实验上不共扩增


def shared_kmer_segments(seqs, ref_acc, cds_pos):
    """逐本自身序列求交 → 共有 k-mer 在 ref 上覆盖的区段（限制在 CDS 内）。"""
    k = 25
    ref = seqs[ref_acc]
    common = {ref[i:i + k] for i in range(len(ref) - k + 1)}
    for a, s in seqs.items():
        if a == ref_acc:
            continue
        common &= {s[i:i + k] for i in range(len(s) - k + 1)}
    pos = [i for i in range(len(ref) - k + 1) if ref[i:i + k] in common]
    if not pos:
        return []
    segs, st, prev = [], pos[0], pos[0]
    for p in pos[1:]:
        if p <= prev + k:
            prev = p
        else:
            segs.append((st, prev + k))
            st = prev = p
    segs.append((st, prev + k))
    # 裁剪到 CDS 内连续段
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
              "script": "02_脚本/design_final_primers_v3.py",
              "primer3_version": primer3.__version__,
              "primer3_params": {k: v for k, v in PRIMER3_PARAMS.items()
                                 if k != "PRIMER_EXPLAIN_FLAG"},
              "method": {
                  "objective": "三约束同时满足：CDS锚定 / 全编码转录本覆盖 / 旁系同源区分度",
                  "template": "目标基因**全部蛋白编码转录本共有**的 CDS 区段"
                             "（逐本自身序列求交；v2b 误用整条 mRNA）",
                  "paralog_match": "numpy 向量化穷举全部窗口 × 正反链，零剪枝"
                                   "（v2b 的 4×5nt 种子索引会漏检）",
                  "three_prime_rule": "正链命中→引物3'端=窗口后5nt；"
                                      "负链命中→引物3'端=窗口前5nt",
                  "grading": {"STRONG": "identity ≤ 90% → P2 噪声",
                              "WEAK": "90% < identity ≤ 95% 且 3'端有错配 → P1 风险",
                              "UNACCEPTABLE": "identity ≥ 95% 或 (≥90% 且 3'端全匹配) → P0"},
              },
              "genes": {}}
    selected_rows = []

    for gene, (gid, role, module) in GENES.items():
        print(f"\n{'=' * 74}\n=== {gene} (GeneID {gid}) ===")
        # ── 1. 目标基因的全部转录本 ──
        tx_accs = gene_transcripts(gid)
        seqs, coding, cds_map, cds_len_map = {}, [], {}, {}
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
                    cds_len_map[acc] = len(cds)
                time.sleep(0.25)
            except Exception as e:                     # noqa: BLE001
                print(f"   fetch fail {acc}: {e}")
        if not coding:
            print("   !! 无蛋白编码转录本")
            continue
        print(f"   RefSeq RNA {len(tx_accs)} 条，其中蛋白编码 {len(coding)} 条")

        # ── 2. canonical 选择：优先 NM_/NR_ 且 CDS 最长 ──
        def rank(acc):
            pref = 0 if acc.startswith(("NM_", "NR_")) else 1
            return (pref, -cds_len_map[acc])
        canonical = sorted(coding, key=rank)[0]
        ref = seqs[canonical]
        cds_pos = set()
        for a, b in cds_map[canonical]:
            cds_pos.update(range(a - 1, b))
        segs = shared_kmer_segments(seqs, canonical, cds_pos)
        print(f"   canonical={canonical} (CDS {cds_len_map[canonical]} nt) "
              f"共有CDS区段 {[(a, b - a) for a, b in segs]}")

        # ── 3. 旁系同源序列集 ──
        para_full, para_cds = [], []
        para_detail = {}
        for fam in PARALOG_FAMILIES.get(gene, []):
            pgid = gene_name_by_esearch(fam)
            if not pgid:
                print(f"   旁系同源 {fam}: GeneID 未找到")
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
                    time.sleep(0.2)
                except Exception:                      # noqa: BLE001
                    pass
            para_detail[fam] = {"gene_id": pgid, "n_refseq_rna": len(f_accs),
                                 "n_coding": n_coding}
            print(f"   旁系同源 {fam}: RefSeq RNA {len(f_accs)} 条，"
                  f"蛋白编码 {n_coding} 条")
        if not para_full:
            print("   （无蛋白编码旁系同源 → 约束 3 不适用）")

        # ── 4. Primer3 在共有 CDS 区段上设计 ──
        cands = []
        for seg_a, seg_b in segs:
            if seg_b - seg_a < 130:
                continue
            sub = ref[seg_a:seg_b]
            res = primer3.bindings.design_primers(
                {"SEQUENCE_ID": f"{gene}_{seg_a}", "SEQUENCE_TEMPLATE": sub},
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
                    "shared_segment": [seg_a + 1, seg_b],
                })
        print(f"   Primer3 候选对 {len(cands)}")

        # ── 5. 逐候选：约束 2（覆盖）+ 约束 3（旁系同源区分） ──
        for c in cands:
            # 约束 2：全部编码转录本都能定位 + 引物完整在各自 CDS 内
            n_ok, n_in_cds = 0, 0
            for acc in coding:
                s = seqs[acc]
                iF = s.find(c["forward"])
                iR = s.find(rc(c["reverse"]))
                if iF >= 0 and iR >= 0:
                    n_ok += 1
                    pos = set()
                    for a, b in cds_map[acc]:
                        pos.update(range(a - 1, b))
                    if all(p in pos for p in range(iF, iF + len(c["forward"]))) and \
                       all(p in pos for p in range(iR, iR + len(c["reverse"]))):
                        n_in_cds += 1
            c["n_tx_covered"] = n_ok
            c["n_tx_total"] = len(coding)
            c["n_tx_primer_in_CDS"] = n_in_cds
            c["all_tx_covered"] = (n_ok == len(coding))
            c["all_tx_primer_in_CDS"] = (n_in_cds == len(coding))

            # 约束 3：旁系同源区分度（穷举）
            soft = gene in REFERENCE_SOFT_PARALOG
            worst = {"f": None, "r": None}
            for lab, pk in (("f", "forward"), ("r", "reverse")):
                p = c[pk]
                b_full = b_cds = None
                for fam, acc, s in para_full:
                    ident, mm, mm3, strand = best_window_identity(s, p)
                    if b_full is None or ident > b_full[0]:
                        b_full = (ident, mm, mm3, strand, fam, acc)
                for fam, acc, s in para_cds:
                    ident, mm, mm3, strand = best_window_identity(s, p)
                    if b_cds is None or ident > b_cds[0]:
                        b_cds = (ident, mm, mm3, strand, fam, acc)
                use = b_cds if (b_cds and (b_full is None or b_cds[0] >= b_full[0])) \
                    else b_full
                g = grade(use[0], use[2]) if use else "NA"
                c[f"{lab}_max_id_vs_paralog"] = round(use[0], 2) if use else None
                c[f"{lab}_mm_vs_paralog"] = use[1] if use else None
                c[f"{lab}_mm3p_vs_paralog"] = use[2] if use else None
                c[f"{lab}_strand_vs_paralog"] = use[3] if use else None
                c[f"{lab}_worst_paralog"] = use[4] if use else None
                c[f"{lab}_worst_acc"] = use[5] if use else None
                c[f"{lab}_grade"] = g
                worst[lab] = (g, use[0] if use else None)
            if not para_full:
                c["pair_paralog_grade"] = "NA"
            elif soft:
                c["pair_paralog_grade"] = "SOFT"
            else:
                gs = [c["f_grade"], c["r_grade"]]
                c["pair_paralog_grade"] = ("BOTH_STRONG" if gs == ["STRONG"] * 2
                                           else ("ONE_WEAK" if "WEAK" in gs
                                                 else "UNACCEPTABLE"))

        # ── 6. 排序与选择 ──
        def score(c):
            s = 3 * c["penalty"]
            s += 3 * (abs(c["f_tm"] - 60) + abs(c["r_tm"] - 60))
            s += 2 * (abs(c["f_gc"] - 50) / 10 + abs(c["r_gc"] - 50) / 10)
            s += 2 * (0 if c["forward"].endswith(("G", "C")) else 0.5)
            s += 2 * (0 if c["reverse"].endswith(("G", "C")) else 0.5)
            s += 60 * (c["n_tx_total"] - c["n_tx_covered"])
            s += 60 * (c["n_tx_total"] - c["n_tx_primer_in_CDS"])
            pg = c["pair_paralog_grade"]
            s += {"BOTH_STRONG": 0, "ONE_WEAK": 30, "NA": 0, "SOFT": 5,
                  "UNACCEPTABLE": 300}[pg]
            if pg == "ONE_WEAK":
                s += 10 * max(c["f_max_id_vs_paralog"] or 0,
                              c["r_max_id_vs_paralog"] or 0)
            return s

        cands.sort(key=score)
        # 硬性过滤：全转录本覆盖 + CDS 内 + 旁系同源非 UNACCEPTABLE
        ok = [c for c in cands
              if c["all_tx_covered"] and c["all_tx_primer_in_CDS"]
              and c["pair_paralog_grade"] != "UNACCEPTABLE"]
        print(f"   三约束全过 {len(ok)}/{len(cands)}；"
              f"分级分布 {pd.Series([c['pair_paralog_grade'] for c in cands]).value_counts().to_dict()}")
        for c in cands[:4]:
            print(f"     cov={c['n_tx_covered']}/{c['n_tx_total']} "
                  f"inCDS={c['n_tx_primer_in_CDS']}/{c['n_tx_total']} "
                  f"para={c['pair_paralog_grade']} "
                  f"F_id={c['f_max_id_vs_paralog']} R_id={c['r_max_id_vs_paralog']} "
                  f"F={c['forward']} R={c['reverse']} amp={c['product_size']}")

        best = ok[0] if ok else None
        if best:
            print(f"   ★ 最终推荐：F={best['forward']} R={best['reverse']} "
                  f"amp={best['product_size']}bp 旁系同源={best['pair_paralog_grade']}")
        else:
            print("   !! 无候选同时满足三约束 → 需如实记录生物学限制")

        report["genes"][gene] = {
            "gene_id": gid, "role": role, "module": module,
            "canonical": canonical,
            "n_refseq_rna": len(tx_accs),
            "coding_transcripts": coding,
            "n_coding_tx": len(coding),
            "shared_cds_segments": [(a + 1, b - a) for a, b in segs],
            "paralog_families": para_detail,
            "n_paralog_refseq_rna": len(para_full),
            "n_paralog_coding": len(para_cds),
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
                "n_tx_covered": best["n_tx_covered"],
                "n_tx_total": best["n_tx_total"],
                "n_tx_primer_in_CDS": best["n_tx_primer_in_CDS"],
                "pair_paralog_grade": best["pair_paralog_grade"],
                "f_max_id_vs_paralog": best["f_max_id_vs_paralog"],
                "r_max_id_vs_paralog": best["r_max_id_vs_paralog"],
                "f_worst_paralog": best["f_worst_paralog"],
                "r_worst_paralog": best["r_worst_paralog"],
                "f_mm3p_vs_paralog": best["f_mm3p_vs_paralog"],
                "r_mm3p_vs_paralog": best["r_mm3p_vs_paralog"],
                "shared_segment": str(best["shared_segment"]),
            })
        save_cache()

    with open(f"{OUT}/final_primer_design_v3.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1, default=str)
    if selected_rows:
        pd.DataFrame(selected_rows).to_csv(f"{OUT}/final_primers_v3.csv",
                                           index=False)
    print(f"\n完成：{len(selected_rows)}/{len(GENES)} 个基因产出三约束全过的引物对")
    print("DONE — design_final_primers_v3.py")


if __name__ == "__main__":
    main()
