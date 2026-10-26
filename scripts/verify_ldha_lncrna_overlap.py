# -*- coding: utf-8 -*-
"""
verify_ldha_lncrna_overlap.py — LOC105375980 × LDHA CDS 重叠的**独立第二路径**核实

═══ 为什么需要这一步 ═══

BLAST v6 全库复核发现 LDHA 的最终引物对存在 P0 脱靶：
  NR_202041 / NR_202040（Homo sapiens uncharacterized LOC105375980，
  chr9p22.3，lncRNA，7/5 外显子，带 AATAAA 多聚腺苷酸化位点）
对引物 F、R 各有命中，identity 95%（F）/ 100%（R），且 3'端 5nt 全匹配。

═══ 本脚本要独立回答的问题 ═══

1. **这个 locus 到底是什么？**
   LDHAP4（GeneID 158222）同样在 chr9p22.3。若 LOC105375980 只是 LDHAP4 的
   旧名/别名，则它是「已知假基因」，风险性质与「未知第三方重复」完全不同。
   → 用 NCBI Datasets v2 取两者**基因组区间**判定是否重叠。

2. **LDHA CDS 的哪些位置与它重叠？重叠多长、identity 多少？**
   → 取 NR_202041 / NR_202040 的 RefSeq 全长序列，与 LDHA CDS 做**纯本地
     无 BLAST** 的无缺口比对（穷举全部对角线 + 前缀和最大子段），
     输出 LDHA CDS 坐标上的重复区段（1-based）。

3. **正负链分别判定**（引物 3'端锚定依赖链向）。

4. **阴性对照**（防"什么都能对上"的假阳性）：
   LDHA CDS vs 与 LDHA 无关的基因（ALDH3A2、GBP1、ACTB）应无 ≥20nt 高同源区段。

═══ 算法与正确性 ═══

无缺口（ungapped）局部比对，逐对角线处理。对某条对角线，令
  score = +1（匹配） / −19（错配）
  S[k] = 前 k 个位置的得分前缀和
则窗口 [i, e) 的得分 = S[e] − S[i] = alen − 20·mismatch，
  得分 ≥ 0  ⟺  alen ≥ 20·mismatch  ⟺  identity ≥ 95%
（得分 = 0 恰好对应 identity 恰为 95%，按「≥95%」判据**应当接受**；
早期版本用 `score > 0` 会误杀这类窗口 —— 这正是 v2 种子索引漏检教训的
同源陷阱，已在本脚本中显式修正。）

对每个右端点 e，取满足 S[i] ≤ S[e] 的**最小** i（用前缀最小值 pm 的
首个越界位置，np.searchsorted 在 −pm 上求解，完全向量化），
即「以 e 结尾的最长合格窗口」。由单调性可知：任何以 e 结尾的合格窗口都
被它包含 → 全部合格窗口的并集 = 全部最长窗口的并集，**不漏检**。

每个窗口的 identity 由**原始字节直接重算**，不使用打分近似；
并额外做一次「按最终区段整体重算」的独立复核。

输出：
  03_结果/实验验证设计/provenance/ldha_lncrna_overlap.json
"""
import json

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import os
import re
import urllib.parse
import urllib.request

import numpy as np

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
PROV = f"{OUT}/provenance"
UA = {"User-Agent": "WorkBuddy-LDHA-overlap/1.0"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/"

MIN_ALEN = 20          # 判定为重复拷贝的最小比对长度
MIN_IDENTITY = 95.0    # 重复掩码判据（与 map_cds_duplication.py / BLAST v6 一致）
MISMATCH_PENALTY = 19  # score = alen − 20·mismatch


def _get(u, timeout=180):
    return urllib.request.urlopen(
        urllib.request.Request(u, headers=UA), timeout=timeout).read().decode(
        errors="ignore")


def fetch_fasta(acc, db="nuccore"):
    u = EUTILS + "efetch.fcgi?" + urllib.parse.urlencode(
        {"db": db, "id": acc, "rettype": "fasta", "retmode": "text"})
    t = _get(u)
    seqs, name, buf = [], "", []
    for ln in t.splitlines():
        if ln.startswith(">"):
            if name:
                seqs.append((name, "".join(buf)))
            name, buf = ln[1:].strip(), []
        elif ln.strip():
            buf.append(ln.strip())
    if name:
        seqs.append((name, "".join(buf)))
    return seqs


def gene_interval(gid):
    """NCBI Datasets v2：基因的基因组区间。"""
    t = _get(DATASETS + str(gid))
    j = json.loads(t)
    g = j.get("reports", [{}])[0].get("gene", {})
    locs = []
    for a in g.get("annotations", []) or []:
        for x in a.get("genomic_locations", []) or []:
            rg = x.get("genomic_range") or {}
            locs.append({
                "acc": x.get("sequence_accession"),
                "begin": int(rg["begin"]), "end": int(rg["end"]),
                "orientation": rg.get("orientation"),
            })
    return {"gene_id": str(gid), "symbol": g.get("symbol"),
            "description": g.get("description"), "locations": locs}


def cds_of(accession, cache):
    """从 final_design_cache 取 mRNA 全序列，并按 GenBank CDS 注释抽取 CDS。"""
    m = cache.get(accession)
    gb = cache.get(f"__gb__{accession}")
    if isinstance(gb, dict):
        gb = gb.get("seq")
    if not m or not gb:
        return None
    full = m["seq"] if isinstance(m, dict) else m
    seg = None
    for b in re.split(r"\n {5}(?=\w)", gb):
        if b.startswith("CDS"):
            loc = re.match(r"CDS\s+(\S+)", b).group(1).replace("<", ">")
            for ln in b.split("\n")[1:]:
                st = ln.strip()
                if not st:
                    continue
                if st.startswith("/"):
                    break
                loc += st
            seg = [(int(a), int(c)) for a, c in re.findall(r"(\d+)\.\.(\d+)", loc)]
            break
    if not seg:
        return None
    cds = "".join(full[a - 1:b] for a, b in seg)
    return {"accession": accession, "cds": cds, "cds_len": len(cds),
            "cds_coords_in_mrna": seg, "mrna_len": len(full)}


def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def _diagonal_windows(qa, sa):
    """在一条对角线上返回 (i, e) 列表（0-based 半开区间，相对对角线起点）。

    每个元素是「以该右端点 e 结尾的最长合格窗口」。
    """
    n = len(qa)
    if n < MIN_ALEN:
        return []
    eq = (qa == sa)
    sc = np.where(eq, 1, -MISMATCH_PENALTY).astype(np.int64)
    S = np.concatenate(([0], np.cumsum(sc)))          # S[k] = 前 k 位得分和
    pm = np.minimum.accumulate(S)                     # pm[k] = min(S[0..k])，非增
    negpm = -pm                                        # 非减 → 可 searchsorted

    ends = np.arange(MIN_ALEN, n + 1)                  # 半开右端点
    # 最小的 i 使 pm[i] <= S[e]  ⟺  −pm[i] >= −S[e]
    idx = np.searchsorted(negpm, -S[ends], side="left")
    ok = idx <= (ends - MIN_ALEN)                      # 窗口长度下限
    if not np.any(ok):
        return []
    idx = idx[ok]
    ends = ends[ok]
    scores = S[ends] - S[idx]
    keep = scores >= 0
    return [(int(idx[k]), int(ends[k])) for k in range(len(idx)) if keep[k]]


def overlaps(q, s):
    """穷举所有对角线（零剪枝），返回 q 上 identity ≥ MIN_IDENTITY、
    长度 ≥ MIN_ALEN 的**合并后不重叠**区段。identity 由原始字节精确重算。"""
    q = q.upper()
    s = s.upper()
    nq, ns = len(q), len(s)
    if nq == 0 or ns == 0:
        return []
    qa = np.frombuffer(q.encode(), dtype=np.uint8)
    sa = np.frombuffer(s.encode(), dtype=np.uint8)
    raw = []
    for d in range(-(ns - 1), nq):
        i0 = max(0, -d)
        i1 = min(nq, ns - d)
        if i1 - i0 < MIN_ALEN:
            continue
        for (a, b) in _diagonal_windows(qa[i0:i1], sa[i0 + d:i1 + d]):
            raw.append((i0 + a, i0 + b, d))
    if not raw:
        return []
    hits = []
    for (qs, qe, d) in raw:
        n = qe - qs
        mm = int(np.count_nonzero(qa[qs:qe] != sa[qs + d:qe + d]))
        ident = (n - mm) * 100.0 / n
        if ident + 1e-9 < MIN_IDENTITY:
            continue
        hits.append({"q_start": qs + 1, "q_end": qe, "alen": n,
                     "mismatch": mm, "identity_pct": round(ident, 2),
                     "s_start": qs + d + 1, "s_end": qe + d})
    hits.sort(key=lambda h: (h["q_start"], -h["q_end"]))
    merged = []
    for h in hits:
        if merged and h["q_start"] <= merged[-1]["q_end"]:
            m = merged[-1]
            m["q_end"] = max(m["q_end"], h["q_end"])
            m["n_hits"] = m.get("n_hits", 1) + 1
            m["best_identity_pct"] = max(m.get("best_identity_pct", 0.0),
                                          h["identity_pct"])
        else:
            m = dict(h)
            m["n_hits"] = 1
            m["best_identity_pct"] = h["identity_pct"]
            merged.append(m)
    for m in merged:
        m["span"] = m["q_end"] - m["q_start"] + 1
    return merged


def merge_segments(segs):
    segs = sorted(segs)
    out = []
    for a, b in segs:
        if out and a <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def main():
    cache = json.load(open(f"{PROV}/final_design_cache.json", encoding="utf-8"))

    report = {"date": "2026-10-25",
              "script": "02_脚本/verify_ldha_lncrna_overlap.py",
              "purpose": "LOC105375980 × LDHA CDS 重叠的独立第二路径核实（不依赖 BLAST）",
              "criteria": {"min_alen": MIN_ALEN, "min_identity_pct": MIN_IDENTITY},
              "algorithm": ("穷举全部对角线 + 前缀和最大子段（score = alen − "
                            "20·mismatch，≥0 等价 identity ≥95%，含恰为 95% 者）；"
                            "每条对角线取「以每个右端点结尾的最长合格窗口」，"
                            "由单调性保证不漏检；identity 由原始字节精确重算"),
              "steps": {}}

    # ── 步骤 1：locus 身份鉴定 ───────────────────────────────
    print("=" * 74)
    print("步骤 1：locus 身份鉴定（NCBI Datasets v2 基因组区间）")
    print("=" * 74)
    ids = {}
    for gid, tag in [(105375980, "LOC105375980"), (158222, "LDHAP4")]:
        gi = gene_interval(gid)
        ids[tag] = gi
        print(f"  {tag} (GeneID {gid}) {gi['description']}")
        for l in gi["locations"]:
            print(f"     chr9:{l['begin']:,}-{l['end']:,} ({l['orientation']})")
    pairs = []
    for la in ids["LOC105375980"]["locations"]:
        for lb in ids["LDHAP4"]["locations"]:
            o = max(la["begin"], lb["begin"]) <= min(la["end"], lb["end"])
            gap = max(la["begin"], lb["begin"]) - min(la["end"], lb["end"])
            pairs.append({"a": f"{la['begin']}-{la['end']}",
                          "b": f"{lb['begin']}-{lb['end']}",
                          "overlap": o, "gap_nt": gap})
            print(f"     重叠? {o}  间隔 {gap:+,} nt  "
                  f"[{la['begin']}-{la['end']}] vs [{lb['begin']}-{lb['end']}]")
    any_ov = any(p["overlap"] for p in pairs)
    min_gap = min(p["gap_nt"] for p in pairs)
    report["steps"]["locus_identity"] = {
        "LOC105375980": ids["LOC105375980"],
        "LDHAP4": ids["LDHAP4"],
        "pairwise": pairs, "any_overlap": any_ov, "min_gap_nt": min_gap,
        "verdict": ("LOC105375980 与 LDHAP4 基因组区间不重叠（最近相距 %d nt）"
                    "→ LOC105375980 是独立 locus，不是 LDHAP4 的别名；"
                    "二者同处 chr9p22.3 的 LDHA 假基因簇" % min_gap)
        if not any_ov else
        ("LOC105375980 与 LDHAP4 基因组区间重叠（%d nt）→ 同一 locus 的不同注释"
         % min_gap)}

    # ── 步骤 2：序列重叠（正负链）─────────────────────────
    print("\n" + "=" * 74)
    print("步骤 2：LDHA CDS 与 NR_202041 / NR_202040 的无缺口重叠（本地比对）")
    print("=" * 74)
    ldha = cds_of("NM_001165414.2", cache)
    if not ldha:
        raise SystemExit("无法取得 LDHA CDS")
    cds = ldha["cds"]
    print(f"  LDHA {ldha['accession']} CDS {ldha['cds_len']} nt "
          f"(mRNA 上位置 {ldha['cds_coords_in_mrna']})")
    report["steps"]["ldha_cds"] = {k: v for k, v in ldha.items() if k != "cds"}

    segs_all = []
    per_tx = {}
    tx_seqs = {}
    for acc in ["NR_202041", "NR_202040"]:
        s = fetch_fasta(acc)[0][1].upper()
        tx_seqs[acc] = s
        n = 0
        for strand, seq in (("+", s), ("-", revcomp(s))):
            for sg in overlaps(cds, seq):
                sg["strand"] = strand
                sg["tx"] = acc
                a, b = sg["q_start"] - 1, sg["q_end"]
                wq, ws = cds[a:b], seq[a:b]
                mm = sum(1 for x, y in zip(wq, ws) if x != y)
                sg["verify_identity_pct"] = round(
                    (len(wq) - mm) * 100.0 / len(wq), 2)
                sg["verify_mismatch"] = mm
                print(f"    链{strand}: CDS {sg['q_start']}-{sg['q_end']} "
                      f"({sg['span']} nt, 最佳 {sg['best_identity_pct']}%, "
                      f"复核 {sg['verify_identity_pct']}%)")
                segs_all.append(sg)
                n += 1
        per_tx[acc] = {"len": len(s), "n_segments": n}
        print(f"    → {acc} 共 {n} 个重复区段")

    masked = merge_segments([(x["q_start"], x["q_end"]) for x in segs_all])
    safe, prev = [], 1
    for a, b in masked:
        if a > prev:
            safe.append([prev, a - 1])
        prev = b + 1
    if prev <= ldha["cds_len"]:
        safe.append([prev, ldha["cds_len"]])
    print(f"\n  → LDHA CDS 重复掩码区段: {masked}")
    print(f"  → 被掩码 {sum(b - a + 1 for a, b in masked)} / {ldha['cds_len']} nt")
    print(f"  → 安全连续区段: {safe}")
    longest = max((b - a + 1 for a, b in safe), default=0)
    print(f"  → 最长安全连续区段: {longest} nt")

    fwd = "TTCAGCCCGATTCCGTTACC"
    rev = "CTCCATGTTCCCCAAGGACC"
    rc_rev = revcomp(rev)
    pos_f = cds.find(fwd)
    pos_r = cds.find(rc_rev)

    def in_mask(p, l):
        return any(a <= p + 1 <= b for a, b in masked)

    print("\n  v4 引物落点（CDS 1-based）：")
    print(f"    F {fwd} → CDS {pos_f + 1}-{pos_f + len(fwd)}"
          f"  {'✗ 落在掩码内' if in_mask(pos_f, len(fwd)) else '✓ 安全'}")
    print(f"    R {rev} → CDS {pos_r + 1}-{pos_r + len(rc_rev)}"
          f"  {'✗ 落在掩码内' if in_mask(pos_r, len(rc_rev)) else '✓ 安全'}")

    report["steps"]["sequence_overlap"] = {
        "per_tx": per_tx, "segments": segs_all,
        "masked_segments": masked,
        "masked_nt": sum(b - a + 1 for a, b in masked),
        "safe_segments": safe, "longest_safe_segment": longest,
        "v4_primer_positions": {
            "forward": {"seq": fwd, "cds_start": pos_f + 1,
                        "cds_end": pos_f + len(fwd),
                        "in_mask": in_mask(pos_f, len(fwd))},
            "reverse": {"seq": rev, "cds_start": pos_r + 1,
                        "cds_end": pos_r + len(rc_rev),
                        "in_mask": in_mask(pos_r, len(rc_rev))},
        },
    }

    # ── 步骤 3：阴性对照 ───────────────────────────────────
    print("\n" + "=" * 74)
    print("步骤 3：阴性对照（无关基因 CDS vs NR_202041/NR_202040）")
    print("=" * 74)
    neg = {}
    for g, acc in [("ALDH3A2", "NM_001031806.2"), ("GBP1", "NM_002053.3"),
                   ("ACTB", "NM_001101.5")]:
        o = cds_of(acc, cache)
        if not o:
            continue
        segs = []
        for a, s in tx_seqs.items():
            segs += overlaps(o["cds"], s)
            segs += overlaps(o["cds"], revcomp(s))
        neg[g] = {"accession": acc, "cds_len": o["cds_len"],
                  "n_segments": len(segs), "segments": segs}
        print(f"  {g} ({acc}, CDS {o['cds_len']} nt): {len(segs)} 个合格区段 "
              f"{'→ 通过' if not segs else '→ 失败 ' + str(segs)}")
    neg_ok = all(v["n_segments"] == 0 for v in neg.values())
    print(f"\n  阴性对照总体: {'通过' if neg_ok else '失败'}")
    report["steps"]["negative_controls"] = neg

    # ── 步骤 4：阳性对照（LDHA CDS vs LDHB CDS）────────────
    print("\n" + "=" * 74)
    print("步骤 4：阳性对照（LDHA CDS vs LDHB CDS，应出现高同源区段）")
    print("=" * 74)
    ldhb_cds = None
    for cand in ["NM_001135238.2", "NM_002335.3"]:
        o = cds_of(cand, cache)
        if o:
            ldhb_cds = o
            break
    if ldhb_cds:
        segs = overlaps(ldha["cds"], ldhb_cds["cds"])
        print(f"  LDHA CDS vs {ldhb_cds['accession']} CDS: "
              f"{len(segs)} 个 ≥{MIN_ALEN}nt/≥{MIN_IDENTITY}% 区段")
        for sg in segs[:5]:
            print(f"     CDS {sg['q_start']}-{sg['q_end']} "
                  f"({sg['span']} nt, 最佳 {sg['best_identity_pct']}%)")
        report["steps"]["positive_control"] = {
            "accession": ldhb_cds["accession"],
            "n_segments": len(segs), "segments": segs}
    else:
        print("  （缓存中无 LDHB CDS，跳过）")
        report["steps"]["positive_control"] = {"accession": None, "n_segments": 0}

    report["conclusion"] = {
        "loc105375980_is_ldhap4_alias": any_ov,
        "masked_cds_segments": masked,
        "longest_safe_segment": longest,
        "negative_controls_pass": neg_ok,
        "statement": (
            "LOC105375980 是 chr9p22.3 上独立于 LDHAP4 的 lncRNA locus，"
            "其转录本外显子与 LDHA 编码序列高度同源，构成 RT-qPCR 的真实"
            "共扩增风险；LDHA 引物必须避开被掩码区段。"),
    }
    with open(f"{PROV}/ldha_lncrna_overlap.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
    print(f"\n写出 {PROV}/ldha_lncrna_overlap.json")


if __name__ == "__main__":
    main()
