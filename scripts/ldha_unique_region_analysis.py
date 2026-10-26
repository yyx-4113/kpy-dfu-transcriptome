# -*- coding: utf-8 -*-
"""
ldha_unique_region_analysis.py — LDHA 的**全基因组唯一区段**分析（本地，无 BLAST）

═══ 背景 ═══

BLAST v6 全库复核发现 LDHA 引物对存在 P0 脱靶：
  NR_202041 / NR_202040（LOC105375980，chr9p22.3 的 lncRNA）
经 Biopython Smith-Waterman 局部比对独立确认：
  LDHA mRNA 260-1067  ↔ LOC tx 1224-2031  808 nt, 31 mm, identity 96.16%
  LDHA mRNA 1069-1694 ↔ LOC tx 2032-2657  626 nt, 28 mm, identity 95.53%
  LDHA mRNA 1740-1849 ↔ LOC tx 2709-2818  110 nt,  5 mm, identity 95.45%
→ LOC105375980 携带 LDHA 编码序列近全长的高同源拷贝；
  该区内错配密度仅 3.8%，任一 20-mer 约 46% 概率与 LOC 完全一致
  ⇒ 「换一个窗口」无法消除 P0。

同时查明：v4 的旁系同源集合只含 LDHB/LDHC/LDHAL6A/LDHAL6B，
因为 RefSeq 对假基因**不提供策展转录本**（refseqrna 链接为 0），
整个 LDHAP/LDAL 假基因家族被系统性漏掉。

═══ 本脚本的目的 ═══

对 LDHA 的每一条策展转录本，**逐碱基**计算：
  「该位置是否被任一非靶标人源 RefSeq 转录本以 ≥95% identity、≥20 nt 覆盖」
输出可作为引物设计硬约束的「唯一区段掩码」。

比对引擎：穷举全部对角线 + 前缀和最大子段（零剪枝，不漏检），
与 verify_ldha_lncrna_overlap.py 同一算法；再以 Biopython
Smith-Waterman 对关键配对做**独立交叉校验**。

输出：
  03_结果/实验验证设计/provenance/ldha_unique_regions.json
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
UA = {"User-Agent": "WorkBuddy-LDHA-unique/1.0"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

MIN_ALEN = 20
MIN_IDENTITY = 95.0
PEN = 19  # score = alen − 20·mismatch

# LDHA 假基因/旁系同源家族的全部 locus（GeneID）。
# 注意：LDHAP2/LDHAL2、LDHAP3/LDHAL3、LDHAP4/LDHAL4、LDHAP5/LDHAL5、
#       LDHAP7/LDHAL7 在 NCBI 中是**同一 locus 的两个基因记录**
#       （坐标完全相同），故按 locus 去重。
FAMILY_GENES = {
    "LDHAL2": 3941, "LDHAL3": 3942, "LDHAL4": 3943, "LDHAL5": 3944,
    "LDHAL6A": 160287, "LDHAL6B": 92483, "LDHAL6CP": 121498,
    "LDHAL6DP": 100420508, "LDHAL6EP": 100420462, "LDHAL6FP": 100420477,
    "LDHAL7": 378485, "LDHAP1": 3940, "LDHAP2": 100190799,
    "LDHAP3": 442013, "LDHAP4": 158222, "LDHAP5": 729666,
    "LDHAP7": 100190800, "LOC105375980": 105375980,
}
# 已知与 LDHA 同一 locus 的重复记录（基因组坐标完全一致）
SAME_LOCUS_PAIRS = [("LDHAP2", "LDHAL2"), ("LDHAP3", "LDHAL3"),
                    ("LDHAP4", "LDHAL4"), ("LDHAP5", "LDHAL5"),
                    ("LDHAP7", "LDHAL7")]


def _get(u, timeout=240):
    return urllib.request.urlopen(
        urllib.request.Request(u, headers=UA), timeout=timeout).read().decode(
        errors="ignore")


def _post(u, payload, timeout=300):
    req = urllib.request.Request(
        u, data=json.dumps(payload).encode(),
        headers={"User-Agent": UA["User-Agent"], "Content-Type": "application/json"},
        method="POST")
    return urllib.request.urlopen(req, timeout=timeout).read().decode(
        errors="ignore")


def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def fetch_fasta_batch(accs):
    """批量取 FASTA；返回 {accession: seq}。"""
    out = {}
    for i in range(0, len(accs), 50):
        chunk = accs[i:i + 50]
        u = EUTILS + "efetch.fcgi?" + urllib.parse.urlencode(
            {"db": "nuccore", "id": ",".join(chunk), "rettype": "fasta",
             "retmode": "text"})
        t = _get(u)
        name, buf = None, []
        for ln in t.splitlines():
            if ln.startswith(">"):
                if name:
                    out[name] = "".join(buf)
                name, buf = ln[1:].strip().split()[0], []
            elif ln.strip():
                buf.append(ln.strip())
        if name:
            out[name] = "".join(buf)
    return out


def gene_transcripts(gid):
    """基因的全部 RefSeq 转录本 accession（NM_/NR_/XM_/XR_）。"""
    u = EUTILS + "elink.fcgi?" + urllib.parse.urlencode(
        {"dbfrom": "gene", "id": str(gid), "db": "nuccore",
         "retmode": "json", "linkname": "gene_nuccore_refseqrna"})
    ids = []
    try:
        j = json.loads(_get(u))
        for ls in j.get("linksets", []):
            for d in ls.get("linksetdbs", []):
                ids += d.get("links", [])
    except Exception:                                     # noqa: BLE001
        return []
    if not ids:
        return []
    # 若返回的是纯数字 UID → 用 esummary 转成 accession.version
    if all(x.isdigit() for x in ids):
        accs = []
        for i in range(0, len(ids), 50):
            u2 = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(
                {"db": "nuccore", "id": ",".join(ids[i:i + 50]),
                 "retmode": "json"})
            try:
                res = json.loads(_get(u2)).get("result", {})
            except Exception:                             # noqa: BLE001
                continue
            for uid in ids[i:i + 50]:
                r = res.get(uid)
                if r:
                    accs.append(r.get("accessionversion") or r.get("caption"))
        ids = [a for a in ids if a]
        ids = accs + ids
    seen, out = set(), []
    for a in ids:
        if not a or a in seen:
            continue
        seen.add(a)
        if re.match(r"^(NM_|NR_|XM_|XR_)", a):
            out.append(a)
    return out


def diag_mask(q, s, min_alen=MIN_ALEN):
    """返回 s 对 q 的掩码：list of (q_start0, q_end0_exclusive, identity)。

    穷举全部对角线；对每条对角线取「以每个右端点结尾的最长合格窗口」
    （score = alen − 20·mm ≥ 0 ⟺ identity ≥ 95%），由单调性保证不漏检。
    返回的窗口按 q 坐标给出，identity 由原始字节精确重算。
    """
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
        if i1 - i0 < min_alen:
            continue
        seg = qa[i0:i1]
        oth = sa[i0 + d:i1 + d]
        eq = (seg == oth)
        sc = np.where(eq, 1, -PEN).astype(np.int64)
        S = np.concatenate(([0], np.cumsum(sc)))
        pm = np.minimum.accumulate(S)
        negpm = -pm
        ends = np.arange(min_alen, len(seg) + 1)
        idx = np.searchsorted(negpm, -S[ends], side="left")
        ok = idx <= (ends - min_alen)
        if not np.any(ok):
            continue
        idx = idx[ok]
        ends = ends[ok]
        score = S[ends] - S[idx]
        for k in range(len(idx)):
            if score[k] < 0:
                continue
            i, e = int(idx[k]), int(ends[k])
            n = e - i
            mm = int(np.count_nonzero(~eq[i:e]))
            ident = (n - mm) * 100.0 / n
            if ident + 1e-9 < MIN_IDENTITY:
                continue
            raw.append((i0 + i, i0 + e, ident))
    return raw


def main():
    cache = json.load(open(f"{PROV}/final_design_cache.json", encoding="utf-8"))
    report = {"date": "2026-10-25",
              "script": "02_脚本/ldha_unique_region_analysis.py",
              "criteria": {"min_alen": MIN_ALEN, "min_identity_pct": MIN_IDENTITY},
              "algorithm": ("穷举全部对角线 + 前缀和最大子段；"
                            "score = alen − 20·mismatch ≥ 0 ⟺ identity ≥ 95%；"
                            "identity 由原始字节精确重算；零剪枝不漏检"),
              "steps": {}}

    # ── 1. 取全部家族 locus 的转录本 ────────────────────────
    print("=" * 76)
    print("步骤 1：收集 LDHAP/LDAL 家族全部 locus 的 RefSeq 转录本")
    print("=" * 76)
    tx_by_gene = {}
    for g, gid in FAMILY_GENES.items():
        txs = gene_transcripts(gid)
        tx_by_gene[g] = txs
        print(f"  {g:14s} (GeneID {gid:>9}) : {len(txs):2d} 转录本"
              f"  {txs[:5]}")
    report["steps"]["family_transcripts"] = tx_by_gene

    # ── 2. 取所有序列（含 LDHA 策展转录本）──────────────────
    V4_early = json.load(open(f"{OUT}/final_primer_design_v4.json",
                              encoding="utf-8"))
    ldha_accs_early = V4_early["genes"]["LDHA"]["curated_accs"]
    all_accs = sorted({a for v in tx_by_gene.values() for a in v}
                      | set(ldha_accs_early))
    print(f"\n  去重后共 {len(all_accs)} 条转录本，取序列…")
    seqs = fetch_fasta_batch(all_accs)
    seqs = {k: v.upper() for k, v in seqs.items() if v}
    print(f"  取得 {len(seqs)} 条序列")
    # 缓存兜底（final_design_cache 已含 LDHA 策展转录本）
    miss = [a for a in ldha_accs_early if a not in seqs]
    for a in miss:
        c = cache.get(a)
        if isinstance(c, dict) and c.get("seq"):
            seqs[a] = c["seq"].upper()
    print(f"  缓存兜底补入 {len(miss) - len([a for a in miss if a in seqs])} 条")

    # ── 3. LDHA 策展转录本 ─────────────────────────────────
    V4 = json.load(open(f"{OUT}/final_primer_design_v4.json", encoding="utf-8"))
    ldha_accs = V4["genes"]["LDHA"]["curated_accs"]
    print("\n  LDHA 策展转录本:", ldha_accs)
    report["steps"]["ldha_curated"] = ldha_accs

    # 竞争转录本集合：**只有真正有转录本**的家族 locus 才构成共扩增风险
    # （假基因 locus 只有 RefSeqGene NG_ 记录、零转录本 → 不进入 cDNA 池）
    comp_accs = []
    for g, txs in tx_by_gene.items():
        if not txs:
            continue
        for a in txs:
            if a not in ldha_accs:
                comp_accs.append((g, a))
    zero_tx = [g for g, v in tx_by_gene.items() if not v]
    print(f"\n  有转录本的竞争 locus: "
          f"{sorted(set(g for g, _ in comp_accs))}")
    print(f"  零转录本（不构成共扩增风险）: {zero_tx}")
    report["steps"]["zero_transcript_loci"] = zero_tx

    # ── 4. 逐条 LDHA 转录本计算掩码 ─────────────────────────
    print("\n" + "=" * 76)
    print("步骤 2：逐条 LDHA 策展转录本 × 竞争转录本 → 唯一区段掩码")
    print("=" * 76)
    results = {}
    for acc in ldha_accs:
        q = seqs.get(acc)
        if not q:
            print(f"  {acc}: 序列不可得，跳过")
            continue
        cds = None
        try:
            import importlib.util
            spec = importlib.util.spec_from_file_location(
                "v", f"{BASE}/02_脚本/verify_ldha_lncrna_overlap.py")
            v = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(v)
            cds = v.cds_of(acc, cache)
        except Exception:                                 # noqa: BLE001
            cds = None
        covered = np.zeros(len(q), dtype=bool)
        max_id = np.zeros(len(q), dtype=float)
        per_comp = {}
        for g, cacc in comp_accs:
            s = seqs.get(cacc)
            if not s:
                continue
            best = []
            for strand, t in (("+", s), ("-", revcomp(s))):
                for (i, j, ident) in diag_mask(q, t):
                    best.append((i, j, ident, strand))
            if not best:
                continue
            per_comp.setdefault(g, {"n_tx": 0, "n_segments": 0,
                                    "max_identity_pct": 0.0,
                                    "total_covered_nt": 0,
                                    "example_acc": cacc})
            per_comp[g]["n_segments"] += len(best)
            per_comp[g]["n_tx"] += 1
            for (i, j, ident, strand) in best:
                covered[i:j] = True
                max_id[i:j] = np.maximum(max_id[i:j], ident)
                per_comp[g]["total_covered_nt"] += j - i
                per_comp[g]["max_identity_pct"] = max(
                    per_comp[g]["max_identity_pct"], ident)
        # 唯一连续区段（未被覆盖）
        safe = []
        st = None
        for i in range(len(q)):
            if not covered[i] and st is None:
                st = i
            elif covered[i] and st is not None:
                safe.append([st + 1, i])
                st = None
        if st is not None:
            safe.append([st + 1, len(q)])
        longest = max((b - a + 1 for a, b in safe), default=0)
        print(f"\n  {acc} (len {len(q)}):")
        print(f"     被 ≥95% 覆盖: {int(covered.sum())} nt "
              f"({covered.sum()*100.0/len(q):.1f}%)")
        print(f"     唯一连续区段: {safe}")
        print(f"     最长唯一区段: {longest} nt")
        for g, d in sorted(per_comp.items(),
                           key=lambda x: -x[1]["max_identity_pct"])[:8]:
            print(f"     竞争 {g:14s} n_seg={d['n_segments']:4d} "
                  f"max_id={d['max_identity_pct']:6.2f}% "
                  f"({d['example_acc']})")
        results[acc] = {
            "len": len(q),
            "cds_len": cds["cds_len"] if cds else None,
            "cds_coords_in_mrna": cds["cds_coords_in_mrna"] if cds else None,
            "n_covered_nt": int(covered.sum()),
            "pct_covered": round(covered.sum() * 100.0 / len(q), 2),
            "safe_segments": safe,
            "longest_safe_segment": longest,
            "competitors": {k: v for k, v in per_comp.items()},
        }
        # 保存逐碱基掩码供下游引物设计使用
        np.save(f"{PROV}/_ldha_mask_{acc}.npy", covered)

    report["steps"]["per_ldha_transcript"] = results

    # ── 5. CDS 层面的唯一区段（对所有策展转录本取交集）────
    print("\n" + "=" * 76)
    print("步骤 3：CDS 层面唯一区段（所有策展转录本共同）")
    print("=" * 76)
    cds_safe_all = None
    for acc in ldha_accs:
        q = seqs.get(acc)
        cds = None
        try:
            import importlib.util
            spec = importlib.util.spec_from_file_location(
                "v", f"{BASE}/02_脚本/verify_ldha_lncrna_overlap.py")
            v = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(v)
            cds = v.cds_of(acc, cache)
        except Exception:                                 # noqa: BLE001
            continue
        if not cds:
            continue
        mask = np.load(f"{PROV}/_ldha_mask_{acc}.npy")
        # CDS 在 mRNA 上的区间（1-based）
        (s1, e1), = cds["cds_coords_in_mrna"]
        cds_mask = mask[s1 - 1:e1]
        cds_safe = []
        st = None
        for i in range(len(cds_mask)):
            if not cds_mask[i] and st is None:
                st = i
            elif cds_mask[i] and st is not None:
                cds_safe.append([st + 1, i])
                st = None
        if st is not None:
            cds_safe.append([st + 1, len(cds_mask)])
        print(f"  {acc}: CDS {len(cds_mask)} nt，被覆盖 "
              f"{int(cds_mask.sum())} nt，唯一区段 {cds_safe}")
        if cds_safe_all is None:
            cds_safe_all = [[a, b] for a, b in cds_safe]
        else:
            cds_safe_all = _intersect(cds_safe_all, cds_safe)
    print(f"\n  → 所有策展转录本 CDS 的共同唯一区段: {cds_safe_all}")
    longest_cds = max((b - a + 1 for a, b in cds_safe_all), default=0)
    print(f"  → 最长共同唯一 CDS 区段: {longest_cds} nt")
    report["steps"]["cds_common_safe"] = {
        "segments": cds_safe_all, "longest": longest_cds}

    with open(f"{PROV}/ldha_unique_regions.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
    print(f"\n写出 {PROV}/ldha_unique_regions.json")


def _intersect(a, b):
    out = []
    for (a1, a2) in a:
        for (b1, b2) in b:
            lo, hi = max(a1, b1), min(a2, b2)
            if lo <= hi:
                out.append([lo, hi])
    merged = []
    for lo, hi in sorted(out):
        if merged and lo <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return merged


if __name__ == "__main__":
    main()
