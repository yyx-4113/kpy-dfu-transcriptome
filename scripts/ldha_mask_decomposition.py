# -*- coding: utf-8 -*-
"""
ldha_mask_decomposition.py — LDHA 掩码的**竞争者分解**与 3'UTR 可行性评估

═══ 目的 ═══

`ldha_unique_region_analysis.py` 结论：LDHA 的 **CDS 共同唯一区段 = 0 nt**，
即常规「CDS 锚定 + 旁系同源区分」在 LDHA 上**数学上不可满足**。
本脚本进一步回答：

1. **谁造成了掩码？** 逐竞争者（LOC105375980 / LDHAL6A / LDHAL6B）给出
   各自覆盖的 LDHA mRNA 区间；区分「策展转录本」（真会进 cDNA 池）与
   「算法预测转录本」（XM_，生物存在性未证实）。
2. **单引物特异性是否可救？** 共扩增要求 **两条引物都** 能结合竞争者；
   故只要**一条**引物落在竞争者不覆盖的区段即可特异性扩增 →
   逐区间计算「覆盖该区段的竞争者数量」。
3. **3'UTR 唯一段是否可用？** GC 含量、重复元件、跨策展转录本共有性。
4. **跨转录本共有**的可用区段（在全部 5 条策展转录本中都存在且唯一）。

输出：
  03_结果/实验验证设计/provenance/ldha_mask_decomposition.json
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
UA = {"User-Agent": "WorkBuddy-LDHA-decomp/1.0"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

MIN_ALEN = 20
MIN_IDENTITY = 95.0
PEN = 19

COMPETITORS = {
    "LOC105375980": ["NR_202040.1", "NR_202041.1", "NR_202034.1",
                     "NR_202033.1", "NR_202035.1", "NR_202036.1",
                     "NR_202037.1", "NR_202038.1", "NR_202039.1",
                     "NR_202042.1", "NR_202043.1", "NR_202044.1",
                     "NR_202045.1"],
    "LDHAL6A": [f"XM_05436786{i}.1" for i in range(3, 9)],
    "LDHAL6B": ["NM_033195.3"],
}
CURATED_COMPETITORS = {"LOC105375980", "LDHAL6B"}   # 有策展(NM_/NR_)转录本


def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def _get(u, timeout=240):
    return urllib.request.urlopen(
        urllib.request.Request(u, headers=UA), timeout=timeout).read().decode(
        errors="ignore")


def fetch_fasta(accs):
    out = {}
    for i in range(0, len(accs), 40):
        chunk = accs[i:i + 40]
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
    return {k: v.upper() for k, v in out.items() if v}


def diag_hits(q, s, min_alen=MIN_ALEN):
    """q 上被 s 以 ≥95% identity、≥min_alen 覆盖的区段（精确 identity 重算）。"""
    qa = np.frombuffer(q.encode(), dtype=np.uint8)
    sa = np.frombuffer(s.encode(), dtype=np.uint8)
    nq, ns = len(qa), len(sa)
    raw = []
    for d in range(-(ns - 1), nq):
        i0 = max(0, -d)
        i1 = min(nq, ns - d)
        if i1 - i0 < min_alen:
            continue
        eq = qa[i0:i1] == sa[i0 + d:i1 + d]
        sc = np.where(eq, 1, -PEN).astype(np.int64)
        S = np.concatenate(([0], np.cumsum(sc)))
        pm = np.minimum.accumulate(S)
        negpm = -pm
        ends = np.arange(min_alen, len(eq) + 1)
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


def merge(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def safe_segments(covered_mask, off=0):
    """covered_mask: True = 被覆盖（坏）。返回未被覆盖（安全）的连续区段。"""
    out, st = [], None
    for i in range(len(covered_mask)):
        if not covered_mask[i] and st is None:
            st = i
        elif covered_mask[i] and st is not None:
            out.append([st + 1 + off, i + off])
            st = None
    if st is not None:
        out.append([st + 1 + off, len(covered_mask) + off])
    return out


def cds_of(acc, cache):
    m = cache.get(acc)
    gb = cache.get(f"__gb__{acc}")
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
            seg = [(int(a), int(c)) for a, c in
                   re.findall(r"(\d+)\.\.(\d+)", loc)]
            break
    if not seg:
        return None
    return {"cds": "".join(full[a - 1:b] for a, b in seg).upper(),
            "seg": seg, "mrna": full.upper(), "len": len(full)}


def gc(s):
    if not s:
        return 0.0
    return round((s.count("G") + s.count("C")) * 100.0 / len(s), 2)


def max_run(s, base):
    """最长同碱基连续长度（检测 poly-X）。"""
    run = mx = 0
    for ch in s:
        if ch == base:
            run += 1
            mx = max(mx, run)
        else:
            run = 0
    return mx


def main():
    cache = json.load(open(f"{PROV}/final_design_cache.json", encoding="utf-8"))
    R = json.load(open(f"{PROV}/ldha_unique_regions.json", encoding="utf-8"))
    ldha_accs = R["steps"]["ldha_curated"]

    all_comp = sorted({a for v in COMPETITORS.values() for a in v})
    print(f"取 {len(all_comp)} 条竞争转录本序列…")
    cseq = fetch_fasta(all_comp)
    print(f"  取得 {len(cseq)} 条\n")

    report = {"date": "2026-10-25",
              "script": "02_脚本/ldha_mask_decomposition.py",
              "criteria": {"min_alen": MIN_ALEN, "min_identity_pct": MIN_IDENTITY},
              "curated_competitors": sorted(CURATED_COMPETITORS),
              "steps": {}}

    per_acc = {}
    for acc in ldha_accs:
        info = cds_of(acc, cache)
        if not info:
            continue
        q = info["mrna"]
        print("=" * 76)
        print(f"{acc}  (mRNA {info['len']} nt, CDS {info['seg']})")
        print("=" * 76)
        n = len(q)
        # 逐竞争者掩码
        comp_mask = {}
        for g, accs in COMPETITORS.items():
            m = np.zeros(n, dtype=bool)
            per_tx = {}
            for a in accs:
                s = cseq.get(a)
                if not s:
                    continue
                ivs = []
                for strand, t in (("+", s), ("-", revcomp(s))):
                    for (i, j, ident) in diag_hits(q, t):
                        ivs.append((i, j, ident))
                if ivs:
                    for (i, j, ident) in ivs:
                        m[i:j] = True
                    per_tx[a] = {"n_seg": len(ivs),
                                 "max_identity_pct": round(
                                     max(x[2] for x in ivs), 2),
                                 "covered_nt": int(sum(j - i for i, j, _ in ivs))}
            comp_mask[g] = m
            cov = int(m.sum())
            print(f"  {g:15s} 覆盖 {cov:5d} nt ({cov*100.0/n:5.1f}%)  "
                  f"策展={'是' if g in CURATED_COMPETITORS else '否(仅XM_预测)'}")
            for a, d in sorted(per_tx.items(),
                               key=lambda x: -x[1]["max_identity_pct"])[:4]:
                print(f"       {a}  n_seg={d['n_seg']:4d} "
                      f"max_id={d['max_identity_pct']:6.2f}%")
        # 策展竞争者并集
        curated_union = np.zeros(n, dtype=bool)
        for g in CURATED_COMPETITORS:
            curated_union |= comp_mask[g]
        all_union = curated_union.copy()
        for g in COMPETITORS:
            all_union |= comp_mask[g]
        # 「仅策展竞争者不覆盖」区段（真可用）
        s_cur = safe_segments(curated_union)
        s_all = safe_segments(all_union)
        print(f"\n  仅策展竞争者(LOC+LDHAL6B)不覆盖的区段: {s_cur}")
        print(f"  全部竞争者(含 XM_ 预测)都不覆盖的区段: {s_all}")
        # 单引物特异性：区段内「策展竞争者覆盖数」
        for label, segs in (("策展口径", s_cur), ("全口径", s_all)):
            print(f"\n  ── {label} 各区段 GC / 重复检查 ──")
            for (a, b) in segs:
                sub = q[a - 1:b]
                if len(sub) < 20:
                    continue
                gcp = gc(sub)
                px = max(max_run(sub, x) for x in "ACGT")
                # 找 GC 40-60 的窗口
                win = [(i, j) for i in range(len(sub) - 19)
                       for j in [i + 20]
                       if 40 <= gc(sub[i:j]) <= 60]
                print(f"     mRNA {a}-{b} ({b-a+1} nt)  GC={gcp}%  "
                      f"maxPolyX={px}  20-mer GC合格窗口={len(win)}")
        per_acc[acc] = {
            "len": n, "cds_seg": info["seg"],
            "competitor_coverage_nt": {
                g: int(comp_mask[g].sum()) for g in COMPETITORS},
            "safe_curated_only": s_cur,
            "safe_all": s_all,
            "n_curated_competitors_covering": int(curated_union.sum()),
            "n_all_competitors_covering": int(all_union.sum()),
        }
        np.save(f"{PROV}/_ldha_curated_mask_{acc}.npy", curated_union)
        np.save(f"{PROV}/_ldha_allmask_{acc}.npy", all_union)

    report["steps"]["per_ldha_transcript"] = per_acc

    # ── 3'UTR 跨策展转录本共有性 ────────────────────────────
    print("\n" + "=" * 76)
    print("3'UTR 唯一段的跨策展转录本共有性（以 GC 富集段为探针）")
    print("=" * 76)
    probe_acc = "NM_001165414.2"
    info = cds_of(probe_acc, cache)
    q = info["mrna"]
    cmask = np.load(f"{PROV}/_ldha_curated_mask_{probe_acc}.npy")
    cds_end = info["seg"][-1][1]
    s3 = safe_segments(cmask[cds_end:], cds_end)
    print(f"  {probe_acc} 3'UTR（策展口径）唯一段: {s3}")
    common = None
    for (a, b) in s3:
        sub = q[a - 1:b]
        # 在其他策展转录本中的存在性
        present = {}
        for other in ldha_accs:
            o = cds_of(other, cache)
            if not o:
                continue
            s = o["mrna"]
            hit = s.find(sub)
            present[other] = hit + 1 if hit >= 0 else None
        print(f"\n  探针 mRNA {a}-{b} ({len(sub)} nt), GC={gc(sub)}%")
        print(f"    序列: {sub[:70]}")
        for k, v in present.items():
            print(f"    {k}: {'位置 '+str(v) if v else '不含'}")
        if all(v is not None for v in present.values()):
            common = [a, b] if common is None else [
                max(common[0], a), min(common[1], b)]
    print(f"\n  → 全部策展转录本共有且策展唯一的 3'UTR 区间: {common}")
    if common:
        sub = q[common[0] - 1:common[1]]
        print(f"     长度 {len(sub)} nt, GC {gc(sub)}%, "
              f"maxPolyX {max(max_run(sub, x) for x in 'ACGT')}")
    report["steps"]["utr3_common"] = {
        "segments_3utr": s3, "common_across_curated": common,
        "gc_pct": gc(q[common[0] - 1:common[1]]) if common else None,
    }

    with open(f"{PROV}/ldha_mask_decomposition.json", "w",
              encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
    print(f"\n写出 {PROV}/ldha_mask_decomposition.json")


if __name__ == "__main__":
    main()
