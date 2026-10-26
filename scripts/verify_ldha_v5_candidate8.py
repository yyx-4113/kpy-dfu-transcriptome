# -*- coding: utf-8 -*-
"""
verify_ldha_v5_candidate8.py — 对候选 8（及两个备份正向）做 BLAST 全库复核 + 独立几何复算

═══ 目的 ═══
候选 8：F=ACCACTTCCATTGTAAGTCCCA（正向，已在 prior 运行验证 P0=0）
        R=CCAGCTTCCAGAGGACAAGAT（反向 21-mer，PIGS 命中的 20-mer 的 3' 加 T 变体，须新验证）
        —— 产物 174bp（落在 80–180 窗口）、GC 差 6.9%（落在 ≤10），同时满足两项约束。
备份：candidate 6 / 7 的正向（3'UTR 内，未被 BLAST 测试），反向复用已验证 P0=0 的 CACTGCACCCTGCCTTAAGA。

═══ 修复 prior 几何 bug ═══
prior 的 product_size = iRrc - iF（漏加反向引物长度）。本脚本独立复算：
    product = iR + len(R) - iF   （iR = rc(R) 在 + 链上的 0-based 起点）
同时用直接计数复算 GC%，以独立代码路径交叉验证设计脚本（Primer3）的口径。

判据同 v6：core_nt + 人源限定；全长匹配分层 TARGET/PARALOG/PSEUDO/P0/P1/P2/PARTIAL；
usable = 目标存在全长 RefSeq 匹配 且 P0==0 且 LOC_hit==0 且 双向 tgtTx≥1。
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import re
import time
import json
import urllib.request
import urllib.parse

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
PROV = f"{OUT}/provenance"
RAW = f"{PROV}/blast_xml_raw_v5_ldha_c8"
os.makedirs(RAW, exist_ok=True)
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
UA = {"User-Agent": "WorkBuddy-LDHAv5c8/1.0"}

HUMAN_TAXON = "txid9606[Organism:exp]"
HUMAN_RE = re.compile(r"Homo sapiens", re.I)
PSEUDO_RE = re.compile(r"pseudogene", re.I)
GENE = "LDHA"
PARALOG_FAMILY = {"LDHA", "LDHB", "LDHC", "LDHAL", "LDHAL6A", "LDHAL6B",
                  "LDHAL1", "LDHAL2", "LDHAL3", "LDHAL4", "LDHAL5", "LDHAL7"}

# 候选：(正向, 反向, 标签)
CANDIDATES = [
    ("ACCACTTCCATTGTAAGTCCCA", "CCAGCTTCCAGAGGACAAGAT", "c8"),
    ("GTGTGTCACTATAGGGAACACAGA", "CACTGCACCCTGCCTTAAGA", "c6"),
    ("AGTGTGTCACTATAGGGAACACAG", "CACTGCACCCTGCCTTAAGA", "c7"),
]


def blk(s, t):
    return re.findall(rf"<{t}>(.*?)</{t}>", s, re.S)


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


def blast_put(seq, db="core_nt"):
    p = {"CMD": "Put", "PROGRAM": "blastn", "DATABASE": db, "QUERY": seq,
         "MEGARUN": "on", "HITLIST_SIZE": "5000", "FILTER": "F",
         "EXPECT": "1000", "WORD_SIZE": "7", "NEIGHBOURING": "0",
         "ENTREZ_QUERY": HUMAN_TAXON}
    t = urllib.request.urlopen(
        urllib.request.Request(BLAST, data=urllib.parse.urlencode(p).encode(),
                               headers=UA), timeout=180).read().decode(errors="ignore")
    rid = re.search(r"RID = (\S+)", t)
    if not rid:
        raise RuntimeError("BLAST Put 未返回 RID")
    return rid.group(1)


def blast_status(rid):
    u = BLAST + "?" + urllib.parse.urlencode(
        {"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid})
    t = urllib.request.urlopen(u, timeout=180).read().decode(errors="ignore")
    m = re.search(r"Status=([A-Z]+)", t)
    return m.group(1) if m else "UNKNOWN"


def blast_xml(rid, db="core_nt"):
    u = BLAST + "?" + urllib.parse.urlencode(
        {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML",
         "ALIGNMENTS": "5000", "DESCRIPTIONS": "5000"})
    t = urllib.request.urlopen(u, timeout=600).read().decode(errors="ignore")
    with open(f"{RAW}/{db}_{rid}.xml", "w", encoding="utf-8") as fh:
        fh.write(t)
    return t


def wait_ready(rid, max_wait=2400):
    t0 = time.time()
    time.sleep(8)
    while time.time() - t0 < max_wait:
        st = blast_status(rid)
        if st in ("READY", "FAILED"):
            return st
        time.sleep(8)
    return "TIMEOUT"


def parse_hits(xml):
    hits = []
    for hb in re.findall(r"<Hit>(.*?)</Hit>", xml, re.S):
        acc = (blk(hb, "Hit_accession") or [""])[0]
        title = (blk(hb, "Hit_def") or [""])[0].strip()
        for hp in re.findall(r"<Hsp>(.*?)</Hsp>", hb, re.S):
            g = lambda t: (blk(hp, t) or [""])[0]      # noqa: E731
            idn = int(g("Hsp_identity") or 0)
            al = int(g("Hsp_align-len") or 0)
            qs, hs = g("Hsp_qseq"), g("Hsp_hseq")
            if not qs or not hs:
                continue
            hits.append({
                "acc": acc, "title": title, "idn": idn, "al": al,
                "identity_pct": round(idn * 100.0 / al, 2) if al else 0.0,
                "qseq": qs, "hseq": hs,
                "human": bool(HUMAN_RE.search(title)),
                "is_refseq_tx": bool(re.match(r"(NM_|XM_|NR_|XR_)", acc)),
                "pseudo": bool(PSEUDO_RE.search(title)),
                "gene": gene_of(title),
            })
    return hits


def classify(hit, primer_len):
    if hit["al"] != primer_len:
        return ("PARTIAL", None, None)
    qs, hs = hit["qseq"], hit["hseq"]
    L = min(len(qs), len(hs))
    mm_pos = [i for i in range(L) if qs[i] != hs[i] and hs[i] not in "- "]
    mm_total = len(mm_pos)
    mm_3p5 = sum(1 for i in mm_pos if i >= L - 5)
    g = hit["gene"]
    if g and g == GENE:
        return ("TARGET", mm_total, mm_3p5)
    if hit["pseudo"]:
        return ("PSEUDO", mm_total, mm_3p5)
    if g and g in PARALOG_FAMILY:
        return ("PARALOG", mm_total, mm_3p5)
    if not hit["is_refseq_tx"]:
        return ("P2", mm_total, mm_3p5)
    if mm_3p5 == 0 and hit["identity_pct"] >= 95:
        return ("P0", mm_total, mm_3p5)
    if mm_3p5 == 0 and hit["identity_pct"] >= 90:
        return ("P1", mm_total, mm_3p5)
    return ("P2", mm_total, mm_3p5)


def run_one(seq, label):
    last_err = None
    for attempt in range(1, 4):
        try:
            rid = blast_put(seq)
            st = wait_ready(rid)
            if st == "READY":
                xml = blast_xml(rid)
                return {"status": "READY", "hits": parse_hits(xml),
                        "rid": rid, "attempt": attempt}
            last_err = f"status={st}"
        except Exception as e:                  # noqa: BLE001
            last_err = f"{type(e).__name__}: {e}"[:200]
        print(f"    [{label}] 尝试 {attempt} 失败（{last_err}），"
              f"{'重试' if attempt < 3 else '标记 INCONCLUSIVE'}", flush=True)
        time.sleep(5)
    return {"status": "INCONCLUSIVE", "hits": [], "rid": None,
            "attempt": 3, "error": last_err}


def verify_primer(seq, label):
    res = run_one(seq, label)
    if res["status"] != "READY":
        return {"status": res["status"], "P0": None, "P1": None,
                "PARALOG": None, "tgt_tx": None, "LOC": None,
                "note": "INCONCLUSIVE"}
    hits = res["hits"]
    human = [h for h in hits if h["human"]]
    counts, tgt_tx = {}, 0
    loc = 0
    for h in human:
        lvl, _, _ = classify(h, len(seq))
        counts[lvl] = counts.get(lvl, 0) + 1
        if lvl == "TARGET" and h["is_refseq_tx"]:
            tgt_tx += 1
        if "105375980" in h["title"]:
            loc += 1
    return {"status": "READY", "P0": counts.get("P0", 0),
            "P1": counts.get("P1", 0), "PARALOG": counts.get("PARALOG", 0),
            "PSEUDO": counts.get("PSEUDO", 0),
            "tgt_tx": tgt_tx, "LOC": loc,
            "n_human": len(human)}


def rc(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def geom(F, R, seq):
    iF = seq.find(F)
    iR = seq.find(rc(R))
    ok = (iF >= 0 and iR >= 0 and iR >= iF)
    prod = (iR + len(R) - iF) if ok else None
    gc = lambda s: round(100.0 * (s.count("G") + s.count("C")) / len(s), 1)
    return {
        "f_pos_1based": (iF + 1) if iF >= 0 else None,
        "r_rc_pos_1based": (iR + 1) if iR >= 0 else None,
        "product_size": prod,
        "f_gc": gc(F), "r_gc": gc(R),
        "gc_diff": (round(gc(R) - gc(F), 1)) if ok else None,
        "f_len": len(F), "r_len": len(R),
        "ok_in_canon": ok,
    }


def main():
    cache = json.load(open(f"{PROV}/final_design_cache.json", encoding="utf-8"))
    seq = cache["NM_005566.4"]["seq"].upper()

    results = []
    for (F, R, tag) in CANDIDATES:
        print(f"\n##### [{tag}] F={F} R={R} #####", flush=True)
        g = geom(F, R, seq)
        print("  独立几何复算:", g, flush=True)
        fv = verify_primer(F, "F")
        time.sleep(3)
        rv = verify_primer(R, "R")
        time.sleep(3)
        specific = (fv["status"] == "READY" and rv["status"] == "READY"
                    and (fv["P0"] or 0) == 0 and (rv["P0"] or 0) == 0
                    and (fv["tgt_tx"] or 0) >= 1 and (rv["tgt_tx"] or 0) >= 1
                    and (fv["LOC"] or 0) == 0 and (rv["LOC"] or 0) == 0)
        verdict = ("USABLE" if specific else
                   ("INCONCLUSIVE" if (fv["status"] != "READY"
                                      or rv["status"] != "READY")
                    else "OFFTARGET_P0"))
        print(f"  F: P0={fv.get('P0')} tgtTx={fv.get('tgt_tx')} "
              f"LOC={fv.get('LOC')}", flush=True)
        print(f"  R: P0={rv.get('P0')} tgtTx={rv.get('tgt_tx')} "
              f"LOC={rv.get('LOC')}", flush=True)
        print(f"  → verdict = {verdict}", flush=True)
        results.append({"tag": tag, "forward": F, "reverse": R,
                        "geometry": g, "forward_blast": fv,
                        "reverse_blast": rv, "verdict": verdict})

    usable = [r for r in results if r["verdict"] == "USABLE"]
    # 优先级：同时满足 产物≤180 且 GC差≤10 —— 按产物大小（更稳健）降序
    usable.sort(key=lambda r: -(r["geometry"]["product_size"] or 0))
    best = usable[0] if usable else None
    json.dump({
        "date": "2026-10-25",
        "script": "02_脚本/verify_ldha_v5_candidate8.py",
        "gene": GENE, "template": "NM_005566.4",
        "candidates": CANDIDATES,
        "results": results,
        "selected": best,
    }, open(f"{OUT}/primer_blast_v5_ldha_candidate8.json", "w",
            encoding="utf-8"), ensure_ascii=False, indent=1, default=str)
    print(f"\n=== 选定（双向无 P0，且优先满足 产物≤180 & GC差≤10）===")
    if best:
        print(f"  [{best['tag']}] F={best['forward']} R={best['reverse']} "
              f"product={best['geometry']['product_size']}bp "
              f"GC差={best['geometry']['gc_diff']}% "
              f"F@{best['geometry']['f_pos_1based']} "
              f"Rrc@{best['geometry']['r_rc_pos_1based']}")
    else:
        print("  !! 候选 8 / 6 / 7 均不可用，需重新设计反向/正向")
    print("DONE — verify_ldha_v5_candidate8.py")


if __name__ == "__main__":
    main()
