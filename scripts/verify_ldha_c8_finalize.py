# -*- coding: utf-8 -*-
"""
verify_ldha_c8_finalize.py — 仅验证候选 8 的新反向引物，产出 finalize 所需的 JSON

═══ 设计 ═══
- 正向 ACCACTTCCATTGTAAGTCCCA 已在 prior 工件 primer_blast_v5_ldha.json
  （verify_ldha_v5_specificity.py）中实测 BLAST：P0=0 / tgt_tx=6 / LOC=0 /
  n_human=3667 —— 序列未变，直接复用该权威结果，避免冗余 BLAST。
- 反向 CCAGCTTCCAGAGGACAAGAT（21-mer，3' 端加 T 逃逸 PIGS/NM_033198 的
  3' 端 5 nt 全匹配脱靶）是**唯一未验证的新片段**，本脚本仅对它做 BLAST。
- 几何独立复算（iR + len(R) - iF，修复 prior 漏加 len(R) 的 bug）。

判据同 v6：core_nt + txid9606；全长匹配分层 TARGET/PARALOG/PSEUDO/P0/P1/P2/PARTIAL；
USABLE = 正向 P0==0 且 tgt_tx>=1 且 LOC==0，且反向 P0==0 且 tgt_tx>=1 且 LOC==0。

输出 primer_blast_v5_ldha_candidate8.json（finalize_ldha_v5.py 读取的 schema）。
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
RAW = f"{PROV}/blast_xml_raw_c8rev"
os.makedirs(RAW, exist_ok=True)
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
UA = {"User-Agent": "WorkBuddy-LDHAc8rev/1.0"}

HUMAN_TAXON = "txid9606[Organism:exp]"
HUMAN_RE = re.compile(r"Homo sapiens", re.I)
PSEUDO_RE = re.compile(r"pseudogene", re.I)
GENE = "LDHA"
PARALOG_FAMILY = {"LDHA", "LDHB", "LDHC", "LDHAL", "LDHAL6A", "LDHAL6B",
                  "LDHAL1", "LDHAL2", "LDHAL3", "LDHAL4", "LDHAL5", "LDHAL7"}

# 候选 8（最终选定）
F = "ACCACTTCCATTGTAAGTCCCA"
R = "CCAGCTTCCAGAGGACAAGAT"


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


def wait_ready(rid, max_wait=1200):
    t0 = time.time()
    time.sleep(10)
    while time.time() - t0 < max_wait:
        st = blast_status(rid)
        if st in ("READY", "FAILED"):
            return st
        time.sleep(10)
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


def verify_reverse(seq, label):
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


def geom(Fp, Rp, seq):
    iF = seq.find(Fp)
    iR = seq.find(rc(Rp))
    ok = (iF >= 0 and iR >= 0 and iR >= iF)
    prod = (iR + len(Rp) - iF) if ok else None
    gc = lambda s: round(100.0 * (s.count("G") + s.count("C")) / len(s), 1)
    return {
        "f_pos_1based": (iF + 1) if iF >= 0 else None,
        "r_rc_pos_1based": (iR + 1) if iR >= 0 else None,
        "product_size": prod,
        "f_gc": gc(Fp), "r_gc": gc(Rp),
        "gc_diff": (round(gc(Rp) - gc(Fp), 1)) if ok else None,
        "f_len": len(Fp), "r_len": len(Rp),
        "ok_in_canon": ok,
    }


def main():
    # 复用 prior 已验证正向
    prior = json.load(open(f"{OUT}/primer_blast_v5_ldha.json", encoding="utf-8"))
    prior_fwd = prior["results"][0]["forward_blast"]
    print(f"复用 prior 正向 BLAST：P0={prior_fwd['P0']} tgtTx={prior_fwd['tgt_tx']} "
          f"LOC={prior_fwd['LOC']} n_human={prior_fwd['n_human']}", flush=True)

    cache = json.load(open(f"{PROV}/final_design_cache.json", encoding="utf-8"))
    seq = cache["NM_005566.4"]["seq"].upper()

    g = geom(F, R, seq)
    print("独立几何复算:", g, flush=True)

    print(f"\n##### 新反向 R={R} #####", flush=True)
    rv = verify_reverse(R, "R")
    print(f"  R: status={rv['status']} P0={rv.get('P0')} tgtTx={rv.get('tgt_tx')} "
          f"LOC={rv.get('LOC')}", flush=True)

    f_ok = (prior_fwd["status"] == "READY" and prior_fwd["P0"] == 0
            and prior_fwd["tgt_tx"] >= 1 and prior_fwd["LOC"] == 0)
    r_ok = (rv["status"] == "READY" and (rv["P0"] or 0) == 0
            and (rv["tgt_tx"] or 0) >= 1 and (rv["LOC"] or 0) == 0)
    if f_ok and r_ok:
        verdict = "USABLE"
    elif rv["status"] != "READY":
        verdict = "INCONCLUSIVE"
    else:
        verdict = "OFFTARGET_P0"
    print(f"  → verdict = {verdict}", flush=True)

    entry = {"tag": "c8", "forward": F, "reverse": R, "geometry": g,
             "forward_blast": prior_fwd, "reverse_blast": rv, "verdict": verdict}
    c8 = {
        "date": "2026-10-25",
        "script": "02_脚本/verify_ldha_c8_finalize.py",
        "gene": GENE, "template": "NM_005566.4",
        "reuse_forward_from": "primer_blast_v5_ldha.json (verify_ldha_v5_specificity.py): "
                              "P0=0, tgt_tx=6, LOC=0, n_human=3667",
        "new_reverse_blasted": R,
        "candidates": [[F, R]],
        "results": [entry],
        "selected": entry,
    }
    json.dump(c8, open(f"{OUT}/primer_blast_v5_ldha_candidate8.json", "w",
                       encoding="utf-8"), ensure_ascii=False, indent=1, default=str)
    print(f"\n=== 候选 8 写回：F={F} R={R} product={g['product_size']}bp "
          f"GC差={g['gc_diff']}% verdict={verdict} ===")
    print("DONE — verify_ldha_c8_finalize.py")


if __name__ == "__main__":
    main()
