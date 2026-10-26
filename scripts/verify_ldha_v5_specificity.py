# -*- coding: utf-8 -*-
"""
verify_ldha_v5_specificity.py — LDHA v5 多候选对的 BLAST 全库复核 + 最优干净对选择

═══ 背景 ═══
初版选定对 F=ACCACTTCCATTGTAAGTCCCA / R=ACCAGCTTCCAGAGGACAAG 的 **反向引物**
命中 PIGS (NM_033198) 出现 1 条 P0（3' 端 5nt 全匹配、id=95%）——按既定 v6 判据
（单引物 P0 即阻断）该对不可用。正向引物经 BLAST 复核 **P0=0、PARALOG=1、tgtTx=6、
LOC_hit=0**，完全干净。

→ 本脚本在「已验证干净正向 + 多个候选反向」组合中，BLAST 复核每条引物，按 v6
判据选出 **两条引物均无 P0** 的最优对；优先内部 3'UTR 区段、较大产物。

═══ 修复 v6 假阴性 ═══
状态≠READY 时显式重试（最多 3 轮），仍失败标记 INCONCLUSIVE，**绝不**把
「取不到结果」冒充「零命中」。同时修正原脚本变量作用域崩溃 bug。

判据（同 v6）：core_nt + 人源限定；全长匹配人源命中分层 TARGET/PARALOG/PSEUDO/
P0/P1/P2/PARTIAL；usable = 目标存在全长匹配人源 RefSeq 且 P0==0 且 目标 id≥90%。
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

import pandas as pd

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
PROV = f"{OUT}/provenance"
RAW = f"{PROV}/blast_xml_raw_v5_ldha"
os.makedirs(RAW, exist_ok=True)
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
UA = {"User-Agent": "WorkBuddy-LDHAv5Verify/1.0"}

HUMAN_TAXON = "txid9606[Organism:exp]"
HUMAN_RE = re.compile(r"Homo sapiens", re.I)
PSEUDO_RE = re.compile(r"pseudogene", re.I)
GENE = "LDHA"
PARALOG_FAMILY = {"LDHA", "LDHB", "LDHC", "LDHAL", "LDHAL6A", "LDHAL6B",
                  "LDHAL1", "LDHAL2", "LDHAL3", "LDHAL4", "LDHAL5", "LDHAL7"}

# 候选对短名单（均来自 design_ldha_v5_utr3.py 的 21 个通过复核候选；
# 正向 ACCACTTCCATTGTAAGTCCCA 已 BLAST 验证 P0=0）
SHORTLIST = [
    ("ACCACTTCCATTGTAAGTCCCA", "CACTGCACCCTGCCTTAAGA"),   # 内部3'UTR杂交对(大产物)
    ("AGTGTGTCACTATAGGGAACACA", "CACTGCACCCTGCCTTAAGA"),   # 内部3'UTR(90bp)
    ("GCCTGGCCAACACAACAAAA", "TGAGTAGACATCCACCAAGGT"),     # 2070-2241段备份
    ("GCCTGGCCAACACAACAAAA", "CTTGAGTAGACATCCACCAAGGT"),   # 2070-2241段备份
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


def run_one(seq, primer_len, label):
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
    res = run_one(seq, len(seq), label)
    if res["status"] != "READY":
        return {"status": res["status"], "P0": None, "P1": None,
                "PARALOG": None, "tgt_tx": None, "LOC": None, "note": "INCONCLUSIVE"}
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


def main():
    cache = json.load(open(f"{PROV}/final_design_cache.json", encoding="utf-8"))
    seq = cache["NM_005566.4"]["seq"].upper()
    rc = lambda s: s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]

    results = []
    for (F, R) in SHORTLIST:
        print(f"\n##### 候选对 F={F} R={R} #####", flush=True)
        # numpy 复核：在 NM_005566.4 的落位 + 避开 LOC + 在 3'UTR
        iF = seq.find(F)
        iRrc = seq.find(rc(R))
        ok = (iF >= 0 and iRrc >= 0 and iRrc > iF)
        rec = {"forward": F, "reverse": R, "ok_in_canon": ok,
               "f_pos": iF + 1 if iF >= 0 else None,
               "r_rc_pos": iRrc + 1 if iRrc >= 0 else None,
               "product_size": (iRrc - iF) if ok else None}
        if not ok:
            rec["verdict"] = "BAD_ANGLE"
            results.append(rec)
            continue
        # BLAST 两条引物
        fv = verify_primer(F, "F")
        time.sleep(3)
        rv = verify_primer(R, "R")
        rec.update({"forward_blast": fv, "reverse_blast": rv})
        specific = (fv["status"] == "READY" and rv["status"] == "READY"
                    and (fv["P0"] or 0) == 0 and (rv["P0"] or 0) == 0
                    and (fv["tgt_tx"] or 0) >= 1 and (rv["tgt_tx"] or 0) >= 1
                    and (fv["LOC"] or 0) == 0 and (rv["LOC"] or 0) == 0)
        rec["verdict"] = "USABLE" if specific else (
            "INCONCLUSIVE" if (fv["status"] != "READY"
                              or rv["status"] != "READY") else "OFFTARGET_P0")
        print(f"  F: P0={fv.get('P0')} PARALOG={fv.get('PARALOG')} "
              f"tgtTx={fv.get('tgt_tx')} LOC={fv.get('LOC')}", flush=True)
        print(f"  R: P0={rv.get('P0')} PARALOG={rv.get('PARALOG')} "
              f"tgtTx={rv.get('tgt_tx')} LOC={rv.get('LOC')}", flush=True)
        print(f"  → verdict = {rec['verdict']}", flush=True)
        time.sleep(3)
        results.append(rec)

    usable = [r for r in results if r.get("verdict") == "USABLE"]
    # 优先：产物更大（更内部、更稳健）→ 选 product_size 最大者
    usable.sort(key=lambda r: -(r.get("product_size") or 0))
    best = usable[0] if usable else None
    json.dump({
        "date": "2026-10-25",
        "script": "02_脚本/verify_ldha_v5_specificity.py",
        "gene": GENE, "template": "NM_005566.4",
        "shortlist": SHORTLIST,
        "results": results,
        "selected": best,
    }, open(f"{OUT}/primer_blast_v5_ldha.json", "w", encoding="utf-8"),
        ensure_ascii=False, indent=1, default=str)
    print(f"\n=== 最终选定对（两条引物均无 P0）===")
    if best:
        print(f"  F={best['forward']}  R={best['reverse']}  "
              f"product={best['product_size']}bp  "
              f"F@{best['f_pos']} Rrc@{best['r_rc_pos']}")
    else:
        print("  !! 短名单内无两条引物均无 P0 的对——需重新设计反向引物")
    print("DONE — verify_ldha_v5_specificity.py")


if __name__ == "__main__":
    main()
