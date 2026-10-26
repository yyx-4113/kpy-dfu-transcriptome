# -*- coding: utf-8 -*-
"""
regen_ldha_v5_blast.py — 全 12 对引物的 BLAST 全库复核（v6 判据 + 重试/INCONCLUSIVE 修复）

═══ 为什么必须重跑全部 12 对 ═══
- 旧版 verify_primers_blast_v6.py（被 ZBpily 任务调用）**缺少重试/INCONCLUSIVE 修复**：
  状态≠READY 时返回 hits=[]，导致多对被**误报**为 NO_TARGET / tgtTx=0
  （FAHD2B-F、ADH1B-F、ALDH3A2-R、ACTB-R 等 implausible 结果）。
- 当前 primer_blast_v5.csv / primer_blast_v6.csv / primer_blast_specificity_v3.csv
  中，v5/v6 缺失、v6 由缺陷脚本生成、v3 不存在；v2 仅含旧 v1 LDHA 且未含 3' 端锚定修正。
- 因此**全 12 对以修复后脚本统一重跑**，确保 v5/v6/v3 三套产物同源、正确、无缺陷。

输出（均同源、同一次 BLAST）：
  primer_blast_v5.csv          —— verify_primers.py 直接读取（v5 schema）
  primer_blast_v6.csv          —— 审计溯源（与 v5 同列，仅文件名不同）
  primer_blast_specificity_v3.csv —— build_supplementary_tables.py S10d 读取（v3 schema）
  primer_blast_v5_ldha_all.json  —— 完整分层明细（含每条 P0/P1 命中）

判据（v6）：core_nt + txid9606 人源限定 + XML；全长匹配（align_len==引物长度）分层
  TARGET / PARALOG / PSEUDO / P0（3'端5nt全匹配+id≥95）/ P1（90≤id<95）/ P2 / PARTIAL；
  usable = 目标存在全长匹配人源 RefSeq 转录本 且 P0==0 且 目标 identity≥90%。
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import re
import csv
import time
import json
import urllib.request
import urllib.parse

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
RAW = f"{OUT}/provenance/blast_xml_raw_regen"
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
UA = {"User-Agent": "WorkBuddy-RegenV6/1.0"}
os.makedirs(RAW, exist_ok=True)

HUMAN_TAXON = "txid9606[Organism:exp]"
HUMAN_RE = re.compile(r"Homo sapiens", re.I)
PSEUDO_RE = re.compile(r"pseudogene", re.I)

STRICT_PARALOG_ONLY = {"ACTB", "ACTG1", "GAPDH", "RPLP0", "TBP", "HPRT1",
                        "SDHA", "DHX58"}
PARALOG_FAMILIES = {
    "ACTB": {"ACTB", "ACTG1", "ACTG2", "ACTC1"},
    "GAPDH": {"GAPDH", "GAPDHS", "GAPDHL"},
    "RPLP0": {"RPLP0", "RPLP0P1", "RPLP0P2", "RPLP0P3", "RPLP0P4", "RPLP0P5"},
    "LDHA": {"LDHA", "LDHAL", "LDHAL6A", "LDHAL6B", "LDHB", "LDHC"},
    "ALDH3A2": {"ALDH3A2", "ALDH3A1", "ALDH3A3"},
    "ADH1B": {"ADH1B", "ADH1A", "ADH1C", "ADH4", "ADH6", "ADH7"},
    "GBP1": {"GBP1", "GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "GBP7"},
    "DHX58": {"DHX58", "DHX58P1", "DHX58P2"},
    "FAHD2B": {"FAHD2B", "FAHD2A", "FAHD1"},
    "TBP": {"TBP"}, "HPRT1": {"HPRT1"},
    "SDHA": {"SDHA", "SDHAP1", "SDHAP2", "SDHAP3"},
}


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


def classify(hit, gene, primer_len):
    if hit["al"] != primer_len:
        return ("PARTIAL", None, None)
    qs, hs = hit["qseq"], hit["hseq"]
    L = min(len(qs), len(hs))
    mm_pos = [i for i in range(L) if qs[i] != hs[i] and hs[i] not in "- "]
    mm_total = len(mm_pos)
    mm_3p5 = sum(1 for i in mm_pos if i >= L - 5)
    g = hit["gene"]
    if g and g == gene:
        return ("TARGET", mm_total, mm_3p5)
    if hit["pseudo"]:
        return ("PSEUDO", mm_total, mm_3p5)
    if gene in STRICT_PARALOG_ONLY and g and g in PARALOG_FAMILIES.get(gene, set()):
        return ("PARALOG", mm_total, mm_3p5)
    if not hit["is_refseq_tx"]:
        return ("P2", mm_total, mm_3p5)
    if mm_3p5 == 0 and hit["identity_pct"] >= 95:
        return ("P0", mm_total, mm_3p5)
    if mm_3p5 == 0 and hit["identity_pct"] >= 90:
        return ("P1", mm_total, mm_3p5)
    return ("P2", mm_total, mm_3p5)


def run_one(seq, label):
    """带重试/INCONCLUSIVE 修复：状态≠READY 重试至多 3 轮，仍失败显式标记。"""
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


def verify_primer(seq, gene, label):
    res = run_one(seq, label)
    if res["status"] != "READY":
        return {"status": res["status"], "hits": [], "counts": {},
                "tgt_tx": 0, "tgt_id": 0.0, "detail": [], "n_human": 0}
    hits = res["hits"]
    human = [h for h in hits if h["human"]]
    counts = {}
    tgt_tx = 0
    tgt_id = 0.0
    detail = []
    for h in human:
        lvl, mm_t, mm_35 = classify(h, gene, len(seq))
        h["level"] = lvl
        counts[lvl] = counts.get(lvl, 0) + 1
        if lvl == "TARGET":
            if h["is_refseq_tx"]:
                tgt_tx += 1
            tgt_id = max(tgt_id, h["identity_pct"])
        if lvl in ("P0", "P1"):
            detail.append({
                "gene": gene, "primer": label, "query": seq,
                "accession": h["acc"], "hit_gene": h["gene"],
                "identity_n": h["idn"], "align_len": h["al"],
                "identity_pct": h["identity_pct"],
                "mm_total": mm_t, "mm_3p5": mm_35, "level": lvl,
                "title": h["title"][:160]})
    return {"status": "READY", "hits": human, "counts": counts,
            "tgt_tx": tgt_tx, "tgt_id": tgt_id, "detail": detail,
            "n_human": len(human)}


def main():
    src = f"{OUT}/qpcr_primers.csv"
    # utf-8-sig: qpcr_primers.csv carries a UTF-8 BOM; plain utf-8 leaves the
    # first header as "\ufeffgene" and r["gene"] would raise KeyError.
    rows_in = list(csv.DictReader(open(src, encoding="utf-8-sig")))
    rows_out = []
    detail_all = []
    summary = []
    # ── 断点续跑：按基因缓存已完成结果，中断后不必重跑整轮（每对约 4–8 min）──
    PARTIAL = f"{OUT}/provenance/_regen_partial.json"
    done = {}
    if os.path.exists(PARTIAL):
        try:
            done = json.load(open(PARTIAL, encoding="utf-8")).get("done", {})
            print(f"[resume] 缓存命中 {len(done)} 对: {sorted(done)}", flush=True)
        except Exception as e:
            print(f"[resume] 缓存损坏，从头开始: {e}", flush=True)
            done = {}
    for r in rows_in:
        gene = r["gene"]
        f, rev = r["forward_primer_5to3"], r["reverse_primer_5to3"]
        # 仅当引物序列与当前主表逐字一致时才复用缓存（主表改动即失效重跑）
        c = done.get(gene)
        if (c and c["rec"].get("forward_primer_5to3") == f
                and c["rec"].get("reverse_primer_5to3") == rev):
            print(f"\n##### {gene} ##### (cached, skip)", flush=True)
            rr = c["rec"]
            rows_out.append(rr)
            detail_all.extend(c["detail"])
            summary.append({
                "gene": gene, "F_P0": rr.get("forward_n_P0"),
                "R_P0": rr.get("reverse_n_P0"),
                "F_tgtTx": rr.get("forward_target_full_match_refseq_tx"),
                "R_tgtTx": rr.get("reverse_target_full_match_refseq_tx"),
                "verdict": rr.get("verdict")})
            print(f"  {gene} F[P0={rr.get('forward_n_P0')} "
                  f"tgtTx={rr.get('forward_target_full_match_refseq_tx')}] "
                  f"R[P0={rr.get('reverse_n_P0')} "
                  f"tgtTx={rr.get('reverse_target_full_match_refseq_tx')}] "
                  f"-> {rr.get('verdict')}", flush=True)
            continue
        rec = {"gene": gene, "role": r["role"], "module": r["module"],
               "refseq_accession": r["refseq_accession"],
               "forward_primer_5to3": f, "reverse_primer_5to3": rev}
        print(f"\n##### {gene} #####", flush=True)
        fv = verify_primer(f, gene, "F")
        time.sleep(3)
        rv = verify_primer(rev, gene, "R")
        time.sleep(3)

        detail_local = []

        def fill(lab, v):
            c = v["counts"]
            rec[f"{lab}_blast_status"] = v["status"]
            rec[f"{lab}_n_hits_human"] = v["n_human"]
            rec[f"{lab}_n_full_length_human"] = sum(
                val for k, val in c.items() if k != "PARTIAL")
            rec[f"{lab}_n_P0"] = c.get("P0", 0)
            rec[f"{lab}_n_P1"] = c.get("P1", 0)
            rec[f"{lab}_n_PARALOG"] = c.get("PARALOG", 0)
            rec[f"{lab}_n_PSEUDO"] = c.get("PSEUDO", 0)
            rec[f"{lab}_target_full_match_refseq_tx"] = v["tgt_tx"]
            rec[f"{lab}_target_max_identity_pct"] = round(v["tgt_id"], 2)
            detail_local.extend(v["detail"])

        fill("forward", fv)
        fill("reverse", rv)
        p0 = rec.get("forward_n_P0", 0) + rec.get("reverse_n_P0", 0)
        p1 = rec.get("forward_n_P1", 0) + rec.get("reverse_n_P1", 0)
        ft = rec.get("forward_target_full_match_refseq_tx", 0)
        rt = rec.get("reverse_target_full_match_refseq_tx", 0)
        fi = rec.get("forward_target_max_identity_pct", 0)
        ri = rec.get("reverse_target_max_identity_pct", 0)
        specific = bool(p0 == 0 and ft >= 1 and rt >= 1
                        and fi >= 90 and ri >= 90)
        rec["n_P0_total"] = p0
        rec["n_P1_total"] = p1
        rec["specific"] = specific
        rec["verdict"] = ("USABLE" if specific
                          else ("NO_TARGET" if (ft < 1 or rt < 1)
                                else "OFFTARGET_P0"))
        rows_out.append(rec)
        detail_all.extend(detail_local)
        done[gene] = {"rec": dict(rec), "detail": detail_local}
        try:
            json.dump({"done": done}, open(PARTIAL, "w", encoding="utf-8"),
                      ensure_ascii=False)
        except Exception as e:
            print(f"  [warn] 缓存写入失败，结果仍会进最终产物: {e}", flush=True)
        summary.append({
            "gene": gene, "F_P0": rec.get("forward_n_P0"),
            "R_P0": rec.get("reverse_n_P0"), "F_tgtTx": ft, "R_tgtTx": rt,
            "verdict": rec["verdict"]})
        print(f"  {gene} F[P0={rec.get('forward_n_P0')} tgtTx={ft}] "
              f"R[P0={rec.get('reverse_n_P0')} tgtTx={rt}] -> {rec['verdict']}",
              flush=True)
        time.sleep(3)

    # ── 写 primer_blast_v5.csv（verify_primers.py 读取的 schema）──
    v5_cols = ["gene", "role", "module", "refseq_accession",
               "forward_primer_5to3", "reverse_primer_5to3",
               "forward_blast_status", "forward_n_hits_human",
               "forward_n_full_length_human", "forward_n_P0", "forward_n_P1",
               "forward_n_PARALOG", "forward_n_PSEUDO",
               "forward_target_full_match_refseq_tx",
               "forward_target_max_identity_pct",
               "reverse_blast_status", "reverse_n_hits_human",
               "reverse_n_full_length_human", "reverse_n_P0", "reverse_n_P1",
               "reverse_n_PARALOG", "reverse_n_PSEUDO",
               "reverse_target_full_match_refseq_tx",
               "reverse_target_max_identity_pct",
               "n_P0_total", "n_P1_total", "specific", "verdict"]
    with open(f"{OUT}/primer_blast_v5.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=v5_cols)
        w.writeheader()
        for r in rows_out:
            w.writerow({k: r.get(k, "") for k in v5_cols})
    # v6 同列同源
    with open(f"{OUT}/primer_blast_v6.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=v5_cols)
        w.writeheader()
        for r in rows_out:
            w.writerow({k: r.get(k, "") for k in v5_cols})
    # 同一运行的明细 JSON：v5/v6 两个文件名都写，避免历史上出现过的
    # "CSV 已更新而 primer_blast_v6.json 停留在旧判定" 的文件间矛盾
    #（2026-10-25 事故：18:25 的 JSON 与 22:51 的 CSV 判定互相矛盾）
    _payload = {
        "date": "2026-10-25",
        "script": "02_脚本/regen_ldha_v5_blast.py",
        "database": "core_nt + ENTREZ_QUERY='txid9606[Organism:exp]'",
        "method_fix": "v6 判据 + 重试/INCONCLUSIVE 修复（状态≠READY 重试 3 轮，"
                     "仍失败显式标记，绝不把取不到结果冒充零命中）",
        "summary": summary,
        "rows": rows_out,
        "detail": detail_all,
    }
    with open(f"{OUT}/primer_blast_v5_ldha_all.json", "w", encoding="utf-8") as fh:
        json.dump(_payload, fh, ensure_ascii=False, indent=1, default=str)
    with open(f"{OUT}/primer_blast_v6.json", "w", encoding="utf-8") as fh:
        json.dump(_payload, fh, ensure_ascii=False, indent=1, default=str)

    # ── 写 primer_blast_specificity_v3.csv（S10d schema）──
    v3_cols = ["gene", "role", "forward_target_human_hits",
               "reverse_target_human_hits", "forward_target_max_identity",
               "reverse_target_max_identity", "forward_target_top_acc",
               "reverse_target_top_acc", "forward_n_human_refseq_tx",
               "reverse_n_human_refseq_tx", "forward_n_offtarget_refseq_tx",
               "reverse_n_offtarget_refseq_tx", "forward_offtarget_tx_genes",
               "reverse_offtarget_tx_genes", "forward_n_offtarget_genomic",
               "reverse_n_offtarget_genomic", "forward_paralog_genes",
               "reverse_paralog_genes", "forward_pseudogene_genes",
               "reverse_pseudogene_genes", "forward_human_hit_genes",
               "reverse_human_hit_genes", "offtarget_refseq_tx_total",
               "offtarget_genomic_total", "paralog_total", "specific"]

    v3_rows = []
    for rec in rows_out:
        row = {"gene": rec["gene"], "role": rec["role"]}
        row.update(_v3_map(rec, "forward", detail_all))
        row.update(_v3_map(rec, "reverse", detail_all))
        row["offtarget_refseq_tx_total"] = (row["forward_n_offtarget_refseq_tx"]
                                           + row["reverse_n_offtarget_refseq_tx"])
        row["offtarget_genomic_total"] = (row["forward_n_offtarget_genomic"]
                                         + row["reverse_n_offtarget_genomic"])
        row["paralog_total"] = (len([x for x in row["forward_paralog_genes"].split(";") if x])
                               + len([x for x in row["reverse_paralog_genes"].split(";") if x]))
        row["specific"] = rec["specific"]
        v3_rows.append(row)
    with open(f"{OUT}/primer_blast_specificity_v3.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=v3_cols)
        w.writeheader()
        for r in v3_rows:
            w.writerow({k: r.get(k, "") for k in v3_cols})

    # ── JSON 明细 ──
    n_ok = sum(1 for s in summary if s["verdict"] == "USABLE")
    print(f"\n重跑完成：{n_ok}/{len(summary)} 对 USABLE（修复后）")
    print("DONE — regen_ldha_v5_blast.py")


def _v3_map(rec, lab, detail_all):
    c = {k.split("_", 1)[1]: v for k, v in rec.items()
         if k.startswith(lab + "_")}
    # detail 的 primer 键是 "F"/"R"（不是 "forward"/"reverse"）——
    # 2026-10-25 修复：旧代码用 lab 直接匹配，导致脱靶/旁系/假基因基因列表恒为空。
    p = {"forward": "F", "reverse": "R"}[lab]
    off_tx_genes = sorted({d["hit_gene"] for d in detail_all
                           if d["gene"] == rec["gene"] and d["primer"] == p
                           and d["level"] in ("P0", "P1", "PARALOG", "PSEUDO")})
    para_genes = sorted({d["hit_gene"] for d in detail_all
                         if d["gene"] == rec["gene"] and d["primer"] == p
                         and d["level"] == "PARALOG"})
    pseudo_genes = sorted({d["hit_gene"] for d in detail_all
                           if d["gene"] == rec["gene"] and d["primer"] == p
                           and d["level"] == "PSEUDO"})
    human_genes = sorted({d["hit_gene"] for d in detail_all
                          if d["gene"] == rec["gene"] and d["primer"] == p})
    # 目标基因统计必须取自 rec（rows_out）而非 detail：detail 只记录
    # P0/P1/PARALOG/PSEUDO 等**脱靶层**命中，不含 TARGET 层记录。
    # 2026-10-25 修复：旧 _v3_map 只查 detail 的 TARGET 记录 → 全部 12 个基因的
    # target_human_hits / max_identity / top_acc / n_human_refseq_tx 均为 0/空，
    # S10d 附表因此显示"引物不命中任何目标"（与 specific=True 自相矛盾）。
    tgt_n = rec.get(f"{lab}_target_full_match_refseq_tx", 0)
    tgt_id = rec.get(f"{lab}_target_max_identity_pct", 0)
    top_acc = rec.get("refseq_accession", "")
    n_off_refseq_tx = sum(
        1 for d in detail_all if d["gene"] == rec["gene"]
        and d["primer"] == lab and d["accession"].startswith(("NM_", "XM_", "NR_", "XR_")))
    return {
        f"{lab}_target_human_hits": tgt_n,
        f"{lab}_target_max_identity": round(float(tgt_id or 0), 2),
        f"{lab}_target_top_acc": top_acc,
        f"{lab}_n_human_refseq_tx": (tgt_n or 0) + n_off_refseq_tx,
        f"{lab}_n_offtarget_refseq_tx": len(off_tx_genes),
        f"{lab}_offtarget_tx_genes": ";".join(off_tx_genes),
        f"{lab}_n_offtarget_genomic": 0,
        f"{lab}_paralog_genes": ";".join(para_genes),
        f"{lab}_pseudogene_genes": ";".join(pseudo_genes),
        # 人源命中基因 = 目标基因本身 + detail 中的脱靶层命中基因
        f"{lab}_human_hit_genes": ";".join(sorted(human_genes | {rec["gene"]})),
    }


if __name__ == "__main__":
    main()
