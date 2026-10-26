# -*- coding: utf-8 -*-
"""
verify_primers_blast.py — 引物单靶特异性 BLAST 验证（任务 #11d）

方法学要点（本脚本的设计依据，均经实测）：
  1. NCBI **esearch 不支持真正的序列检索**（"[Sequence]" 字段被降级为引号短语
     "[All Fields]" 匹配，实测 count=0）→ 特异性验证必须走 BLAST API。
  2. BLAST API 的 FORMAT_TYPE 仅 **Text** 可用（JSON/CSV/Tabular 均返回空或
     "Invalid FORMAT_TYPE"）→ 解析 Text 格式，命中列表行以
     "Sequences producing significant alignments" 段落后逐行解析。
  3. refseq_rna 数据库含**所有物种**（实测 primate/plant/insect 大量命中）。
     对"人源组织 qPCR"而言，非人源序列**不构成脱靶**（人体样本中不存在该模板），
     故判据分两层：
       - **真脱靶（off-target）**：人源 Homo sapiens 转录本中命中非目标基因 → 计入脱靶
       - **跨物种命中（cross-species）**：非人源 → 单独记录，不计入脱靶
  4. 判为特异的条件：人源命中基因集合 == {目标基因}，且正反向引物最高比对
     identity 均 ≥ 90%（人源目标命中）。

输出：03_结果/实验验证设计/primer_blast_specificity.csv / .json
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
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
UA = {"User-Agent": "WorkBuddy-PrimerSpecificity/1.0"}
os.makedirs(OUT, exist_ok=True)

HUMAN_RE = re.compile(r"Homo sapiens", re.I)
GENE_RE = re.compile(r"\(([A-Za-z0-9\-\.@_]+)\)")


def blast_put(seq):
    p = {"CMD": "Put", "PROGRAM": "blastn", "DATABASE": "refseq_rna",
         "QUERY": seq, "MEGARUN": "on", "HITLIST_SIZE": "1000",
         "FILTER": "F", "EXPECT": "1000", "WORD_SIZE": 7,
         "NEIGHBOURING": "0"}
    req = urllib.request.Request(BLAST, data=urllib.parse.urlencode(p).encode(),
                                 headers=UA)
    t = urllib.request.urlopen(req, timeout=150).read().decode(errors="ignore")
    rid = re.search(r"RID = (\S+)", t)
    if not rid:
        raise RuntimeError("BLAST Put 未返回 RID")
    return rid.group(1)


def blast_status(rid):
    u = BLAST + "?" + urllib.parse.urlencode(
        {"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid})
    t = urllib.request.urlopen(u, timeout=150).read().decode(errors="ignore")
    m = re.search(r"Status=([A-Z]+)", t)
    return m.group(1) if m else "UNKNOWN"


def blast_text(rid):
    u = BLAST + "?" + urllib.parse.urlencode(
        {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "Text",
         "ALIGNMENTS": "1000", "DESCRIPTIONS": "1000"})
    t = urllib.request.urlopen(u, timeout=240).read().decode(errors="ignore")
    m = re.search(r"<PRE>(.*?)</PRE>", t, re.S)
    return m.group(1) if m else t


def wait_ready(rid, max_wait=900):
    t0 = time.time()
    time.sleep(12)
    while time.time() - t0 < max_wait:
        st = blast_status(rid)
        if st in ("READY", "FAILED"):
            return st
        time.sleep(10)
    return "TIMEOUT"


HIT_LINE = re.compile(
    r"^(\S+)\s+(.{10,}?)\s+([\d.]+)\s+([\deE\.\-]+)\s+(\d+)%")


def parse_hits(body):
    """解析 BLAST Text 命中列表段"""
    lines = body.split("\n")
    try:
        start = next(i for i, l in enumerate(lines)
                     if "Sequences producing significant alignments" in l)
    except StopIteration:
        return []
    hits = []
    for l in lines[start + 2:]:
        if not l.strip():
            if hits:
                break
            continue
        m = HIT_LINE.match(l.rstrip())
        if not m:
            break
        acc, title, bits, ev, ident = m.groups()
        g = GENE_RE.search(title)
        hits.append({
            "acc": acc, "title": title.strip()[:160],
            "gene": g.group(1) if g else "",
            "bits": float(bits), "evalue": float(ev),
            "identity_pct": int(ident),
            "human": bool(HUMAN_RE.search(title)),
        })
    return hits


def run_one(seq):
    rid = blast_put(seq)
    st = wait_ready(rid)
    if st != "READY":
        return {"status": st, "hits": []}
    body = blast_text(rid)
    hits = parse_hits(body)
    return {"status": st, "hits": hits, "rid": rid}


def main():
    df = pd.read_csv(f"{OUT}/qpcr_primers.csv")
    rows = []
    for _, r in df.iterrows():
        gene = r["gene"]
        rec = {"gene": gene, "role": r["role"],
               "forward_primer_5to3": r["forward_primer_5to3"],
               "reverse_primer_5to3": r["reverse_primer_5to3"]}
        for label, seq in [("forward", r["forward_primer_5to3"]),
                           ("reverse", r["reverse_primer_5to3"])]:
            try:
                res = run_one(seq)
                hits = res["hits"]
                rec[f"{label}_blast_status"] = res["status"]
                rec[f"{label}_n_hits_total"] = len(hits)
                human = [h for h in hits if h["human"]]
                rec[f"{label}_n_hits_human"] = len(human)
                rec[f"{label}_n_hits_nonshuman"] = len(hits) - len(human)
                hg = sorted({h["gene"] for h in human if h["gene"]})
                rec[f"{label}_human_hit_genes"] = ";".join(hg)
                rec[f"{label}_human_hit_genes_n"] = len(hg)
                off = [g for g in hg if g != gene]
                rec[f"{label}_human_offtarget_genes"] = ";".join(off)
                rec[f"{label}_human_offtarget_n"] = len(off)
                # 目标基因命中的一致性（identity / accession）
                tgt = [h for h in human if h["gene"] == gene]
                rec[f"{label}_target_human_hits"] = len(tgt)
                rec[f"{label}_target_max_identity"] = max(
                    [h["identity_pct"] for h in tgt], default=0)
                rec[f"{label}_target_top_acc"] = tgt[0]["acc"] if tgt else ""
                # 跨物种命中仅记录数量
                rec[f"{label}_cross_species_hits"] = len(hits) - len(human)
            except Exception as e:                       # noqa: BLE001
                rec[f"{label}_error"] = str(e)[:220]
            time.sleep(4)          # 尊重 NCBI 速率限制
        off_n = rec.get("forward_human_offtarget_n", 0) + \
            rec.get("reverse_human_offtarget_n", 0)
        rec["human_offtarget_total"] = off_n
        rec["specific"] = bool(
            off_n == 0
            and rec.get("forward_human_hit_genes_n", 0) >= 1
            and rec.get("reverse_human_hit_genes_n", 0) >= 1
            and rec.get("forward_target_max_identity", 0) >= 90
            and rec.get("reverse_target_max_identity", 0) >= 90)
        rows.append(rec)
        print(f"{gene:8s} F[human {rec.get('forward_n_hits_human','?')} hits, "
              f"genes={rec.get('forward_human_hit_genes','?')[:36]}, "
              f"off={rec.get('forward_human_offtarget_n','?')}] "
              f"R[human {rec.get('reverse_n_hits_human','?')}, "
              f"genes={rec.get('reverse_human_hit_genes','?')[:36]}, "
              f"off={rec.get('reverse_human_offtarget_n','?')}] "
              f"xsp={rec.get('forward_cross_species_hits','?')}/"
              f"{rec.get('reverse_cross_species_hits','?')} "
              f"-> specific={rec['specific']}", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(f"{OUT}/primer_blast_specificity.csv", index=False)
    with open(f"{OUT}/primer_blast_specificity.json", "w", encoding="utf-8") as fh:
        json.dump({"date": "2026-10-25",
                   "database": "refseq_rna (all species)",
                   "program": "blastn MEGARUN, expect=1000, hitlist=1000, word_size=7",
                   "criteria": ("off-target counted ONLY within Homo sapiens "
                                "transcripts; cross-species hits recorded but NOT "
                                "counted as off-target; target identity >= 90%"),
                   "results": rows}, fh, ensure_ascii=False, indent=1, default=str)
    print(f"\nBLAST 特异性验证: {int(out.specific.sum())}/{len(out)} 对通过"
          f"（人源层面无非目标基因命中）")
    print("DONE — verify_primers_blast.py")


if __name__ == "__main__":
    main()
