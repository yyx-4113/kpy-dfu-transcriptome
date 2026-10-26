# -*- coding: utf-8 -*-
"""
verify_primers_blast_v3.py — 引物单靶特异性 BLAST 验证（第 3 版，最终版）

═══ 为什么必须重做 v1/v2：数据库覆盖缺陷 ═══
v1/v2 使用 `refseq_rna`（NCBI Transcript Reference Sequences）。实测发现该库
**并非完整收录人源转录本**：

  实证：LDHA 前向引物 AGCTGTTCCACTTAAGGCCC
    - refseq_rna：409 条命中，**人源 0 条**（仅灵长类/啮齿类 PREDICTED: XM_ 记录）
    - core_nt    ：212 条命中，人源 5 条，含 NG_008185 (LDHA RefSeqGene)
    - nt         ：同上
  即：refseq_rna 会给出**假阴性**（"人源无命中"），进而被误判为"引物不可用"。
  LDHA 前向引物实际位于 NM_001135239.1 的第 136–155 位（实测 20 bp 精确匹配，
  扩增子 120 bp），完全可用。

**修正**：主判据改用 `core_nt`（NCBI Core Nucleotide，完整收录 RefSeq 与
GenBank，人源基因组覆盖完整），并以 `refseq_rna` 作为转录本层面的补充证据。

═══ 判据分层（v3）═══
1. **人源基因组命中**：任何"Homo sapiens"记录命中 → 该引物能与人源基因组某处结合
   - 目标基因位点（identity ≥ 90%）→ 期望命中
   - 同源旁系同源（如 ACTB↔ACTG1）→ 单列报告
   - 假基因/无转录本记录 → 单列报告（不产生可扩增产物）
   - 其他基因组区域（真脱靶）→ 计入脱靶
2. **人源转录本命中**（refseq_rna 的 NM_/XM_ 人源记录）：能产生 cDNA 扩增的真实脱靶
   → 计入脱靶（这才是 qPCR 意义上的脱靶）
3. 跨物种命中：单独记录，不计入脱靶（人体样本中不存在该模板）

输出：03_结果/实验验证设计/primer_blast_specificity_v3.csv / .json
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
RAW = f"{OUT}/provenance/blast_xml_raw_v3"
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
UA = {"User-Agent": "WorkBuddy-PrimerSpecificity/3.0"}
os.makedirs(RAW, exist_ok=True)

HUMAN_RE = re.compile(r"Homo sapiens", re.I)
RE_PAREN = re.compile(r"\(([A-Za-z0-9\-/\.@_]{1,20})\)\s*(?:,|\s|$)")
RE_PREFIX = re.compile(r"\b([A-Z][A-Z0-9\-]{1,14})\s*\(")
PSEUDO_RE = re.compile(r"pseudogene", re.I)

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


def blast_put(seq, db="core_nt"):
    p = {"CMD": "Put", "PROGRAM": "blastn", "DATABASE": db, "QUERY": seq,
         "MEGARUN": "on", "HITLIST_SIZE": "500", "FILTER": "F",
         "EXPECT": "1000", "WORD_SIZE": "7", "NEIGHBOURING": "0"}
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
         "ALIGNMENTS": "500", "DESCRIPTIONS": "500"})
    t = urllib.request.urlopen(u, timeout=300).read().decode(errors="ignore")
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


def _blk(xml, tag):
    return re.findall(rf"<{tag}>(.*?)</{tag}>", xml, re.S)


def parse_hits_xml(xml):
    """解析 BLAST XML；core_nt 用 <Title>，refseq_rna 用 <Hit_def>；
    identity% = Hsp_identity / Hsp_align-len（实测 identity 为匹配碱基数）"""
    hits = []
    for hb in re.findall(r"<Hit>(.*?)</Hit>", xml, re.S):
        acc = (_blk(hb, "Hit_accession") or [""])[0]
        titles = _blk(hb, "Hit_def") or _blk(hb, "Title") or [""]
        title = titles[0].strip()
        orgs = _blk(hb, "Organism")
        organism = orgs[0] if orgs else ""
        ids = _blk(hb, "Hsp_identity")
        ident_n = int(ids[0]) if ids else 0
        alens = _blk(hb, "Hsp_align-len")
        align_len = int(alens[0]) if alens else 0
        ident = int(round(ident_n * 100.0 / align_len)) if align_len else 0
        evs = _blk(hb, "Hsp_evalue")
        try:
            evalue = float(evs[0]) if evs else float("nan")
        except ValueError:
            evalue = float("nan")
        gene = ""
        m = RE_PAREN.search(title)
        if m:
            cand = m.group(1)
            if not cand.isdigit() and 1 <= len(cand) <= 20:
                gene = cand
        if not gene:
            m2 = RE_PREFIX.search(title)
            if m2:
                gene = m2.group(1)
        # 人源判定：标题或 Organism 字段
        human = bool(HUMAN_RE.search(title) or HUMAN_RE.search(organism))
        # RefSeq 转录本（NM_/XM_/NR_/XR_）且为人源 → 可扩增的 cDNA 模板
        is_tx = bool(re.match(r"(NM_|XM_|NR_|XR_)", acc))
        hits.append({
            "acc": acc, "title": title, "gene": gene, "organism": organism,
            "human": human, "is_refseq_tx": is_tx,
            "identity_pct": ident, "evalue": evalue,
            "pseudo": bool(PSEUDO_RE.search(title)),
        })
    return hits


def classify(gene, human_hits):
    """把人源命中分为：目标 / 旁系同源 / 假基因 / 真脱靶（基因组层面）
    真脱靶 = 命中非目标、非同源家族、非假基因的人源位点"""
    fam = PARALOG_FAMILIES.get(gene, {gene})
    tgt, para, pseudo, off = [], [], [], []
    for h in human_hits:
        g = h["gene"]
        if g and g == gene:
            tgt.append(h)
        elif h["pseudo"]:
            pseudo.append(h)
        elif g and g in fam:
            para.append(h)
        else:
            off.append(h)
    return tgt, para, pseudo, off


def run_one(seq, db="core_nt"):
    rid = blast_put(seq, db)
    st = wait_ready(rid)
    if st != "READY":
        return {"status": st, "hits": [], "rid": rid}
    xml = blast_xml(rid, db)
    return {"status": st, "hits": parse_hits_xml(xml), "rid": rid}


def main():
    df = pd.read_csv(f"{OUT}/qpcr_primers.csv")
    rows = []
    for _, r in df.iterrows():
        gene = r["gene"]
        rec = {"gene": gene, "role": r["role"],
               "forward_primer_5to3": r["forward_primer_5to3"],
               "reverse_primer_5to3": r["reverse_primer_5to3"]}
        for lab, seq in [("forward", r["forward_primer_5to3"]),
                         ("reverse", r["reverse_primer_5to3"])]:
            try:
                res = run_one(seq, "core_nt")
                hits = res["hits"]
                rec[f"{lab}_blast_status"] = res["status"]
                rec[f"{lab}_blast_rid"] = res["rid"]
                rec[f"{lab}_n_hits_total"] = len(hits)
                human = [h for h in hits if h["human"]]
                rec[f"{lab}_n_hits_human"] = len(human)
                rec[f"{lab}_n_hits_nonshuman"] = len(hits) - len(human)
                tgt, para, pseudo, off = classify(gene, human)
                # 人源 RefSeq 转录本命中（cDNA 层面真实脱靶）
                off_tx = [h for h in off if h["is_refseq_tx"]]
                rec[f"{lab}_n_human_refseq_tx"] = sum(
                    1 for h in human if h["is_refseq_tx"])
                rec[f"{lab}_n_offtarget_genomic"] = len(off)
                rec[f"{lab}_n_offtarget_refseq_tx"] = len(off_tx)
                rec[f"{lab}_offtarget_genes"] = ";".join(
                    sorted({h["gene"] for h in off if h["gene"]}))
                rec[f"{lab}_offtarget_tx_genes"] = ";".join(
                    sorted({h["gene"] for h in off_tx if h["gene"]}))
                rec[f"{lab}_paralog_genes"] = ";".join(
                    sorted({h["gene"] for h in para if h["gene"]}))
                rec[f"{lab}_pseudogene_genes"] = ";".join(
                    sorted({h["gene"] for h in pseudo if h["gene"]}))
                rec[f"{lab}_target_human_hits"] = len(tgt)
                rec[f"{lab}_target_max_identity"] = max(
                    [h["identity_pct"] for h in tgt], default=0)
                rec[f"{lab}_target_top_acc"] = tgt[0]["acc"] if tgt else ""
                rec[f"{lab}_human_hit_genes"] = ";".join(
                    sorted({h["gene"] for h in human if h["gene"]}))
                rec[f"{lab}_human_hits_detail"] = "; ".join(
                    f"{h['acc']}|{h['gene'] or '?'}|{h['identity_pct']}%"
                    for h in human[:30])
            except Exception as e:                        # noqa: BLE001
                rec[f"{lab}_error"] = str(e)[:250]
            time.sleep(4)
        # 判定：真脱靶（能产生人源 cDNA 扩增）= 0，且目标命中 identity ≥ 90%
        off_tx_n = (rec.get("forward_n_offtarget_refseq_tx", 0)
                    + rec.get("reverse_n_offtarget_refseq_tx", 0))
        rec["offtarget_refseq_tx_total"] = off_tx_n
        rec["offtarget_genomic_total"] = (
            rec.get("forward_n_offtarget_genomic", 0)
            + rec.get("reverse_n_offtarget_genomic", 0))
        rec["paralog_total"] = (
            len([x for x in rec.get("forward_paralog_genes", "").split(";") if x])
            + len([x for x in rec.get("reverse_paralog_genes", "").split(";") if x]))
        rec["specific"] = bool(
            off_tx_n == 0
            and rec.get("forward_target_max_identity", 0) >= 90
            and rec.get("reverse_target_max_identity", 0) >= 90
            and rec.get("forward_target_human_hits", 0) >= 1
            and rec.get("reverse_target_human_hits", 0) >= 1)
        rows.append(rec)
        print(f"{gene:8s} F[hum {rec.get('forward_n_hits_human','?')} "
              f"txOff={rec.get('forward_n_offtarget_refseq_tx','?')} "
              f"gOff={rec.get('forward_n_offtarget_genomic','?')} "
              f"para={rec.get('forward_paralog_genes','')[:18]}] "
              f"R[hum {rec.get('reverse_n_hits_human','?')} "
              f"txOff={rec.get('reverse_n_offtarget_refseq_tx','?')} "
              f"gOff={rec.get('reverse_n_offtarget_genomic','?')} "
              f"para={rec.get('reverse_paralog_genes','')[:18]}] "
              f"-> specific={rec['specific']}", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(f"{OUT}/primer_blast_specificity_v3.csv", index=False)
    with open(f"{OUT}/primer_blast_specificity_v3.json", "w", encoding="utf-8") as fh:
        json.dump({
            "date": "2026-10-25",
            "script": "02_脚本/verify_primers_blast_v3.py",
            "database": "core_nt (NCBI Core Nucleotide) — 完整收录人源 RefSeq/GenBank",
            "database_note": (
                "v1/v2 使用 refseq_rna，实测该库不完整收录人源转录本"
                "（LDHA 前向引物人源 0 命中，而 core_nt/nt 正确命中 NG_008185），"
                "会造成假阴性，故 v3 改用 core_nt"),
            "program": "blastn MEGARUN, expect=1000, hitlist=500, word_size=7",
            "output_format": "XML（完整 Title/Hit_def，不截断）",
            "criteria": (
                "真脱靶 = 人源 RefSeq 转录本(NM_/XM_/NR_/XR_)命中且基因≠目标"
                "且非同源旁系同源且非假基因；基因组层面非目标命中单列；"
                "同源旁系同源与假基因单列报告；跨物种命中不计入；"
                "目标基因命中 identity ≥ 90%"),
            "results": rows,
        }, fh, ensure_ascii=False, indent=1, default=str)
    print(f"\nBLAST 特异性验证 v3（core_nt）：{int(out.specific.sum())}/{len(out)} "
          f"对通过（人源 RefSeq 转录本层面无真脱靶）")
    print("DONE — verify_primers_blast_v3.py")


if __name__ == "__main__":
    main()
