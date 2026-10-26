# -*- coding: utf-8 -*-
"""
verify_primers_blast_v4.py — 引物单靶特异性 BLAST 验证（第 4 版，最终版）

═══ 为什么必须重做 v1/v2/v3：三处工具层缺陷 ═══

1. **v1 Text 格式截断**：`FORMAT_TYPE=Text` 的命中列表 description 被截断为
   `... uncharacterized... 39.2 1.2 100%`，基因符号落在截断点之后 → genes 为空
   → 误判为「非特异」。v2 改用 XML 修复。

2. **v2 refseq_rna 覆盖缺陷**：该库**未完整收录人源转录本**。实证：LDHA 前向
   引物在 refseq_rna 中人源 0 命中（409 条全为灵长类/啮齿类 XM_ 记录），
   而 core_nt/nt 正确命中 NG_008185 (LDHA RefSeqGene)。→ v3 改用 core_nt。

3. **v3 core_nt 命中列表截断（本版修复的核心问题）**：
   core_nt 查询的命中列表被 `HITLIST_SIZE=500` 截断，而排序按 bit-score 优先，
   **进化上更近的非人源物种会挤占前 500 名**。实证（v3 留档 XML）：
   - 第 3 个查询起，命中数**全部等于 500**（恰好等于上限 → 确证截断）；
   - ACTB 反向引物查询：前 500 名中 `Pseudopipra pipra` 24 条、
     `Neopelma chrysocephalum` 22 条、`Corapipo altera` 22 条占主导，
     人源仅剩 7 条 → 被误判为「人源 0 命中」；
   - 该查询最后一个 XML（`core_nt_DV33BJXZ014.xml`）中
     `Vibrio cholerae` 独占 250/500 条。
   → **v4 使用 NCBI 的物种限定检索**，只查人源记录，从根本上排除截断。

**v4 判据**
- 数据库：`core_nt` + `ENTREZ_QUERY="txid9606[Organism:exp]"`（人源 taxonomy id 9606）
- 格式：XML（`<Hit_def>` 完整不截断）
- identity% = `Hsp_identity / Hsp_align-len * 100`（原始 identity 是匹配碱基数）
- 真脱靶 = 人源 RefSeq 转录本(NM_/XM_/NR_/XR_)命中 且 基因≠目标 且 非同源旁系同源
- 同源旁系同源 / 假基因 / 基因组层面命中 分层单列报告，不计入真脱靶

输出：03_结果/实验验证设计/primer_blast_specificity_v4.csv / .json
      03_结果/实验验证设计/provenance/blast_xml_raw_v4/<RID>.xml
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
RAW = f"{OUT}/provenance/blast_xml_raw_v4"
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
UA = {"User-Agent": "WorkBuddy-PrimerSpecificity/4.0"}
os.makedirs(RAW, exist_ok=True)

HUMAN_RE = re.compile(r"Homo sapiens", re.I)
RE_PAREN = re.compile(r"\(([A-Za-z0-9\-/\.@_]{1,20})\)\s*(?:,|\s|$)")
RE_PREFIX = re.compile(r"\b([A-Z][A-Z0-9\-]{1,14})\s*\(")
PSEUDO_RE = re.compile(r"pseudogene", re.I)

# 人源物种限定：txid9606 = Homo sapiens
HUMAN_TAXON = "txid9606[Organism:exp]"

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


def blast_put(seq, db="core_nt", entrez=None):
    p = {"CMD": "Put", "PROGRAM": "blastn", "DATABASE": db, "QUERY": seq,
         "MEGARUN": "on", "HITLIST_SIZE": "5000", "FILTER": "F",
         "EXPECT": "1000", "WORD_SIZE": "7", "NEIGHBOURING": "0"}
    if entrez:
        p["ENTREZ_QUERY"] = entrez
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


def wait_ready(rid, max_wait=1800):
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
    """解析 BLAST XML。identity% = Hsp_identity / Hsp_align-len。"""
    hits = []
    for hb in re.findall(r"<Hit>(.*?)</Hit>", xml, re.S):
        acc = (_blk(hb, "Hit_accession") or [""])[0]
        titles = _blk(hb, "Hit_def") or _blk(hb, "Title") or [""]
        title = titles[0].strip()
        ids = _blk(hb, "Hsp_identity")
        ident_n = int(ids[0]) if ids else 0
        alens = _blk(hb, "Hsp_align-len")
        align_len = int(alens[0]) if alens else 0
        ident = int(round(ident_n * 100.0 / align_len)) if align_len else 0
        # 关键：BLAST 允许局部比对（align-len < 引物长度），
        # 短于引物的部分匹配不构成可扩增的脱靶（例如 LDHA 前向引物在
        # 2,261 条人源命中中，1,081 条仅 13 bp 匹配，且多为基因组克隆记录）。
        # 故必须同时记录 align_len 与 identity_n，由判定层按「全长匹配」筛选。
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
        human = bool(HUMAN_RE.search(title))
        is_tx = bool(re.match(r"(NM_|XM_|NR_|XR_)", acc))
        hits.append({
            "acc": acc, "title": title, "gene": gene,
            "human": human, "is_refseq_tx": is_tx,
            "identity_pct": ident, "evalue": evalue,
            "align_len": align_len, "identity_n": ident_n,
            "pseudo": bool(PSEUDO_RE.search(title)),
        })
    return hits


def classify(gene, human_hits, primer_len):
    """把命中分层。**全长匹配**（align_len == 引物长度）才可能真正扩增；
    部分匹配（align_len < 引物长度）无法与引物 3' 端完全退火，判为噪声。
    真脱靶 = 全长匹配 且 人源 RefSeq 转录本 且 基因≠目标 且 非旁系同源 且 非假基因。
    """
    fam = PARALOG_FAMILIES.get(gene, {gene})
    tgt, para, pseudo, off = [], [], [], []
    for h in human_hits:
        full = h.get("align_len", 0) == primer_len
        g = h["gene"]
        if not full:
            continue
        if g and g == gene:
            tgt.append(h)
        elif h["pseudo"]:
            pseudo.append(h)
        elif g and g in fam:
            para.append(h)
        else:
            off.append(h)
    return tgt, para, pseudo, off


def run_one(seq, db="core_nt", entrez=None):
    rid = blast_put(seq, db, entrez)
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
                res = run_one(seq, "core_nt", HUMAN_TAXON)
                hits = res["hits"]
                rec[f"{lab}_blast_status"] = res["status"]
                rec[f"{lab}_blast_rid"] = res["rid"]
                rec[f"{lab}_n_hits_total"] = len(hits)
                human = [h for h in hits if h["human"]]
                rec[f"{lab}_n_hits_human"] = len(human)
                _plen = len(seq)
                rec[f"{lab}_primer_len"] = _plen
                rec[f"{lab}_n_hits_full_length"] = sum(
                    1 for h in human if h.get("align_len", 0) == _plen)
                rec[f"{lab}_n_hits_partial"] = len(human) - rec[
                    f"{lab}_n_hits_full_length"]
                tgt, para, pseudo, off = classify(gene, human, _plen)
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
                rec[f"{lab}_target_full_match_tx"] = sum(
                    1 for h in tgt if h.get("is_refseq_tx"))
                rec[f"{lab}_target_top_acc"] = tgt[0]["acc"] if tgt else ""
                rec[f"{lab}_human_hit_genes"] = ";".join(
                    sorted({h["gene"] for h in human if h["gene"]}))
                rec[f"{lab}_human_hits_detail"] = "; ".join(
                    f"{h['acc']}|{h['gene'] or '?'}|{h['identity_pct']}%"
                    for h in human[:30])
            except Exception as e:                        # noqa: BLE001
                rec[f"{lab}_error"] = str(e)[:250]
            time.sleep(4)

        off_tx_n = (rec.get("forward_n_offtarget_refseq_tx", 0)
                    + rec.get("reverse_n_offtarget_refseq_tx", 0))
        rec["offtarget_refseq_tx_total"] = off_tx_n
        rec["offtarget_genomic_total"] = (
            rec.get("forward_n_offtarget_genomic", 0)
            + rec.get("reverse_n_offtarget_genomic", 0))
        rec["paralog_total"] = (
            len([x for x in rec.get("forward_paralog_genes", "").split(";") if x])
            + len([x for x in rec.get("reverse_paralog_genes", "").split(";") if x]))
        # 判定：①人源 RefSeq 转录本层面无全长真脱靶；
        #       ②目标基因存在**全长匹配的人源 RefSeq 转录本**（否则引物不可用）；
        #       ③目标位点 identity ≥ 90%。
        rec["specific"] = bool(
            off_tx_n == 0
            and rec.get("forward_target_max_identity", 0) >= 90
            and rec.get("reverse_target_max_identity", 0) >= 90
            and rec.get("forward_target_full_match_tx", 0) >= 1
            and rec.get("reverse_target_full_match_tx", 0) >= 1)
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
    out.to_csv(f"{OUT}/primer_blast_specificity_v4.csv", index=False)
    with open(f"{OUT}/primer_blast_specificity_v4.json", "w", encoding="utf-8") as fh:
        json.dump({
            "date": "2026-10-25",
            "script": "02_脚本/verify_primers_blast_v4.py",
            "database": "core_nt (NCBI Core Nucleotide) + ENTREZ_QUERY='txid9606[Organism:exp]'",
            "database_note": (
                "v3 使用 core_nt 但未限定物种，命中列表被 HITLIST_SIZE=500 截断，"
                "进化上更近的非人源物种挤占前 500 名（实测 ACTB 反向引物："
                "Pseudopipra pipra 24 条、Neopelma chrysocephalum 22 条、"
                "Corapipo altera 22 条占主导，人源仅剩 7 条 → 误判为人源 0 命中）。"
                "v4 加物种限定，只查人源记录，从根本上排除截断"),
            "program": "blastn MEGARUN, expect=1000, hitlist=5000, word_size=7",
            "organism_filter": HUMAN_TAXON,
            "output_format": "XML（完整 Hit_def，不截断）",
            "criteria": (
                "★ 全长匹配优先：BLAST 允许局部比对，align_len < 引物长度的部分匹配"
                "无法与引物 3' 端完全退火，判为噪声并单列计数；"
                "真脱靶 = 全长匹配 且 人源 RefSeq 转录本(NM_/XM_/NR_/XR_) "
                "且 基因≠目标 且 非同源旁系同源 且 非假基因；"
                "基因组层面命中（NG_/AC/AL 等）单列；"
                "同源旁系同源与假基因单列报告；"
                "目标基因须存在全长匹配的人源 RefSeq 转录本且 identity ≥ 90%"),
            "partial_match_note": (
                "实测 LDHA 前向引物在 2,261 条人源命中中，"
                "1,081 条仅 13 bp 匹配、342 条 16 bp，多为基因组克隆记录；"
                "全长 20 bp 匹配仅 27 条，其中 100% identity 者全部为 LDHA 位点"
                "及其邻近增强子（NG_008185/NG_122038/NG_122039/AC084117），"
                "人源 RefSeq 转录本层面无 100% 全长非目标命中"),
            "results": rows,
        }, fh, ensure_ascii=False, indent=1, default=str)
    print(f"\nBLAST 特异性验证 v4（core_nt + 人源限定）：{int(out.specific.sum())}/{len(out)} "
          f"对通过（人源 RefSeq 转录本层面无真脱靶）")
    print("DONE — verify_primers_blast_v4.py")


if __name__ == "__main__":
    main()
