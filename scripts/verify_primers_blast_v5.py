# -*- coding: utf-8 -*-
"""
verify_primers_blast_v5.py — 新引物集的 BLAST 全库复核 + 严格脱靶判读
（一次完成，取代 v1–v4 四轮）

═══ 为什么做 v5 ═══

v1–v4 走了四轮，每轮都在**判据**上打补丁，而非在**流程**上补全：
  v1  Text 格式 description 截断 → 基因符号提取失败
  v2  refseq_rna 库未完整收录人源转录本 → LDHA 人源 0 命中
  v3  core_nt 命中列表被 HITLIST_SIZE=500 截断，非人源挤占 → ACTB 人源 0 命中
  v4  加全长匹配筛选（align_len == 引物长度）+ 物种限定

v4 已解决"截断"与"部分匹配"两个工具层缺陷，但**判据本身仍然粗糙**：
把「全长匹配的人源非目标 RefSeq 转录本」一律计为真脱靶。实测：

- ALDH3A2 反向引物：4,713 条全长人源命中，散布 60+ 无关基因座
  （DYM 48、A1CF 20、LRCH3 16…），identity 全为 95%，其中绝大多数
  来自 **3'UTR 的 Alu 重复元件** → 属重复序列噪声，非"会扩增的脱靶"；
- ADH1B 前向：ADH1C 100% 全长匹配（4 条）→ **真实脱靶**；
- LDHA 反向：LYST 95%、3' 端全匹配（16 条）→ **真实脱靶**。

v5 把三件事合并为一次执行，并引入 **3' 端锚定**这一 PCR 决定性判据：
PCR 延伸依赖引物 3' 端的完美退火。3' 端末 5 nt 存在错配者无法有效起始
延伸，即使全长 identity 95% 也不会共扩增。

═══ v5 判据 ═══

数据库：core_nt + ENTREZ_QUERY='txid9606[Organism:exp]'（人源限定）
格式：XML（Hit_def 不截断）；HITLIST_SIZE=5000；EXPECT=1000；FILTER=F

命中分层（对**全长**匹配 align_len == 引物长度的人源命中）：
  TARGET   基因 == 目标基因
  PARALOG  严格旁系同源白名单内的同家族成员（同基因座/经典可接受）
  PSEUDO   假基因
  P0       真实脱靶（阻断级）：非目标 + 非旁系同源 + 非假基因 +
           人源 RefSeq 转录本 + **3' 端 5 nt 全匹配** + identity ≥ 95%
  P1       风险脱靶：同上但 90% ≤ identity < 95%
  P2       噪声：其余（3' 端有错配，或 identity < 90%）
  PARTIAL  align_len < 引物长度 → 无法完全退火，判为噪声

最终 usable 判据：
  ① 目标基因存在**全长匹配的人源 RefSeq 转录本**（否则引物不可用）
  ② P0 == 0
  ③ 目标位点 identity ≥ 90%

输出：
  03_结果/实验验证设计/primer_blast_v5.csv / .json
  03_结果/实验验证设计/provenance/blast_xml_raw_v5/<RID>.xml
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
RAW = f"{OUT}/provenance/blast_xml_raw_v5"
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
UA = {"User-Agent": "WorkBuddy-PrimerSpecificityV5/1.0"}
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
                "acc": acc, "title": title,
                "idn": idn, "al": al,
                "identity_pct": round(idn * 100.0 / al, 2) if al else 0.0,
                "qseq": qs, "hseq": hs,
                "human": bool(HUMAN_RE.search(title)),
                "is_refseq_tx": bool(re.match(r"(NM_|XM_|NR_|XR_)", acc)),
                "pseudo": bool(PSEUDO_RE.search(title)),
                "gene": gene_of(title),
            })
    return hits


def classify(hit, gene, primer_len):
    """返回 (level, mm_total, mm_3p5)"""
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


def run_one(seq):
    rid = blast_put(seq)
    st = wait_ready(rid)
    if st != "READY":
        return {"status": st, "hits": [], "rid": rid}
    xml = blast_xml(rid)
    return {"status": st, "hits": parse_hits(xml), "rid": rid}


def main():
    src = f"{OUT}/primer_redesign_v2_selected.csv"
    if not os.path.exists(src):
        raise FileNotFoundError(f"缺少重设计结果：{src}")
    df = pd.read_csv(src)
    rows, detail = [], []
    for _, r in df.iterrows():
        gene = r["gene"]
        rec = {"gene": gene, "role": r["role"], "module": r["module"],
               "refseq_accession": r["refseq_accession"],
               "forward_primer_5to3": r["forward_primer_5to3"],
               "reverse_primer_5to3": r["reverse_primer_5to3"]}
        for lab, seq in [("forward", r["forward_primer_5to3"]),
                         ("reverse", r["reverse_primer_5to3"])]:
            try:
                res = run_one(seq)
                hits = res["hits"]
                plen = len(seq)
                rec[f"{lab}_blast_status"] = res["status"]
                rec[f"{lab}_blast_rid"] = res["rid"]
                rec[f"{lab}_primer_len"] = plen
                rec[f"{lab}_n_hits_total"] = len(hits)
                human = [h for h in hits if h["human"]]
                rec[f"{lab}_n_hits_human"] = len(human)
                counts = {}
                off_gene = {}
                tgt_tx = 0
                tgt_id = 0.0
                for h in human:
                    lvl, mm_t, mm_35 = classify(h, gene, plen)
                    h["level"] = lvl
                    counts[lvl] = counts.get(lvl, 0) + 1
                    if lvl == "TARGET":
                        if h["is_refseq_tx"]:
                            tgt_tx += 1
                        tgt_id = max(tgt_id, h["identity_pct"])
                    if lvl in ("P0", "P1"):
                        off_gene.setdefault(lvl, {}).setdefault(
                            h["gene"] or "?", 0)
                        off_gene[lvl][h["gene"] or "?"] += 1
                        detail.append({
                            "gene": gene, "primer": lab, "query": seq,
                            "accession": h["acc"], "hit_gene": h["gene"],
                            "identity_n": h["idn"], "align_len": h["al"],
                            "identity_pct": h["identity_pct"],
                            "mm_total": mm_t, "mm_3p5": mm_35,
                            "level": lvl, "title": h["title"][:160]})
                rec[f"{lab}_n_full_length_human"] = sum(
                    v for k, v in counts.items() if k != "PARTIAL")
                rec[f"{lab}_n_partial_human"] = counts.get("PARTIAL", 0)
                rec[f"{lab}_n_P0"] = counts.get("P0", 0)
                rec[f"{lab}_n_P1"] = counts.get("P1", 0)
                rec[f"{lab}_n_PARALOG"] = counts.get("PARALOG", 0)
                rec[f"{lab}_n_PSEUDO"] = counts.get("PSEUDO", 0)
                rec[f"{lab}_P0_genes"] = ";".join(
                    f"{k}:{v}" for k, v in sorted(off_gene.get("P0", {}).items()))
                rec[f"{lab}_P1_genes"] = ";".join(
                    f"{k}:{v}" for k, v in sorted(off_gene.get("P1", {}).items()))
                rec[f"{lab}_target_full_match_refseq_tx"] = tgt_tx
                rec[f"{lab}_target_max_identity_pct"] = tgt_id
            except Exception as e:                        # noqa: BLE001
                rec[f"{lab}_error"] = str(e)[:250]
            time.sleep(4)

        p0 = rec.get("forward_n_P0", 0) + rec.get("reverse_n_P0", 0)
        p1 = rec.get("forward_n_P1", 0) + rec.get("reverse_n_P1", 0)
        ft = rec.get("forward_target_full_match_refseq_tx", 0)
        rt = rec.get("reverse_target_full_match_refseq_tx", 0)
        fi = rec.get("forward_target_max_identity_pct", 0)
        ri = rec.get("reverse_target_max_identity_pct", 0)
        rec["n_P0_total"] = p0
        rec["n_P1_total"] = p1
        rec["specific"] = bool(p0 == 0 and ft >= 1 and rt >= 1
                               and fi >= 90 and ri >= 90)
        rec["verdict"] = ("USABLE" if rec["specific"]
                          else ("NO_TARGET" if (ft < 1 or rt < 1)
                                else "OFFTARGET_P0"))
        rows.append(rec)
        print(f"{gene:9s} F[full={rec.get('forward_n_full_length_human','?')} "
              f"P0={rec.get('forward_n_P0','?')} P1={rec.get('forward_n_P1','?')} "
              f"tgtTx={ft}] R[full={rec.get('reverse_n_full_length_human','?')} "
              f"P0={rec.get('reverse_n_P0','?')} P1={rec.get('reverse_n_P1','?')} "
              f"tgtTx={rt}] -> {rec['verdict']}", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(f"{OUT}/primer_blast_v5.csv", index=False)
    with open(f"{OUT}/primer_blast_v5.json", "w", encoding="utf-8") as fh:
        json.dump({
            "date": "2026-10-25",
            "script": "02_脚本/verify_primers_blast_v5.py",
            "database": "core_nt + ENTREZ_QUERY='txid9606[Organism:exp]'",
            "program": "blastn MEGARUN, expect=1000, hitlist=5000, word_size=7",
            "output_format": "XML（完整 Hit_def）",
            "criteria": {
                "TARGET": "基因 == 目标基因",
                "PARALOG": "严格旁系同源白名单内同家族成员",
                "PSEUDO": "假基因",
                "P0": "真实脱靶（阻断级）：全长匹配 + 人源 RefSeq 转录本 + 非目标"
                      " + 非旁系同源 + 非假基因 + 3'端末5nt全匹配 + identity ≥ 95%",
                "P1": "风险脱靶：同上但 90% ≤ identity < 95%",
                "P2": "噪声：3'端有错配或 identity < 90%",
                "PARTIAL": "align_len < 引物长度（局部比对，无法完全退火）",
                "usable": "目标存在全长匹配的人源 RefSeq 转录本 且 P0 == 0 "
                          "且目标位点 identity ≥ 90%",
            },
            "method_note": ("v5 引入 3' 端锚定判据：PCR 延伸依赖引物 3' 端"
                            "完美退火，3' 端末 5 nt 有错配者无法有效起始延伸，"
                            "故不计入真实脱靶（v4 曾将其一律误判为脱靶，"
                            "导致 ALDH3A2 因 3'UTR Alu 重复元件产生 4,713 条"
                            "假阳性命中）。"),
            "results": rows,
        }, fh, ensure_ascii=False, indent=1, default=str)
    if detail:
        pd.DataFrame(detail).to_csv(
            f"{OUT}/blast_offtarget_adjudication_v5.csv", index=False)
    n_ok = int(out.specific.sum())
    print(f"\nBLAST v5 复核：{n_ok}/{len(out)} 对可用")
    print("DONE — verify_primers_blast_v5.py")


if __name__ == "__main__":
    main()
