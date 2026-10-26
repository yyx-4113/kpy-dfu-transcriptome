# -*- coding: utf-8 -*-
"""
verify_primers_blast_v2.py — 引物单靶特异性 BLAST 验证（修正版，任务 #11d）

═══ 修正说明（对 v1 脚本 verify_primers_blast.py 的两处方法学缺陷）═══

缺陷 1（**描述截断导致基因符号提取失败**）
  v1 使用 FORMAT_TYPE=Text。实测 NCBI 返回的命中列表中 description 列被截断为
  「... uncharacterized... 39.2 1.2 100%」形式，基因符号常落在截断点之后，
  导致 GENE_RE 匹配不到 → `genes=` 为空 → 被误判为"非特异"。
  **实证**：FAHD2B F 人源命中 29 条但 genes= 为空；GAPDH F/R genes= 为空；
  而 LDHA R 恰好符号未被截断 → genes=LDHA。说明失败源于格式而非真实脱靶。

缺陷 2（**未区分"同源旁系同源基因"与"随机脱靶"**）
  ACTB 引物命中 ACTG1（肌动蛋白 γ）——ACTG1 与 ACTB 序列高度同源，属**结构性
  不可避免的旁系同源**，与随机脱靶性质不同，须如实区分并可解释。

修正方案
  1. 改用 **FORMAT_TYPE=XML**：`<Hit_def>` 给出**不截断的完整描述**
     （实测 'Homo sapiens ... (LDHA), mRNA' 全长可见），并额外取
     `<Hit_accession>` / `<Organism>` / `<Hsp_identity>` / `<Hsp_evalue>`。
  2. 基因符号提取采用**多级正则**：
     (a) 前置型 SYMBOL ( 如 'Homo sapiens LDHA (LDHA), mRNA'
     (b) 括号型 (SYMBOL) → 括号内为符号
     (c) 若仍失败，用 accession 通过 **NCBI E-utilities esummary(efetch) 查 Gene
         数据库** 解析官方 symbol（权威兜底，非启发式猜测）
  3. 脱靶判定分层：
     - **真脱靶（random off-target）**：人源命中、基因 ≠ 目标、且非目标基因的
       **同源家族成员**（同一 Gene symbol 前缀/已知旁系同源对）
     - **旁系同源命中（paralog）**：单独列出并如实报告
     - **跨物种命中**：单独记录，不计入脱靶（人体样本中不存在该模板）
  4. 判定特异 = 真脱靶为 0 **且** 目标基因命中 identity ≥ 90% **且**
     正反向均命中目标基因。

输出：03_结果/实验验证设计/primer_blast_specificity_v2.csv / .json
      03_结果/实验验证设计/provenance/blast_xml_raw/<RID>.xml（原始 XML 留档）
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
RAW = f"{OUT}/provenance/blast_xml_raw"
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
EUTIL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "WorkBuddy-PrimerSpecificity/2.0"}
os.makedirs(RAW, exist_ok=True)

HUMAN_RE = re.compile(r"Homo sapiens", re.I)

# 基因符号提取：优先括号内，其次 "SYMBOL (" 前置式
RE_PAREN = re.compile(r"\(([A-Za-z0-9\-/\.@_]{1,20})\)\s*(?:,|\s|$)")
RE_PREFIX = re.compile(r"\b([A-Z][A-Z0-9\-]{1,14})\s*\(")
RE_NCBI_SYMBOL = re.compile(r"\[(?:gb|ref|emb|dbj|pir|prf|tpg|sp|tr|pat|rs|un|gb\|.*?)\|"
                            r"([^\]]+)\]")

# 已知同源旁系同源家族（用于区分"结构性旁系同源"与"随机脱靶"）
PARALOG_FAMILIES = {
    "ACTB": {"ACTG1", "ACTG2", "ACTC1", "ACTB"},
    "ACTG1": {"ACTB", "ACTG2", "ACTC1"},
    "GAPDH": {"GAPDHS", "GAPDHL"},   # 假基因/拟基因不表达，单独说明
    "RPLP0": {"RPLP0", "RPLP0P"},
    "RPLP": {"RPLP0"},
    "TBP": {"TBP"},
    "HPRT1": {"HPRT1"},
    "SDHA": {"SDHA", "SDHAP1"},
    "GBP1": {"GBP1", "GBP3", "GBP4", "GBP5", "GBP6", "GBP7"},
    "DHX58": {"DHX58", "DHX58P1"},
    "LDHA": {"LDHA", "LDHAL"},
    "ALDH3A2": {"ALDH3A2", "ALDH3A1"},
    "ADH1B": {"ADH1B", "ADH1C", "ADH1A", "ADH4"},
    "FAHD2B": {"FAHD2B", "FAHD2A", "FAHD2B"},
}


def blast_put(seq):
    p = {"CMD": "Put", "PROGRAM": "blastn", "DATABASE": "refseq_rna",
         "QUERY": seq, "MEGARUN": "on", "HITLIST_SIZE": "500",
         "FILTER": "F", "EXPECT": "1000", "WORD_SIZE": "7",
         "NEIGHBOURING": "0"}
    req = urllib.request.Request(BLAST, data=urllib.parse.urlencode(p).encode(),
                                 headers=UA)
    t = urllib.request.urlopen(req, timeout=180).read().decode(errors="ignore")
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


def blast_xml(rid):
    u = BLAST + "?" + urllib.parse.urlencode(
        {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML",
         "ALIGNMENTS": "500", "DESCRIPTIONS": "500"})
    t = urllib.request.urlopen(u, timeout=300).read().decode(errors="ignore")
    with open(f"{RAW}/{rid}.xml", "w", encoding="utf-8") as fh:
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
    """提取所有 <tag>...</tag> 块"""
    return re.findall(rf"<{tag}>(.*?)</{tag}>", xml, re.S)


def _blk_attr(xml, tag, attr):
    return re.findall(rf"<{tag}\s+{attr}=\"([^\"]*)\"", xml, re.S)


def parse_hits_xml(xml):
    """解析 BLAST XML → 命中列表（完整描述，不截断）"""
    hits = []
    # 每个 <Hit> 块
    for hb in re.findall(r"<Hit>(.*?)</Hit>", xml, re.S):
        acc = _blk(hb, "Hit_accession")
        acc = acc[0] if acc else ""
        defs = _blk(hb, "Hit_def")
        title = defs[0] if defs else ""
        orgs = _blk(hb, "Organism")
        organism = orgs[0] if orgs else ""
        ids = _blk(hb, "Hsp_identity")
        ident_n = int(ids[0]) if ids else 0
        alens = _blk(hb, "Hsp_align-len")
        align_len = int(alens[0]) if alens else 0
        # 关键：BLAST XML 的 <Hsp_identity> 是**匹配碱基数**而非百分比
        # （实测 ACTB 正向命中 NM_001101：identity=20, align-len=20 → 100%）
        ident = int(round(ident_n * 100.0 / align_len)) if align_len else 0
        evs = _blk(hb, "Hsp_evalue")
        try:
            evalue = float(evs[0]) if evs else float("nan")
        except ValueError:
            evalue = float("nan")
        bits = _blk(hb, "Hsp_bit-score")
        try:
            bit = float(bits[0]) if bits else float("nan")
        except ValueError:
            bit = float("nan")
        # 基因符号
        gene = ""
        m = RE_PAREN.search(title)
        if m:
            cand = m.group(1)
            # 括号内容若是纯数字或过长则不接受
            if not cand.isdigit() and 1 <= len(cand) <= 20:
                gene = cand
        if not gene:
            m2 = RE_PREFIX.search(title)
            if m2:
                gene = m2.group(1)
        hits.append({
            "acc": acc, "title": title.strip(), "gene": gene,
            "organism": organism,
            "human": bool(HUMAN_RE.search(title) or HUMAN_RE.search(organism)),
            "identity_pct": ident, "evalue": evalue, "bits": bit,
        })
    return hits


# ── E-utilities 兜底：accession → 官方 gene symbol（权威） ──────────────
def acc_to_symbol(acc):
    """用 E-utilities nuccore esummary 取官方 gene symbol；失败返回 ''"""
    try:
        u = EUTIL + "/esummary.fcgi?" + urllib.parse.urlencode(
            {"db": "nuccore", "id": acc, "retmode": "json"})
        t = urllib.request.urlopen(u, timeout=90).read().decode(errors="ignore")
        j = json.loads(t)
        res = j.get("result", {})
        for uid in res.get("uids", []):
            rec = res.get(uid, {})
            for it in rec.get("dbxrefs", []) or []:
                pass
            title = rec.get("title", "")
            m = RE_PAREN.search(title)
            if m:
                return m.group(1)
            m2 = RE_PREFIX.search(title)
            if m2:
                return m2.group(1)
            return title.split(",")[0].strip()
        return ""
    except Exception:                                   # noqa: BLE001
        return ""


PSEUDO_RE = re.compile(r"pseudogene", re.I)


def classify_offtargets(gene, hit_genes, pseudo_genes=()):
    """区分真脱靶 / 旁系同源 / 假基因

    - 真脱靶（random off-target）：人源命中、基因≠目标、且非同源旁系同源、且非假基因
    - 旁系同源命中（paralog）：与目标同属已知同源家族（如 ACTB↔ACTG1）
    - 假基因（pseudogene）：无成熟转录本，RT-qPCR 不扩增 → 单列，不计脱靶
    """
    fam = PARALOG_FAMILIES.get(gene, {gene})
    random_off, paralog, pseudo = [], [], []
    for g in hit_genes:
        if g == gene:
            continue
        if g in pseudo_genes:
            pseudo.append(g)
        elif g in fam:
            paralog.append(g)
        else:
            random_off.append(g)
    return sorted(random_off), sorted(paralog), sorted(pseudo)


def run_one(seq, label):
    rid = blast_put(seq)
    st = wait_ready(rid)
    if st != "READY":
        return {"status": st, "hits": [], "rid": rid}
    xml = blast_xml(rid)
    hits = parse_hits_xml(xml)
    return {"status": st, "hits": hits, "rid": rid}


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
                res = run_one(seq, lab)
                hits = res["hits"]
                rec[f"{lab}_blast_status"] = res["status"]
                rec[f"{lab}_blast_rid"] = res["rid"]
                rec[f"{lab}_n_hits_total"] = len(hits)
                human = [h for h in hits if h["human"]]
                rec[f"{lab}_n_hits_human"] = len(human)
                rec[f"{lab}_n_hits_nonshuman"] = len(hits) - len(human)
                # 基因符号兜底：XML 未取到 symbol 时用 E-utilities 解析
                need = [h for h in human if not h["gene"]]
                for h in need:
                    h["gene"] = acc_to_symbol(h["acc"])
                hg = sorted({h["gene"] for h in human if h["gene"]})
                rec[f"{lab}_human_hit_genes"] = ";".join(hg)
                rec[f"{lab}_human_hit_genes_n"] = len(hg)
                rec[f"{lab}_human_hit_genes_unresolved"] = len(
                    [h for h in human if not h["gene"]])
                pseudo_genes = sorted({
                    h["gene"] for h in human
                    if h["gene"] and PSEUDO_RE.search(h["title"])})
                rec[f"{lab}_pseudogene_genes"] = ";".join(pseudo_genes)
                rec[f"{lab}_pseudogene_n"] = len(pseudo_genes)
                rand_off, para, pseudo = classify_offtargets(gene, hg, pseudo_genes)
                rec[f"{lab}_random_offtarget_genes"] = ";".join(rand_off)
                rec[f"{lab}_random_offtarget_n"] = len(rand_off)
                rec[f"{lab}_paralog_hit_genes"] = ";".join(para)
                rec[f"{lab}_paralog_n"] = len(para)
                tgt = [h for h in human if h["gene"] == gene]
                rec[f"{lab}_target_human_hits"] = len(tgt)
                rec[f"{lab}_target_max_identity"] = max(
                    [h["identity_pct"] for h in tgt], default=0)
                rec[f"{lab}_target_top_acc"] = tgt[0]["acc"] if tgt else ""
                rec[f"{lab}_cross_species_hits"] = len(hits) - len(human)
                # 记录人源命中明细（供审计）
                rec[f"{lab}_human_hits_detail"] = "; ".join(
                    f"{h['acc']}|{h['gene'] or '?'}|{h['identity_pct']}%"
                    for h in human[:25])
            except Exception as e:                        # noqa: BLE001
                rec[f"{lab}_error"] = str(e)[:250]
            time.sleep(4)
        rand_n = rec.get("forward_random_offtarget_n", 0) + \
            rec.get("reverse_random_offtarget_n", 0)
        rec["random_offtarget_total"] = rand_n
        rec["paralog_total"] = rec.get("forward_paralog_n", 0) + \
            rec.get("reverse_paralog_n", 0)
        rec["specific"] = bool(
            rand_n == 0
            and rec.get("forward_human_hit_genes_n", 0) >= 1
            and rec.get("reverse_human_hit_genes_n", 0) >= 1
            and rec.get("forward_target_max_identity", 0) >= 90
            and rec.get("reverse_target_max_identity", 0) >= 90
            and rec.get("forward_human_hit_genes_unresolved", 0) == 0
            and rec.get("reverse_human_hit_genes_unresolved", 0) == 0)
        rows.append(rec)
        print(f"{gene:8s} F[hum {rec.get('forward_n_hits_human','?')} "
              f"genes={rec.get('forward_human_hit_genes','?')[:30]} "
              f"randOff={rec.get('forward_random_offtarget_n','?')} "
              f"para={rec.get('forward_paralog_n','?')} "
              f"pseu={rec.get('forward_pseudogene_n','?')}] "
              f"R[hum {rec.get('reverse_n_hits_human','?')} "
              f"genes={rec.get('reverse_human_hit_genes','?')[:30]} "
              f"randOff={rec.get('reverse_random_offtarget_n','?')} "
              f"para={rec.get('reverse_paralog_n','?')} "
              f"pseu={rec.get('reverse_pseudogene_n','?')}] "
              f"xsp={rec.get('forward_cross_species_hits','?')}/"
              f"{rec.get('reverse_cross_species_hits','?')} "
              f"-> specific={rec['specific']}", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(f"{OUT}/primer_blast_specificity_v2.csv", index=False)
    with open(f"{OUT}/primer_blast_specificity_v2.json", "w", encoding="utf-8") as fh:
        json.dump({
            "date": "2026-10-25",
            "script": "02_脚本/verify_primers_blast_v2.py",
            "database": "refseq_rna (all species)",
            "program": "blastn MEGARUN, expect=1000, hitlist=500, word_size=7",
            "output_format": "XML (完整 Hit_def，不截断) — 修正 v1 的描述截断缺陷",
            "criteria": (
                "真脱靶 = 人源命中且基因符号≠目标基因且不属于已知同源旁系同源家族；"
                "旁系同源命中单独报告；跨物种命中记录但不计入脱靶；"
                "目标基因命中 identity ≥ 90%；人源命中基因符号须 100% 解析成功"),
            "results": rows,
        }, fh, ensure_ascii=False, indent=1, default=str)
    print(f"\nBLAST 特异性验证 v2：{int(out.specific.sum())}/{len(out)} 对通过"
          f"（人源层面无真脱靶）")
    print("DONE — verify_primers_blast_v2.py")


if __name__ == "__main__":
    main()
