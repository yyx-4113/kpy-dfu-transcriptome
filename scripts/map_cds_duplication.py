# -*- coding: utf-8 -*-
"""
map_cds_duplication.py — 目标基因 CDS 的**全基因组重复图谱**

═══ 为什么需要这一步 ═══

BLAST v6 全库复核发现：LDHA 新引物命中 **LOC105375980**（9p22.3 的
lncRNA，8 外显子），其外显子序列与 LDHA 编码序列 95% 同源、3'端全匹配
→ 判为 P0 真实共扩增风险。

这不是「声明的旁系同源家族」能预见的：LOC105375980 与 LDHA **既不同源
家族名、也不是蛋白编码基因**（lncRNA），却是基因组上一段真实的重复拷贝。
它揭示了 v4 设计期的结构性盲区：
    **只在「声明家族」内比对，无法发现基因组上任何位置的第三方重复拷贝。**

═══ 本脚本的解法 ═══

把每个目标基因的**整条 CDS**（而不是单条引物）作为一次 BLAST 查询，
用一张图回答："CDS 的哪些位置在人类基因组别处有 ≥95% 同源、且两端
完美退火的拷贝？"

输出 = CDS 逐位置的「重复掩码」：
- 对 CDS 每一碱基，标记该位置是否被某个**非目标、非假基因**的人源
  RefSeq 转录本以 ≥95% identity（且 3'端锚定干净）覆盖
- 设计引物时**避开被掩码的位置** → 从源头消除该类脱靶

这比"设计完再 BLAST、发现 P0 再换引物"的迭代快得多，且判据一致。

═══ 参数与判据（与 verify_primers_blast_v6.py 对齐）═══

数据库：core_nt + ENTREZ_QUERY='txid9606[Organism:exp]'
程序：blastn MEGARUN, expect=1000, hitlist=5000, word_size=7
格式：XML

重复掩码判据（对命中 Hsp 的**查询侧**坐标，即 CDS 位置）：
  - 该 Hsp 的 identity ≥ 95%
  - 该 Hsp 覆盖 CDS 的位置 ≥ 20 nt（一个引物长度）
  - 命中对象：人源、RefSeq 转录本（NM_/XM_/NR_/XR_）、非目标基因、
    非假基因
  → 该 Hsp 覆盖的 CDS 区间标记为「高风险重复区」

输出：
  03_结果/实验验证设计/cds_duplication_map.json
  03_结果/实验验证设计/provenance/blast_xml_raw_cdsmap/*.xml
"""
import json

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import os
import re
import time
import urllib.request
import urllib.parse

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
PROV = f"{OUT}/provenance"
RAW = f"{PROV}/blast_xml_raw_cdsmap"
os.makedirs(RAW, exist_ok=True)
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
UA = {"User-Agent": "WorkBuddy-CDDuplicationMap/1.0"}
HUMAN_TAXON = "txid9606[Organism:exp]"
HUMAN_RE = re.compile(r"Homo sapiens", re.I)
PSEUDO_RE = re.compile(r"pseudogene", re.I)

# CDS 重复判据
REPEAT_MIN_IDENTITY = 95.0
REPEAT_MIN_ALEN = 20

GENES = ["LDHA", "FAHD2B", "ADH1B", "ALDH3A2", "DHX58", "GBP1",
         "ACTB", "GAPDH", "TBP", "HPRT1", "RPLP0", "SDHA"]


def _get(u, data=None, timeout=180):
    if data is None:
        return urllib.request.urlopen(
            urllib.request.Request(u, headers=UA), timeout=timeout).read().decode(
                errors="ignore")
    return urllib.request.urlopen(
        urllib.request.Request(u, data=data, headers=UA),
        timeout=timeout).read().decode(errors="ignore")


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


def blast_put(seq):
    p = {"CMD": "Put", "PROGRAM": "blastn", "DATABASE": "core_nt",
         "QUERY": seq, "MEGARUN": "on", "HITLIST_SIZE": "5000",
         "FILTER": "F", "EXPECT": "1000", "WORD_SIZE": "7",
         "NEIGHBOURING": "0", "ENTREZ_QUERY": HUMAN_TAXON}
    t = _get(BLAST, urllib.parse.urlencode(p).encode())
    rid = re.search(r"RID = (\S+)", t)
    if not rid:
        raise RuntimeError("BLAST Put 未返回 RID")
    return rid.group(1)


def blast_status(rid):
    u = BLAST + "?" + urllib.parse.urlencode(
        {"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid})
    t = _get(u)
    m = re.search(r"Status=([A-Z]+)", t)
    return m.group(1) if m else "UNKNOWN"


def blast_xml(rid):
    u = BLAST + "?" + urllib.parse.urlencode(
        {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML",
         "ALIGNMENTS": "5000", "DESCRIPTIONS": "5000"})
    t = _get(u, timeout=900)
    with open(f"{RAW}/core_nt_{rid}.xml", "w", encoding="utf-8") as fh:
        fh.write(t)
    return t


def parse_hsps(xml):
    out = []
    for hb in re.findall(r"<Hit>(.*?)</Hit>", xml, re.S):
        acc = (blk(hb, "Hit_accession") or [""])[0]
        title = (blk(hb, "Hit_def") or [""])[0].strip()
        for hp in re.findall(r"<Hsp>(.*?)</Hsp>", hb, re.S):
            g = lambda t: (blk(hp, t) or [""])[0]      # noqa: E731
            try:
                idn = int(g("Hsp_identity") or 0)
                al = int(g("Hsp_align-len") or 0)
                qfrom = int(g("Hsp_query-from") or 0)
                qto = int(g("Hsp_query-to") or 0)
                hfrom = int(g("Hsp_hit-from") or 0)
                hto = int(g("Hsp_hit-to") or 0)
            except ValueError:
                continue
            qs, hs = g("Hsp_qseq"), g("Hsp_hseq")
            if not qs or not hs:
                continue
            out.append({
                "acc": acc, "title": title, "idn": idn, "al": al,
                "qfrom": qfrom, "qto": qto, "hfrom": hfrom, "hto": hto,
                "qs": qs, "hs": hs,
                "identity_pct": round(idn * 100.0 / al, 2) if al else 0.0,
                "human": bool(HUMAN_RE.search(title)),
                "is_refseq_tx": bool(re.match(r"(NM_|XM_|NR_|XR_)", acc)),
                "pseudo": bool(PSEUDO_RE.search(title)),
                "gene": gene_of(title),
            })
    return out


def main():
    src = f"{OUT}/final_primers_v4.csv"
    import pandas as pd
    df = pd.read_csv(src)
    cds_src = json.load(open(f"{PROV}/final_design_cache.json", encoding="utf-8"))

    # 取每个基因的 canonical CDS（来自 v4 报告）
    V4 = json.load(open(f"{OUT}/final_primer_design_v4.json", encoding="utf-8"))

    queries = {}
    for gene in GENES:
        acc = V4["genes"][gene]["canonical"]
        s = cds_src.get(acc, {}).get("seq")
        gb = cds_src.get(f"__gb__{acc}", {}).get("seq")
        if not s or not gb:
            print(f"  {gene}: CDS 不可得")
            continue
        import re as _re
        seg = None
        for b in _re.split(r"\n {5}(?=\w)", gb):
            if b.startswith("CDS"):
                loc = _re.match(r"CDS\s+(\S+)", b).group(1).replace("<", ">")
                for ln in b.split("\n")[1:]:
                    st = ln.strip()
                    if not st:
                        continue
                    if st.startswith("/"):
                        break
                    loc += st
                seg = [(int(a), int(c)) for a, c in
                       _re.findall(r"(\d+)\.\.(\d+)", loc)]
                break
        if not seg:
            continue
        cds = "".join(s[a - 1:b] for a, b in seg)
        queries[gene] = {"accession": acc, "cds": cds,
                         "cds_len": len(cds),
                         "cds_coords": seg}
        print(f"  {gene}: {acc} CDS {len(cds)} nt")

    # ── 全部提交（并行），再统一轮询 ──
    rids = {}
    for gene, q in queries.items():
        try:
            rid = blast_put(q["cds"])
            rids[gene] = rid
            print(f"  submitted {gene}: RID={rid}", flush=True)
            time.sleep(3)
        except Exception as e:                        # noqa: BLE001
            print(f"  submit fail {gene}: {e}")
    print(f"\n已提交 {len(rids)} 个查询，等待完成…")

    t0 = time.time()
    done = set()
    while len(done) < len(rids) and time.time() - t0 < 5400:
        for gene, rid in list(rids.items()):
            if gene in done:
                continue
            try:
                st = blast_status(rid)
            except Exception:                          # noqa: BLE001
                continue
            if st == "READY":
                done.add(gene)
                print(f"  {gene} READY ({time.time()-t0:.0f}s)", flush=True)
            elif st == "FAILED":
                done.add(gene)
                print(f"  {gene} FAILED", flush=True)
        if len(done) < len(rids):
            time.sleep(20)

    # ── 解析并生成掩码 ──
    report = {"date": "2026-10-25",
              "script": "02_脚本/map_cds_duplication.py",
              "database": "core_nt + ENTREZ_QUERY='txid9606[Organism:exp]'",
              "criteria": {
                  "repeat_mask": f"非目标/非假基因人源 RefSeq 转录本，identity ≥ "
                                 f"{REPEAT_MIN_IDENTITY}%，Hsp 覆盖 CDS ≥ "
                                 f"{REPEAT_MIN_ALEN} nt",
              },
              "genes": {}}
    for gene, rid in rids.items():
        p = f"{RAW}/core_nt_{rid}.xml"
        if not os.path.exists(p):
            try:
                blast_xml(rid)
            except Exception as e:                    # noqa: BLE001
                print(f"  xml fail {gene}: {e}")
                continue
        xml = open(p, encoding="utf-8", errors="ignore").read()
        hsps = parse_hsps(xml)
        q = queries[gene]
        cds_len = q["cds_len"]
        # query 侧是 CDS（1-based），Hsp_query-from/to 是 CDS 坐标
        mask = bytearray(cds_len)
        detail = {}
        n_full = 0
        for h in hsps:
            if not h["human"] or not h["is_refseq_tx"]:
                continue
            if h["pseudo"]:
                continue
            if h["gene"] == gene:
                continue
            if h["identity_pct"] < REPEAT_MIN_IDENTITY:
                continue
            q1 = min(h["qfrom"], h["qto"])
            q2 = max(h["qfrom"], h["qto"])
            span = q2 - q1 + 1
            if span < REPEAT_MIN_ALEN:
                continue
            n_full += 1
            for i in range(q1 - 1, min(q2, cds_len)):
                mask[i] = 1
            key = h["gene"] or "?"
            detail.setdefault(key, {"n_hsp": 0, "max_identity_pct": 0.0,
                                     "accs": set()})
            detail[key]["n_hsp"] += 1
            detail[key]["max_identity_pct"] = max(
                detail[key]["max_identity_pct"], h["identity_pct"])
            detail[key]["accs"].add(h["acc"])
        # 连续掩码区段
        segs, st = [], None
        for i, v in enumerate(mask):
            if v and st is None:
                st = i
            elif not v and st is not None:
                segs.append([st + 1, i])
                st = None
        if st is not None:
            segs.append([st + 1, cds_len])
        # 最长安全连续区段
        safe, cur = [], None
        for i, v in enumerate(mask):
            if not v and cur is None:
                cur = i
            elif v and cur is not None:
                safe.append([cur + 1, i])
                cur = None
        if cur is not None:
            safe.append([cur + 1, cds_len])
        detail_out = {k: {"n_hsp": v["n_hsp"],
                          "max_identity_pct": v["max_identity_pct"],
                          "n_accessions": len(v["accs"]),
                          "example_accs": sorted(v["accs"])[:6]}
                      for k, v in detail.items()}
        print(f"\n  {gene}: 重复命中 Hsp {n_full}，掩码区段 {len(segs)}，"
              f"安全区段 {len(safe)}")
        for k, v in sorted(detail_out.items(), key=lambda x: -x[1]["n_hsp"])[:6]:
            print(f"     {k:14s} n_hsp={v['n_hsp']:3d} "
                  f"max_id={v['max_identity_pct']} n_acc={v['n_accessions']}")
        report["genes"][gene] = {
            "accession": q["accession"], "cds_len": cds_len,
            "cds_coords": q["cds_coords"],
            "rid": rid,
            "n_repeat_hsp": n_full,
            "repeat_masked_segments": segs,
            "repeat_masked_nt": sum(b - a + 1 for a, b in segs),
            "safe_segments": safe,
            "longest_safe_segment": max((b - a + 1 for a, b in safe),
                                        default=0),
            "offtarget_genes": detail_out,
        }
    with open(f"{OUT}/cds_duplication_map.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1, default=str)
    print("\nDONE — map_cds_duplication.py")


if __name__ == "__main__":
    main()
