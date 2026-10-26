# -*- coding: utf-8 -*-
"""
audit_refseq_version_drift.py — 全部 12 个引物基因的 RefSeq 版本漂移与
转录本覆盖核查（2026-10-25）

═══ 触发事件 ═══

BLAST v4 独立复核中发现 LDHA 前向引物在 **全部 6 条当前 RefSeq 转录本**中
匹配位置均为 -1 → 该引物对无法扩增 LDHA cDNA。根因：设计时所用模板
`NM_001135239.1` 的 5'UTR 在 RefSeq 版本更新中被重排（`.1`→`.2`
位置一致性仅 0.2524）。

单点修复 LDHA 属于"只审包装"式处理——必须做**全量同类排查**：
其余 11 个基因同样存在"设计版本 ≠ 当前版本"的可能（本次已实测：
ADH1B .5→.6、ACTB .3→.5、GAPDH .6→.7 均发生版本更新），
且每个基因的真实 RefSeq 转录本数 > 1（单引物对未必覆盖全部变体）。

═══ 核查项 ═══

对每个基因（6 靶基因 + 6 内参候选）：
1. **当前版本 vs 设计版本**：逐条取当前 RefSeq（efetch 不带版本号 → 取最新），
   记录版本号、长度、序列一致性；
2. **设计引物在当前版本中的落位**：正向引物须在正义链出现，
   反向引物须以反向互补出现；
3. **转录本覆盖**：枚举该基因全部 RefSeq 转录本（gene→nuccore elink
   `gene_nuccore_refseqrna` + esummary），逐条核查引物对落位，
   报告覆盖数 / 总数；
4. **CDS 落位**：引物须位于 CDS 内（RT-qPCR 测的是 mRNA 丰度，
   5'UTR 变体高度变异，UTR 内引物易漏检变体）；
5. **分类**：NONCODING（无 CDS）转录本单列，不计入覆盖分母但如实报告。

输出：
  03_结果/实验验证设计/refseq_version_drift_audit.csv   逐基因汇总
  03_结果/实验验证设计/refseq_version_drift_audit.json   逐转录本明细
  03_结果/实验验证设计/provenance/refseq_audit_fetch.log 抓取留档
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import json
import re
import time
import urllib.request

import pandas as pd

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
PROV = f"{OUT}/provenance"
os.makedirs(PROV, exist_ok=True)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "WorkBuddy-RefSeqAudit/1.0"}
RC_MAP = str.maketrans("ACGT", "TGCA")

# gene_id 与 design_rtqpcr_primers.py 一致
GENES = {
    "LDHA": (3939, "target"), "FAHD2B": (151313, "target"),
    "ADH1B": (125, "target"), "DHX58": (79132, "target"),
    "ALDH3A2": (224, "target"), "GBP1": (2633, "target"),
    "ACTB": (60, "reference_candidate"), "GAPDH": (2597, "reference_candidate"),
    "TBP": (6908, "reference_candidate"), "HPRT1": (3251, "reference_candidate"),
    "RPLP0": (6175, "reference_candidate"), "SDHA": (6389, "reference_candidate"),
}
LOG = []


def rc(s):
    return s.translate(RC_MAP)[::-1]


def _get(url, tries=4):
    last = None
    for k in range(tries):
        try:
            t = urllib.request.urlopen(
                urllib.request.Request(url, headers=UA), timeout=120).read().decode()
            LOG.append(f"GET {url}\n  -> {len(t)} chars")
            return t
        except Exception as e:                        # noqa: BLE001
            last = e
            time.sleep(2 + 3 * k)
    raise RuntimeError(f"efetch failed: {last}")


def efetch_fasta(acc):
    raw = _get(f"{EUTILS}/efetch.fcgi?db=nuccore&id={acc}"
               f"&rettype=fasta&retmode=text")
    lines = raw.strip().split("\n")
    header = lines[0].lstrip(">")
    seq = "".join(l.strip() for l in lines[1:]
                   if re.fullmatch(r"[A-Za-z]+", l.strip())
                   and not l.startswith("BASE"))
    return header, seq.upper()


def efetch_gb(acc):
    return _get(f"{EUTILS}/efetch.fcgi?db=nuccore&id={acc}"
                f"&rettype=gb&retmode=text")


def parse_cds(gb):
    m = re.search(r"^\s{5}CDS\s+(\S+)", gb, re.M)
    if not m:
        return []
    return [(int(a), int(b))
            for a, b in re.findall(r"(\d+)\.\.(\d+)", m.group(1).replace("<", ">"))]


def list_transcripts(gid):
    d = json.loads(_get(f"{EUTILS}/elink.fcgi?dbfrom=gene&db=nuccore&id={gid}"
                        f"&linkname=gene_nuccore_refseqrna&retmode=json"))
    ids = [l for ls in d.get("linksets", [])
           for dbk in ls.get("linksetdbs", []) for l in dbk.get("links", [])]
    if not ids:
        return []
    s = json.loads(_get(f"{EUTILS}/esummary.fcgi?db=nuccore&id={','.join(ids)}"
                        f"&retmode=json"))
    out = []
    for u in s["result"]["uids"]:
        v = s["result"][u]
        out.append((v.get("accessionversion"), v.get("title", ""),
                    v.get("slen")))
    return sorted(out)


def in_cds(pos0, cds):
    if pos0 < 0:
        return False
    return any(a - 1 <= pos0 < b for a, b in cds)


def main():
    primers = pd.read_csv(f"{OUT}/qpcr_primers.csv")
    pmap = {r["gene"]: r for _, r in primers.iterrows()}
    summary, detail = [], {}

    for gene, (gid, role) in GENES.items():
        r = pmap[gene]
        F, R = r["forward_primer_5to3"], r["reverse_primer_5to3"]
        design_acc = r["refseq_accession"]
        design_base = design_acc.split(".")[0]

        # --- 设计版本序列（用于与当前版本比对）---
        try:
            _, s_design = efetch_fasta(design_acc)
        except Exception as e:                        # noqa: BLE001
            s_design = ""
            LOG.append(f"  design version fetch failed: {e}")
        time.sleep(0.35)

        tx_list = list_transcripts(gid)
        time.sleep(0.35)

        recs = []
        for acc, title, slen in tx_list:
            try:
                header, seq = efetch_fasta(acc)
                time.sleep(0.35)
                gb = efetch_gb(acc)
                time.sleep(0.35)
                cds = parse_cds(gb)
                iF = seq.find(F)
                iRrc = seq.find(rc(R))
                recs.append({
                    "acc": acc, "title": title, "len": len(seq),
                    "cds": cds, "is_coding": bool(cds),
                    "F_pos": iF, "R_rc_pos": iRrc,
                    "F_in_CDS": in_cds(iF, cds),
                    "R_in_CDS": in_cds(iRrc, cds),
                    "primers_present": (iF >= 0 and iRrc >= 0),
                })
            except Exception as e:                    # noqa: BLE001
                recs.append({"acc": acc, "title": title, "error": str(e)[:120]})
        detail[gene] = recs

        coding = [x for x in recs if x.get("is_coding")]
        noncoding = [x for x in recs if not x.get("is_coding")]
        covered = [x for x in coding if x.get("primers_present")]
        cur_ver = next((x["acc"] for x in recs if x["acc"].startswith(design_base)
                        and x["acc"].split(".")[-1] != design_acc.split(".")[-1]),
                       "")
        # 当前版本 = 枚举结果中该 base 的最高版本号
        cands = [x["acc"] for x in recs if x["acc"].split(".")[0] == design_base]
        cur_ver = max(cands, key=lambda a: int(a.split(".")[1])) if cands else ""

        seq_ident = None
        if s_design and cur_ver:
            try:
                _, s_cur = efetch_fasta(cur_ver)
                time.sleep(0.35)
                n = min(len(s_design), len(s_cur))
                seq_ident = round(sum(1 for i in range(n)
                                      if s_design[i] == s_cur[i]) / n, 4)
            except Exception:                         # noqa: BLE001
                pass

        row = {
            "gene": gene, "role": role, "gene_id": gid,
            "design_accession": design_acc, "current_accession": cur_ver,
            "version_changed": (cur_ver != design_acc),
            "design_len": len(s_design),
            "current_len": next((x["len"] for x in recs
                                 if x["acc"] == cur_ver), None),
            "seq_positional_identity": seq_ident,
            "n_refseq_tx": len(recs),
            "n_coding_tx": len(coding),
            "n_noncoding_tx": len(noncoding),
            "n_coding_tx_covered": len(covered),
            "coverage_coding": (f"{len(covered)}/{len(coding)}" if coding else "n/a"),
            "all_coding_covered": (len(covered) == len(coding) and bool(coding)),
            "F_in_CDS_all_coding": all(x.get("F_in_CDS") for x in covered)
                                    if covered else False,
            "R_in_CDS_all_coding": all(x.get("R_in_CDS") for x in covered)
                                    if covered else False,
        }
        # ── 分级判定（区分三类不同严重度）──
        #   OK          : 模板版本未变（或变更但序列完全一致），
        #                 引物覆盖全部编码转录本，且完整位于 CDS 内
        #   VERSION_ONLY: 仅 RefSeq 版本号变更、序列 100% 一致，
        #                 引物本身仍可用 → 只需更新溯源中的版本号
        #   NEEDS_REDESIGN: 引物实际失效（落位丢失 / 覆盖不全 / 移出 CDS）
        primers_still_valid = (len(covered) == len(coding) and bool(coding)
                               and all(x.get("F_in_CDS") for x in covered)
                               and all(x.get("R_in_CDS") for x in covered))
        if primers_still_valid:
            row["status"] = ("OK" if cur_ver == design_acc
                             else "VERSION_ONLY")
            row["status_note"] = (
                "模板版本未变，引物覆盖全部编码转录本且位于 CDS 内" if cur_ver == design_acc
                else f"RefSeq 版本 {design_acc}→{cur_ver}，但序列 100% 一致，"
                     f"引物仍可原样使用（仅需更新溯源版本号）")
        else:
            row["status"] = "NEEDS_REDESIGN"
            why = []
            if len(covered) != len(coding):
                why.append(f"仅覆盖 {len(covered)}/{len(coding)} 条编码转录本")
            if covered and not all(x.get("F_in_CDS") for x in covered):
                why.append("前向引物未完全位于 CDS 内")
            if covered and not all(x.get("R_in_CDS") for x in covered):
                why.append("反向引物未完全位于 CDS 内")
            row["status_note"] = "；".join(why) or "引物落位失效"
        summary.append(row)
        print(f"{gene:9s} design={design_acc:16s} cur={cur_ver:16s} "
              f"tx={len(recs)}(编码{len(coding)}) cov={row['coverage_coding']:6s} "
              f"id={seq_ident} → {row['status']}")

    df = pd.DataFrame(summary)
    df.to_csv(f"{OUT}/refseq_version_drift_audit.csv", index=False)
    with open(f"{OUT}/refseq_version_drift_audit.json", "w", encoding="utf-8") as fh:
        json.dump({"date": "2026-10-25",
                   "script": "02_脚本/audit_refseq_version_drift.py",
                   "summary": summary, "per_transcript": detail},
                  fh, ensure_ascii=False, indent=1, default=str)
    with open(f"{PROV}/refseq_audit_fetch.log", "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG))

    print("\n=== 分级汇总 ===")
    for st in ("OK", "VERSION_ONLY", "NEEDS_REDESIGN"):
        g = [r["gene"] for r in summary if r["status"] == st]
        print(f"  {st:16s} ({len(g)}): {g}")
    need = [r["gene"] for r in summary if r["status"] == "NEEDS_REDESIGN"]
    print(f"\n需重新设计：{len(need)}/{len(summary)} → {need}")
    print("DONE — audit_refseq_version_drift.py")


if __name__ == "__main__":
    main()
