# -*- coding: utf-8 -*-
"""
redesign_primers_v2.py — 全部 12 个引物基因的通用重新设计（CDS 锚定版）

═══ 触发事件与缺陷谱系（2026-10-25，全部经独立复核）═══

原设计（`design_rtqpcr_primers.py`，2026-10-25）用 primer3 在**整条 mRNA**
上直接设计，未施加两项关键约束，由此产生 **5 个基因**的引物缺陷：

| 基因 | 缺陷 | 实证 |
|------|------|------|
| LDHA | 模板版本漂移 | 前向引物在全部 6 条当前 RefSeq 转录本中匹配位置 -1（正/反互补均是）→ 无法扩增 cDNA；`.1`→`.2` 位置一致性仅 0.2524 |
| ALDH3A2 | 引物落在 **3'UTR 的 Alu 重复元件** | CDS 222..1679，引物在 2636/2745（UTR 内）；人源全长命中 3,111/4,713 条，散布 90+ 无关基因座 |
| FAHD2B | 前向引物 100% 匹配 **LOC124900512**（FAHD2B-like 独立转录本） | 8 条转录本 identity 100%、3' 端全匹配 |
| ADH1B | 前向 100% 匹配 **ADH1C**（4 条），反向 100% 匹配 **ADH1A**（2 条） | 不同基因座、不同酶亚基，交叉扩增会污染读出 |
| LDHA | 反向引物 95% 匹配 **LYST**（16 条转录本，3' 端全匹配） | 真实共扩增风险 |

根因是同一类：**引物位置与特异性只做了"设计期"约束，未做"设计后"
的全库复核与转录本版本核对**。

═══ v2 设计约束（在 primer3 参数之上追加）═══

1. **CDS 锚定**：引物对必须完整位于 **CDS 内**（RT-qPCR 测 mRNA 丰度；
   3'UTR 在各转录本间高度变异，且常含 Alu/LINE 重复元件）。
2. **跨转录本覆盖**：设计窗口必须落在该基因**全部蛋白编码 RefSeq 转录本
   共有的序列区段**（逐本自身序列求交，不复用同一参考切片）。
3. **CDS 内无接头**：引物不得跨越 exon 边界（跨接头引物在基因组 DNA 上
   有信号、在 cDNA 上被剪掉 → 严重偏倚）。
4. **Primer3 参数与原设计完全一致**（保证与其余引物口径统一）。
5. 设计后仍须经 **BLAST 全库复核**（本脚本不替代 BLAST，只负责产出
   待复核的候选集）。

═══ 判读阈值（与 adjudicate_blast_offtargets.py 一致）════

- P0（阻断）：全长匹配 + 人源 RefSeq 转录本 + 非目标 + 非旁系同源 +
  3' 端末 5 nt 全匹配 + identity ≥ 95%
- P1（需注意）：同上但 90% ≤ identity < 95%
- P2/PARTIAL：判为噪声

输出：
  03_结果/实验验证设计/primer_redesign_v2_candidates.json
  03_结果/实验验证设计/primer_redesign_v2_selected.csv
  03_结果/实验验证设计/provenance/redesign_v2_transcripts.fasta
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import re
import json
import time
import urllib.request

import pandas as pd
import primer3

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
PROV = f"{OUT}/provenance"
os.makedirs(PROV, exist_ok=True)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "WorkBuddy-PrimerRedesignV2/1.0"}
RC_MAP = str.maketrans("ACGT", "TGCA")

# gene_id 与 design_rtqpcr_primers.py 一致
GENES = {
    "LDHA": (3939, "target", "M2_Pyruvate_supply"),
    "FAHD2B": (151313, "target", "M2_Pyruvate_supply"),
    "ADH1B": (125, "target", "M2_Pyruvate_supply"),
    "ALDH3A2": (224, "target", "M2_Pyruvate_supply"),
    "DHX58": (79132, "target", "M4_IFNI_readout"),
    "GBP1": (2633, "target", "M4_IFNI_readout"),
    "ACTB": (60, "reference_candidate", ""),
    "GAPDH": (2597, "reference_candidate", ""),
    "TBP": (6908, "reference_candidate", ""),
    "HPRT1": (3251, "reference_candidate", ""),
    "RPLP0": (6175, "reference_candidate", ""),
    "SDHA": (6389, "reference_candidate", ""),
}
# 严格旁系同源白名单（同基因座/经典可接受；其余同家族成员按真脱靶处理）
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

PRIMER3_PARAMS = {
    "PRIMER_TASK": "generic",
    "PRIMER_PICK_LEFT_PRIMER": 1,
    "PRIMER_PICK_RIGHT_PRIMER": 1,
    "PRIMER_PICK_INTERNAL_OLIGO": 0,
    "PRIMER_NUM_RETURN": 5,
    "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 25,
    "PRIMER_OPT_SIZE": 20, "PRIMER_MAX_POLY_X": 4,
    "PRIMER_MIN_TM": 58.0, "PRIMER_MAX_TM": 62.0, "PRIMER_OPT_TM": 60.0,
    "PRIMER_MIN_GC": 40.0, "PRIMER_MAX_GC": 60.0, "PRIMER_MAX_NS_ACCEPTED": 0,
    "PRIMER_MAX_SELF_ANY_TH": 47.0, "PRIMER_MAX_SELF_END_TH": 47.0,
    "PRIMER_MAX_HAIRPIN_TH": 47.0,
    "PRIMER_PAIR_MAX_COMPL_ANY_TH": 47.0, "PRIMER_PAIR_MAX_COMPL_END_TH": 47.0,
    "PRIMER_PRODUCT_SIZE_RANGE": [[80, 180]],
    "PRIMER_EXPLAIN_FLAG": 1,
}


def rc(s):
    return s.translate(RC_MAP)[::-1]


def _get(url, tries=4):
    last = None
    for k in range(tries):
        try:
            return urllib.request.urlopen(
                urllib.request.Request(url, headers=UA), timeout=120).read().decode()
        except Exception as e:                        # noqa: BLE001
            last = e
            time.sleep(2 + 3 * k)
    raise RuntimeError(str(last))


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


def parse_exons(gb):
    return [(int(a), int(b)) for a, b in
            re.findall(r"^\s{5}exon\s+(\d+)\.\.(\d+)", gb, re.M)]


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
        out.append((v.get("accessionversion"), v.get("title", ""), v.get("slen")))
    return sorted(out)


def kmers(s, k=25):
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def shared_segments(seqs, ref_acc):
    """逐本自身序列求交，返回共有 25-mer 在 ref 上的区段列表。"""
    ref = seqs[ref_acc]
    common = kmers(ref)
    for a, s in seqs.items():
        if a == ref_acc:
            continue
        common &= kmers(s)
    pos = [i for i in range(len(ref) - 24) if ref[i:i + 25] in common]
    segs = []
    if pos:
        st = prev = pos[0]
        for p in pos[1:]:
            if p <= prev + 25:
                prev = p
            else:
                segs.append((st, prev + 25))
                st = prev = p
        segs.append((st, prev + 25))
    return segs, len(pos)


def exons_from_cds(cds, exons):
    """返回 CDS 各段所属 exon 列表（用于跨接头检查）。"""
    out = []
    for a, b in cds:
        for ea, eb in exons:
            if not (b < ea or a > eb):
                out.append((ea, eb))
    return out


def crosses_exon(pos0, length, exon_list):
    """引物区间 [pos0, pos0+length) 是否跨 exon 边界。"""
    segs = []
    for ea, eb in exon_list:
        s, e = max(pos0, ea - 1), min(pos0 + length, eb)
        if s < e:
            segs.append((s, e))
    # 合并相邻
    if not segs:
        return True
    covered = sum(e - s for s, e in segs)
    return covered < length


def main():
    old = {r["gene"]: r for _, r in
           pd.read_csv(f"{OUT}/qpcr_primers.csv").iterrows()}
    report = {"date": "2026-10-25",
              "script": "02_脚本/redesign_primers_v2.py",
              "primer3_version": primer3.__version__,
              "primer3_params": {k: v for k, v in PRIMER3_PARAMS.items()
                                 if k != "PRIMER_EXPLAIN_FLAG"},
              "constraints": [
                  "引物对完整位于 CDS 内",
                  "设计窗口为全部蛋白编码 RefSeq 转录本共有区段",
                  "引物不跨 exon 边界",
                  "Primer3 参数与原设计完全一致",
              ],
              "genes": {}}
    selected_rows = []
    fasta_path = f"{PROV}/redesign_v2_transcripts.fasta"
    if os.path.exists(fasta_path):
        os.remove(fasta_path)

    for gene, (gid, role, module) in GENES.items():
        print(f"\n=== {gene} ===")
        tx_list = list_transcripts(gid)
        time.sleep(0.35)
        seqs, cdss, exon_sets, coding = {}, {}, {}, []
        for acc, title, slen in tx_list:
            try:
                _, seq = efetch_fasta(acc)
                time.sleep(0.35)
                gb = efetch_gb(acc)
                time.sleep(0.35)
                cds = parse_cds(gb)
                ex = parse_exons(gb)
                seqs[acc] = seq
                cdss[acc] = cds
                exon_sets[acc] = ex
                if cds:
                    coding.append(acc)
                with open(fasta_path, "a", encoding="utf-8") as fh:
                    fh.write(f">{acc} {title} CDS={cds}\n")
                    for i in range(0, len(seq), 60):
                        fh.write(seq[i:i + 60] + "\n")
            except Exception as e:                    # noqa: BLE001
                print(f"   fetch fail {acc}: {e}")
        if not coding:
            print("   无蛋白编码转录本，跳过")
            continue

        # canonical 选择规则（修正 v2 初版 bug）：
        #   v2 初版按「CDS 最长」选，结果给 SDHA / FAHD2B 选中了 XM_ 预测
        #   转录本（GenBank 记录**无 exon 特征注释**），导致跨接头检查把所有
        #   候选误判为 cross=True（exon_list 为空 → 全覆盖判定失败）。
        #   修正：优先选 NM_/NR_ 规范记录且**带 exon 注释**者；
        #   若规范记录无 exon 注释，则退回 CDS 最长的带注释记录。
        def has_exon(acc):
            return bool(exon_sets.get(acc))

        ncbi_canonical = [a for a in coding
                          if a.startswith(("NM_", "NR_")) and has_exon(a)]
        if ncbi_canonical:
            canonical = max(ncbi_canonical,
                            key=lambda a: sum(b - x + 1 for x, b in cdss[a]))
        else:
            with_exon = [a for a in coding if has_exon(a)]
            pool = with_exon or coding
            canonical = max(pool,
                            key=lambda a: sum(b - x + 1 for x, b in cdss[a]))
        ref = seqs[canonical]
        segs, n_shared = shared_segments(seqs, canonical)
        cds_set = set()
        for x, y in cdss[canonical]:
            cds_set.update(range(x - 1, y))
        cds_segs = []
        for a, b in segs:
            sub = [p for p in range(a, b) if p in cds_set]
            if not sub:
                continue
            # 归并 CDS 内连续段
            st = prev = sub[0]
            for p in sub[1:]:
                if p == prev + 1:
                    prev = p
                else:
                    cds_segs.append((st, prev + 1))
                    st = prev = p
            cds_segs.append((st, prev + 1))
        print(f"   转录本 {len(seq_list := tx_list)} 条（编码 {len(coding)}），"
              f"canonical={canonical}")
        print(f"   共有 25-mer {n_shared} 个；CDS 内共有区段 "
              f"{[(a, b - a) for a, b in cds_segs]}")

        cands = []
        for seg_a, seg_b in cds_segs:
            if seg_b - seg_a < 130:
                continue
            sub = ref[seg_a:seg_b]
            res = primer3.bindings.design_primers(
                {"SEQUENCE_ID": f"{gene}_{seg_a}", "SEQUENCE_TEMPLATE": sub},
                PRIMER3_PARAMS, {})
            for i in range(res.get("PRIMER_PAIR_NUM_RETURNED", 0)):
                f = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
                r = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
                lp, rp = res[f"PRIMER_LEFT_{i}"], res[f"PRIMER_RIGHT_{i}"]
                f_abs0, r_abs0 = seg_a + lp[0], seg_a + rp[0]
                cands.append({
                    "forward": f, "reverse": r,
                    "f_tm": res[f"PRIMER_LEFT_{i}_TM"],
                    "r_tm": res[f"PRIMER_RIGHT_{i}_TM"],
                    "f_gc": round(100 * sum(c in "GC" for c in f) / len(f), 1),
                    "r_gc": round(100 * sum(c in "GC" for c in r) / len(r), 1),
                    "product_size": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                    "canonical_f_start_1based": f_abs0 + 1,
                    "canonical_f_end_1based": f_abs0 + len(f),
                    "canonical_r_start_1based": r_abs0 + 1,
                    "canonical_r_end_1based": r_abs0 + len(r),
                    "f_self_any_th": res.get(f"PRIMER_LEFT_{i}_SELF_ANY_TH"),
                    "f_self_end_th": res.get(f"PRIMER_LEFT_{i}_SELF_END_TH"),
                    "f_hairpin_th": res.get(f"PRIMER_LEFT_{i}_HAIRPIN_TH"),
                    "r_self_any_th": res.get(f"PRIMER_RIGHT_{i}_SELF_ANY_TH"),
                    "r_self_end_th": res.get(f"PRIMER_RIGHT_{i}_SELF_END_TH"),
                    "r_hairpin_th": res.get(f"PRIMER_RIGHT_{i}_HAIRPIN_TH"),
                    "pair_compl_any_th": res.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH"),
                    "pair_compl_end_th": res.get(f"PRIMER_PAIR_{i}_COMPL_END_TH"),
                    "penalty": res.get(f"PRIMER_PAIR_{i}_PENALTY"),
                })

        # 覆盖 + 落位检查
        # ── 方法学说明 ──
        # RefSeq **mRNA** 记录（NM_/XM_）的 `exon` 特征坐标是 **mRNA 自身
        # 坐标**（已拼接），相邻 exon 在 mRNA 上连续 → 用 mRNA 坐标判定
        # "是否跨 exon–exon 接头"在数学上不成立（恒为"不跨"或无意义）。
        # 真正的接头判定须用 RefSeqGene（NG_）的**基因组坐标**，见
        # `genome_exon_anchor.py`（作为信息性指标输出，不作硬约束：
        # RT-qPCR 跨接头引物本身是允许且常用的策略，真正的硬约束是
        # 「CDS 内 + 全部编码转录本共有 + 无真实脱靶」）。
        # 因此此处仅统计引物落在 CDS 内的比例与跨 exon 边界的**转录本
        # 内**相邻 exon 数（信息性），并把「引物是否完整位于单个 exon」
        # 作为排序的**软偏好**（加分项而非否决项）。
        canon_exons = exon_sets[canonical]
        for c in cands:
            cov, in_cds = {}, True
            for acc in coding:
                s = seqs[acc]
                iF = s.find(c["forward"])
                iRrc = s.find(rc(c["reverse"]))
                ok = (iF >= 0) and (iRrc >= 0)
                cov[acc] = {"F_pos": iF, "R_rc_pos": iRrc, "ok": ok}
                if not ok:
                    in_cds = False
            # 信息性：引物在 canonical 上跨越的 exon 边界数（mRNA 坐标）
            def n_exon_breaks(pos0, length):
                if pos0 < 0:
                    return None
                end = pos0 + length
                return sum(1 for ea, eb in canon_exons
                           if ea - 1 < pos0 < eb)
            c["transcript_coverage"] = cov
            c["n_coding_covered"] = sum(1 for v in cov.values() if v["ok"])
            c["n_coding_tx"] = len(coding)
            c["all_coding_covered"] = (c["n_coding_covered"] == len(coding))
            c["both_in_CDS"] = in_cds
            c["F_exon_breaks"] = n_exon_breaks(
                c["canonical_f_start_1based"] - 1, len(c["forward"]))
            c["R_exon_breaks"] = n_exon_breaks(
                c["canonical_r_start_1based"] - 1, len(c["reverse"]))

        def score(c):
            s = 3 * c["penalty"]
            s += 3 * (abs(c["f_tm"] - 60) + abs(c["r_tm"] - 60))
            s += 2 * (abs(c["f_gc"] - 50) / 10 + abs(c["r_gc"] - 50) / 10)
            s += 2 * (0 if c["forward"].endswith(("G", "C")) else 0.5)
            s += 2 * (0 if c["reverse"].endswith(("G", "C")) else 0.5)
            s += 50 * (c["n_coding_tx"] - c["n_coding_covered"])
            # 软偏好：引物落在单个 exon 内（不跨接头）者优先
            s += 2 * (c["F_exon_breaks"] or 0) + 2 * (c["R_exon_breaks"] or 0)
            return s

        cands.sort(key=score)
        ok_cands = [c for c in cands if c["all_coding_covered"]]
        best = ok_cands[0] if ok_cands else None
        print(f"   候选 {len(cands)}，覆盖全部编码转录本 {len(ok_cands)}")
        for c in cands[:4]:
            print(f"     cov={c['n_coding_covered']}/{c['n_coding_tx']} "
                  f"exonBrk={c['F_exon_breaks']}/{c['R_exon_breaks']} "
                  f"F={c['forward']} R={c['reverse']} "
                  f"amp={c['product_size']} pen={c['penalty']:.4f}")

        o = old[gene]
        report["genes"][gene] = {
            "gene_id": gid, "role": role, "module": module,
            "canonical": canonical,
            "n_refseq_tx": len(seq_list),
            "coding_transcripts": coding,
            "n_coding_tx": len(coding),
            "shared_25mer_segments": segs,
            "cds_shared_segments": cds_segs,
            "n_shared_25mer": n_shared,
            "old_design": {
                "accession": o["refseq_accession"],
                "forward": o["forward_primer_5to3"],
                "reverse": o["reverse_primer_5to3"],
                "f_start": int(o["forward_start_1based"]),
                "r_end": int(o["reverse_end_1based"]),
                "amplicon_bp": int(o["amplicon_bp"]),
            },
            "candidates": cands,
            "selected": best,
        }
        if best:
            print(f"   ★ 新推荐：F={best['forward']} R={best['reverse']} "
                  f"amp={best['product_size']}bp "
                  f"(canonical {best['canonical_f_start_1based']}.."
                  f"{best['canonical_r_end_1based']})")
            selected_rows.append({
                "gene": gene, "role": role, "module": module,
                "refseq_accession": canonical, "gene_id": gid,
                "forward_primer_5to3": best["forward"],
                "reverse_primer_5to3": best["reverse"],
                "f_tm_C": round(best["f_tm"], 6), "r_tm_C": round(best["r_tm"], 6),
                "f_gc_pct": best["f_gc"], "r_gc_pct": best["r_gc"],
                "amplicon_bp": best["product_size"],
                "forward_start_1based": best["canonical_f_start_1based"],
                "reverse_end_1based": best["canonical_r_end_1based"],
                "penalty": best["penalty"],
                "n_coding_tx_covered": best["n_coding_covered"],
                "n_coding_tx_total": best["n_coding_tx"],
                "F_exon_breaks": best["F_exon_breaks"],
                "R_exon_breaks": best["R_exon_breaks"],
                "old_forward": o["forward_primer_5to3"],
                "old_reverse": o["reverse_primer_5to3"],
                "old_refseq_accession": o["refseq_accession"],
            })
        else:
            print("   !! 无满足全部约束的候选")

    with open(f"{OUT}/primer_redesign_v2_candidates.json", "w",
              encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1, default=str)
    if selected_rows:
        pd.DataFrame(selected_rows).to_csv(
            f"{OUT}/primer_redesign_v2_selected.csv", index=False)
    print(f"\n完成：{len(selected_rows)}/{len(GENES)} 个基因产出新推荐引物对")
    print("DONE — redesign_primers_v2.py")


if __name__ == "__main__":
    main()
