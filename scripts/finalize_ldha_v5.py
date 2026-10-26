# -*- coding: utf-8 -*-
"""
finalize_ldha_v5.py — 将经 BLAST 确认的 LDHA v5 引物写回主表与快照

读取 primer_blast_v5_ldha_candidate8.json 中判定为 USABLE 的对（优先候选 8），
用 Primer3 重算 Tm/penalty（与 11 对其他对口径一致），写回：
  - 03_结果/实验验证设计/qpcr_primers.csv  （LDHA 行替换；其余 11 行与列全部保留）
  - 03_结果/实验验证设计/final_primers_v4.csv（同步 LDHA 行，保持快照一致）
  - 03_结果/实验验证设计/primer_redesign_v2_selected.csv（同上）

偏差说明列（scope_note）一并写入 qpcr_primers.csv，记录：
  “3'UTR-anchored on NM_005566.4 only; CDS-anchoring infeasible (LOC105375980
   carries 96.16% of LDHA CDS); quantifies NM_005566.4 transcript only.”
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import csv
import json
import importlib.util

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"

# Primer3（envs/default venv 已装）
import primer3  # noqa: E402


def tm(seq):
    return round(primer3.bindings.calc_tm(seq), 6)


def gc(seq):
    return round(100.0 * sum(c in "GC" for c in seq) / len(seq), 1)


def main():
    c8 = json.load(open(f"{OUT}/primer_blast_v5_ldha_candidate8.json", encoding="utf-8"))
    chosen = c8.get("selected") or next(
        (r for r in c8["results"] if r["verdict"] == "USABLE"), None)
    if not chosen or chosen["verdict"] != "USABLE":
        raise SystemExit("候选 8 / 6 / 7 均无 USABLE 对，终止写回，需重新设计")
    F = chosen["forward"]
    R = chosen["reverse"]
    geom = chosen["geometry"]
    print(f"选中 LDHA v5 对：F={F}  R={R}  product={geom['product_size']} "
          f"GC差={geom['gc_diff']}%  verdict={chosen['verdict']}")

    f_tm = tm(F)
    r_tm = tm(R)
    f_gc = gc(F)
    r_gc = gc(R)
    amplicon = geom["product_size"]
    f_start = geom["f_pos_1based"]
    r_end = f_start + amplicon - 1
    # ★ 2026-10-25 更正：先前写作 "quantifies NM_005566.4 only"，与 BLAST 实测矛盾。
    #   forward/reverse 的 target_full_match_refseq_tx = 6，即两端引物各自与 6 条
    #   LDHA RefSeq 转录本全长匹配（3'UTR 该区段为多条转录本共有；独立复核亦显示
    #   NM_001135239.1 上 F@1832 / amp=174，与 NM_005566.4 相同）。
    #   故正确表述是"覆盖多条 LDHA RefSeq 转录本、报告总 LDHA 转录本丰度"，
    #   偏离点只在于"不与 LOC105375980 共扩增"，而非"只测单一转录本"。
    note = ("3'UTR-anchored: both primers lie in the LDHA 3'UTR and fully match 6 "
            "LDHA RefSeq transcripts (incl. NM_005566.4), so the assay reports "
            "aggregate LDHA transcript abundance rather than a single isoform. "
            "CDS-anchoring was infeasible because LOC105375980 (NR_202040.1) carries "
            "96.16% of the LDHA CDS and would co-amplify; this pair does not "
            "co-amplify LOC105375980 (BLAST LOC hits = 0).")

    # ── qpcr_primers.csv ──
    qp = list(csv.DictReader(open(f"{OUT}/qpcr_primers.csv", encoding="utf-8")))
    cols = list(qp[0].keys())
    if "scope_note" not in cols:
        cols.append("scope_note")
    for r in qp:
        if r["gene"] == "LDHA":
            r["refseq_accession"] = "NM_005566.4"
            r["gene_id"] = "3939"
            r["template_len"] = "2241"
            r["forward_primer_5to3"] = F
            r["reverse_primer_5to3"] = R
            r["f_tm_C"] = f_tm
            r["r_tm_C"] = r_tm
            r["f_gc_pct"] = f_gc
            r["r_gc_pct"] = r_gc
            r["amplicon_bp"] = amplicon
            r["forward_start_1based"] = f_start
            r["reverse_end_1based"] = r_end
            r["penalty"] = ""
            r["scope_note"] = note
    with open(f"{OUT}/qpcr_primers.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(qp)
    print(f"已写回 qpcr_primers.csv：LDHA -> NM_005566.4 / {F} / {R} / "
          f"{amplicon}bp / Tm {f_tm:.2f}/{r_tm:.2f} / GC {f_gc}/{r_gc}")

    # ── 同步 final_primers_v4.csv 与 primer_redesign_v2_selected.csv（若存在）──
    for snap in ("final_primers_v4.csv", "primer_redesign_v2_selected.csv"):
        p = f"{OUT}/{snap}"
        if not os.path.exists(p):
            continue
        rows = list(csv.DictReader(open(p, encoding="utf-8")))
        scols = list(rows[0].keys())
        if "scope_note" not in scols:
            scols.append("scope_note")
        for r in rows:
            if r.get("gene") == "LDHA":
                r["refseq_accession"] = "NM_005566.4"
                r["forward_primer_5to3"] = F
                r["reverse_primer_5to3"] = R
                r["f_tm_C"] = f_tm
                r["r_tm_C"] = r_tm
                r["f_gc_pct"] = f_gc
                r["r_gc_pct"] = r_gc
                r["amplicon_bp"] = amplicon
                r["scope_note"] = note
                for k in ("forward_start_1based", "reverse_end_1based",
                          "template_len", "gene_id", "penalty"):
                    if k in scols:
                        r[k] = {"forward_start_1based": f_start,
                                "reverse_end_1based": r_end,
                                "template_len": "2241", "gene_id": "3939",
                                "penalty": ""}.get(k, r.get(k))
        with open(p, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=scols)
            w.writeheader()
            w.writerows(rows)
        print(f"已同步 {snap}")

    # 落位实测校验（独立代码路径）：F 与 rc(R) 必须精确命中 NM_005566.4
    cache = json.load(open(f"{OUT}/provenance/final_design_cache.json", encoding="utf-8"))
    seq = cache["NM_005566.4"]["seq"].upper()
    rc = lambda s: s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]
    assert seq.find(F) >= 0, "F 未在 NM_005566.4 命中"
    assert seq.find(rc(R)) >= 0, "rc(R) 未在 NM_005566.4 命中"
    print("落位校验通过：F 与 rc(R) 均精确命中 NM_005566.4")
    print("DONE — finalize_ldha_v5.py")


if __name__ == "__main__":
    main()
