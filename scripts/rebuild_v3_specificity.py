#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从 primer_blast_v5_ldha_all.json 重建 primer_blast_specificity_v3.csv（S10d 数据源）

2026-10-25 修复记录：
  regen_ldha_v5_blast.py 的 _v3_map 原先只从 detail 中找 level=='TARGET' 的记录，
  而 detail 只记录 P0/P1/PARALOG/PSEUDO 等脱靶层命中、不含 TARGET 层 →
  重建出的 v3 CSV 中全部 12 个基因的 target_human_hits / target_max_identity /
  target_top_acc / n_human_refseq_tx 均为 0 或空，与 specific=True 自相矛盾
  （"引物不命中任何目标" 却 "特异"），S10d 附表因此不可信。
  regen 脚本的 _v3_map 已改为从 rec（rows_out）取目标统计；
  本脚本用同一逻辑从既有产物重建 CSV，无需重跑 BLAST（无网络请求）。
"""
import csv

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import json

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"

SRC = f"{OUT}/primer_blast_v5_ldha_all.json"
DST = f"{OUT}/primer_blast_specificity_v3.csv"

V3_COLS = ["gene", "role", "forward_target_human_hits",
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


def _v3_map(rec, lab, detail_all):
    """与修复后 regen_ldha_v5_blast.py::_v3_map 同一逻辑。"""
    # detail 的 primer 键是 "F"/"R"（不是 "forward"/"reverse"）——
    # 2026-10-25 修复：旧代码用 lab 直接匹配，导致脱靶基因列表恒为空。
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
        f"{lab}_human_hit_genes": ";".join(sorted(set(human_genes) | {rec["gene"]})),
    }


def main():
    src = json.load(open(SRC, encoding="utf-8"))
    rows_out = src["rows"]
    detail_all = src["detail"]

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

    with open(DST, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=V3_COLS)
        w.writeheader()
        for r in v3_rows:
            w.writerow({k: r.get(k, "") for k in V3_COLS})

    # ── 自检 ──
    # 1) specific 列与 v6 CSV 一致（CSV 读出为字符串，统一转 str 比较）
    v6 = {r["gene"]: r for r in csv.DictReader(
        open(f"{OUT}/primer_blast_v6.csv", encoding="utf-8"))}
    bad = [g for g, r in v6.items()
           if str(next(x for x in v3_rows if x["gene"] == g)["specific"]) != str(r["specific"])]
    # 2) 每个基因至少一个方向 target_human_hits >= 1（引物必须命中自身目标）
    zero = [r["gene"] for r in v3_rows
            if (r["forward_target_human_hits"] or 0) < 1
            and (r["reverse_target_human_hits"] or 0) < 1]
    # 3) target_top_acc 非空且为设计模板 accession
    noacc = [r["gene"] for r in v3_rows
             if not r["forward_target_top_acc"] or not r["reverse_target_top_acc"]]
    # 4) 有 P0 脱靶的基因，其 P0 基因必须出现在 offtarget_tx_genes 中
    import json as _json
    det = _json.load(open(SRC, encoding="utf-8"))["detail"]
    missing = []
    for r in v3_rows:
        for lab, p in (("forward", "F"), ("reverse", "R")):
            p0 = {d["hit_gene"] for d in det
                  if d["gene"] == r["gene"] and d["primer"] == p and d["level"] == "P0"}
            listed = set(x for x in r[f"{lab}_offtarget_tx_genes"].split(";") if x)
            if p0 - listed:
                missing.append((r["gene"], lab, sorted(p0 - listed)))

    print(f"重建 {DST}（{len(v3_rows)} 行）")
    print("自检1 specific↔v6：", "一致" if not bad else f"失败 {bad}")
    print("自检2 每基因至少一端命中目标：", "通过" if not zero else f"失败 {zero}")
    print("自检3 top_acc 非空：", "通过" if not noacc else f"失败 {noacc}")
    print("自检4 P0 基因出现在 offtarget_tx_genes：", "通过" if not missing else f"失败 {missing}")
    for r in v3_rows:
        print(f"  {r['gene']:8s} tgt F/R={r['forward_target_human_hits']}/{r['reverse_target_human_hits']} "
              f"maxid={r['forward_target_max_identity']}/{r['reverse_target_max_identity']} "
              f"offtx={r['offtarget_refseq_tx_total']} offtx_genes(F)={r['forward_offtarget_tx_genes'][:40]} "
              f"specific={r['specific']}")


if __name__ == "__main__":
    main()
