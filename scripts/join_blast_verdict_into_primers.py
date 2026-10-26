#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 BLAST v6 重跑判定写回 qpcr_primers.csv（新增 blast_v6_verdict 列）

2026-10-25：regen_ldha_v5_blast.py 完成 12 对重跑（6 USABLE / 6 OFFTARGET_P0），
但 qpcr_primers.csv 只有设计信息、没有验证判定列，读者无法从引物表直接看出
哪些对不可用。本脚本以 primer_blast_v6.csv 为权威源做 join，写回判定列与
P0 脱靶基因列表，并在 scope_note 中追加判定摘要（保留原有全部文字）。
"""
import csv

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))


BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"

P = f"{OUT}/qpcr_primers.csv"
B = f"{OUT}/primer_blast_v6.csv"

blast = {r["gene"]: r for r in csv.DictReader(open(B, encoding="utf-8"))}
rows = list(csv.DictReader(open(P, encoding="utf-8-sig")))

FIELDS = list(rows[0].keys()) + ["blast_v6_verdict", "blast_v6_P0_genes"]

for r in rows:
    b = blast[r["gene"]]
    r["blast_v6_verdict"] = b["verdict"]
    # P0 脱靶基因列表来自 v5/v6 明细 JSON（CSV 只有计数）
    import json
    det = json.load(open(f"{OUT}/primer_blast_v5_ldha_all.json", encoding="utf-8"))["detail"]
    p0 = sorted({d["hit_gene"] for d in det
                 if d["gene"] == r["gene"] and d["level"] == "P0"})
    r["blast_v6_P0_genes"] = ";".join(p0)
    if b["verdict"] != "USABLE":
        r["scope_note"] = (
            (r["scope_note"] + " " if r["scope_note"] else "")
            + f"BLAST v6 verdict = {b['verdict']} "
              f"(F_P0={b['forward_n_P0']}, R_P0={b['reverse_n_P0']}; "
              f"P0 off-target genes: {';'.join(p0) or 'none'}); "
              "disqualified for wet-lab use — retained as design record only.")

with open(P, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    w.writeheader()
    for r in rows:
        w.writerow({k: r.get(k, "") for k in FIELDS})

# 自检：判定列与 blast CSV 完全一致
rows2 = {r["gene"]: r for r in csv.DictReader(open(P, encoding="utf-8"))}
bad = [g for g, b in blast.items() if rows2[g]["blast_v6_verdict"] != b["verdict"]]
print("写回完成：qpcr_primers.csv 新增 blast_v6_verdict / blast_v6_P0_genes 列")
print("自检：判定列与 primer_blast_v6.csv 一致" if not bad else f"自检失败：{bad}")
for g, b in blast.items():
    print(f"  {g:8s} {b['verdict']:14s} P0_genes={rows2[g]['blast_v6_P0_genes']}")
