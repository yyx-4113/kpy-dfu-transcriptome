#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""同步 primer_blast_v6.json 与最新重跑产物（修复陈旧 JSON 与新鲜 CSV 互相矛盾的隐患）

背景（2026-10-25）：
  regen_ldha_v5_blast.py 重跑后只写 primer_blast_v5.csv / primer_blast_v6.csv /
  primer_blast_specificity_v3.csv / primer_blast_v5_ldha_all.json，
  而 primer_blast_v6.json 由旧的 verify_primers_blast_v6.py 在 18:25 写入，
  内容是修复前的判定（LDHA 用旧引物判 OFFTARGET_P0、GAPDH 误判 USABLE），
  与 22:51 的新 CSV 直接矛盾 → 任何读者拿到的结论取决于读哪个文件。

本脚本不做任何网络请求，仅从权威产物 primer_blast_v5_ldha_all.json
重建 primer_blast_v6.json（保持旧扁平 schema，键名不变），使文件名继续可用
且内容与 CSV 一致。reglen 脚本已同步打补丁：今后重跑会直接写该 JSON。
"""
import json

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))


BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"

SRC = f"{OUT}/primer_blast_v5_ldha_all.json"   # 权威：重跑全量明细
DST = f"{OUT}/primer_blast_v6.json"            # 目标：与 CSV 同步的 v6 JSON


def _genes(detail, gene, primer, levels):
    return sorted({d["hit_gene"] for d in detail
                   if d["gene"] == gene and d["primer"] == primer
                   and d["level"] in levels})


def main():
    src = json.load(open(SRC, encoding="utf-8"))
    results = []
    for r in src["rows"]:
        g = r["gene"]
        rec = dict(r)
        for lab in ("forward", "reverse"):
            rec[f"{lab}_P0_genes"] = _genes(src["detail"], g, lab[0], ("P0",))
            rec[f"{lab}_P1_genes"] = _genes(src["detail"], g, lab[0], ("P1",))
            # 旧 schema 中无这些键的置空，保持键集合一致
            for k in (f"{lab}_n_hits_total", f"{lab}_n_partial_human",
                      f"{lab}_blast_rid", f"{lab}_primer_len"):
                rec.setdefault(k, "")
        results.append(rec)

    payload = {
        "date": src.get("date", "2026-10-25"),
        "script": "regen_ldha_v5_blast.py（重跑）→ sync_blast_v6_json.py（本脚本同步）",
        "database": src.get("database"),
        "program": "blastn",
        "output_format": "XML",
        "criteria": "v6 分层判据：TARGET/PARALOG/PSEUDO/P0(3'端5nt全匹配且id>=95%)/"
                    "P1(90<=id<95%)；usable = 目标全长 RefSeq 转录本>=1 且 P0==0 且 "
                    "目标 identity>=90%",
        "method_fix": src.get("method_fix"),
        "method_note": "本文件由 sync_blast_v6_json.py 从 primer_blast_v5_ldha_all.json "
                       "重建；此前 18:25 版本为修复前判定（含旧 LDHA 引物），已废弃。",
        "summary": src["summary"],
        "results": results,
        "detail": src["detail"],
    }
    with open(DST, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1, default=str)

    # ── 自检：JSON 与 CSV 判定必须逐基因一致 ──
    import csv
    csv_rows = {r["gene"]: r["verdict"] for r in csv.DictReader(
        open(f"{OUT}/primer_blast_v6.csv", encoding="utf-8"))}
    bad = [g for g, v in csv_rows.items()
           if next(x for x in results if x["gene"] == g)["verdict"] != v]
    print(f"写入 {DST}")
    print(f"自检：{len(csv_rows)} 个基因 CSV↔JSON 判定一致"
          if not bad else f"自检失败：{bad}")
    for s in src["summary"]:
        print(f"  {s['gene']:8s} F_P0={s['F_P0']:>3} R_P0={s['R_P0']:>3} "
              f"tgtTx={s['F_tgtTx']}/{s['R_tgtTx']} -> {s['verdict']}")


if __name__ == "__main__":
    main()
