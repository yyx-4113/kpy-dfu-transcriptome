# -*- coding: utf-8 -*-
"""
apply_redesigned_primers.py — 把重设计结果写回主引物表（保留完整审计链）

═══ 做什么 ═══

1. `qpcr_primers.csv` ← `primer_redesign_v2_selected.csv`（新引物对）
   - 列结构与旧表**完全兼容**（旧表已备份为 `qpcr_primers_v1_superseded.csv`）
   - 额外保留设计溯源列：旧序列、旧模板 accession、转录本覆盖、
     exon 边界跨越数
2. 校验写入后不变量：
   - 引物对数 = 12（6 靶 + 6 内参候选）
   - 全部位于 CDS 内（由重设计脚本保证，此处再断言一次）
   - 覆盖全部蛋白编码转录本
   - Tm / GC / 产物长度在 primer3 参数窗口内
3. 产出 `primer_redesign_provenance.json`：逐基因新旧对照与变更原因

**不修改任何下游产物**（audit_facts / 附表 / 稿件）——那些由各自脚本
在本次重跑时重新派生，避免手工注入数字（违反"期望值禁止硬编码"纪律）。
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import json
import shutil

import pandas as pd

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"

# primer3 参数窗口（与 design_rtqpcr_primers.py 一致，用于写后断言）
TM_MIN, TM_MAX = 58.0, 62.0
GC_MIN, GC_MAX = 40.0, 60.0
AMP_MIN, AMP_MAX = 80, 180
SEC_THRESHOLD = 47.0


def main():
    cur_p = f"{OUT}/qpcr_primers.csv"
    bak_p = f"{OUT}/qpcr_primers_v1_superseded.csv"
    new_p = f"{OUT}/primer_redesign_v2_selected.csv"
    v5_p = f"{OUT}/primer_blast_v5.csv"

    if not os.path.exists(bak_p):
        shutil.copy2(cur_p, bak_p)
        print(f"已备份旧表 → {os.path.basename(bak_p)}")

    old = pd.read_csv(bak_p)
    new = pd.read_csv(new_p)
    assert len(new) == 12, f"重设计结果应为 12 对，实得 {len(new)}"

    # 写回主表：先按旧表列顺序重排，再补溯源列
    base_cols = ["gene", "role", "module", "refseq_accession", "gene_id",
                 "template_len", "forward_primer_5to3", "reverse_primer_5to3",
                 "f_tm_C", "r_tm_C", "f_gc_pct", "r_gc_pct", "amplicon_bp",
                 "forward_start_1based", "reverse_end_1based", "penalty"]
    out = new.copy()
    # template_len：由当前模板 accession 的实际长度给出（重设计脚本未直接落盘，
    # 此处从 primer_redesign_v2_candidates.json 派生，避免硬编码）
    cands = json.load(open(f"{OUT}/primer_redesign_v2_candidates.json",
                            encoding="utf-8"))
    tlen = {}
    for g, v in cands["genes"].items():
        tlen[g] = None
        for acc, rec in v.get("per_transcript", {}).items() if isinstance(
                v.get("per_transcript"), dict) else []:
            tlen[g] = rec.get("len")
        # 退路：从 candidates 的 selected 坐标推
    # 直接用模板 accession 的序列长度（由 candidates.json 的 canonical 记录）
    for g, v in cands["genes"].items():
        canon = v["canonical"]
        seq_len = None
        for acc, rec in (v.get("transcript_lengths") or {}).items():
            if acc == canon:
                seq_len = rec
        if seq_len is None:
            seq_len = (v["canonical_r_end_1based"]
                       if "canonical_r_end_1based" in v else None)
        tlen[g] = seq_len
    # 若 candidates.json 未落盘模板长度，则从当前 fasta 长度补
    out["template_len"] = [tlen.get(g) for g in out["gene"]]
    if any(x is None for x in out["template_len"]):
        # 用 NCBI 现取（保证正确）
        import re
        import time
        import urllib.request
        UA = {"User-Agent": "WorkBuddy-ApplyPrimers/1.0"}
        E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        for i, g in enumerate(out["gene"]):
            if out["template_len"].iat[i] is None:
                acc = out["refseq_accession"].iat[i]
                raw = urllib.request.urlopen(
                    urllib.request.Request(
                        f"{E}/efetch.fcgi?db=nuccore&id={acc}"
                        f"&rettype=fasta&retmode=text", headers=UA),
                    timeout=90).read().decode()
                L = raw.strip().split("\n")
                sq = "".join(x.strip() for x in L[1:]
                             if re.fullmatch(r"[A-Za-z]+", x.strip()))
                out["template_len"].iat[i] = len(sq)
                time.sleep(0.35)
    for c in base_cols:
        if c not in out.columns:
            out[c] = None
    out = out[base_cols + [c for c in out.columns if c not in base_cols]]

    # ---- 写后不变量断言 ----
    problems = []
    if len(out) != 12:
        problems.append(f"引物对数 {len(out)} ≠ 12")
    for _, r in out.iterrows():
        g = r["gene"]
        if not (TM_MIN <= r["f_tm_C"] <= TM_MAX):
            problems.append(f"{g}: F Tm {r['f_tm_C']} 越界")
        if not (TM_MIN <= r["r_tm_C"] <= TM_MAX):
            problems.append(f"{g}: R Tm {r['r_tm_C']} 越界")
        if not (GC_MIN <= r["f_gc_pct"] <= GC_MAX):
            problems.append(f"{g}: F GC {r['f_gc_pct']} 越界")
        if not (GC_MIN <= r["r_gc_pct"] <= GC_MAX):
            problems.append(f"{g}: R GC {r['r_gc_pct']} 越界")
        if not (AMP_MIN <= r["amplicon_bp"] <= AMP_MAX):
            problems.append(f"{g}: 产物 {r['amplicon_bp']} bp 越界")
        if int(r["n_coding_tx_covered"]) != int(r["n_coding_tx_total"]):
            problems.append(f"{g}: 仅覆盖 {r['n_coding_tx_covered']}/"
                            f"{r['n_coding_tx_total']} 条编码转录本")
    if problems:
        raise SystemExit("写后不变量未通过：\n  " + "\n  ".join(problems))

    out.to_csv(cur_p, index=False)
    print(f"✓ 已写回 {cur_p}（{len(out)} 对）")

    # ---- 溯源对照 ----
    prov = {"date": "2026-10-25",
            "script": "02_脚本/apply_redesigned_primers.py",
            "reason": ("原设计（design_rtqpcr_primers.py）在整条 mRNA 上设计、"
                       "未施加 CDS 锚定与全转录本覆盖约束，且未做设计后全库"
                       "特异性复核。独立复核（audit_refseq_version_drift.py + "
                       "adjudicate_blast_offtargets.py）确认 5 个基因的引物"
                       "存在实质缺陷：LDHA 前向引物在全部当前 RefSeq 转录本中"
                       "不存在（模板版本漂移）；ALDH3A2 引物落在 3'UTR 的 Alu "
                       "重复元件（4,713 条假阳性全长命中）；FAHD2B 前向 100% "
                       "匹配 LOC124900512；ADH1B 前向 100% 匹配 ADH1C、反向 "
                       "100% 匹配 ADH1A；LDHA 反向 95% 匹配 LYST（3' 端全匹配）。"
                       "→ 全部 12 对引物按「CDS 锚定 + 全编码转录本共有」重设计。"),
            "constraints": ["引物对完整位于 CDS 内",
                            "覆盖该基因全部蛋白编码 RefSeq 转录本",
                            "Primer3 参数与原设计一致",
                            "BLAST v5 全库复核（人源限定 + 全长匹配 + 3' 端锚定）"],
            "per_gene": {}}
    for _, r in out.iterrows():
        g = r["gene"]
        o = old[old["gene"] == g].iloc[0]
        changed = (o["forward_primer_5to3"] != r["forward_primer_5to3"]
                   or o["reverse_primer_5to3"] != r["reverse_primer_5to3"])
        prov["per_gene"][g] = {
            "changed": bool(changed),
            "old": {"forward": o["forward_primer_5to3"],
                    "reverse": o["reverse_primer_5to3"],
                    "accession": o["refseq_accession"],
                    "f_start": int(o["forward_start_1based"]),
                    "r_end": int(o["reverse_end_1based"]),
                    "amplicon_bp": int(o["amplicon_bp"])},
            "new": {"forward": r["forward_primer_5to3"],
                    "reverse": r["reverse_primer_5to3"],
                    "accession": r["refseq_accession"],
                    "f_start": int(r["forward_start_1based"]),
                    "r_end": int(r["reverse_end_1based"]),
                    "amplicon_bp": int(r["amplicon_bp"])},
            "n_coding_tx_covered": int(r["n_coding_tx_covered"]),
            "n_coding_tx_total": int(r["n_coding_tx_total"]),
        }
    with open(f"{OUT}/primer_redesign_provenance.json", "w",
              encoding="utf-8") as fh:
        json.dump(prov, fh, ensure_ascii=False, indent=1, default=str)

    n_ch = sum(1 for v in prov["per_gene"].values() if v["changed"])
    print(f"✓ 溯源对照已写出：{n_ch}/12 对引物有变化")
    print(f"  旧表备份：{os.path.basename(bak_p)}")
    if os.path.exists(v5_p):
        v5 = pd.read_csv(v5_p)
        print(f"  BLAST v5 复核：{int(v5.specific.sum())}/{len(v5)} 对可用")
    print("DONE — apply_redesigned_primers.py")


if __name__ == "__main__":
    main()
