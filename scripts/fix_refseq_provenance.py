# -*- coding: utf-8 -*-
"""
fix_refseq_provenance.py — 修复 RT-qPCR 引物模板的 RefSeq 溯源文件（任务 #11d）

═══ 缺陷说明 ═══
原 `refseq_mrna_sequences.json` 存在两处**溯源错误**（非引物设计错误）：
  1. LDHA 存入的是 NM_001165414.2（transcript variant 3），而设计时实际使用的
     模板是 NM_001135239.1（transcript variant 2）；
  2. ACTB 存入的是 **XM_054343706.1（POTE ankyrin domain family member F / POTEF）**
     ——这是 NCBI Gene→nuccore elink 返回列表中排在前面的一个**预测型**记录，
     被错误地当成 ACTB 的模板，而设计时实际使用的是 NM_001101.3。
  后果：用该溯源文件做"引物是否落在模板上"的复核会得到假阴性
  （LDHA 与 ACTB 显示 F/R 均不命中）。

本脚本按**设计脚本实际使用的 accession** 重新 efetch 并写入溯源文件，
同时逐对实测引物在模板上的落位坐标与扩增子长度（独立于 primer3 的输出）。

输出：
  03_结果/实验验证设计/provenance/refseq_mrna_sequences.json  （修正版）
  03_结果/实验验证设计/provenance/primer_template_mapping.csv（逐对落位实测）
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import json
import time
import urllib.request

import pandas as pd

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
PROV = f"{OUT}/provenance"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
os.makedirs(PROV, exist_ok=True)

# 与 02_脚本/design_rtqpcr_primers.py 的 TARGETS/REFERENCES 完全一致
DESIGN_ACC = {
    "LDHA": ("NM_001135239.1", 3939, "target"),
    "FAHD2B": ("NM_001320848.2", 151313, "target"),
    "ADH1B": ("NM_000668.5", 125, "target"),
    "DHX58": ("NM_024119.3", 79132, "target"),
    "ALDH3A2": ("NM_000382.3", 224, "target"),
    "GBP1": ("NM_002053.3", 2633, "target"),
    "ACTB": ("NM_001101.3", 60, "reference"),
    "GAPDH": ("NM_002046.6", 2597, "reference"),
    "TBP": ("NM_003194.5", 6908, "reference"),
    "HPRT1": ("NM_000194.3", 3251, "reference"),
    "RPLP0": ("NM_001002.4", 6175, "reference"),
    "SDHA": ("NM_004168.4", 6389, "reference"),
}
RC = str.maketrans("ACGT", "TGCA")


def efetch(acc):
    u = f"{EUTILS}/efetch.fcgi?db=nuccore&id={acc}&rettype=fasta&retmode=text"
    t = urllib.request.urlopen(u, timeout=120).read().decode()
    lines = t.split("\n")
    header = lines[0].lstrip(">")
    seq = "".join(l.strip() for l in lines[1:]
                  if l.strip() and not l.startswith(">")).upper()
    return header, seq


def main():
    qp = pd.read_csv(f"{OUT}/qpcr_primers.csv")
    out_seq = {}
    rows = []
    for _, r in qp.iterrows():
        g = r["gene"]
        # ★ accession 以主表 qpcr_primers.csv 为准：引物重设计后主表会换 accession
        #   （如 LDHA 由 NM_001135239.1 改为 NM_005566.4 的 3'UTR 锚定），
        #   若继续用下方硬编码表 efetch，会在**错误模板**上落位并给出假 OK。
        #   硬编码表降级为兜底（主表缺该列/该单元格为空时使用）。
        acc_fb, gid_fb, role_fb = DESIGN_ACC[g]
        acc = str(r.get("refseq_accession") or acc_fb).strip() or acc_fb
        gid_raw = str(r.get("gene_id", "")).strip()
        gid = int(float(gid_raw)) if gid_raw.replace(".", "", 1).isdigit() else gid_fb
        role = str(r.get("role") or role_fb).strip() or role_fb
        header, seq = efetch(acc)
        F = r["forward_primer_5to3"]
        R = r["reverse_primer_5to3"]
        rcs = R.translate(RC)[::-1]
        fp = seq.find(F)
        rp = seq.find(rcs)
        amp = (rp + len(R) - fp) if (fp >= 0 and rp >= 0) else None
        out_seq[g] = {
            "acc": header,
            "seq": seq,
            "header": header,
            "design_accession": acc,
            "gene_id": gid,
            "role": role,
        }
        rows.append({
            "gene": g, "role": role, "design_accession": acc,
            "gene_id": gid, "header": header,
            "mRNA_length": len(seq),
            "forward_primer_5to3": F, "forward_start_1based": fp + 1 if fp >= 0 else -1,
            "reverse_primer_5to3": R,
            "reverse_rc_start_1based": rp + 1 if rp >= 0 else -1,
            "amplicon_length": amp,
            "both_primers_on_template": bool(fp >= 0 and rp >= 0),
        })
        print(f"{g:8s} {acc:15s} F@{fp+1 if fp>=0 else -1:5d} "
              f"Rrc@{rp+1 if rp>=0 else -1:5d} amp={amp}  "
              f"{'OK' if fp>=0 and rp>=0 else 'FAIL'}")
        time.sleep(0.4)

    with open(f"{PROV}/refseq_mrna_sequences.json", "w", encoding="utf-8") as fh:
        json.dump(out_seq, fh, ensure_ascii=False, indent=1)
    df = pd.DataFrame(rows)
    df.to_csv(f"{PROV}/primer_template_mapping.csv", index=False)

    n_ok = int(df.both_primers_on_template.sum())
    print(f"\n溯源文件已修正；{n_ok}/{len(df)} 对引物实测落在设计模板上")
    print("-> provenance/refseq_mrna_sequences.json（修正版）")
    print("-> provenance/primer_template_mapping.csv")
    print("DONE — fix_refseq_provenance.py")


if __name__ == "__main__":
    main()
