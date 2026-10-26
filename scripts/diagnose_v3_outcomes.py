# -*- coding: utf-8 -*-
"""
diagnose_v3_outcomes.py — v3 结果的三项专项诊断（全部基于已有缓存，不联网）

═══ 诊断项 ═══

D1 **FAHD2B 0 通过的根因**：共有 CDS 区段（315..776）内是否存在可区分
   FAHD2A 的窗口？→ 逐窗口计算对 FAHD2A 全长/编码区的最佳 identity，
   给出「可区分窗口数」与「3'端有错配的窗口数」。

D2 **curated vs predicted 覆盖分离**：RefSeq 转录本分 NM_/NR_（已策展，
   真实剪接变体）与 XM_（算法预测）。要求「全部蛋白编码转录本」共有
   在 29 条 FAHD2B 编码转录本下可能过严（多数为 XM_）。→ 分别统计
   curated-only 与 all-transcript 的覆盖与共有区段长度。

D3 **重算 mm3 与重分级**：v3 首跑使用了修复前的 mm3 实现（跨窗口求和，
   产生 67/80/84 这类不可能的值）与旧版 grade()。→ 以缓存重算全部候选的
   真实 mm3 与新分级，判定各基因最终通过数是否改变。

输出：
  03_结果/实验验证设计/provenance/v3_outcome_diagnostic.json
"""
import json

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import os
import re

import numpy as np

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
PROV = f"{OUT}/provenance"
RC_MAP = str.maketrans("ACGT", "TGCA")

ENC = np.full(256, 4, dtype=np.uint8)
for _i, _ch in enumerate(b"ACGT"):
    ENC[_ch] = _i


def rc(s):
    return s.translate(RC_MAP)[::-1]


def enc(s):
    return ENC[np.frombuffer(s.encode("ascii"), dtype=np.uint8)]


def best_window(subject: str, primer: str):
    """穷举全部窗口 × 正反链 → (identity%, mm_total, mm_3p5, strand)

    3' 端锚定按**引物自身方向**：
      - 正链（P 出现于 S）：引物 3' 端 = 窗口后 5 nt
      - 负链（rc(P) 出现于 S）：引物 3' 端 = 窗口前 5 nt
    """
    L = len(primer)
    S = enc(subject)
    if len(S) < L:
        return (0.0, L, L, "none")
    best = (0.0, L, L, "none")
    n = len(S) - L + 1
    for strand, P in (("+", enc(primer)), ("-", enc(rc(primer)))):
        mm = np.zeros(n, dtype=np.int16)
        for k in range(L):
            mm += (S[k:k + n] != P[k])
        i = int(np.argmin(mm))
        m = int(mm[i])
        if m >= L:
            continue
        ident = 100.0 * (L - m) / L
        if ident <= best[0]:
            continue
        pos_mm = (S[i:i + L] != P)
        mm3 = int(pos_mm[L - 5:L].sum()) if strand == "+" else int(pos_mm[0:5].sum())
        best = (ident, m, mm3, strand)
    return best


def grade(ident, mm3):
    """与 BLAST v5 判据严格对齐（修复后版本）。"""
    if mm3 == 0 and ident >= 95.0:
        return "UNACCEPTABLE"
    if mm3 == 0 and ident >= 90.0:
        return "UNACCEPTABLE"
    if ident >= 90.0:
        return "WEAK"
    return "STRONG"


def parse_cds_location(gb):
    for b in re.split(r"\n {5}(?=\w)", gb):
        if not b.startswith("CDS"):
            continue
        loc = re.match(r"CDS\s+(\S+)", b).group(1).replace("<", ">")
        for ln in b.split("\n")[1:]:
            st = ln.strip()
            if not st:
                continue
            if st.startswith("/"):
                break
            loc += st
        return [(int(a), int(c)) for a, c in
                re.findall(r"(\d+)\.\.(\d+)", loc)]
    return []


def main():
    C = json.load(open(f"{PROV}/final_design_cache.json", encoding="utf-8"))
    V = json.load(open(f"{OUT}/final_primer_design_v3.json", encoding="utf-8"))
    rep = {"date": "2026-10-25",
           "script": "02_脚本/diagnose_v3_outcomes.py"}

    # ═══ D1：FAHD2B 共有 CDS 区段的可区分性 ═══
    print("═══ D1 FAHD2B 共有 CDS 区段 vs FAHD2A 逐窗口可区分性 ═══")
    g = V["genes"]["FAHD2B"]
    canon = g["canonical"]
    ref = C[canon]["seq"]
    segs = g["shared_cds_segments"]
    seg_a, seg_b = segs[0][0] - 1, segs[0][0] - 1 + segs[0][1]
    print(f"   canonical={canon} len={len(ref)} 共有CDS区段={segs}")

    # FAHD2A / FAHD1 序列
    # ═══ 只取**蛋白编码**转录本 ═══
    # 理由：RT-qPCR 的底物是 cDNA —— 非编码 RNA（XR_/NR_ non-coding）
    # 与预测转录本（XM_）若不产生可扩增的转录本，不构成共扩增风险。
    # （v3 首跑把 51 条 FAHD2A RefSeq RNA 全算入 → 混入了 XR_ 非编码记录，
    #   夸大了脱靶面。）
    fam_seqs = {}
    for fam, gid_key in (("FAHD2A", "51011"), ("FAHD1", "81889")):
        accs = C.get(f"__tx__{gid_key}", {}).get("seq", [])
        seqs = []
        for a in accs:
            s = C.get(a, {}).get("seq")
            g = C.get(f"__gb__{a}", {}).get("seq")
            if not s:
                continue
            has_cds = bool(g and parse_cds_location(g))
            seqs.append((a, s, has_cds))
        fam_seqs[fam] = seqs
        n_cod = sum(1 for _, _, c in seqs if c)
        print(f"   {fam}: {len(seqs)} 条 RefSeq RNA，其中含 CDS {n_cod} 条")

    # 逐 20-mer 窗口计算对 FAHD2A 的最佳 identity
    # ═══ 修正说明（2026-10-25）═══
    # 旧版把 `worst_id` 初始化为 100.0 且仅在新值**更大**时更新
    # → 由于 identity 上界恰为 100.0，该变量恒为 100.0 → 所有窗口都被
    #    误报为「完全不可区分」，产生**假阴性**结论。
    # 正确：初始化为 -1.0（低于任何合法 identity），并取所有旁系同源
    # 转录本中的**最大** identity（即最保守的脱靶风险估计）。
    win_stats = {"n_win": 0, "n_le90": 0, "n_le90_3pmm": 0,
                 "n_le95": 0, "hardest": None, "easiest": None}
    L = 20
    for i in range(seg_a, seg_b - L + 1):
        w = ref[i:i + L]
        worst_id, worst_mm3 = -1.0, 0
        for fam, accs in fam_seqs.items():
            for acc, s, has_cds in accs:
                if not has_cds:
                    continue          # 仅蛋白编码转录本构成共扩增风险
                ident, mm, mm3, strand = best_window(s, w)
                if ident > worst_id:
                    worst_id, worst_mm3 = ident, mm3
        if worst_id < 0:
            continue
        win_stats["n_win"] += 1
        if worst_id <= 90.0:
            win_stats["n_le90"] += 1
            if worst_mm3 >= 1:
                win_stats["n_le90_3pmm"] += 1
        if worst_id <= 95.0:
            win_stats["n_le95"] += 1
        if win_stats["hardest"] is None or worst_id < win_stats["hardest"][0]:
            win_stats["hardest"] = (worst_id, worst_mm3, i + 1, w)
        if win_stats["easiest"] is None or worst_id > win_stats["easiest"][0]:
            win_stats["easiest"] = (worst_id, worst_mm3, i + 1, w)
    print(f"   20-mer 窗口总数 {win_stats['n_win']}")
    print(f"   对 FAHD2A/FAHD1（仅编码转录本）identity ≤90% 的窗口："
          f"{win_stats['n_le90']}（其中 3'端有错配 {win_stats['n_le90_3pmm']}）")
    print(f"   identity ≤95% 的窗口：{win_stats['n_le95']}")
    print(f"   最难区分窗口 id={win_stats['hardest'][0]:.1f}% mm3={win_stats['hardest'][1]} "
          f"pos={win_stats['hardest'][2]} {win_stats['hardest'][3]}")
    print(f"   最易区分窗口 id={win_stats['easiest'][0]:.1f}% mm3={win_stats['easiest'][1]} "
          f"pos={win_stats['easiest'][2]} {win_stats['easiest'][3]}")
    rep["D1_fahd2b_window_distinguishability"] = {
        "canonical": canon,
        "shared_cds_segment": segs,
        "comparison_face": "仅蛋白编码转录本（RT-qPCR 底物为 cDNA）",
        "n_windows": win_stats["n_win"],
        "n_windows_identity_le90": win_stats["n_le90"],
        "n_windows_identity_le90_and_3p_mismatch": win_stats["n_le90_3pmm"],
        "n_windows_identity_le95": win_stats["n_le95"],
        "hardest_to_distinguish_window": {
            "identity_pct": win_stats["hardest"][0],
            "mm_3p5": win_stats["hardest"][1],
            "pos_1based": win_stats["hardest"][2],
            "seq": win_stats["hardest"][3]},
        "easiest_to_distinguish_window": {
            "identity_pct": win_stats["easiest"][0],
            "mm_3p5": win_stats["easiest"][1],
            "pos_1based": win_stats["easiest"][2],
            "seq": win_stats["easiest"][3]},
    }

    # ═══ D2：curated vs predicted 覆盖分离 ═══
    print("\n═══ D2 curated (NM_/NR_) vs predicted (XM_) 覆盖分离 ═══")
    d2 = {}
    for gene, gv in V["genes"].items():
        canon = gv["canonical"]
        coding = gv["coding_transcripts"]
        curated = [a for a in coding if a.startswith(("NM_", "NR_"))]
        predicted = [a for a in coding if a.startswith("XM_")]
        # curated-only 共有 25-mer 区段
        ref = C[canon]["seq"]
        k = 25

        def shared_for(subset):
            common = {ref[i:i + k] for i in range(len(ref) - k + 1)}
            for a in subset:
                s = C.get(a, {}).get("seq")
                if not s:
                    continue
                common &= {s[i:i + k] for i in range(len(s) - k + 1)}
            pos = [i for i in range(len(ref) - k + 1) if ref[i:i + k] in common]
            if not pos:
                return []
            out, st, pv = [], pos[0], pos[0]
            for p in pos[1:]:
                if p <= pv + k:
                    pv = p
                else:
                    out.append((st + 1, pv + k - st))
                    st = pv = p
            out.append((st + 1, pv + k - st))
            return out

        seg_all = shared_for(coding)
        seg_cur = shared_for(curated) if curated else []
        print(f"   {gene:8s} coding={len(coding):2d} (curated {len(curated):2d} / "
              f"predicted {len(predicted):2d})  all-seg={[(a,b) for a,b in seg_all]}  "
              f"curated-seg={[(a,b) for a,b in seg_cur]}")
        d2[gene] = {"n_coding": len(coding), "n_curated": len(curated),
                    "n_predicted": len(predicted),
                    "shared_segments_all": seg_all,
                    "shared_segments_curated": seg_cur,
                    "curated_accs": curated}
    rep["D2_curated_vs_predicted"] = d2

    # ═══ D3：重算选中引物的真实 mm3 与新分级 ═══
    print("\n═══ D3 重算选中引物 mm3（修复后实现）与新分级 ═══")
    d3 = {}
    for gene, gv in V["genes"].items():
        sel = gv.get("selected")
        if not sel:
            continue
        canon = gv["canonical"]
        ref = C[canon]["seq"]
        # 该基因的旁系同源序列
        para = []
        for fam, gv2 in gv["paralog_families"].items():
            gid = None
            for k2, v2 in C.items():
                if k2.startswith("__gid__") and v2["seq"] and \
                        C.get(f"__tx__{v2['seq']}", {}).get("seq") and \
                        len(C[f"__tx__{v2['seq']}"]["seq"]) == gv2["n_refseq_rna"]:
                    pass
            gid = gv2["gene_id"]
            for acc in C.get(f"__tx__{gid}", {}).get("seq", []):
                s = C.get(acc, {}).get("seq")
                if s:
                    para.append((fam, acc, s))
        row = {"v3_selected": {"F": sel["forward"], "R": sel["reverse"]},
               "v3_reported_mm3": {"F": sel["f_mm3p_vs_paralog"],
                                   "R": sel["r_mm3p_vs_paralog"]},
               "v3_reported_id": {"F": sel["f_max_id_vs_paralog"],
                                  "R": sel["r_max_id_vs_paralog"]}}
        for lab, pk in (("F", "forward"), ("R", "reverse")):
            p = sel[pk]
            best = None
            for fam, acc, s in para:
                ident, mm, mm3, strand = best_window(s, p)
                if best is None or ident > best[0]:
                    best = (ident, mm, mm3, strand, fam, acc)
            # 目标自比对
            iF = ref.find(p)
            iR = ref.find(rc(p))
            row[lab] = {"primer": p, "target_pos_1based": (iF + 1) if iF >= 0
                        else ((iR + 1) if iR >= 0 else None),
                        "vs_paralog": ({"identity_pct": round(best[0], 2),
                                        "mm_total": best[1], "mm_3p5": best[2],
                                        "strand": best[3], "gene": best[4],
                                        "acc": best[5],
                                        "new_grade": grade(best[0], best[2])}
                                       if best else None)}
            print(f"   {gene:8s} {lab} {p}  v3报 mm3={row['v3_reported_mm3'][lab]} "
                  f"id={row['v3_reported_id'][lab]}  →  重算 mm3={best[2] if best else None} "
                  f"id={round(best[0],2) if best else None} "
                  f"grade={row[lab]['vs_paralog']['new_grade'] if best else 'NA'}")
        d3[gene] = row
    rep["D3_recomputed_mm3"] = d3

    with open(f"{PROV}/v3_outcome_diagnostic.json", "w", encoding="utf-8") as fh:
        json.dump(rep, fh, ensure_ascii=False, indent=1, default=str)
    print("\nDONE — diagnose_v3_outcomes.py")


if __name__ == "__main__":
    main()
