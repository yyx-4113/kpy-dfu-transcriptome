# -*- coding: utf-8 -*-
"""
design_ldha_v5_utr3.py — LDHA 引物 v5 重设计（3'UTR 锚定，靶向规范转录本 NM_005566.4）

═══ 决策依据（阻塞 #22/#23 的用户决策，2026-10-25）═══

LDHA 无法用单一 CDS 锚定引物对做到跨转录本位点特异：
  - CDS 共同唯一区段 = 0 nt（ldha_unique_region_analysis.py 实证）
  - LOC105375980（9p22.3 lncRNA）携带 LDHA CDS 808 nt / 96.16% 拷贝，
    任一 20-mer 约 46% 概率与 LOC 完全一致（verify_ldha_lncrna_overlap.py
    双路径交叉验证：Biopython SW + 自写 numpy 对角线精确评分）

→ 用户选择方案 ①：**3'UTR 锚定、靶向规范转录本 NM_005566.4**
  - 避开 LOC105375980 覆盖的 74.1% 区域（5'UTR + CDS + 近侧 3'UTR）；
  - 仅在 NM_005566.4 的 **3'UTR 唯一区段**（由已验证掩码动态派生）内设计；
  - 透明披露：本对仅定量规范转录本 NM_005566.4，偏离"CDS 锚定 +
    全编码转录本共有"范式（其余 11 对保持该范式）。

═══ 设计约束（与 v4 其余 11 对口径统一）═══

  - Primer3 参数（Tm 58–62、GC 40–60、产物 80–180、二级结构 ≤47°C）
    与 design_final_primers_v4.py / redesign_ldha_primers.py 一致；
  - 1667–2047 段整体 GC 31.76%（偏低），该段允许 GC 放宽至 30–60；
    2070–2241 段 GC 43%（达标），优先用标准 40–60；
  - 扩增子**完整落在** 3'UTR 唯一区段内（两条引物 3' 端均在唯一区段，
    LOC 不携带该序列 → 不可被 LOC 共扩增）；
  - 设计后独立 numpy 复核：① 引物在 NM_005566.4 的落位与方向；
    ② 两条引物均**不出现**于 LOC105375980（含正/反链）→ 确证避开 LOC；
    ③ 两条引物坐标均 > CDS 终止（确证在 3'UTR）。

输出：
  03_结果/实验验证设计/ldha_v5_utr3_candidates.json
"""
import json

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import os
import numpy as np
import primer3

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
PROV = f"{OUT}/provenance"
CANON = "NM_005566.4"
LOC = "LOC105375980"

RC_MAP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(RC_MAP)[::-1]


def safe_segments(covered_mask, off=0):
    """covered_mask: True=被覆盖（坏）。返回未被覆盖（安全）的连续区段(1-based)。"""
    out, st = [], None
    for i in range(len(covered_mask)):
        if not covered_mask[i] and st is None:
            st = i
        elif covered_mask[i] and st is not None:
            out.append([st + 1 + off, i + off])
            st = None
    if st is not None:
        out.append([st + 1 + off, len(covered_mask) + off])
    return out


def gc(s):
    if not s:
        return 0.0
    return round((s.count("G") + s.count("C")) * 100.0 / len(s), 2)


def main():
    # ---------- 1. 取 NM_005566.4 序列 + CDS + 已验证掩码 ----------
    cache = json.load(open(f"{PROV}/final_design_cache.json", encoding="utf-8"))
    seq = cache[CANON]["seq"].upper()
    # CDS 来自已验证的掩码分解结果（per_ldha_transcript NM_005566.4 cds_seg）
    decomp = json.load(open(f"{PROV}/ldha_mask_decomposition.json",
                            encoding="utf-8"))
    cds_seg = decomp["steps"]["per_ldha_transcript"][CANON]["cds_seg"]
    cds_end = cds_seg[-1][1]            # 1097（1-based）
    cmask = np.load(f"{PROV}/_ldha_curated_mask_{CANON}.npy")  # True=被LOC+LDHAL6B覆盖
    n = len(seq)
    assert len(cmask) == n

    safe = safe_segments(cmask)         # 动态派生，1-based
    # 3'UTR 唯一区段（起点 > CDS 终止）
    utr3_safe = [a for a in safe if a[0] > cds_end]
    print(f"NM_005566.4: len={n}, CDS={cds_seg}, CDS_end={cds_end}")
    print(f"动态派生的安全区段: {safe}")
    print(f"3'UTR 唯一区段: {utr3_safe}")

    # LOC 序列（已下载，避免网络）：gene_105375980.xml 含 NR_ 转录本
    # 为独立复核，我们从 NUCCORE 现取 LOC 的 NR_202040.1（最完整注释）
    loc_seq = None
    try:
        import urllib.request, urllib.parse, time
        E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        UA = {"User-Agent": "WorkBuddy-LDHAv5/1.0"}
        raw = urllib.request.urlopen(
            urllib.request.Request(
                f"{E}/efetch.fcgi?db=nuccore&id=NR_202040.1"
                f"&rettype=fasta&retmode=text", headers=UA),
            timeout=120).read().decode()
        loc_seq = "".join(l.strip() for l in raw.splitlines()[1:]
                          if l.strip() and not l.startswith("BASE")).upper()
        print(f"LOC105375980 (NR_202040.1) 现取长度 = {len(loc_seq)}")
        time.sleep(0.3)
    except Exception as e:                  # noqa: BLE001
        print(f"  ! LOC 序列现取失败（{e}），改用本地 gene_105375980.xml 解析")
        loc_seq = _loc_from_xml()

    # ---------- 2. Primer3 设计（在 3'UTR 唯一区段内）----------
    base_params = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_NUM_RETURN": 8,
        "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 25,
        "PRIMER_OPT_SIZE": 20, "PRIMER_MAX_POLY_X": 4,
        "PRIMER_MIN_TM": 58.0, "PRIMER_MAX_TM": 62.0, "PRIMER_OPT_TM": 60.0,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_MAX_SELF_ANY_TH": 47.0, "PRIMER_MAX_SELF_END_TH": 47.0,
        "PRIMER_MAX_HAIRPIN_TH": 47.0,
        "PRIMER_PAIR_MAX_COMPL_ANY_TH": 47.0, "PRIMER_PAIR_MAX_COMPL_END_TH": 47.0,
        "PRIMER_EXPLAIN_FLAG": 1,
    }
    policies = [("GC40_60", 40.0, 60.0, 80, 180),
                ("GC30_60", 30.0, 60.0, 80, 180)]

    candidates = []
    for (a, b) in utr3_safe:
        sub = seq[a - 1:b]               # 1-based→0-based 切片
        if len(sub) < 60:
            print(f"  区段 {a}-{b} 太短({len(sub)}nt)，跳过")
            continue
        sub_gc = gc(sub)
        for pname, gmin, gmax, pmin, pmax in policies:
            if pname == "GC40_60" and sub_gc < 36:
                continue                 # 整体 GC 过低，标准策略必败，直接跳
            params = dict(base_params)
            params["PRIMER_MIN_GC"] = gmin
            params["PRIMER_MAX_GC"] = gmax
            params["PRIMER_PRODUCT_SIZE_RANGE"] = [[pmin, pmax]]
            try:
                res = primer3.bindings.design_primers(
                    {"SEQUENCE_ID": f"NM_005566.4_UTR3_{a}_{b}_{pname}",
                     "SEQUENCE_TEMPLATE": sub}, params, {})
            except Exception as e:        # noqa: BLE001
                print(f"  Primer3 异常 {a}-{b}/{pname}: {e}")
                continue
            npairs = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
            for i in range(npairs):
                f = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
                r = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
                lp = res[f"PRIMER_LEFT_{i}"]
                rp = res[f"PRIMER_RIGHT_{i}"]
                f0 = a - 1 + lp[0]        # 0-based 绝对
                r0 = a - 1 + rp[0]
                cand = {
                    "policy": pname, "segment": [a, b], "sub_gc": sub_gc,
                    "forward": f, "reverse": r,
                    "f_tm": res[f"PRIMER_LEFT_{i}_TM"],
                    "r_tm": res[f"PRIMER_RIGHT_{i}_TM"],
                    "f_gc": round(100 * sum(c in "GC" for c in f) / len(f), 1),
                    "r_gc": round(100 * sum(c in "GC" for c in r) / len(r), 1),
                    "product_size": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                    "f_abs_1based": f0 + 1, "f_end_1based": f0 + len(f),
                    "r_abs_1based": r0 + 1, "r_end_1based": r0 + len(r),
                    "f_self_any_th": res.get(f"PRIMER_LEFT_{i}_SELF_ANY_TH"),
                    "f_self_end_th": res.get(f"PRIMER_LEFT_{i}_SELF_END_TH"),
                    "f_hairpin_th": res.get(f"PRIMER_LEFT_{i}_HAIRPIN_TH"),
                    "r_self_any_th": res.get(f"PRIMER_RIGHT_{i}_SELF_ANY_TH"),
                    "r_self_end_th": res.get(f"PRIMER_RIGHT_{i}_SELF_END_TH"),
                    "r_hairpin_th": res.get(f"PRIMER_RIGHT_{i}_HAIRPIN_TH"),
                    "pair_compl_any_th": res.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH"),
                    "pair_compl_end_th": res.get(f"PRIMER_PAIR_{i}_COMPL_END_TH"),
                    "penalty": res.get(f"PRIMER_PAIR_{i}_PENALTY"),
                }
                candidates.append(cand)
            if npairs:
                print(f"  区段 {a}-{b}/{pname}: {npairs} 对候选 "
                      f"(GC{sub_gc}%)")
        if candidates and any(c["segment"] == [a, b]
                              for c in candidates):
            pass

    print(f"\n候选对总数 = {len(candidates)}")

    # ---------- 3. 独立复核 + 打分排序 ----------
    def in_loc(s):
        return (loc_seq.find(s) >= 0) or (loc_seq.find(rc(s)) >= 0)

    def score(c):
        s = c["penalty"]
        s += abs(c["f_tm"] - 60) + abs(c["r_tm"] - 60)
        s += abs(c["f_gc"] - 50) / 10 + abs(c["r_gc"] - 50) / 10
        s += 0 if c["forward"].endswith(("G", "C")) else 0.5
        s += 0 if c["reverse"].endswith(("G", "C")) else 0.5
        # 标准 GC 策略优先于放宽策略（1667-2047 段虽整体 GC 低，但其 primers
        # 本身仍在 40-60 标准窗口内，仅因区段整体 GC 偏低才放宽 floor）
        s += 3 if c["policy"] == "GC40_60" else 0
        # 偏好更「内部」的 3'UTR 扩增子（距 polyA 越远越不易受不完全逆转录影响）：
        # 段长越长（1667-2047=380nt）越内部，给段长加负分（优先）
        s -= (c["segment"][1] - c["segment"][0]) * 0.02
        return s

    for c in candidates:
        # 独立复核 1：NM_005566.4 落位与方向
        iF = seq.find(c["forward"])
        iRrc = seq.find(rc(c["reverse"]))
        c["ok_in_canon"] = (iF >= 0 and iRrc >= 0)
        c["f_found"] = iF + 1
        c["r_rc_found"] = iRrc + 1
        # 独立复核 2：避开 LOC（正/反链均不出现）
        c["f_in_loc"] = in_loc(c["forward"])
        c["r_in_loc"] = in_loc(c["reverse"])
        c["avoids_loc"] = (not c["f_in_loc"]) and (not c["r_in_loc"])
        # 独立复核 3：在 3'UTR（坐标 > CDS 终止）
        c["in_utr3"] = (c["f_abs_1based"] > cds_end
                         and c["r_end_1based"] > cds_end)
        # 扩增子完整落在唯一区段（两条引物均在 safe 段内）
        c["amplicon_in_safe"] = (c["f_abs_1based"] >= c["segment"][0]
                                 and c["r_end_1based"] <= c["segment"][1])

    # 过滤：必须通过三项独立复核
    valid = [c for c in candidates if c["ok_in_canon"] and c["avoids_loc"]
             and c["in_utr3"] and c["amplicon_in_safe"]]
    print(f"通过三项独立复核的候选 = {len(valid)}")
    for c in valid[:6]:
        print(f"  [{c['policy']}] F={c['forward']} R={c['reverse']} "
              f"amp={c['product_size']} Tm={c['f_tm']:.1f}/{c['r_tm']:.1f} "
              f"GC={c['f_gc']}/{c['r_gc']} seg={c['segment']} "
              f"F@loc={c['f_in_loc']} R@loc={c['r_in_loc']}")

    valid.sort(key=score)
    best = valid[0] if valid else None

    out = {
        "date": "2026-10-25",
        "script": "02_脚本/design_ldha_v5_utr3.py",
        "decision": "方案①：3'UTR 锚定，靶向规范转录本 NM_005566.4",
        "template": CANON,
        "cds_seg": cds_seg, "cds_end_1based": cds_end,
        "template_len": n,
        "curated_mask_covered_pct": round(100 * int(cmask.sum()) / n, 1),
        "dynamic_safe_segments_1based": safe,
        "utr3_unique_segments_1based": utr3_safe,
        "loc_accession_used": "NR_202040.1",
        "loc_len": len(loc_seq) if loc_seq else None,
        "primer3_version": primer3.__version__,
        "n_candidates": len(candidates),
        "n_valid_pass_audit": len(valid),
        "candidates": candidates,
        "selected": best,
    }
    with open(f"{OUT}/ldha_v5_utr3_candidates.json", "w",
              encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1, default=str)

    print("\n=== 选定对（v5, 3'UTR 锚定 NM_005566.4）===")
    if best:
        print(f"  F = {best['forward']}  ({best['f_abs_1based']}.."
              f"{best['f_end_1based']}, Tm {best['f_tm']:.2f}, GC {best['f_gc']}%)")
        print(f"  R = {best['reverse']}  ({best['r_abs_1based']}.."
              f"{best['r_end_1based']}, Tm {best['r_tm']:.2f}, GC {best['r_gc']}%)")
        print(f"  产物 = {best['product_size']} bp  (模板段 {best['segment']})")
        print(f"  避开 LOC: F_in_LOC={best['f_in_loc']} R_in_LOC="
              f"{best['r_in_loc']} → avoids_LOC={best['avoids_loc']}")
        print(f"  在 3'UTR: {best['in_utr3']}  扩增子完整在唯一段: "
              f"{best['amplicon_in_safe']}")
    else:
        print("  !! 未找到通过复核的候选——需进一步放宽或换策略")
    print("\nDONE — design_ldha_v5_utr3.py")


def _loc_from_xml():
    """从已下载的 gene_105375980.xml 解析 NR_202040.1 序列（兜底）。"""
    import re
    p = f"{PROV}/gene_105375980.xml"
    if not os.path.exists(p):
        return None
    t = open(p, encoding="utf-8", errors="ignore").read()
    # 找 NR_202040.1 的 GBSeq_sequence
    m = re.search(r"NR_202040\.1.*?(?:GBSeq_sequence|>)([ACGTN]+)",
                  t, re.S)
    if m:
        return m.group(1).upper()
    return None


if __name__ == "__main__":
    main()
