# -*- coding: utf-8 -*-
"""
audit_v2_paralog_claims.py — 对 design_paralog_aware_primers_v2.py 自报结论的
独立裁决（穷举滑动窗口，零种子剪枝）

═══ 为什么必须独立裁决 ═══

v2 脚本自报：除 GBP1 外全部 12 条引物的
`max_identity_vs_paralog = 0.0`、`min_mismatch = 20`（恰等于引物长度）、
`worst_paralog = ''`。字面含义是「对所有旁系同源序列未找到任何可比的窗口」。

这与生物学事实直接冲突：
- FAHD2A / FAHD2B 的 CDS 仅 15 个 SNP（一致性 98.41%），蛋白 319 aa 逐字
  相同 → 任何 20-mer 窗口在 FAHD2A 上都应能找到 identity 85–100% 的匹配；
- ADH1B / ADH1A / ADH1C 的 CDS 一致性 0.9495–0.9566 → 同理。

所以「0.0 / mm=20」几乎不可能是真实测量结果，而是**算法漏检**。

═══ v2 算法的两处具体缺陷（本脚本逐一实证）═══

**缺陷 A — 分段种子索引漏检（漏检率取决于错配分布）**
`best_match` 把引物切成 4 段 × 5 nt，只有当某一段 5-mer 在旁系同源序列中
**精确出现**时，才把该命中扩展为候选窗口。若最佳窗口的错配**均匀分散到
全部 4 段**（每段 ≥1 个错配），则 4 个种子全部失配 → 整个最佳窗口被漏掉，
该旁系同源被记为「完全无匹配」。
   - 检出条件：错配数 ≤ 3 且未覆盖全部 4 段（鸽巢原理）
   - 漏检条件：错配数 ≥ 4 且覆盖全部 4 段

**缺陷 B — 反向引物在负链命中时 3' 端锚定测错端**
v2 对反向引物搜 `rc(primer)`，却把 3' 端错配数算在 `s[L-5:L]`。
但 `rc(primer)[0]` 才是引物 3' 端碱基的互补 → 引物自身的 3' 端对应
**窗口的前 5 nt**，不是后 5 nt。于是：
   - 方向性偏差：真正 3' 端有错配者被误判为「3' 端干净」→ 漏报 P0
     （危险方向）
   - 反向偏差：3' 端干净者被误判为「有错配」→ 误杀可用引物

**缺陷 C — 模板用整条 mRNA 而非 CDS**
v2 的 `SEQUENCE_TEMPLATE = best_t["seq"]`（全长 mRNA），违反
`redesign_primers_v2.py` 已确立的 **CDS 锚定**原则 → 选出的引物可能落在
3'UTR（ALDH3A2 的旧缺陷正是 3'UTR Alu 重复元件）。

═══ 本脚本的裁决方法（与 v2 零代码共享）═══

- **穷举滑动窗口**：对每个引物 × 每条旁系同源序列 × 正反链，枚举**全部**
  起点位置，用 numpy 向量化计算每个窗口的错配数，不做任何种子剪枝 →
  得到的 identity 一定是真实最优值，不存在漏检。
- **链向正确的 3' 端锚定**：
    - 正链命中（P 直接出现在 S 中）→ 引物 3' 端 = 窗口**后** 5 nt
    - 负链命中（rc(P) 出现在 S 中）→ 引物 3' 端 = 窗口**前** 5 nt
- **两个比较面**：
    - `full`：旁系同源**全长转录本**（与 BLAST 行为一致，保守）
    - `cds`：旁系同源**编码区**（RT-qPCR 语境下真正会共扩增的面）
- **落位核查**：从 GenBank CDS 特征取真实坐标，判定引物是否在 CDS 内。

输出：
  03_结果/实验验证设计/paralog_claim_audit.csv   逐引物裁决明细
  03_结果/实验验证设计/paralog_claim_audit.json   汇总 + v2 声明对照
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
CACHE_P = f"{PROV}/paralog_seq_cache.json"

RC_MAP = str.maketrans("ACGT", "TGCA")
B = np.frombuffer(b"ACGT", dtype=np.uint8)
ENC = np.full(256, 4, dtype=np.uint8)
for i, ch in enumerate(b"ACGT"):
    ENC[ch] = i


def rc(s):
    return s.translate(RC_MAP)[::-1]


def enc(s):
    return ENC[np.frombuffer(s.encode("ascii"), dtype=np.uint8)]


def windows_mm(subject: np.ndarray, primer: np.ndarray):
    """返回 (mm 数组, 窗口数)。mm[i] = subject[i:i+L] vs primer 的错配数。"""
    M, L = len(subject), len(primer)
    if M < L:
        return None, 0
    n = M - L + 1
    mm = np.zeros(n, dtype=np.int16)
    for k in range(L):
        mm += (subject[k:k + n] != primer[k])
    return mm, n


def exhaustive_best(primer: str, subject: str):
    """穷举全部窗口 × 正反链 → (max_identity%, mm_total, mm_3p5, strand)

    3' 端锚定按**引物自身方向**判定：
      - 正链（P 出现在 S 中）：引物 3' 端 = 窗口后 5 nt
      - 负链（rc(P) 出现在 S 中）：引物 3' 端 = 窗口前 5 nt
    """
    L = len(primer)
    best = (0.0, L, L, "none")
    S = enc(subject)
    for strand, P in (("+", enc(primer)), ("-", enc(rc(primer)))):
        mm, n = windows_mm(S, P)
        if n == 0:
            continue
        i = int(np.argmin(mm))
        m = int(mm[i])
        if m >= L:
            continue
        ident = 100.0 * (L - m) / L
        if ident <= best[0]:
            continue
        # 逐位置错配（必须在**最佳窗口内部**计算，不能跨窗口求和）
        pos_mm = (S[i:i + L] != P)
        if strand == "+":
            mm3 = int(pos_mm[L - 5:L].sum())
        else:
            mm3 = int(pos_mm[0:5].sum())
        best = (ident, m, mm3, strand)
    return best


def cds_of(gb_text):
    """从 GenBank 记录取 CDS 坐标（1-based inclusive）。

    ═══ 修正说明 ═══
    旧实现把 **qualifier 续行**（如 `/note="...67kDa..."`）拼进 location
    字符串 → `69..1847` + `67` → 误解析为 `69..184767`（超出 mRNA 长度
    2883）。正确做法：location 只取特征块首行 + 紧跟其后、不以 `/` 开头的
    续行；一旦遇到 qualifier（`/` 开头）即停止拼接。
    """
    for b in re.split(r"\n {5}(?=\w)", gb_text):
        if not b.startswith("CDS"):
            continue
        loc = re.match(r"CDS\s+(\S+)", b).group(1).replace("<", ">")
        for ln in b.split("\n")[1:]:
            st = ln.strip()
            if not st:
                continue
            if st.startswith("/"):
                break                     # 进入 qualifier 区 → 停止
            loc += st
        return [(int(a), int(c)) for a, c in
                re.findall(r"(\d+)\.\.(\d+)", loc)]
    return []


def in_cds(pos0, length, cds):
    """引物区间 [pos0, pos0+length)（0-based）是否完整位于 CDS 内。"""
    if not cds:
        return None
    covered = set()
    for a, b in cds:
        covered.update(range(a - 1, b))
    return all(p in covered for p in range(pos0, pos0 + length))


def main():
    C = json.load(open(CACHE_P, encoding="utf-8"))
    v2 = json.load(open(f"{OUT}/paralog_aware_design_v2.json", encoding="utf-8"))
    rows = []
    summary = {}

    for gene, gv in v2["genes"].items():
        sel = gv.get("selected")
        if not sel:
            continue
        print(f"\n=== {gene} (v2 目标 {gv['target_accession']}, "
              f"CDS {gv['target_cds_len']} nt) ===")

        # 目标落位
        tgt_acc = gv["target_accession"]
        tgt_seq = C.get(tgt_acc, {}).get("seq", "")
        gb = C.get(f"__gb__{tgt_acc}", {}).get("seq", "")
        t_cds = cds_of(gb) if gb else []

        # 旁系同源序列（全部转录本）
        para_full = []       # (gene, acc, seq)
        para_cds = []
        for fam in gv["paralog_family"]:
            for pgid, acc, title in C.get(f"__tx__{fam}", {}).get("seq", []):
                s = C.get(acc, {}).get("seq")
                if not s:
                    continue
                para_full.append((fam, acc, s))
                g = C.get(f"__gb__{acc}", {}).get("seq")
                if g:
                    seg = cds_of(g)
                    if seg:
                        cds = "".join(s[a - 1:b] for a, b in seg)
                        para_cds.append((fam, acc, cds))
        print(f"   旁系同源：全长转录本 {len(para_full)} 条，"
              f"含 CDS 定义 {len(para_cds)} 条")

        entry = {"gene": gene, "target_accession": tgt_acc,
                 "v2_claim": {}, "audit": {}}

        for lab, key in (("forward", "forward"), ("reverse", "reverse")):
            p = sel[key]
            L = len(p)
            # v2 自报
            claim_id = sel[f"{lab[0]}_max_identity_vs_paralog"]
            claim_mm = sel[f"{lab[0]}_min_mismatch"]
            claim_mm3 = sel[f"{lab[0]}_mm_3p5"]
            claim_worst = sel[f"{lab[0]}_worst_paralog"]

            # 独立裁决
            res_full = {"best": None}
            for fam, acc, s in para_full:
                ident, mm, mm3, strand = exhaustive_best(p, s)
                if res_full["best"] is None or ident > res_full["best"][0]:
                    res_full["best"] = (ident, mm, mm3, strand, fam, acc)
            res_cds = {"best": None}
            for fam, acc, s in para_cds:
                ident, mm, mm3, strand = exhaustive_best(p, s)
                if res_cds["best"] is None or ident > res_cds["best"][0]:
                    res_cds["best"] = (ident, mm, mm3, strand, fam, acc)

            bf = res_full["best"]
            bc = res_cds["best"]

            # 目标自比对（引物是否精确出现在目标转录本中、落在何处）
            self_pos = tgt_seq.find(p)
            self_rc = tgt_seq.find(rc(p))
            inside = in_cds(self_pos if self_pos >= 0 else self_rc, L, t_cds)

            def grade(ident, mm3):
                if ident <= 90.0 and mm3 >= 1:
                    return "STRONG"
                if ident <= 95.0 and mm3 >= 1:
                    return "WEAK"
                if ident <= 90.0:
                    return "STRONG"
                return "UNACCEPTABLE"

            g_full = grade(bf[0], bf[2]) if bf else "NA"
            g_cds = grade(bc[0], bc[2]) if bc else "NA"

            print(f"   {lab:8s} {p}")
            print(f"      v2  自报 : id={claim_id:.1f}% mm={claim_mm} "
                  f"mm3'={claim_mm3} worst={claim_worst!r}")
            if bf:
                print(f"      独立(full): id={bf[0]:.1f}% mm={bf[1]} "
                      f"mm3'={bf[2]} strand={bf[3]} worst={bf[4]}({bf[5]}) "
                      f"→ {g_full}")
            if bc:
                print(f"      独立(cds) : id={bc[0]:.1f}% mm={bc[1]} "
                      f"mm3'={bc[2]} strand={bc[3]} worst={bc[4]}({bc[5]}) "
                      f"→ {g_cds}")
            print(f"      目标落位 : self_pos={self_pos} rc_pos={self_rc} "
                  f"in_CDS={inside}")

            rows.append({
                "gene": gene, "primer": lab, "query": p, "length": L,
                "v2_claim_identity_pct": claim_id,
                "v2_claim_min_mismatch": claim_mm,
                "v2_claim_mm_3p5": claim_mm3,
                "v2_claim_worst_paralog": claim_worst,
                "audit_full_identity_pct": round(bf[0], 2) if bf else None,
                "audit_full_mm_total": bf[1] if bf else None,
                "audit_full_mm_3p5": bf[2] if bf else None,
                "audit_full_strand": bf[3] if bf else None,
                "audit_full_worst_gene": bf[4] if bf else None,
                "audit_full_worst_acc": bf[5] if bf else None,
                "audit_full_grade": g_full,
                "audit_cds_identity_pct": round(bc[0], 2) if bc else None,
                "audit_cds_mm_total": bc[1] if bc else None,
                "audit_cds_mm_3p5": bc[2] if bc else None,
                "audit_cds_worst_gene": bc[4] if bc else None,
                "audit_cds_grade": g_cds,
                "target_self_pos_1based": (self_pos + 1) if self_pos >= 0 else None,
                "target_in_CDS": inside,
            })
            entry["v2_claim"][lab] = {
                "identity_pct": claim_id, "min_mismatch": claim_mm,
                "mm_3p5": claim_mm3, "worst_paralog": claim_worst}
            entry["audit"][lab] = {
                "full": {"identity_pct": round(bf[0], 2) if bf else None,
                         "mm_total": bf[1] if bf else None,
                         "mm_3p5": bf[2] if bf else None,
                         "strand": bf[3] if bf else None,
                         "worst_gene": bf[4] if bf else None,
                         "worst_acc": bf[5] if bf else None,
                         "grade": g_full},
                "cds": {"identity_pct": round(bc[0], 2) if bc else None,
                        "mm_total": bc[1] if bc else None,
                        "mm_3p5": bc[2] if bc else None,
                        "worst_gene": bc[4] if bc else None,
                        "grade": g_cds},
                "target_in_CDS": inside,
                "target_self_pos_1based": (self_pos + 1) if self_pos >= 0 else None,
            }
        summary[gene] = entry

    import csv
    with open(f"{OUT}/paralog_claim_audit.csv", "w", encoding="utf-8",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    with open(f"{OUT}/paralog_claim_audit.json", "w", encoding="utf-8") as fh:
        json.dump({"date": "2026-10-25",
                   "script": "02_脚本/audit_v2_paralog_claims.py",
                   "method": {
                       "match": "穷举滑动窗口（numpy 向量化，零种子剪枝）",
                       "strand_rule": ("正链命中→引物3'端=窗口后5nt；"
                                       "负链命中(rc(P)出现于S)→引物3'端=窗口前5nt"),
                       "comparison_faces": "full=旁系同源全长转录本; cds=旁系同源编码区",
                   },
                   "genes": summary},
                  fh, ensure_ascii=False, indent=1, default=str)
    print(f"\nDONE — {len(rows)} 条引物裁决 → {OUT}/paralog_claim_audit.csv")


if __name__ == "__main__":
    main()
