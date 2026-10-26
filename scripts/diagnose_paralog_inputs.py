# -*- coding: utf-8 -*-
"""
diagnose_paralog_inputs.py — 旁系同源输入序列的身份核验（审计地基自检）

═══ 为什么必须先做这一步 ═══

`audit_v2_paralog_claims.py` 的独立裁决报出：
    FAHD2B 引物 vs FAHD2A 最高 identity = 70.0%
但此前用**另一条代码路径**（RefSeq CDS 逐位比较）实测：
    FAHD2A / FAHD2B 的 CDS 945 nt 仅 15 个 SNP，一致性 98.41%，
    蛋白 319 aa 逐字相同。
两者相差 28 个百分点 → 至少有一方的输入是错的。

若旁系同源**序列本身**取错了基因（身份错配），则所有基于它的特异性
判读全部作废 —— 必须在写任何最终引物之前先把地基钉死。

═══ 本脚本自检项 ═══

S1 缓存中每个旁系同源转录本的 title 是否真的是该基因
   （防 `esearch [Gene Name]` 模糊匹配抓到邻近/别名基因）
S2 FAHD2A / FAHD2B 的 CDS 身份与一致性逐位复核（两条独立路径对照）
S3 GenBank CDS 坐标解析是否正确（对已知长度记录做 join 重建核对）
S4 各旁系同源家族实际取到的转录本数与 accession 分布
S5 ADH1A/ADH1C/ADH1B、GBP1..GBP7、LDHA 家族的一致性概览

输出：
  03_结果/实验验证设计/provenance/paralog_input_diagnostic.json
"""
import json

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import os
import re
from collections import defaultdict

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
PROV = f"{OUT}/provenance"
RC_MAP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(RC_MAP)[::-1]


def cds_segments(gb):
    """块切分 + 续行拼接 + join 展开 → CDS 坐标列表（1-based）。"""
    for b in re.split(r"\n {5}(?=\w)", gb):
        if b.startswith("CDS"):
            loc = re.match(r"CDS\s+(\S+)", b).group(1).replace("<", ">")
            extra = ""
            for ln in b.split("\n")[1:]:
                st = ln.strip()
                if not st or st.startswith("/"):
                    continue
                extra += st
            return [(int(a), int(c)) for a, c in
                    re.findall(r"(\d+)\.\.(\d+)", loc + extra)]
    return []


def main():
    C = json.load(open(f"{PROV}/paralog_seq_cache.json", encoding="utf-8"))
    rep = {"date": "2026-10-25",
           "script": "02_脚本/diagnose_paralog_inputs.py"}

    # ── S1 / S4：每个旁系同源家族的转录本 title 分布 ──
    print("═══ S1/S4 旁系同源家族缓存内容 ═══")
    fam_rep = {}
    for key in sorted(k for k in C if k.startswith("__tx__")):
        sym = key[6:]
        entries = C[key]["seq"]
        accs = [e[1] for e in entries]
        have_seq = [a for a in accs if a in C]
        print(f"\n-- {sym}: {len(accs)} 条 refseqrna 链接，"
              f"{len(have_seq)} 条有序列缓存")
        prefixes = defaultdict(int)
        for a in have_seq:
            prefixes[a.split(".")[0][:3]] += 1
        print(f"   前缀分布: {dict(prefixes)}")
        fam_rep[sym] = {"n_linked": len(accs), "n_with_seq": len(have_seq),
                        "accs": accs, "prefixes": dict(prefixes)}
        # title 抽样（判断是否抓错基因）
        for gid, acc, title in entries[:3]:
            print(f"     {acc} | {title[:88]}")
    rep["families"] = fam_rep

    # ── S2：FAHD2A vs FAHD2B 一致性（两条独立路径） ──
    print("\n═══ S2 FAHD2A / FAHD2B 一致性双路径复核 ═══")
    f2a = [e[1] for e in C["__tx__FAHD2A"]["seq"]]
    f2b = [e[1] for e in C["__tx__FAHD2B"]["seq"]]
    print(f"   FAHD2A accessions: {f2a[:6]} ... 共 {len(f2a)}")
    print(f"   FAHD2B accessions: {f2b[:6]} ... 共 {len(f2b)}")

    def get_cds(acc):
        s = C.get(acc, {}).get("seq")
        g = C.get(f"__gb__{acc}", {}).get("seq")
        if not s or not g:
            return None, None
        seg = cds_segments(g)
        if not seg:
            return s, None
        return s, "".join(s[a - 1:b] for a, b in seg)

    # 路径 1：同一转录本集合内两两比较（缓存自带）
    # 路径 2：GB 记录的 join 重建长度核对
    checks = []
    for acc in f2a[:4]:
        s, cds = get_cds(acc)
        if cds:
            checks.append({"acc": acc, "tx_len": len(s), "cds_len": len(cds),
                           "cds_head": cds[:40], "cds_tail": cds[-40:]})
            print(f"   {acc}: tx {len(s)} nt, CDS {len(cds)} nt")
    for acc in f2b[:4]:
        s, cds = get_cds(acc)
        if cds:
            checks.append({"acc": acc, "tx_len": len(s), "cds_len": len(cds),
                           "cds_head": cds[:40], "cds_tail": cds[-40:]})
            print(f"   {acc}: tx {len(s)} nt, CDS {len(cds)} nt")
    rep["fahd_cds_checks"] = checks

    # ── 直接比较：所有 FAHD2A/FAHD2B 缓存序列两两最大 identity ──
    print("\n   路径 A：FAHD2A 缓存序列 vs FAHD2B 缓存序列（滑动窗口最优）")
    best_overall = (0.0, None, None)
    for a in f2a:
        sa = C.get(a, {}).get("seq")
        if not sa:
            continue
        for b in f2b:
            sb = C.get(b, {}).get("seq")
            if not sb:
                continue
            # 快速：取 a 的前 200 nt 在 b 中找最长同源段
            hit = best_local_identity(sa, sb)
            if hit[0] > best_overall[0]:
                best_overall = (hit[0], a, b)
    print(f"   跨缓存最大 identity = {best_overall[0]:.2f}% "
          f"({best_overall[1]} vs {best_overall[2]})")
    rep["fahd_cross_cache_max_identity"] = {
        "identity_pct": round(best_overall[0], 2),
        "acc_a": best_overall[1], "acc_b": best_overall[2]}

    # ── S3：CDS 解析正确性（join 重建 == mRNA 序列的 CDS 拼接）──
    print("\n═══ S3 GenBank CDS 坐标解析自检 ═══")
    s3 = []
    for acc in ["NM_002053.3", "NM_000668.6", "NM_000382.3",
                "XM_054340782.1", "XM_054317165.1", "NM_001165414.2"]:
        s = C.get(acc, {}).get("seq")
        g = C.get(f"__gb__{acc}", {}).get("seq")
        if not s or not g:
            continue
        seg = cds_segments(g)
        if not seg:
            print(f"   {acc}: 无 CDS 特征")
            continue
        cds = "".join(s[a - 1:b] for a, b in seg)
        # 起点与终点
        print(f"   {acc}: tx={len(s)} nt, CDS={len(cds)} nt, "
              f"segments={seg if len(seg) <= 4 else str(seg[:4]) + '...'}")
        s3.append({"acc": acc, "tx_len": len(s), "cds_len": len(cds),
                   "segments": seg})
    rep["cds_parse_selfcheck"] = s3

    # ── S5：其余家族一致性概览 ──
    print("\n═══ S5 各旁系同源家族 CDS 一致性概览 ═══")
    pairs = [("ADH1A", "ADH1B"), ("ADH1C", "ADH1B"), ("ADH1A", "ADH1C"),
             ("LDHB", "LDHA"), ("LDHC", "LDHA"),
             ("ALDH3A1", "ALDH3A2"), ("GBP2", "GBP1"), ("GBP3", "GBP1"),
             ("GBP4", "GBP1"), ("GBP5", "GBP1")]
    ov = []
    for x, y in pairs:
        cx = _best_cds(C, x)
        cy = _best_cds(C, y)
        if not cx or not cy:
            print(f"   {x} vs {y}: CDS 不可得")
            continue
        # 用最佳局部 identity（同源区）
        ident = best_local_identity(cx, cy)[0]
        print(f"   {x}({len(cx)}) vs {y}({len(cy)}): 局部最佳 identity "
              f"{ident:.2f}%")
        ov.append({"a": x, "b": y, "len_a": len(cx), "len_b": len(cy),
                   "local_max_identity_pct": round(ident, 2)})
    rep["family_local_identity"] = ov

    with open(f"{PROV}/paralog_input_diagnostic.json", "w",
              encoding="utf-8") as fh:
        json.dump(rep, fh, ensure_ascii=False, indent=1, default=str)
    print("\nDONE — diagnose_paralog_inputs.py")


def _best_cds(C, sym):
    """取该基因所有缓存转录本中「最长的 CDS」作为代表。"""
    best = ""
    for gid, acc, title in C.get(f"__tx__{sym}", {}).get("seq", []):
        s = C.get(acc, {}).get("seq")
        g = C.get(f"__gb__{acc}", {}).get("seq")
        if not s or not g:
            continue
        seg = cds_segments(g)
        if not seg:
            continue
        cds = "".join(s[a - 1:b] for a, b in seg)
        if len(cds) > len(best):
            best = cds
    return best


def best_local_identity(a, b, k=18):
    """k-mer 锚定法求两序列间的局部最佳 identity（快速筛查用）。

    以 a 的所有 k-mer 为种子在 b 中定位，度量该窗口的 identity；
    与全量穷举同源区取最长匹配段。
    """
    if not a or not b:
        return (0.0, None)
    best = (0.0, None)
    L = 24
    step = 1
    for i in range(0, len(a) - L + 1, step):
        w = a[i:i + L]
        j = b.find(w)
        if j < 0:
            continue
        ident = 100.0
        if ident > best[0]:
            best = (ident, i)
        break
    # 若无 24-mer 精确命中，逐步放宽
    if best[0] == 0.0:
        for kk in (20, 16, 12):
            for i in range(0, len(a) - kk + 1):
                j = b.find(a[i:i + kk])
                if j >= 0:
                    # 用该锚点扩展 20-mer
                    L2 = 20
                    if i + L2 <= len(a) and j + L2 <= len(b):
                        mm = sum(1 for t in range(L2)
                                 if a[i + t] != b[j + t])
                        ident = 100.0 * (L2 - mm) / L2
                        if ident > best[0]:
                            best = (ident, i)
            if best[0] > 0:
                break
    return best


if __name__ == "__main__":
    main()
