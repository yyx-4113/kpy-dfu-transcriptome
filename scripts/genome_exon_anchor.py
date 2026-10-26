# -*- coding: utf-8 -*-
"""
genome_exon_anchor.py — 基于**基因组坐标**的引物跨接头（exon–exon junction）判定

═══ 为什么必须用基因组坐标 ═══

RefSeq **mRNA** 记录（NM_/XM_）的 `exon` 特征坐标是 **mRNA 自身坐标**
（已拼接），相邻 exon 在 mRNA 上连续，因此**无法**用 mRNA 坐标判断
"引物是否跨越 exon–exon 接头"。

真正的判定必须基于 **RefSeqGene（NG_）** 记录：
- NG_ 记录是基因组区段，`exon` 特征为**基因组坐标**；
- 其 `mRNA` 特征给出该转录本在基因组上的 `join(...)` 拼接坐标；
- 若引物在基因组上落在两个外显子之间（即跨过内含子），
  则它在 **gDNA** 上有扩增信号、在 **cDNA** 上被剪掉 →
  这是 RT-qPCR 的经典伪影来源（基因组 DNA 残留造成假阳性），
  也是引物不理想的原因。

═══ 判定逻辑 ═══

对每个候选引物：
1. 取该基因的 RefSeqGene（NG_）记录，解析其 `mRNA` 特征的 join 坐标
   （该转录本各外显子在基因组上的位置）；
2. 把引物序列在 **基因组序列**（NG_ 的 ORIGIN 段）中定位（正/反链均查）；
3. 若引物区间 **跨过 join 坐标的断点** → `spans_junction = True`
   （即引物横跨外显子边界）。

注：能横跨接头本身**不是绝对禁用**（接头特异引物是常用策略），
但本项目要求**在全部编码转录本共有区段**内设计，且要求引物完整位于
CDS 内——横跨接头的引物必然只存在于部分转录本，与"共有区段"约束
本身已部分互斥。此处输出该指标作为**独立、可复核**的证据列，
用于在候选排序中优先选择不跨接头的引物，并写入补充材料。

输出：
  03_结果/实验验证设计/provenance/refseqgene_records.json  留档
  03_结果/实验验证设计/provenance/genome_exon_anchor.json   判定结果
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import re
import json
import time
import urllib.request

BASE = _ROOT
OUT = f"{BASE}/03_结果/实验验证设计"
PROV = f"{OUT}/provenance"
os.makedirs(PROV, exist_ok=True)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "WorkBuddy-GenomeAnchor/1.0"}
RC_MAP = str.maketrans("ACGT", "TGCA")

# gene_id → RefSeqGene accession（实测，见 refseqgene_records.json）
NG_BY_GENE = {
    "LDHA": ("NG_008185.1", "LDHA"),
    "ALDH3A2": ("NG_007095.2", "ALDH3A2"),
    "SDHA": ("NG_012339.1", "SDHA"),
}


def rc(s):
    return s.translate(RC_MAP)[::-1]


def _get(url, tries=4):
    last = None
    for k in range(tries):
        try:
            return urllib.request.urlopen(
                urllib.request.Request(url, headers=UA), timeout=180).read().decode()
        except Exception as e:                        # noqa: BLE001
            last = e
            time.sleep(2 + 3 * k)
    raise RuntimeError(str(last))


def efetch_gb_with_origin(acc):
    raw = _get(f"{EUTILS}/efetch.fcgi?db=nuccore&id={acc}"
               f"&rettype=gb&retmode=text")
    i = raw.find("ORIGIN")
    if i < 0:
        return raw, ""
    body = raw[i:]
    seq = []
    for line in body.split("\n"):
        line = line.strip()
        if line.startswith("//"):
            break
        parts = line.split()
        if parts and parts[0].isdigit():
            seq.append("".join(parts[1:]))
    return raw, "".join(seq).upper()


def parse_mrna_join(gb):
    """解析 NG_ 记录中 mRNA 特征的 join 坐标（基因组坐标，1-based 闭区间）。

    ═══ 两个必须处理的格式细节 ═══
    1. `join(...)` **跨行折断**：长转录本的 location 可占多行
       （如 SDHA NM_004168.4 的 join 有 17 段，跨 4 行）。
       → 必须把特征块（feature block）整体拼接后再解析，不能只看首行。
    2. `complement(...)` 表示**负链**：此时 join 区间按基因组正链降序给出，
       重建的转录本序列 = 各段正链序列**反向互补后逆序拼接**。
       → 必须记录 strand 并在重建时正确处理。
    """
    out = []
    # 按 feature block 切分：块首为 5 空格缩进 + 关键字
    blocks = re.split(r"\n {5}(?=\w)", gb)
    for b in blocks:
        if not b.startswith("mRNA"):
            continue
        head = b.split("\n")[0]
        # location 可能跨行：把后续的续行（21 空格缩进、非 qualifier）
        # 拼到 head 上
        loc_text = head[len("mRNA"):].strip()
        for ln in b.split("\n")[1:]:
            stripped = ln.strip()
            if not stripped or stripped.startswith("/"):
                continue
            loc_text += stripped
        tx = ""
        mt = re.search(r'/transcript_id="([^"]+)"', b)
        if mt:
            tx = mt.group(1)
        # 基因名：NG_ 记录覆盖整个基因组区段，其中可能包含**邻近基因**的
        # 转录本（如 NG_007095.2 内含 SLC47A2 的 NM_001099646.3；
        # NG_012339.1 内含 CCDC127 的 NM_145265.3）。必须按目标基因名过滤，
        # 否则会把邻近基因的外显子结构误当作目标基因的。
        gene_name = ""
        mg = re.search(r'/gene="([^"]+)"', b)
        if mg:
            gene_name = mg.group(1)
        strand = "-" if "complement" in loc_text else "+"
        parts = re.findall(r"(\d+)\.\.(\d+)", loc_text.replace("<", ">"))
        segs = [(int(a), int(b)) for a, b in parts]
        # 单段形式 `complement(34567..>35156)` → parts 为空，需单独处理
        ms = re.search(r"(\d+)\.\.>?(\d+)", loc_text)
        if not segs and ms:
            segs = [(int(ms.group(1)), int(ms.group(2)))]
        out.append({"transcript_id": tx, "strand": strand, "gene_name": gene_name,
                    "location_raw": loc_text[:400],
                    "segments": segs})
    return out


def locate_genome(gseq, primer):
    """在基因组序列中定位引物（正链与反链均查），返回 [(start0,end0,strand)]"""
    hits = []
    for s, strand in ((primer, "+"), (rc(primer), "-")):
        st = 0
        while True:
            i = gseq.find(s, st)
            if i < 0:
                break
            hits.append((i, i + len(s), strand))
            st = i + 1
    return hits


def spans_junction(seg_start, seg_end, segments):
    """引物区间 [seg_start, seg_end) 是否跨过任一外显子边界。"""
    for k in range(len(segments) - 1):
        a_end = segments[k][1]           # 1-based 闭区间右端
        b_start = segments[k + 1][0]     # 下一段左端
        # 基因组区间（0-based 半开）与边界的关系：
        # 断点位于 a_end 与 b_start 之间的内含子中
        if seg_start < b_start - 1 and seg_end > a_end:
            return True
    return False


def main():
    recs = {}
    for gene, (acc, gene_sym) in NG_BY_GENE.items():
        gb, gseq = efetch_gb_with_origin(acc)
        time.sleep(0.4)
        mrnas = parse_mrna_join(gb)
        # 只保留**目标基因自身**的转录本（NG_ 记录覆盖整个基因组区段，
        # 内部还含邻近基因如 SLC47A2 / CCDC127 的转录本）
        mrnas = [m for m in mrnas if m.get("gene_name") == gene_sym]
        print(f"  {gene}: {acc} 内 {len(mrnas)} 条属于 {gene_sym} 的 mRNA join")
        recs[gene] = {"refseqgene": acc, "genome_len": len(gseq),
                      "mrna_joins": mrnas,
                      # 落盘基因组序列（供 verify_primers.py 的 V4 判定使用；
                      # NG_ 记录长 3–45 kb，体积可接受）
                      "genome_seq": gseq}
        print(f"  {gene:9s} {acc:14s} genome={len(gseq):7d} bp, "
              f"mRNA join 记录 {len(mrnas)} 条")
        for m in mrnas[:3]:
            print(f"      {m['transcript_id']:16s} segments={len(m['segments'])}")
    with open(f"{PROV}/refseqgene_records.json", "w", encoding="utf-8") as fh:
        json.dump(recs, fh, ensure_ascii=False, indent=1)

    # 独立交叉核对：把每条 join 的转录本序列重建出来，
    # 与该转录本的 mRNA 序列比对 → 验证 join 坐标解析正确
    print("\n  --- join 坐标解析交叉核对 ---")
    for gene, r in recs.items():
        for m in r["mrna_joins"]:
            tx = m["transcript_id"]
            if not tx:
                continue
            try:
                raw = _get(f"{EUTILS}/efetch.fcgi?db=nuccore&id={tx}"
                           f"&rettype=fasta&retmode=text")
                time.sleep(0.35)
                L = raw.strip().split("\n")
                mrna = "".join(x.strip() for x in L[1:]
                               if re.fullmatch(r"[A-Za-z]+", x.strip())).upper()
                segs = m["segments"]
                if m["strand"] == "-":
                    # 负链：各段取反向互补，且按基因组升序逆序拼接
                    built = "".join(
                        r["genome_seq"][a - 1:b].translate(RC_MAP)[::-1]
                        for a, b in sorted(segs))
                else:
                    built = "".join(r["genome_seq"][a - 1:b]
                                    for a, b in sorted(segs))
                ok = (built == mrna)
                print(f"    {gene:9s} {tx:16s} strand={m['strand']} "
                      f"segments={len(segs):2d} join重建 "
                      f"{'== mRNA ✓' if ok else '!= mRNA ✗'} "
                      f"(len {len(built)} vs {len(mrna)})")
                m["join_reconstruction_matches_mrna"] = bool(ok)
            except Exception as e:                    # noqa: BLE001
                print(f"    {gene:9s} {tx:16s} 核对失败: {e}")
                m["join_reconstruction_matches_mrna"] = None
    with open(f"{PROV}/refseqgene_records.json", "w", encoding="utf-8") as fh:
        json.dump(recs, fh, ensure_ascii=False, indent=1)

    print("\nDONE — genome_exon_anchor.py（记录采集 + join 交叉核对）")


if __name__ == "__main__":
    main()
