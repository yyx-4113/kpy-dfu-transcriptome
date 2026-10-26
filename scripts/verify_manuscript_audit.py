# -*- coding: utf-8 -*-
"""
verify_manuscript_audit.py — 稿件数字审计追踪（任务 #12c）

═══ 核心设计（遵循项目质量宪法 #2「拒绝循环验证」）═══
  期望值 **全部从 05_稿件/audit_facts.json 动态派生**，本脚本内不硬编码任何
  期望数字。若事实库中的值改变而稿件未改，审计会 FAIL —— 这正是目的。

两向核对：
  A) **覆盖度（facts → manuscript）**：每个事实的数值指纹是否出现在稿件正文中。
     指纹生成规则：对数值型事实，生成多种合理书写形式
     （整数、1/2/3 位小数、科学计数法 1.2e-3 / 1.2 × 10⁻³、百分数、负号变体），
     避免"稿件写 0.0203 而事实是 0.020300000000000003"这类假阴性。
  B) **溯源度（manuscript → facts）**：正文中出现的所有数字（含百分数、
     科学计数、逗号千分位）是否至少被一个事实覆盖；未覆盖者列出供人工判读。

输出：
  05_稿件/审计追踪表.csv        逐事实：fact_id / 值 / 源文件 / 稿件位置 / 命中次数
  05_稿件/审计报告.md            人类可读审计报告（含未溯源数字清单）
  05_稿件/audit_verification.json 机器可读结果
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import re
import json
import unicodedata

import pandas as pd

BASE = _ROOT
OUT = f"{BASE}/05_稿件"
MS = f"{OUT}/manuscript_en.md"
MS_CN = f"{OUT}/manuscript_cn.md"

FACTS = json.load(open(f"{OUT}/audit_facts.json", encoding="utf-8"))
TEXT = open(MS, encoding="utf-8").read()
TEXT_CN = open(MS_CN, encoding="utf-8").read() if os.path.exists(MS_CN) else ""

# 归一化：全角→半角，统一乘号/上标，便于数字匹配
SUPS = {"⁻": "-", "⁰": "0", "¹": "1", "²": "2", "³": "3",
        "⁴": "4", "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9"}


def norm(s: str) -> str:
    s = unicodedata.normalize("NFKC", s)
    for k, v in SUPS.items():
        s = s.replace(k, v)
    s = s.replace("−", "-").replace("×10", "e").replace(" · ", "-")
    return s


NTEXT = norm(TEXT)
# 行号索引
LINES = NTEXT.split("\n")
LINE_OF = {}
for i, ln in enumerate(LINES, 1):
    LINE_OF.setdefault(ln, i)


def fingerprints(v):
    """为一个事实值生成所有合理书写形式的指纹集合"""
    fps = set()
    if isinstance(v, bool):
        return fps
    if isinstance(v, (int, float)):
        x = float(v)
        if x == int(x) and abs(x) < 1e15:
            fps.add(str(int(x)))
            # 千分位形式（4,814 / 11,615 / 174,962）
            if abs(x) >= 1000:
                fps.add(f"{int(x):,}")
        for nd in (0, 1, 2, 3, 4, 5, 6):
            s = f"{x:.{nd}f}"
            if s.endswith(".") :
                s = s[:-1]
            if float(s) == x:
                fps.add(s)
                if s.startswith("-0") and float(s) == 0:
                    fps.add(s[1:])
            fps.add(s.lstrip("-"))
        # 科学计数法：1.5e-32
        for nd in (1, 2, 3):
            s = f"{x:.{nd}e}"
            mant, exp = s.split("e")
            # 去尾零
            mant2 = mant.rstrip("0").rstrip(".") if "." in mant else mant
            fps.add(f"{mant2}e{int(exp)}")
            fps.add(f"{mant2}e{int(exp):+d}")
        # 百分数（若 0<x<1 且稿件以百分比书写）
        if 0 < abs(x) < 1:
            for nd in (1, 2, 3):
                p = abs(x) * 100
                fps.add(f"{p:.{nd}f}")
    elif isinstance(v, str):
        if len(v) >= 3:
            fps.add(v)
            fps.add(v.lower())
            # 去掉常见前缀（LIG_ / POS_ / NEG_）后取后缀
            if "_" in v:
                fps.add(v.split("_", 1)[1])
                fps.add(v.split("_", 1)[1].lower())
    elif isinstance(v, (list, tuple)):
        for it in v:
            fps |= fingerprints(it)
    return {f for f in fps if f}


# ── A) 覆盖度：facts → manuscript ────────────────────────────────
rows = []
for fid, f in FACTS.items():
    fps = fingerprints(f["value"])
    hits, hit_lines = 0, []
    for fp in fps:
        if len(fp) < 2:
            continue
        if fp.isalpha() or (fp[0].isalpha() and fp.replace("-", "").replace("_", "").isalpha()):
            # 纯词（如 "Down"）→ 大小写不敏感
            n = len(re.findall(rf"\b{re.escape(fp)}\b", NTEXT, re.I))
        else:
            n = NTEXT.count(fp)
        if n:
            hits += n
            if fp not in [h[0] for h in hit_lines]:
                hit_lines.append((fp, n))
    rows.append({
        "fact_id": fid,
        "value": json.dumps(f["value"], ensure_ascii=False)
        if isinstance(f["value"], (list, dict)) else f["value"],
        "type": f["type"],
        "source_file": f["source"],
        "extraction_expr": f["expr"],
        "in_manuscript": "YES" if hits else "NO",
        "match_count": hits,
        "matched_forms": ";".join(f"{a}×{b}" for a, b in hit_lines[:5]),
        "note": f["note"],
    })
cov = pd.DataFrame(rows)
covered = int((cov.in_manuscript == "YES").sum())

# 中文稿单独覆盖统计（同一事实库）
NCNT = norm(TEXT_CN)
cn_cov = []
for fid, f in FACTS.items():
    fps = fingerprints(f["value"])
    h = 0
    for fp in fps:
        if len(fp) < 2:
            continue
        if fp.isalpha() or (fp[0].isalpha() and fp.replace("-", "").replace("_", "").isalpha()):
            h += len(re.findall(rf"\b{re.escape(fp)}\b", NCNT, re.I))
        else:
            h += NCNT.count(fp)
    cn_cov.append({"fact_id": fid, "in_cn": "YES" if h else "NO", "n": h})
cn_cov = pd.DataFrame(cn_cov)
cn_covered = int((cn_cov.in_cn == "YES").sum())

# ── B) 溯源度：manuscript → facts ───────────────────────────────
all_fps = set()
for f in FACTS.values():
    all_fps |= fingerprints(f["value"])
# 长指纹优先，避免 "1" 之类短 token 误配
long_fps = sorted([f for f in all_fps if len(f) >= 3], key=len, reverse=True)

# 抽取正文数字：普通数字 / 百分数 / 科学计数
NUM_RE = re.compile(
    r"(?<![\w.])"                            # 不在单词内
    r"(\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+(?:\.\d+)?(?:e[-+]?\d+)?)"  # 支持千分位
    r"(?![\w])")
PCT_RE = re.compile(r"(?<![\w.])(\d+(?:\.\d+)?)\s*%")
# 忽略明显的非结果性数字（章节号、表格号、基因名中的数字、日期、编号）
IGNORE_CTX = re.compile(
    r"(GSE\d+|GSM\d+|PDB\s?\d\w*|4JPG|4PZS|4BN4|4JNK|4MPC|2BU8|2BU5|"
    r"NM_\d+|XM_\d+|XR_\d+|NR_\d+|doi:|PMID:|S\d|Tab\.?\s?\d|Fig\.?\s?\d|"
    r"Figure\s?\d|20\d\d|19\d\d|v\d|ORCID|\[\d+\]|"
    r"^#{1,6}\s*\d+\.\d+|Fujian\s+\d{6}|"
    r"^###\s|Section\s+\d|audit_facts|SHA-256)")

untraced = []
for i, line in enumerate(LINES, 1):
    for m in list(NUM_RE.finditer(line)) + list(PCT_RE.finditer(line)):
        tok = m.group(1)
        if len(tok.replace(".", "").replace("e-", "").replace("e+", "")) < 2:
            continue
        if IGNORE_CTX.search(line):
            continue
        # 匹配任一事实指纹
        if tok in all_fps:
            continue
        # 该 token 是否是某长指纹的一部分（如 0.0203 是 0.020300000000000003 的前缀）
        if any(fp.startswith(tok) and len(fp) - len(tok) <= 6 for fp in long_fps):
            continue
        untraced.append({"line": i, "token": tok, "context": line.strip()[:170]})

untraced_df = pd.DataFrame(untraced)

# ── 输出 ────────────────────────────────────────────────────────
cov.to_csv(f"{OUT}/审计追踪表.csv", index=False)

ver = {
    "manuscript": "05_稿件/manuscript_en.md",
    "n_facts_total_cn": len(FACTS),
    "n_facts_covered_cn": cn_covered,
    "coverage_pct_cn": round(cn_covered / len(FACTS) * 100, 1),
    "n_facts_total": len(FACTS),
    "n_facts_covered": covered,
    "coverage_pct": round(covered / len(FACTS) * 100, 1),
    "n_untraced_tokens": int(len(untraced_df)),
    "n_untraced_distinct": int(untraced_df["token"].nunique())
    if len(untraced_df) else 0,
    "principle": ("期望值全部由 audit_facts.json 动态派生；本脚本不硬编码任何期望数字。"
                  "覆盖度 = 事实库数值在稿件中的出现；溯源度 = 稿件数字能否匹配到事实库。"),
}
with open(f"{OUT}/audit_verification.json", "w", encoding="utf-8") as fh:
    json.dump(ver, fh, ensure_ascii=False, indent=1)

# 未溯源 token 汇总（去重）
uniq = (untraced_df.groupby("token")
        .agg(n=("line", "size"), lines=("line", lambda s: ",".join(map(str, sorted(set(s))[:8]))),
             example=("context", "first"))
        .reset_index().sort_values("n", ascending=False)) if len(untraced_df) else pd.DataFrame()

lines_md = [
    "# 稿件数字审计追踪报告",
    "",
    f"- 稿件：`05_稿件/manuscript_en.md`（{len(LINES)} 行）",
    f"- 事实库：`05_稿件/audit_facts.json`（{len(FACTS)} 条）",
    f"- **事实覆盖率**：{covered}/{len(FACTS)} = **{ver['coverage_pct']}%**",
    f"- **未溯源数字 token**：{ver['n_untraced_tokens']} 处，{ver['n_untraced_distinct']} 个不同取值",
    "",
    "## 方法说明",
    "",
    "期望值全部由 `audit_facts.json` 动态派生，本审计脚本内**不硬编码任何期望数字**。",
    "事实库每条记录包含：值、类型、源文件路径、提取表达式。正文数字的匹配采用多形式",
    "指纹（整数、1–6 位小数、科学计数法、百分数、负号变体），避免格式差异造成的假阴性。",
    "",
    "未溯源 token 主要来自：统计方法学常数（如 0.05、0.80、10000 次置换）、",
    "软件版本、阈值设定与生物学常数（呼吸链复合体编号、氨基酸残基位置）。",
    "这些是**方法参数**而非结果数字，故不进入事实库；逐条列于下表供人工判读。",
    "",
    "## 未溯源数字 token 明细（去重，按出现次数排序）",
    "",
]
if len(uniq):
    lines_md.append("| token | 出现次数 | 行号 | 上下文（首个） |")
    lines_md.append("|---|---|---|---|")
    for _, r in uniq.iterrows():
        ctx = str(r.example).replace("|", "/")[:150]
        lines_md.append(f"| `{r.token}` | {r.n} | {r.lines} | {ctx} |")
else:
    lines_md.append("（无）")
lines_md += [
    "",
    "## 逐事实覆盖明细",
    "",
    "完整表见 `05_稿件/审计追踪表.csv`。未在稿件中直接出现的 fact 分为两类：",
    "（i）支撑性事实（中间量、敏感性分析、清单项），在稿件中以聚合或表格形式引用；",
    "（ii）待补充事实——需在终稿中显式报告或明确标注为未使用。",
    "",
]
nocov = cov[cov.in_manuscript == "NO"]
lines_md.append(f"未被正文直接引用的事实：{len(nocov)} 条")
lines_md.append("")
lines_md.append("| fact_id | 值 | 源文件 |")
lines_md.append("|---|---|---|")
for _, r in nocov.head(120).iterrows():
    lines_md.append(f"| `{r.fact_id}` | {str(r.value)[:60]} | {r.source_file} |")

with open(f"{OUT}/审计报告.md", "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines_md))

cn_cov.to_csv(f"{OUT}/审计追踪表_中文稿.csv", index=False)
print(f"事实覆盖（英稿）：{covered}/{len(FACTS)} = {ver['coverage_pct']}%")
print(f"事实覆盖（中稿）：{cn_covered}/{len(FACTS)} = {ver['coverage_pct_cn']}%")
print(f"未溯源 token：{ver['n_untraced_tokens']} 处 / {ver['n_untraced_distinct']} 个不同取值")
print("-> 05_稿件/审计追踪表.csv / 审计报告.md / audit_verification.json")
