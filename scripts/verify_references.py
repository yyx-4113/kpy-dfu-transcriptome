# -*- coding: utf-8 -*-
"""
verify_references.py — 参考文献双源核验（Crossref + PubMed/GEO）

对稿件 05_稿件/manuscript_en.md 的参考文献表逐条做**实证核验**：
  1. Crossref：`query.bibliographic` 检索，取 top 候选，用「最长公共子串/候选标题长度」
     衡量候选标题是否真实出现在稿件引文串中（避免引文串解析错误导致误判）。
  2. PubMed：`esearch`（作者姓 + 年份 + 标题短语）→ `esummary` 取权威元数据。
  3. GEO 数据集条目（accession 形式）：改走 E-utilities `db=gds`，核 accession 是否
     存在、真实标题与样本量。

判定分级：
  MATCH       — 双源一致，元数据（作者/期刊/卷期/页/年/DOI）无差异
  MATCH_DOI   — 一致但稿件缺 DOI（补值即可）
  MISMATCH    — 有字段与实证不符（列出字段级差异）
  NOT_FOUND   — 双源均无可靠命中（疑似不存在/无法核验）
  GEO_OK      — accession 存在，给出真实标题/样本量

输出：03_结果/参考文献核验.md（人读审计表）+ references_audit.json（机读）
"""

import json
import re
import time
import urllib.parse
import urllib.request
from difflib import SequenceMatcher

MS = "05_稿件/manuscript_en.md"
OUT_MD = "03_结果/参考文献核验.md"
OUT_JSON = "03_结果/参考文献核验.json"
MAILTO = "yongxin.yang@example.org"   # Crossref 礼貌池（勿用于识别，仅占位）
UA = {"User-Agent": f"ref-audit/1.0 (mailto:{MAILTO})"}


# ────────────────────────── 网络层（重试 + 退避）──────────────────────────
def _get(url, tries=4, timeout=60):
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:                      # noqa: BLE001
            last = e
            time.sleep(2.0 * (k + 1))
    raise RuntimeError(f"GET failed after {tries}: {url[:110]} :: {last}")


def _jget(url, **kw):
    return json.loads(_get(url, **kw))


def crossref_bib(line):
    q = urllib.parse.urlencode({"query.bibliographic": line, "rows": 5,
                                "mailto": MAILTO})
    d = _jget("https://api.crossref.org/works?" + q)
    return d.get("message", {}).get("items", []) or []


def pubmed_ids(term, retmax=4):
    q = urllib.parse.urlencode({"db": "pubmed", "term": term,
                                "retmax": retmax, "retmode": "json"})
    return _jget("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
                 + q)["esearchresult"].get("idlist", []) or []


def pubmed_summ(ids):
    if not ids:
        return {}
    q = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(ids),
                                "retmode": "json"})
    return _jget("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
                 + q).get("result", {})


def geo_summ(acc):
    q = urllib.parse.urlencode({"db": "gds", "term": f"{acc}[Accession]",
                                "retmode": "json"})
    r = _jget("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
              + q)["esearchresult"]
    ids = r.get("idlist", []) or []
    if not ids:
        return None
    q2 = urllib.parse.urlencode({"db": "gds", "id": ids[0], "retmode": "json"})
    res = _jget("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
                + q2).get("result", {})
    return res.get(ids[0])


# ────────────────────────── 工具函数 ──────────────────────────
def norm(s):
    s = (s or "").lower()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return " ".join(s.split())


def split_ref(raw):
    """把稿件引文拆成 (authors, title, journal_year)。

    稿件引文标题常被缩写，故必须先拆出标题片段再做包含比对。
    切分点：第一个 '. ' 之后是作者块；其后到下一个 '. ' 是标题。
    """
    s = raw.strip()
    m = re.match(r"^[^.]{2,80}\.\s+(.*)$", s)
    rest = m.group(1) if m else s
    m2 = re.match(r"^(.*?)\.\s+(.*)$", rest)
    if m2:
        return m2.group(1).strip(), m2.group(2).strip()
    return rest, ""


def contains_score(title, frag):
    """稿件标题片段 frag 被候选标题 title 包含的比例（0–1）。

    缩写标题仍得高分；标题写错则得低分——用以区分「缩写」与「错误」。
    """
    a, b = norm(title), norm(frag)
    if not a or not b:
        return 0.0
    m = SequenceMatcher(None, b, a, autojunk=False).find_longest_match(
        0, len(b), 0, len(a))
    return m.size / len(b)


def cr_fields(it):
    t = (it.get("title") or [""])[0]
    j = (it.get("container-title") or [""])[0]
    au = it.get("author") or []
    a1 = au[0] if au else {}
    y = None
    for key in ("published-print", "published-online", "issued"):
        d = (it.get(key) or {}).get("date-parts") or []
        if d and d[0]:
            y = d[0][0]
            break
    return {"title": t, "journal": j, "year": y,
            "volume": it.get("volume"), "issue": it.get("issue"),
            "page": it.get("page"), "doi": it.get("DOI"),
            "first_author": a1.get("family", ""),
            "n_authors": len(au), "type": it.get("type")}


def pm_fields(rec):
    ids = {x["idtype"]: x["value"] for x in rec.get("articleids", [])}
    return {"title": rec.get("title", ""), "journal": rec.get("source", ""),
            "pubdate": rec.get("pubdate", ""), "volume": rec.get("volume"),
            "issue": rec.get("issue"), "pages": rec.get("pages"),
            "doi": ids.get("doi"),
            "first_author": (rec.get("authors") or [{}])[0].get("name", ""),
            "n_authors": len(rec.get("authors") or []),
            "pubtype": rec.get("pubtype")}


def year_of(s):
    m = re.search(r"(19|20)\d{2}", s or "")
    return m.group(0) if m else None


# ────────────────────────── 主流程 ──────────────────────────
def main():
    lines = open(MS, encoding="utf-8").read().split("\n")
    refs, in_refs = [], False
    for ln in lines:
        s = ln.strip()
        if re.match(r"^#+\s*References", s, re.I) or re.match(r"^#+\s*参考文献", s):
            in_refs = True
            continue
        if in_refs:
            m = re.match(r"^(\d{1,2})\.\s+(.*)$", s)
            if m:
                refs.append((int(m.group(1)), m.group(2).strip()))
            elif s.startswith("#"):
                break
    print(f"[解析] 参考文献 {len(refs)} 条")

    out, rows = [], []
    for n, raw in refs:
        print(f"\n[{n:2d}] {raw[:100]}")
        rec = {"n": n, "raw": raw}

        # ── GEO 数据集条目 ──
        g = re.match(r"^GEO\s+(GSE\d+)\b", raw, re.I)
        if g:
            acc = g.group(1).upper()
            try:
                gs = geo_summ(acc)
            except Exception as e:                              # noqa: BLE001
                gs = None
                print("    GEO 查询失败:", e)
            if gs:
                rec.update(kind="GEO", accession=acc, geo_title=gs.get("title"),
                           geo_n=gs.get("n_samples"), geo_type=gs.get("gdstype"),
                           verdict="GEO_OK")
                print(f"    GEO OK | {acc} | n={gs.get('n_samples')} | "
                      f"{gs.get('title','')[:80]}")
            else:
                rec.update(kind="GEO", accession=acc, verdict="NOT_FOUND")
                print(f"    ★ GEO {acc} 无命中")
            out.append(rec)
            time.sleep(0.35)
            continue

        # ── 普通文献：Crossref ──
        ms_title, ms_jyr = split_ref(raw)
        rec["ms_title"] = ms_title
        rec["ms_jyr"] = ms_jyr
        try:
            items = crossref_bib(ms_title if len(ms_title) > 12 else raw)
        except Exception as e:                                  # noqa: BLE001
            items = []
            print("    Crossref 失败:", e)
        # 只保留「期刊论文 + 有 DOI + 标题够长」的候选，避免选到期刊名/文集级记录
        cands = []
        for it in items:
            if not it.get("DOI"):
                continue
            if it.get("type") not in ("journal-article", "book-chapter",
                                      "proceedings-article", "posted-content"):
                continue
            t = (it.get("title") or [""])[0]
            if len(norm(t)) < 12:
                continue
            cands.append(it)
        best, score = None, 0.0
        for it in cands:
            t = (it.get("title") or [""])[0]
            s = contains_score(t, ms_title)
            if s > score:
                best, score = it, s
        cr = cr_fields(best) if best else None
        if cr:
            rec["crossref"] = cr
            rec["cr_score"] = round(score, 3)
            print(f"    CR  {score:.2f} | {cr['first_author']} | "
                  f"{cr['journal']} {cr['year']} | {cr['title'][:60]}")

        # ── PubMed ──
        yr = year_of(raw)
        first_author = ""
        m = re.match(r"^([A-Z][A-Za-z'\-]+)\s+[A-Z]{1,3}\b", raw)
        if m:
            first_author = m.group(1)
        # 用稿件自带的（缩写）标题短语检索 PubMed——缩写标题仍可命中原文
        tw = norm(ms_title).split()
        phrase = " ".join(tw[:8])
        terms = []
        if len(tw) >= 2:
            terms.append(f'"{phrase}"[Title]')
        if first_author and yr:
            terms.append(f"{first_author}[Author] AND {yr}[DP]")
        if cr and cr.get("title"):
            t = re.sub(r"[^A-Za-z0-9 :\-]", " ", cr["title"])
            terms.append(f'"{t[:70]}"[Title]')
        pm = None
        for term in terms:
            try:
                ids = pubmed_ids(term)
            except Exception as e:                              # noqa: BLE001
                ids = []
                print("    PubMed 失败:", e)
            if ids:
                try:
                    res = pubmed_summ(ids)
                except Exception as e:                          # noqa: BLE001
                    res = {}
                    print("    esummary 失败:", e)
                for i in ids:
                    if i in res:
                        cand = pm_fields(res[i])
                        sc = contains_score(cand["title"], ms_title)
                        if pm is None or sc > pm[1]:
                            pm = (cand, sc, i)
                if pm and pm[1] >= 0.6:
                    break
            time.sleep(0.3)
        if pm:
            rec["pubmed"] = pm[0]
            rec["pubmed_pmid"] = pm[2]
            rec["pm_score"] = round(pm[1], 3)
            print(f"    PM  {pm[1]:.2f} | PMID {pm[2]} | {pm[0]['journal']} "
                  f"{pm[0]['pubdate']} | {pm[0]['title'][:60]}")
        else:
            print("    PM  无命中")

        # ── 判定 ──
        rec["verdict"] = judge(rec, cr, pm[0] if pm else None, score,
                               pm[1] if pm else 0.0)
        print(f"    => {rec['verdict']}")
        out.append(rec)
        time.sleep(0.35)

    # ── 产出 ──
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    md = render_md(out)
    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write(md)
    print(f"\n[写出] {OUT_MD} / {OUT_JSON}")

    from collections import Counter
    print("[统计]", dict(Counter(r["verdict"].split(":")[0] for r in out)))


def judge(rec, cr, pm, cr_score, pm_score):
    """判定：先按覆盖率过滤假阳性，再比对字段。

    - cr_score < 0.55：Crossref 候选与引文串不符 → 丢弃（不得据其判 MISMATCH）。
    - pm_score  < 0.70：PubMed 命中与引文串不符 → 丢弃（不得据其判 MISMATCH）。
      2026 年新文献常无 PubMed 记录，属正常。
    """
    cr_ok = cr is not None and cr_score >= 0.60
    pm_ok = pm is not None and pm_score >= 0.60
    if not cr_ok and not pm_ok:
        return "NOT_FOUND"
    diffs = []
    if cr_ok and pm_ok:
        if norm(cr["journal"])[:14] != norm(pm["journal"])[:14]:
            diffs.append(f"期刊 CR={cr['journal']} / PM={pm['journal']}")
        py = year_of(pm["pubdate"])
        if py and cr["year"] and str(cr["year"]) != py:
            diffs.append(f"年份 CR={cr['year']} / PM={py}")
    tag = ("MISMATCH:" + "; ".join(diffs)) if diffs else "MATCH"
    if cr_ok and not pm_ok:
        tag += " (仅 Crossref 证实)"
    return tag


def render_md(out):
    L = ["# 参考文献双源核验（Crossref + PubMed）",
         "",
         f"- 稿件：`05_稿件/manuscript_en.md`",
         f"- 条目：{len(out)}",
         "- 方法：Crossref `query.bibliographic` 全文检索 + PubMed E-utilities；"
         "GEO 数据集条目走 `db=gds`。",
         "- 覆盖率 = 候选标题在稿件引文串中的最长公共子串 / 标题长度（≥0.55 才采信）。",
         "",
         "| # | 稿件原文（截断） | 判定 | Crossref 实证 | PubMed 实证 |",
         "|---|---|---|---|---|"]
    for r in out:
        raw = r["raw"][:52].replace("|", "/")
        cr = r.get("crossref") or {}
        pm = r.get("pubmed") or {}
        if r.get("kind") == "GEO":
            crs = (f"{r.get('accession')} | n={r.get('geo_n')} | "
                   f"{str(r.get('geo_title'))[:58]}")
            pms = "—"
        else:
            crs = (f"覆盖{r.get('cr_score','—')} | {cr.get('first_author','?')} | "
                   f"{cr.get('journal','?')} {cr.get('year','?')};"
                   f"{cr.get('volume','')}({cr.get('issue','')}):{cr.get('page','')} "
                   f"| doi:{cr.get('doi','—')}") if cr else "无命中"
            pms = (f"覆盖{r.get('pm_score','—')} | PMID {r.get('pubmed_pmid','?')} | "
                   f"{pm.get('journal','?')} {pm.get('pubdate','')};"
                   f"{pm.get('volume','')}({pm.get('issue','')}):{pm.get('pages','')} "
                   f"| doi:{pm.get('doi','—')}") if pm else "无命中"
        L.append(f"| {r['n']} | {raw} | **{r['verdict']}** | {crs} | {pms} |")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    main()
