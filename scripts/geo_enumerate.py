# -*- coding: utf-8 -*-
"""
geo_enumerate.py — DFU 相关 GEO 数据集系统盘点（E-utilities）
项目：丙酮酰化（Kpy）与糖尿病足（DFU）
用法：python geo_enumerate.py
输出：01_数据/元数据/GEO_dfu_datasets.csv + 01_数据/元数据/GEO_dfu_datasets.json
说明：仅做只读元数据拉取；速率限制 <=3 req/s；不下载表达数据。
"""
import csv

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import json
import time
import urllib.parse
import urllib.request

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT_DIR = _ROOT + "/01_数据/元数据"


def get(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "kpy-dfu-project/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "replace")


def esearch(term: str, retmax: int = 120) -> list:
    q = urllib.parse.urlencode({"db": "gds", "term": term, "retmax": retmax, "retmode": "json"})
    d = json.loads(get(f"{BASE}/esearch.fcgi?{q}"))
    return d["esearchresult"]["idlist"], d["esearchresult"]["count"]


def esummary(ids: list) -> dict:
    q = urllib.parse.urlencode({"db": "gds", "id": ",".join(ids), "retmode": "json"})
    d = json.loads(get(f"{BASE}/esummary.fcgi?{q}"))
    return d["result"]


def main():
    # 检索式 1：diabetic foot 相关的全部 GEO Series
    ids, count = esearch('"diabetic foot"[All Fields] AND "gse"[Entry Type]')
    print(f"[esearch] diabetic foot x gse -> count={count}, fetched={len(ids)}")
    rows = []
    for i in range(0, len(ids), 100):
        chunk = ids[i:i + 100]
        res = esummary(chunk)
        for uid in res.get("uids", []):
            it = res[uid]
            rows.append({
                "uid": uid,
                "accession": it.get("accession", ""),
                "title": it.get("title", ""),
                "gdsType": it.get("gdsType", ""),
                "organism": it.get("taxon", ""),
                "n_samples": it.get("n_samples", it.get("nSamples", "")),
                "platform": it.get("platform", ""),
                "pubmed": ";".join(it.get("pubmedids", []) or []),
                "summary": (it.get("summary", "") or "").replace("\n", " ")[:500],
            })
        time.sleep(0.5)

    rows.sort(key=lambda r: r["accession"])
    with open(f"{OUT_DIR}/GEO_dfu_datasets.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(f"{OUT_DIR}/GEO_dfu_datasets.json", "w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=2)

    print(f"[out] {len(rows)} datasets -> {OUT_DIR}/GEO_dfu_datasets.csv")
    for r in rows:
        print(f"  {r['accession']:<12} {r['organism'][:18]:<18} n={str(r['n_samples']):<5} {r['gdsType'][:38]:<38} {r['title'][:95]}")


if __name__ == "__main__":
    main()
