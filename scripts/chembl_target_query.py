# -*- coding: utf-8 -*-
"""chembl_target_query.py — 从 ChEMBL 公共 API 检索枢纽基因对应蛋白靶点的
小分子生物活性（IC50/Ki/Kd），为机制导向候选药物提供数据基础。
输出：03_结果/药物预测/chembl_targets.csv、chembl_bioactivities.csv
"""
import json

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import time
import urllib.parse
import urllib.request

import pandas as pd

BASE = "https://www.ebi.ac.uk/chembl/api/data"
OUT = _ROOT + "/03_结果/药物预测"
import os
os.makedirs(OUT, exist_ok=True)

# 枢纽基因 → ChEMBL 靶点检索词（人源蛋白）
TARGETS = {
    "PKM": ["Pyruvate kinase PKM", "pyruvate kinase muscle"],
    "EP300": ["Histone acetyltransferase p300", "E1A-associated protein p300"],
    "SIRT3": ["NAD-dependent deacetylase sirtuin-3", "Sirtuin-3"],
    "STAT1": ["Signal transducer and activator of transcription 1"],
    "LDHA": ["L-lactate dehydrogenase A", "Lactate dehydrogenase A"],
    "PDK2": ["Pyruvate dehydrogenase kinase isoform 2"],
    "DHX58": ["ATP-dependent RNA helicase DHX58"],
    "GBP1": ["Guanylate-binding protein 1"],
    "HAT1": ["Histone acetyltransferase 1"],
    "TPI1": ["Triosephosphate isomerase"],
    "PGAM1": ["Phosphoglycerate mutase 1"],
}


def api_get(path, params=None, retries=3):
    url = f"{BASE}/{path}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)
        except Exception as e:
            print(f"  retry {i+1}: {e}")
            time.sleep(3)
    return None


rows_targets = []
for gene, queries in TARGETS.items():
    found = None
    for q in queries:
        d = api_get("target/search.json", {"q": q, "limit": 10})
        if not d:
            continue
        for t in d.get("targets", []):
            if t.get("organism") == "Homo sapiens":
                found = t
                break
        if found:
            break
    if found:
        rows_targets.append({"gene": gene, "chembl_target_id": found["target_chembl_id"],
                             "pref_name": found["pref_name"],
                             "target_type": found["target_type"]})
        print(f"{gene:8s} → {found['target_chembl_id']}  {found['pref_name'][:60]}")
    else:
        print(f"{gene:8s} → 未找到人源靶点")
    time.sleep(0.5)

pd.DataFrame(rows_targets).to_csv(f"{OUT}/chembl_targets.csv", index=False)

# ---- 生物活性 ----
rows_bio = []
for t in rows_targets:
    tid = t["chembl_target_id"]
    offset = 0
    while True:
        d = api_get("activity.json", {
            "target_chembl_id": tid, "limit": 500, "offset": offset,
            "pchembl_value__gte": 5,  # IC50/Ki/Kd <= 10 uM
            "standard_type__in": "IC50,Ki,Kd,EC50"})
        if not d:
            break
        acts = d.get("activities", [])
        if not acts:
            break
        for a in acts:
            rows_bio.append({
                "gene": t["gene"], "target": tid,
                "molecule_chembl_id": a.get("molecule_chembl_id"),
                "pref_name": a.get("molecule_pref_name"),
                "standard_type": a.get("standard_type"),
                "standard_value": a.get("standard_value"),
                "standard_units": a.get("standard_units"),
                "pchembl": a.get("pchembl_value"),
                "assay_type": a.get("assay_type"),
                "document_chembl_id": a.get("document_chembl_id")})
        offset += 500
        if offset >= d.get("page_meta", {}).get("total_count", 0):
            break
        time.sleep(0.3)
    print(f"{t['gene']:8s}: 累计 {sum(1 for r in rows_bio if r['gene']==t['gene'])} 条活性")

bio = pd.DataFrame(rows_bio)
bio.to_csv(f"{OUT}/chembl_bioactivities_raw.csv", index=False)
print(f"\n总活性记录: {len(bio)}")
if len(bio):
    # 汇总：每个分子在每个靶点的最佳活性
    bio["pchembl"] = pd.to_numeric(bio["pchembl"], errors="coerce")
    best = (bio.sort_values("pchembl", ascending=False)
              .drop_duplicates(subset=["gene", "molecule_chembl_id"]))
    best.to_csv(f"{OUT}/chembl_bioactivities_best.csv", index=False)
    print("\n每靶点强活性分子数 (pChEMBL>=7, ~<=100nM):")
    print(best[best["pchembl"] >= 7].groupby("gene").size().to_string())
