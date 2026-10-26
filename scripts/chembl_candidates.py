# -*- coding: utf-8 -*-
"""chembl_candidates.py — 修正靶点 ID 后重新查询 ChEMBL 活性，并获取分子 SMILES。
对接面板：PKM2 / EP300(p300) / SIRT3 / LDHA / PDK2 / STAT1（DHX58 无 ChEMBL 人源靶点，
GBP1 无单蛋白靶点 → 两者退出对接面板，稿件中如实说明）。
输出：03_结果/药物预测/chembl_targets.csv、chembl_bioactivities.csv、chembl_candidates.csv
"""
import json

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import os
import time
import urllib.parse
import urllib.request

import pandas as pd

OUT = _ROOT + "/03_结果/药物预测"
os.makedirs(OUT, exist_ok=True)

# 修正后的 ChEMBL 靶点 ID（2026-10-25 核实）
TARGETS = {
    "PKM":   {"chembl_target_id": "CHEMBL1075189", "role": "丙酮酸生成酶（PKM2）"},
    "EP300": {"chembl_target_id": "CHEMBL6066577", "role": "Kpy 写入酶 p300 HAT 域"},
    "SIRT3": {"chembl_target_id": "CHEMBL4461",    "role": "Kpy 擦除酶 SIRT3"},
    "LDHA":  {"chembl_target_id": "CHEMBL4835",    "role": "丙酮酸分流（乳酸脱氢酶 A）"},
    "PDK2":  {"chembl_target_id": "CHEMBL3861",    "role": "丙酮酸进入 TCA 的闸门 PDK2"},
    "STAT1": {"chembl_target_id": "CHEMBL6101",    "role": "IFN-I 读出 STAT1"},
}


def api(path, params=None, retries=3):
    url = f"https://www.ebi.ac.uk/chembl/api/data/{path}"
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


trows = []
for gene, info in TARGETS.items():
    d = api(f"target/{info['chembl_target_id']}.json")
    pref = d["pref_name"] if d else "?"
    trows.append({"gene": gene, "chembl_target_id": info["chembl_target_id"],
                  "pref_name": pref, "role": info["role"]})
    print(f"{gene:6s} {info['chembl_target_id']:18s} {pref[:60]}")
pd.DataFrame(trows).to_csv(f"{OUT}/chembl_targets.csv", index=False)

brows = []
for gene, info in TARGETS.items():
    tid = info["chembl_target_id"]
    offset, total = 0, None
    while True:
        d = api("activity.json", {
            "target_chembl_id": tid, "limit": 500, "offset": offset,
            "pchembl_value__gte": 5,
            "standard_type__in": "IC50,Ki,Kd,EC50,AC50"})
        if not d:
            break
        acts = d.get("activities", [])
        total = d.get("page_meta", {}).get("total_count", 0)
        if not acts:
            break
        for a in acts:
            brows.append({
                "gene": gene, "target": tid,
                "molecule_chembl_id": a.get("molecule_chembl_id"),
                "pref_name": a.get("molecule_pref_name"),
                "standard_type": a.get("standard_type"),
                "standard_value": a.get("standard_value"),
                "standard_units": a.get("standard_units"),
                "pchembl": a.get("pchembl_value"),
                "assay_type": a.get("assay_type"),
                "document_chembl_id": a.get("document_chembl_id")})
        offset += 500
        if offset >= total:
            break
        time.sleep(0.3)
    print(f"  {gene}: {total} 条")
    time.sleep(0.3)

bio = pd.DataFrame(brows)
bio["pchembl"] = pd.to_numeric(bio["pchembl"], errors="coerce")
bio.to_csv(f"{OUT}/chembl_bioactivities.csv", index=False)
print(f"\n总记录: {len(bio)}")

# 每基因取 pChEMBL 最高的前 8 个分子（去重），并抓取 SMILES
cands = []
for gene in TARGETS:
    sub = bio[bio["gene"] == gene].sort_values("pchembl", ascending=False)
    sub = sub.drop_duplicates(subset=["molecule_chembl_id"]).head(8)
    for _, r in sub.iterrows():
        mid = r["molecule_chembl_id"]
        m = api(f"molecule/{mid}.json")
        smi = None
        if m and m.get("molecule_structures"):
            smi = m["molecule_structures"].get("canonical_smiles")
        cands.append({
            "gene": gene, "molecule_chembl_id": mid,
            "pref_name": r["pref_name"] if pd.notna(r["pref_name"])
                        else (m or {}).get("pref_name"),
            "pchembl": r["pchembl"],
            "standard_type": r["standard_type"], "standard_value": r["standard_value"],
            "standard_units": r["standard_units"],
            "max_phase": (m or {}).get("max_phase"),
            "smiles": smi, "document": r["document_chembl_id"]})
        time.sleep(0.25)
    print(f"{gene}: 入选 {len(sub)} 个候选")

cd = pd.DataFrame(cands)
cd.to_csv(f"{OUT}/chembl_candidates.csv", index=False)
print(f"\ncandidates: {len(cd)}; with SMILES: {cd['smiles'].notna().sum()}")
print(cd[["gene", "pref_name", "pchembl", "max_phase"]].to_string(index=False))
