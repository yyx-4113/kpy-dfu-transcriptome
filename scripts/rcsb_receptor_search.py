# -*- coding: utf-8 -*-
"""rcsb_receptor_search.py — 检索 RCSB PDB（人源、含小分子配体），
为 5 个 Kpy 轴靶点挑选共晶结构。输出 03_结果/药物预测/rcsb_structures.csv
"""
import json

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import subprocess
import time

import pandas as pd

OUT = _ROOT + "/03_结果/药物预测"

TARGETS = {
    "PKM2": "pyruvate kinase PKM2",
    "EP300": "histone acetyltransferase p300",
    "SIRT3": "sirtuin 3 mitochondrial",
    "LDHA": "lactate dehydrogenase A",
    "PDK2": "pyruvate dehydrogenase kinase isoform 2",
}


def curl_json(url, payload=None):
    cmd = ["curl", "-s", "--max-time", "90"]
    if payload is not None:
        cmd += ["-X", "POST", "-H", "Content-Type: application/json",
                "-d", json.dumps(payload)]
    cmd.append(url)
    r = subprocess.run(cmd, capture_output=True, text=True)
    try:
        return json.loads(r.stdout)
    except Exception:
        return None


rows = []
for gene, qtext in TARGETS.items():
    q = {
        "query": {
            "type": "group", "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "full_text",
                 "parameters": {"value": qtext}},
                {"type": "terminal", "service": "text",
                 "parameters": {"attribute": "rcsb_entity_source_organism.scientific_name",
                                "operator": "exact_match", "value": "Homo sapiens"}},
            ],
        },
        "request_options": {"paginate": {"start": 0, "rows": 20},
                            "results_content_type": ["experimental"]},
        "return_type": "polymer_entity",
    }
    d = curl_json("https://search.rcsb.org/rcsbsearch/v2/query", q)
    hits = [h["identifier"] for h in (d or {}).get("result_set", [])]
    print(f"{gene}: {len(hits)} hits")
    for ident in hits:
        pdb_id, entity_id = ident.split("_")
        e = curl_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
        if not e:
            continue
        info = e.get("rcsb_entry_info", {})
        p = curl_json(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}")
        desc = (p or {}).get("rcsb_polymer_entity", {}).get("pdbx_description", "")
        rows.append({
            "gene": gene, "pdb_id": pdb_id, "entity_id": entity_id,
            "description": desc,
            "resolution": info.get("resolution_combined", [None])[0],
            "method": info.get("experimental_method"),
            "n_ligands": info.get("nonpolymer_entity_count"),
            "title": e.get("struct", {}).get("title", "")[:90],
        })
        time.sleep(0.2)
    time.sleep(0.3)

df = pd.DataFrame(rows)
df["n_ligands"] = pd.to_numeric(df["n_ligands"], errors="coerce").fillna(0).astype(int)
df = df[df["n_ligands"] > 0].sort_values(["gene", "resolution"])
df.to_csv(f"{OUT}/rcsb_structures.csv", index=False)
print(f"\n人源含配体结构 {len(df)} 个:")
for g, sub in df.groupby("gene"):
    print(f"\n--- {g} ---")
    for _, r in sub.head(8).iterrows():
        print(f"  {r['pdb_id']}  {r['resolution']} Å  lig={r['n_ligands']:2d}  {r['description'][:60]}")
