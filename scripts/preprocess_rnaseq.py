# -*- coding: utf-8 -*-
"""
preprocess_rnaseq.py — RNA-seq 三队列预处理（GSE134431 / GSE199939 / GSE143735）
项目：丙酮酰化（Kpy）与糖尿病足（DFU）
产出：
  01_数据/处理后数据/<GSE>_log2_expr.csv    基因(行)×样本(列)，log2(x+1)
  01_数据/处理后数据/<GSE>_annotation.csv   样本注释（分组/愈合状态/GSM）
  03_结果/QC/<GSE>_pca.png, _corr.png       质控图
原则：不拼接跨平台矩阵；每队列独立标准化；样本映射以 GEO MINiML 为准。
"""
import io

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import os
import re
import tarfile
import urllib.request
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = _ROOT
RAW = f"{BASE}/01_数据/原始数据"
PROC = f"{BASE}/01_数据/处理后数据"
QC = f"{BASE}/03_结果/QC"
os.makedirs(PROC, exist_ok=True)
os.makedirs(QC, exist_ok=True)

# ----------------------------------------------------------------------
# 工具函数
# ----------------------------------------------------------------------
def qc_plots(mat: pd.DataFrame, gse: str, ann: pd.DataFrame, meta_col: str):
    """mat: gene x sample (log2); ann: index=sample, col=group"""
    groups = ann[meta_col] if meta_col in ann.columns else None
    # PCA
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    X = StandardScaler().fit_transform(mat.T)  # samples x genes
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)
    fig, ax = plt.subplots(figsize=(6, 5))
    if groups is not None:
        for g in groups.unique():
            idx = groups[groups == g].index
            ax.scatter(pcs[groups.index.get_indexer(idx), 0], pcs[groups.index.get_indexer(idx), 1], label=str(g), s=45)
        ax.legend()
    else:
        ax.scatter(pcs[:, 0], pcs[:, 1], s=45)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title(f"{gse} PCA (samples x genes)")
    fig.tight_layout()
    fig.savefig(f"{QC}/{gse}_pca.png", dpi=150)
    plt.close(fig)
    # 样本相关性
    corr = mat.corr(method="spearman")
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(corr.values, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(corr))); ax.set_xticklabels(corr.columns, rotation=90, fontsize=6)
    ax.set_yticks(range(len(corr))); ax.set_yticklabels(corr.index, fontsize=6)
    fig.colorbar(im, ax=ax, shrink=0.8)
    ax.set_title(f"{gse} sample Spearman correlation")
    fig.tight_layout()
    fig.savefig(f"{QC}/{gse}_corr.png", dpi=150)
    plt.close(fig)
    print(f"  [QC] PCA/corr saved; median inter-sample r = {np.median(corr.values[np.triu_indices(len(corr),1)]):.3f}")


# ----------------------------------------------------------------------
# 1. GSE134431（RPKM xlsx + Samples sheet + 本地 MINiML）
# ----------------------------------------------------------------------
def preprocess_gse134431():
    gse = "GSE134431"
    print(f"\n=== {gse} ===")
    f = f"{RAW}/{gse}/{gse}_181016_diabetic_rna-seq_results.gene.rpkm.xlsx"
    raw = pd.read_excel(f, sheet_name="181016_diabetic_rna-seq_results", header=None)
    # 定位 'RPKM (Partek)' 标记行；其下一行=样本名，再下一行=s-id，数据再下一行开始
    marker_row = None
    for i in range(min(15, len(raw))):
        if raw.iloc[i].astype(str).str.contains("RPKM").any():
            marker_row = i
            break
    assert marker_row is not None, "RPKM marker row not found"
    name_row = marker_row + 1
    sid_row = marker_row + 2
    data_start = marker_row + 3
    n_samples = 21  # 该数据集固定 21 个样本；其后为 IPA 注释列（Entrez/Location/Type），剔除
    sample_names = raw.iloc[name_row, 1:1 + n_samples].tolist()
    s_ids = raw.iloc[sid_row, 1:1 + n_samples].tolist()
    expr = raw.iloc[data_start:, 0:1 + n_samples].copy()
    expr.columns = ["gene"] + [str(s) for s in s_ids]
    expr = expr.set_index("gene")
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)  # xlsx 空格=0
    expr = expr.loc[expr.sum(axis=1) > 0]
    if expr.index.duplicated().any():
        expr = expr.groupby(level=0).mean()
    log2 = np.log2(expr + 1)

    # Samples sheet
    smp = pd.read_excel(f, sheet_name="Samples")
    smp.columns = ["s_id", "hihg_id", "external_id", "patient_id", "group", "n_align"]
    smp = smp.set_index("s_id")

    # MINiML -> GSM 标题/特征
    tgz = f"{BASE}/01_数据/元数据/GSE134431_family.xml.tgz"
    tf = tarfile.open(tgz, mode="r:gz")
    root = ET.fromstring(tf.extractfile("GSE134431_family.xml").read())
    NS = "{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}"
    gsm_info = {}
    for s in root.iter(NS + "Sample"):
        acc = s.findtext(NS + "Accession")
        title = s.findtext(NS + "Title") or ""
        chars = {}
        for ch in s.findall(".//" + NS + "Characteristics"):
            chars[ch.get("tag", "")] = (ch.text or "").strip()
        gsm_info[acc] = {"title": title, **chars}

    def key_from_external(ext: str):
        m = re.search(r"(DFS\d+|UM#\d+W\d+|P2#\d+)", str(ext))
        return m.group(1) if m else None

    def key_from_title(title: str):
        m = re.search(r"(DFS|Donor\s*#|UM#|P2#)", title)
        if "DFS" in title or "Donor" in title:
            n = re.search(r"#?(\d+)", title)
            return f"DFS{n.group(1)}" if n else None
        m = re.search(r"(UM#\d+W\d+|P2#\d+)", title)
        return m.group(1) if m else None

    # 建立 external-key -> gsm 映射
    title_map = {}
    for acc, info in gsm_info.items():
        k = key_from_title(info["title"])
        if k:
            title_map[k] = acc

    rows = []
    for s_id in log2.columns:
        ext = smp.loc[s_id, "external_id"] if s_id in smp.index else ""
        grp = smp.loc[s_id, "group"] if s_id in smp.index else "NA"
        k = key_from_external(ext)
        gsm = title_map.get(k, "")
        info = gsm_info.get(gsm, {})
        rows.append({
            "sample": s_id, "gsm": gsm, "external_id": ext,
            "group": grp,  # DFS=diabetic foot skin, DFU=ulcer
            "ulcer_or_skin": info.get("ulcer_or_skin", ""),
            "healer_or_nonhealer": info.get("healer_or_nonhealer", ""),
            "gsm_title": info.get("title", ""),
        })
    ann = pd.DataFrame(rows).set_index("sample")
    ann["group2"] = ann.apply(
        lambda r: r["group"] if r["group"] != "NA" else r["ulcer_or_skin"], axis=1)
    # 愈合状态（仅 DFU）
    ann["healing"] = ann["healer_or_nonhealer"]
    log2.to_csv(f"{PROC}/{gse}_log2_expr.csv")
    ann.to_csv(f"{PROC}/{gse}_annotation.csv", encoding="utf-8-sig")
    print(f"  matrix: {log2.shape}; groups: {ann['group2'].value_counts().to_dict()}; "
          f"healing(DFU): {ann[ann['group2']=='DFU']['healing'].value_counts().to_dict()}")
    qc_plots(log2, gse, ann, "group2")


# ----------------------------------------------------------------------
# 2. GSE199939（TPM 矩阵）
# ----------------------------------------------------------------------
def preprocess_gse199939():
    gse = "GSE199939"
    print(f"\n=== {gse} ===")
    f = f"{RAW}/{gse}/{gse}_gene.tpm.matrix.annot_-_21.txt.gz"
    df = pd.read_csv(f, sep="\t")
    meta_cols = ["gene_id", "gene_name", "description", "length"]
    samples = [c for c in df.columns if c not in meta_cols]
    expr = df.set_index("gene_name")[samples]
    expr = expr[~expr.index.astype(str).str.startswith("NA|")]
    expr = expr.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    expr = expr.loc[expr.sum(axis=1) > 0]
    if expr.index.duplicated().any():
        expr = expr.groupby(level=0).mean()
    log2 = np.log2(expr + 1)
    ann = pd.DataFrame({
        "sample": samples,
        "group2": ["DFU" if s.startswith("DW") else "control" for s in samples],
        "orig_label": samples,
    }).set_index("sample")
    log2.to_csv(f"{PROC}/{gse}_log2_expr.csv")
    ann.to_csv(f"{PROC}/{gse}_annotation.csv", encoding="utf-8-sig")
    print(f"  matrix: {log2.shape}; groups: {ann['group2'].value_counts().to_dict()}")
    print("  注：DW=diabetic wound(DFU), N=normal foot skin（MINiML characteristics 标注为 foot Skin，"
          "与系列标题 DFU vs control 一致，已在数据层核实）")
    qc_plots(log2, gse, ann, "group2")


# ----------------------------------------------------------------------
# 3. GSE143735（featureCounts per-sample counts，ENSG→symbol 经 Ensembl REST）
# ----------------------------------------------------------------------
def preprocess_gse143735():
    import json
    gse = "GSE143735"
    print(f"\n=== {gse} ===", flush=True)
    d = f"{RAW}/{gse}"
    files = sorted([x for x in os.listdir(d) if x.endswith(".COUNTS.txt.gz")])
    counts = {}
    for fn in files:
        gsm = fn.split("_")[0]
        df = pd.read_csv(os.path.join(d, fn), sep="\t", skiprows=1)
        counts[gsm] = df.set_index("Geneid").iloc[:, -1]  # 最后一列为 count
    mat = pd.DataFrame(counts)
    mat = mat.loc[mat.sum(axis=1) > 0]
    print(f"  raw counts: {mat.shape}", flush=True)

    # ENSG -> symbol（GENCODE v46 basic 预提取映射，避免在线 API 超时）
    mapping = json.load(open(f"{BASE}/01_数据/元数据/ensg_symbol_map.json", encoding="utf-8"))
    unmapped = [g for g in mat.index if g not in mapping]
    mat.index = [mapping.get(g, g) for g in mat.index]
    n_unmapped = len(unmapped)
    if mat.index.duplicated().any():
        mat = mat.groupby(level=0).sum()
    # CPM + log2
    lib = mat.sum(axis=0)
    cpm = mat.div(lib, axis=1) * 1e6
    log2 = np.log2(cpm + 1)

    ann = pd.DataFrame({
        "sample": list(log2.columns),
        "gsm": list(log2.columns),
    }).set_index("sample")
    # 分组来自样本级 MINiML characteristics（GSE143735_samples.csv）
    st = pd.read_csv(f"{BASE}/01_数据/元数据/样本表/GSE143735_samples.csv", encoding="utf-8-sig")
    st = st.set_index("gsm")
    def grp_of(gsm: str) -> str:
        ch = str(st.loc[gsm, "characteristics"]) if gsm in st.index else ""
        if "No Ulcer" in ch:
            return "DM_no_ulcer"
        if "Non Healer" in ch:
            return "DFU_nonhealer"
        if "Ulcer - Healer" in ch:
            return "DFU_healer"
        return "unknown"
    ann["group2"] = [grp_of(g) for g in log2.columns]
    log2.to_csv(f"{PROC}/{gse}_log2_expr.csv")
    ann.to_csv(f"{PROC}/{gse}_annotation.csv", encoding="utf-8-sig")
    print(f"  matrix: {log2.shape}; groups: {ann['group2'].value_counts().to_dict()}; "
          f"unmapped ENSG kept as-is: {int(n_unmapped)}")
    qc_plots(log2, gse, ann, "group2")


if __name__ == "__main__":
    import json  # noqa
    preprocess_gse134431()
    preprocess_gse199939()
    preprocess_gse143735()
    print("\nALL RNA-SEQ PREPROCESSING DONE")
