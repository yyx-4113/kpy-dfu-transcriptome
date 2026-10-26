# -*- coding: utf-8 -*-
"""
scrna_load_qc.py — GSE165816 单细胞 Stage 1：逐样本加载 + QC（磁盘分阶段版）
项目：丙酮酰化（Kpy）与糖尿病足（DFU）

数据：GSE165816（10x Genomics，54 样本、174,962 细胞（原始）；
     足/前臂皮肤 + PBMC × 非糖尿病/非溃疡糖尿病/DFU-healer/DFU-nonhealer）
     per-sample Cell Ranger CSV（首行=条码头，无标签列；数据行=基因符号,计数...）

本版修订（2026-10-25，相对初版，均已实测验证）：
  1. 【关键 bug 修复】初版 `barcodes = header[1:]` 假设表头有标签列，
     实测全部 54 个文件表头第 1 个字段即 16 位条码（部分带 -1 后缀），
     该写法会丢弃首个条码并使全部细胞计数错位一位。现改为整行作条码，
     并加断言 len(数据行列) == len(条码) + 1 硬校验。
  2. 【按实情记录】基因 ID 列本就是基因符号（ENSG 前缀 0 个），
     ens2sym 映射为空操作；仅对 ENSG 前缀 ID 生效，日志报告命中数。
  3. 【内存约束】整机 8 GB RAM（可用 ~1.4 GB），174,962 细胞全量稀疏矩阵
     无法常驻，且单文件稠密解析峰值可达 ~2 GB。改为：
     (a) 逐样本处理，样本间零驻留；
     (b) 按基因行分块（chunk=2000 基因）读取 + 稀疏 vstack，解析峰值 <500 MB；
     (c) 每样本 QC 后 save_npz 落盘。
  4. 重复符号按行求和（0/1 稀疏聚合矩阵实现），日志报告重复数。

输出：
  01_数据/处理后数据/scRNA/GSE165816/<GSM>_counts.npz   （cells × genes, CSR float32, 已过滤）
  01_数据/处理后数据/scRNA/GSE165816/<GSM>_genes.csv     （基因符号列表）
  03_结果/QC/scRNA_GSE165816_cell_qc.csv                 （逐细胞 QC 元数据，含过滤标记）
  03_结果/QC/scRNA_GSE165816_stage1_summary.json         （汇总 + 断言）
"""
import json

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

BASE = Path(_ROOT)
RAW = BASE / "01_数据/原始数据/GSE165816"
OUT = BASE / "01_数据/处理后数据/scRNA/GSE165816"
META = BASE / "01_数据/元数据/样本表"
QCDIR = BASE / "03_结果/QC"
OUT.mkdir(parents=True, exist_ok=True)

EXPECTED_TOTAL_CELLS = 174962  # 54 文件表头条码合计（已实测）
CHUNK_GENES = 2000             # 分块读取行数（基因数），控内存峰值
MIN_GENES, MAX_GENES, MAX_MITO = 200, 8000, 20.0

t0 = time.time()


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def read_sample_matrix(fp):
    """分块读取 per-sample CSV → (genes×cells CSR float32, barcodes)。
    表头行 = 条码头（无标签列）；数据行 = 基因符号 + 计数。
    """
    barcodes = pd.read_csv(fp, nrows=0).columns.to_numpy()
    n_bc = len(barcodes)
    chunks, genes = [], []
    for chk in pd.read_csv(fp, header=None, skiprows=1, chunksize=CHUNK_GENES,
                           dtype="float32", converters={0: str}):
        assert chk.shape[1] == n_bc + 1, \
            f"{fp.name}: chunk cols {chk.shape[1]} != barcodes+1 ({n_bc + 1})"
        genes.extend(chk.iloc[:, 0].tolist())
        chunks.append(sp.csr_matrix(chk.iloc[:, 1:].to_numpy(), dtype="float32"))
    X = chunks[0] if len(chunks) == 1 else sp.vstack(chunks, format="csr")
    return X, np.array(genes), barcodes


def aggregate_duplicate_symbols(X, genes):
    """重复符号按行求和（稀疏 0/1 聚合矩阵）；无重复时原样返回。"""
    uniq, inverse = np.unique(genes, return_inverse=True)
    if len(uniq) == len(genes):
        return X, genes, 0
    D = sp.csr_matrix((np.ones(len(genes), dtype="float32"),
                       (inverse, np.arange(len(genes)))),
                      shape=(len(uniq), len(genes)))
    return D @ X, uniq, int(len(genes) - len(uniq))


# ---------- Ensembl → symbol（仅对 ENSG 前缀生效；本数据集为符号，预期空操作） ----------
with open(BASE / "01_数据/元数据/ensg_symbol_map.json", encoding="utf-8") as fh:
    ens2sym = json.load(fh)
log(f"ensg->symbol map loaded: {len(ens2sym)} entries")

# ---------- 样本元数据 ----------
ann = pd.read_csv(META / "GSE165816_samples.csv")
ann.columns = [c.strip() for c in ann.columns]
meta = {}
for _, r in ann.iterrows():
    chars = dict(kv.split("=", 1) for kv in str(r["characteristics"]).split(" | "))
    meta[r["gsm"]] = {"tissue": chars.get("tissue", "NA"),
                      "disease": chars.get("disease", "NA")}
log(f"samples: {len(meta)}; groups: "
    f"{pd.Series([tuple(v.values()) for v in meta.values()]).value_counts().to_dict()}")

# ---------- 逐样本处理 ----------
files = sorted(RAW.glob("*counts.csv.gz"))
assert len(files) == 54, f"expect 54 files, got {len(files)}"

qc_rows = []
summary = {"n_files": len(files), "per_sample": {}}
ensg_hits_total = 0
dup_total = 0

for i, fp in enumerate(files):
    gsm = fp.name.split("_")[0]
    # 1) 分块读取（genes × cells）
    X, genes, barcodes = read_sample_matrix(fp)
    # 2) ENSG→symbol（预期空操作，按实情记录）
    is_ensg = np.array([str(g).startswith("ENSG") for g in genes])
    ensg_hits_total += int(is_ensg.sum())
    if is_ensg.any():
        genes = np.array([ens2sym.get(g, g) for g in genes])
    # 3) 重复符号求和
    X, genes, n_dup = aggregate_duplicate_symbols(X, genes)
    dup_total += n_dup
    # 4) QC 指标（按细胞=列）
    n_counts = np.asarray(X.sum(axis=0)).ravel()
    n_genes = X.getnnz(axis=0).astype(np.int32)
    mt_mask = np.array([str(g).upper().startswith("MT-") for g in genes])
    mt_counts = np.asarray(X[mt_mask, :].sum(axis=0)).ravel() if mt_mask.any() \
        else np.zeros(X.shape[1], dtype="float32")
    pct_mt = np.where(n_counts > 0, 100.0 * mt_counts / np.maximum(n_counts, 1), 0.0)
    # 5) 细胞过滤 + 基因过滤（本样本 ≥1 细胞）
    keep = (n_genes >= MIN_GENES) & (n_genes <= MAX_GENES) & (pct_mt <= MAX_MITO)
    gene_ncells = X.getnnz(axis=1)
    gkeep = gene_ncells >= 1
    Xf = X[gkeep][:, keep].T.tocsr()  # cells × genes
    genes_k = genes[gkeep]
    # 6) 落盘
    sp.save_npz(OUT / f"{gsm}_counts.npz", Xf)
    pd.Series(genes_k).to_csv(OUT / f"{gsm}_genes.csv", index=False, header=False)
    # 7) QC 记录（含被过滤细胞，供审计）
    for bc, nc, ng, pm, kp in zip(barcodes, n_counts, n_genes, pct_mt, keep):
        qc_rows.append({"gsm": gsm, "barcode": bc, "tissue": meta[gsm]["tissue"],
                        "disease": meta[gsm]["disease"], "n_counts": float(nc),
                        "n_genes": int(ng), "pct_mt": float(pm), "pass_qc": bool(kp)})
    summary["per_sample"][gsm] = {
        "tissue": meta[gsm]["tissue"], "disease": meta[gsm]["disease"],
        "cells_raw": int(len(barcodes)), "cells_pass": int(keep.sum()),
        "pass_rate": round(float(keep.mean()), 4),
        "median_n_genes": float(np.median(n_genes)),
        "median_pct_mt": round(float(np.median(pct_mt)), 3),
        "genes_kept": int(gkeep.sum()),
    }
    if (i + 1) % 10 == 0 or i == len(files) - 1:
        log(f"processed {i + 1}/{len(files)} samples "
            f"(last: {gsm}, {keep.sum()}/{len(barcodes)} pass; "
            f"{time.time() - t0:.0f}s)")

# ---------- 汇总与断言 ----------
qc_df = pd.DataFrame(qc_rows)
qc_df.to_csv(QCDIR / "scRNA_GSE165816_cell_qc.csv", index=False)

total_raw = int(qc_df.shape[0])
total_pass = int(qc_df["pass_qc"].sum())
assert total_raw == EXPECTED_TOTAL_CELLS, \
    f"total cells {total_raw} != expected {EXPECTED_TOTAL_CELLS}"

summary.update({
    "total_cells_raw": total_raw,
    "total_cells_pass": total_pass,
    "pass_rate_overall": round(total_pass / total_raw, 4),
    "ensg_prefixed_ids": ensg_hits_total,
    "id_system": "gene_symbol" if ensg_hits_total == 0 else "mixed",
    "duplicate_symbol_rows_summed": dup_total,
    "qc_rule": f"{MIN_GENES}<=n_genes<={MAX_GENES} & pct_mt<={MAX_MITO}",
    "runtime_sec": round(time.time() - t0, 1),
    "fixes": ["barcode header off-by-one (header[1:] -> header, col-count assertion)",
              "gene IDs are symbols already (ensg map is no-op)",
              "chunked read + sparse vstack (peak <500MB) + per-sample disk-backed stage"],
})
with open(QCDIR / "scRNA_GSE165816_stage1_summary.json", "w", encoding="utf-8") as fh:
    json.dump(summary, fh, ensure_ascii=False, indent=2)

log(f"DONE: {total_pass}/{total_raw} cells pass QC "
    f"({100 * total_pass / total_raw:.1f}%); ensg ids: {ensg_hits_total}; "
    f"dup rows: {dup_total}; total {time.time() - t0:.0f}s")
