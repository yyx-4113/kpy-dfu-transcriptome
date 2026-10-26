# -*- coding: utf-8 -*-
"""
mr_eqtl_glycemic.py — 枢纽基因 cis-eQTL → 糖代谢结局 的孟德尔随机化（预注册替代方案 a）
项目：丙酮酰化（Kpy）与糖尿病足（DFU）

预注册设计（见 03_结果/MR可行性评估.md，2026-10-25 定稿）：
  暴露：eQTLGen significant cis-eQTL（2018-10-17 release，GRCh37）
        每枢纽基因取 1Mb 窗口内 P 最小的 lead SNP 作为唯一工具变量（Wald ratio / 单 IV MR）
        beta/se 按 eQTLGen 官方公式换算：beta = Z / sqrt(2f(1-f)(N+Z^2))，
        se = 1 / sqrt(2f(1-f)(N+Z^2))（README_cis_full_SMR，f=esi 文件等位基因频率，N=NrSamples）
        等位基因：Z 相对 AssessedAllele；esi 文件 A1=AssessedAllele、A2=OtherAllele（已逐 SNP 核对一致）
  结局：T2D GCST90018926（EUR；Xue 2018 GCST006867 作第二来源/敏感性）
        FI  GCST90002238（Lagou 2021，EUR 151,013）
  HbA1c：GWAS Catalog 无公开全文统计（已核实）→ 不用
  DFU 结局：无公开全文统计（Meng 2017 GCST004725 full_stats=false；FinnGen 无专用端点）
             → 不制造伪 MR，如实报告不可行（见 MR可行性评估.md）
  工具变量准入：F = Z^2 > 10；结局中 rsid/chr:pos 命中；等位基因可协调；频率差 |Δfreq| < 0.2
  多重检验：BH-FDR 跨全部 gene × outcome 检验
输出：03_结果/MR/*.csv、figures/fig_mr_forest.png
"""
import gzip

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import hashlib
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(_ROOT)
MR = BASE / "01_数据/MR"
OUT = BASE / "03_结果/MR"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

t0 = time.time()


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)
    with open(OUT / "mr.log", "a", encoding="utf-8") as fh:
        fh.write(str(msg) + "\n")


def sha256(path, block=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


# 预注册枢纽基因清单（ML 共识 6 基因 + meta-DEG KpyRN 显著 + Kpy 机器核心）
# 注意：ADH1B 与 PKM（PKM2 的基因符号）在 eQTLGen 无显著 cis-eQTL → 按预注册规则如实 SKIPPED
HUB_GENES = ["ADH1B", "ALDH3A2", "BPGM", "DHX58", "EP300", "FAHD2B", "GBP1",
             "HAT1", "LDHA", "PFKP", "PGAM1", "PKM", "SIRT3", "STAT1", "TPI1"]

log("=" * 70)
log("mr_eqtl_glycemic.py — 启动（eQTLGen lead cis-eQTL → T2D/FI，单 IV Wald ratio）")
log("=" * 70)

# ---------------- 0. MANIFEST ----------------
manifest_rows = []
for f in sorted(MR.glob("*.gz")):
    manifest_rows.append({"file": f.name, "bytes": f.stat().st_size,
                          "sha256": sha256(f)})
pd.DataFrame(manifest_rows).to_csv(MR / "MANIFEST_mr.csv", index=False)
log(f"MANIFEST_mr.csv: {len(manifest_rows)} files")

# ---------------- 1. 暴露：eQTLGen lead cis-eQTL ----------------
sig_file = MR / "eqtlgen_cis_significant_20181017.txt.gz"
eqtl = pd.read_csv(sig_file, sep="\t")
log(f"eQTLGen significant rows: {len(eqtl)}")

instruments = []
for g in HUB_GENES:
    sub = eqtl[eqtl["GeneSymbol"] == g]
    if not len(sub):
        log(f"  {g}: eQTLGen 无显著 cis-eQTL → SKIPPED（如实记录，不降级用非显著 eQTL）")
        continue
    lead = sub.sort_values("Pvalue").iloc[0]
    instruments.append({
        "gene": g, "rsid": lead["SNP"], "chr": int(lead["SNPChr"]),
        "pos": int(lead["SNPPos"]), "a1_assessed": lead["AssessedAllele"],
        "a2_other": lead["OtherAllele"], "zscore": float(lead["Zscore"]),
        "n_samples": int(lead["NrSamples"]), "n_cohorts": int(lead["NrCohorts"]),
        "pvalue": float(lead["Pvalue"]), "n_cis_sig": len(sub)})
iv = pd.DataFrame(instruments)
log(f"lead cis-eQTL instruments: {len(iv)}/{len(HUB_GENES)} 基因")

# ---------------- 2. 等位基因频率（SMR .esi 文件） ----------------
esi_file = MR / "eqtlgen_esi_freq.txt.gz"
esi_rows = []
rsid_set = set(iv["rsid"])
with gzip.open(esi_file, "rt") as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 7 and parts[1] in rsid_set:
            esi_rows.append({"rsid": parts[1], "esi_a1": parts[4],
                             "esi_a2": parts[5], "freq_a1": float(parts[6])})
esi = pd.DataFrame(esi_rows)
iv = iv.merge(esi, on="rsid", how="left")
# 等位基因一致性核对（esi A1 应 == AssessedAllele）
mismatch = iv[(iv["esi_a1"] != iv["a1_assessed"]) | (iv["esi_a2"] != iv["a2_other"])]
if len(mismatch):
    log(f"  !! 等位基因标注不一致: {mismatch['rsid'].tolist()} → 剔除")
    iv = iv[~iv["rsid"].isin(mismatch["rsid"])]
else:
    log("  esi 与 significant 文件等位基因标注逐 SNP 一致 ✓")

# 官方公式换算 beta/se（README_cis_full_SMR）
denom = np.sqrt(2 * iv["freq_a1"] * (1 - iv["freq_a1"]) * (iv["n_samples"] + iv["zscore"] ** 2))
iv["eqtl_beta"] = iv["zscore"] / denom
iv["eqtl_se"] = 1.0 / denom
iv["F"] = iv["zscore"] ** 2
weak = iv[iv["F"] < 10]
if len(weak):
    log(f"  !! F<10 剔除: {weak['gene'].tolist()}")
    iv = iv[~iv["gene"].isin(weak["gene"])]
log("  instruments (beta/se per eQTLGen official formula):")
log(iv[["gene", "rsid", "chr", "pos", "a1_assessed", "freq_a1", "zscore",
        "eqtl_beta", "eqtl_se", "F"]].to_string(index=False))
iv.to_csv(OUT / "mr_instruments.csv", index=False)

# ---------------- 3. 结局 GWAS 分块加载（内存安全：仅保留命中 SNP） ----------------
OUTCOMES = {
    "T2D_GCST90018926": "T2D_GCST90018926.tsv.gz",
    "FI_GCST90002238": "FI_GCST90002238.tsv.gz",
}
if (MR / "T2D_Xue2018_GCST006867.tsv.gz").exists():
    OUTCOMES["T2D_Xue2018_GCST006867"] = "T2D_Xue2018_GCST006867.tsv.gz"

pos_keys = set(zip(iv["chr"], iv["pos"]))
rsid_keys = set(iv["rsid"])


def load_outcome(path, name):
    """分块读取，仅保留工具变量命中的行（rsid 或 chr:pos）。
    自动探测分隔符（GWAS Catalog 文件为 TAB，Xue 2018 为空格）。"""
    with gzip.open(path, "rt") as fh:
        first = fh.readline()
    sep = "\t" if first.count("\t") >= first.count(" ") else None  # None → \s+
    keep = []
    n_rows = 0
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t") if sep == "\t" \
            else fh.readline().split()
        if sep == "\t":  # TAB 文件走 C 引擎（快）；python 引擎不支持 low_memory
            read_kw = dict(sep="\t", header=None, names=header, chunksize=500_000,
                           low_memory=False)
        else:             # 空白分隔文件（Xue 2018）走 python 引擎
            read_kw = dict(sep=r"\s+", header=None, names=header, chunksize=500_000,
                           engine="python")
        for chunk in pd.read_csv(fh, **read_kw):
            n_rows += len(chunk)
            cols = {c.lower(): c for c in chunk.columns}
            rsid_col = next((cols[c] for c in ["variant_id", "rsid", "snp", "rs_id"]
                             if c in cols), None)
            chr_col = next((cols[c] for c in ["chromosome", "chr"] if c in cols), None)
            pos_col = next((cols[c] for c in ["base_pair_location", "pos",
                                              "position"] if c in cols), None)
            m = pd.Series(False, index=chunk.index)
            if rsid_col is not None:
                m = m | chunk[rsid_col].astype(str).isin(rsid_keys)
            if chr_col is not None and pos_col is not None:
                key = list(zip(pd.to_numeric(chunk[chr_col], errors="coerce"),
                               pd.to_numeric(chunk[pos_col], errors="coerce")))
                m = m | pd.Series([k in pos_keys for k in key], index=chunk.index)
            if m.any():
                keep.append(chunk[m])
    df = pd.concat(keep) if keep else pd.DataFrame(columns=header)
    log(f"[{name}] 总行 {n_rows:,} → 命中 {len(df)}")
    return df


outcomes = {}
for name, fname in OUTCOMES.items():
    p = MR / fname
    if p.exists():
        outcomes[name] = load_outcome(p, name)
    else:
        log(f"[{name}] 文件缺失: {p}")

# ---------------- 4. 等位基因协调 + Wald ratio ----------------
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
rows, audit = [], []
for oname, odf in outcomes.items():
    if not len(odf):
        continue
    cols = {c.lower(): c for c in odf.columns}
    ea_col = next((cols[c] for c in ["effect_allele", "a1", "ea"] if c in cols), None)
    oa_col = next((cols[c] for c in ["other_allele", "a2", "oa"] if c in cols), None)
    b_col = next((cols[c] for c in ["beta", "effect", "or", "b"] if c in cols), None)
    se_col = next((cols[c] for c in ["standard_error", "se", "std err"] if c in cols), None)
    p_col = next((cols[c] for c in ["p_value", "pval", "p"] if c in cols), None)
    f_col = next((cols[c] for c in ["effect_allele_frequency", "eaf", "freq",
                                    "frq_a1"] if c in cols), None)
    rsid_col = next((cols[c] for c in ["variant_id", "rsid", "snp"] if c in cols), None)
    chr_col = next((cols[c] for c in ["chromosome", "chr"] if c in cols), None)
    pos_col = next((cols[c] for c in ["base_pair_location", "pos"] if c in cols), None)
    log(f"[{oname}] 列识别: rsid={rsid_col} ea={ea_col} oa={oa_col} beta={b_col} "
        f"se={se_col} p={p_col} freq={f_col}")

    for _, r in iv.iterrows():
        # 定位：优先 rsid，其次 chr:pos
        hit = pd.DataFrame()
        if rsid_col is not None:
            hit = odf[odf[rsid_col].astype(str) == str(r["rsid"])]
        if not len(hit) and chr_col is not None and pos_col is not None:
            hit = odf[(pd.to_numeric(odf[chr_col], errors="coerce") == r["chr"]) &
                      (pd.to_numeric(odf[pos_col], errors="coerce") == r["pos"])]
        if not len(hit):
            audit.append({"gene": r["gene"], "rsid": r["rsid"], "outcome": oname,
                          "status": "NOT_IN_OUTCOME"})
            log(f"  {r['gene']} {r['rsid']} 不在 {oname} → SKIPPED")
            continue
        h = hit.iloc[0]
        gwas_ea = str(h[ea_col]).upper() if ea_col else ""
        gwas_oa = str(h[oa_col]).upper() if oa_col else ""
        a1, a2 = r["a1_assessed"], r["a2_other"]

        # 等位基因协调：eQTL Z 相对 a1；GWAS beta 相对 effect_allele
        if gwas_ea == a1 and gwas_oa == a2:
            flip = 1.0; allele_status = "concordant"
        elif gwas_ea == a2 and gwas_oa == a1:
            flip = -1.0; allele_status = "flipped"
        elif gwas_ea == COMPLEMENT.get(a1) and gwas_oa == COMPLEMENT.get(a2):
            flip = 1.0; allele_status = "concordant_strandflip"
        elif gwas_ea == COMPLEMENT.get(a2) and gwas_oa == COMPLEMENT.get(a1):
            flip = -1.0; allele_status = "flipped_strandflip"
        else:
            audit.append({"gene": r["gene"], "rsid": r["rsid"], "outcome": oname,
                          "status": f"ALLELE_MISMATCH ea={gwas_ea} oa={gwas_oa} vs {a1}/{a2}"})
            log(f"  {r['gene']} {r['rsid']} 等位基因不匹配 ({gwas_ea}/{gwas_oa} vs {a1}/{a2}) → SKIPPED")
            continue

        # 频率一致性核查：eQTLGen freq 是对 assessed allele(a1) 报告的。
        # 优先级：先精确匹配 a1/a2，再做互补链匹配（避免回文 SNP 上
        # COMPLEMENT[a1]==a2 导致的误判，如 G/C 情形 comp(G)=C=a2）。
        freq_note = "NA"
        if f_col and not pd.isna(h[f_col]):
            gwas_freq_ea = float(h[f_col])
            if gwas_ea == a1:
                eqtl_freq_ea = r["freq_a1"]
            elif gwas_ea == a2:
                eqtl_freq_ea = 1 - r["freq_a1"]
            elif gwas_ea == COMPLEMENT.get(a1, "") and COMPLEMENT[a1] != a2:
                eqtl_freq_ea = r["freq_a1"]
            elif gwas_ea == COMPLEMENT.get(a2, "") and COMPLEMENT[a2] != a1:
                eqtl_freq_ea = 1 - r["freq_a1"]
            else:
                eqtl_freq_ea = np.nan
            if not np.isnan(eqtl_freq_ea):
                d = abs(gwas_freq_ea - eqtl_freq_ea)
                freq_note = f"{d:.3f}"
                if d > 0.2:
                    audit.append({"gene": r["gene"], "rsid": r["rsid"], "outcome": oname,
                                  "status": f"FREQ_DIFF {d:.3f}"})
                    log(f"  {r['gene']} {r['rsid']} 频率差 {d:.3f}>0.2 → SKIPPED")
                    continue

        gwas_beta = float(h[b_col]) * flip
        gwas_se = float(h[se_col])
        beta_ratio = gwas_beta / r["eqtl_beta"]
        wald_se = np.sqrt(gwas_se ** 2 + beta_ratio ** 2 * r["eqtl_se"] ** 2) / abs(r["eqtl_beta"])
        wald_z = beta_ratio / wald_se
        wald_p = 2 * (1 - stats.norm.cdf(abs(wald_z)))
        rows.append({
            "gene": r["gene"], "rsid": r["rsid"], "outcome": oname,
            "eqtl_beta": r["eqtl_beta"], "eqtl_se": r["eqtl_se"], "F": r["F"],
            "gwas_beta": gwas_beta, "gwas_se": gwas_se, "gwas_p": float(h[p_col]) if p_col else np.nan,
            "wald_beta": beta_ratio, "wald_se": wald_se,
            "wald_z": wald_z, "wald_p": wald_p,
            "allele_status": allele_status, "freq_diff": freq_note})
        audit.append({"gene": r["gene"], "rsid": r["rsid"], "outcome": oname,
                      "status": f"INCLUDED ({allele_status})"})

res = pd.DataFrame(rows)
pd.DataFrame(audit).to_csv(OUT / "mr_harmonization_audit.csv", index=False)

if len(res):
    res["p_fdr"] = multipletests(res["wald_p"], method="fdr_bh")[1]
    # OR 仅对 binary 结局（T2D，logOR 尺度）有意义；FI 为定量性状（ln 单位），不报 OR
    res["OR"] = np.where(res["outcome"].str.startswith("T2D"),
                         np.exp(res["wald_beta"]), np.nan)
    res = res.sort_values("p_fdr")
    res.to_csv(OUT / "mr_results.csv", index=False)
    log(f"\nMR 完成: {len(res)} 个基因-结局对")
    log(res[["gene", "outcome", "wald_beta", "wald_se", "wald_p", "p_fdr", "OR", "F",
             "allele_status"]].to_string(index=False))
    log(f"FDR<0.05: {len(res[res['p_fdr'] < 0.05])} 对; P<0.05: {len(res[res['wald_p'] < 0.05])} 对")
else:
    log("无有效 MR 结果")

# ---------------- 5. 森林图 ----------------
if len(res):
    fig, ax = plt.subplots(figsize=(9, max(3.5, 0.42 * len(res))))
    res_sorted = res.sort_values("wald_beta").reset_index(drop=True)
    colors = ["#C44E52" if p < 0.05 else "#4C72B0" for p in res_sorted["p_fdr"]]
    y = np.arange(len(res_sorted))
    ax.errorbar(res_sorted["wald_beta"], y, xerr=1.96 * res_sorted["wald_se"],
                fmt="o", color="black", ecolor="grey", markersize=4, capsize=2, lw=0)
    for i, (b, c) in enumerate(zip(res_sorted["wald_beta"], colors)):
        ax.scatter([b], [i], color=c, s=30, zorder=3)
    ax.axvline(0, color="grey", lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{g} → {o.replace('_GCST', '')}"
                        for g, o in zip(res_sorted["gene"], res_sorted["outcome"])], fontsize=8)
    ax.set_xlabel("Wald ratio effect (per unit cis-eQTL expression change; logOR scale for T2D)")
    ax.set_title("MR: hub gene cis-eQTL → glycemic outcomes (red: FDR<0.05, blue: FDR>=0.05)")
    fig.tight_layout()
    fig.savefig(FIG / "fig_mr_forest.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
    log("fig_mr_forest.png saved")

log(f"DONE — mr_eqtl_glycemic.py ({time.time() - t0:.0f}s)")
