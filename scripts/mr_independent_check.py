# -*- coding: utf-8 -*-
"""mr_independent_check.py — MR 关键统计量的独立重算（区别于 mr_eqtl_glycemic.py 的代码路径）
纪律要求：验证数字时须区分"与脚本输出一致"和"脚本本身正确"——本脚本从原始 gz 文件
经 subprocess+awk 独立解析，手写 Wald ratio 公式，与脚本产物对比。
"""
import math

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import subprocess

import pandas as pd
from scipy.stats import norm

BASE = _ROOT
MR = BASE + "/01_数据/MR"
OUT = BASE + "/03_结果/MR"


def sh(cmd):
    r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"cmd failed: {cmd}\n{r.stderr}")
    return r.stdout.strip()


def wald_ratio(beta_y, se_y, beta_x, se_x):
    ratio = beta_y / beta_x
    se = math.sqrt(se_y ** 2 + ratio ** 2 * se_x ** 2) / abs(beta_x)
    z = ratio / se
    p = 2 * (1 - norm.cdf(abs(z)))
    return ratio, se, z, p


print("=" * 72)
print("独立重算 1：DHX58 rs2277621 → T2D_Xue2018（FDR 显著对）")
print("=" * 72)
# eQTLGen significant 行
line = sh(f'zcat {MR}/eqtlgen_cis_significant_20181017.txt.gz | '
          f'awk -F"\\t" \'$2=="rs2277621" && $9=="DHX58"\' | sort -k1,1g | head -1')
f = line.split("\t")
A1, A2, Z, NrSamp = f[4], f[5], float(f[6]), float(f[12])
print(f"eQTLGen raw: rs2277621 DHX58 A1={A1} A2={A2} Z={Z} N={int(NrSamp)} P={f[0]}")

# esi 频率
esi = sh(f'zcat {MR}/eqtlgen_esi_freq.txt.gz | awk \'$2=="rs2277621"\'')
p2 = esi.split("\t")
freq_a1, esi_a1, esi_a2 = float(p2[6]), p2[4], p2[5]
print(f"esi raw:     A1={esi_a1} A2={esi_a2} freq_a1={freq_a1}")
assert (esi_a1, esi_a2) == (A1, A2), "allele mismatch eQTLGen vs esi"

# 官方公式换算
denom = math.sqrt(2 * freq_a1 * (1 - freq_a1) * (NrSamp + Z ** 2))
beta_x, se_x = Z / denom, 1 / denom
print(f"eQTL beta={beta_x:.6f} se={se_x:.6f} F={Z ** 2:.1f}")

# Xue 2018 原始行
xue = sh(f'zcat {MR}/T2D_Xue2018_GCST006867.tsv.gz | awk \'$3=="rs2277621"\'')
x = xue.split()
XA1, XA2, FRQ, B, SE, PX, NX = x[3], x[4], float(x[5]), float(x[6]), float(x[7]), x[8], x[9]
print(f"Xue raw:     A1={XA1} A2={XA2} frq_A1={FRQ:.4f} beta={B} se={SE} p={PX} N={NX}")
assert XA1 == A1 and XA2 == A2, "allele mismatch eQTLGen vs Xue"
assert abs(FRQ - freq_a1) < 0.2, "freq mismatch"

ratio, wse, wz, wp = wald_ratio(B, SE, beta_x, se_x)
print(f"INDEP Wald:  beta={ratio:.6f} se={wse:.6f} z={wz:.4f} p={wp:.6g} OR={math.exp(ratio):.4f}")

res = pd.read_csv(OUT + "/mr_results.csv")
r = res[(res["gene"] == "DHX58") & (res["outcome"] == "T2D_Xue2018_GCST006867")].iloc[0]
print(f"SCRIPT Wald: beta={r['wald_beta']:.6f} se={r['wald_se']:.6f} "
      f"z={r['wald_z']:.4f} p={r['wald_p']:.6g} OR={r['OR']:.4f} q={r['p_fdr']:.6g}")
ok = abs(ratio - r["wald_beta"]) < 1e-6 and abs(wp - r["wald_p"]) < 1e-9
print("VERDICT:", "MATCH ✓" if ok else "MISMATCH ✗")

print()
print("=" * 72)
print("独立重算 2：LDHA rs6486423 → FI（flipped 等位基因 + 频率互补）")
print("=" * 72)
line = sh(f'zcat {MR}/eqtlgen_cis_significant_20181017.txt.gz | '
          f'awk -F"\\t" \'$2=="rs6486423" && $9=="LDHA"\' | sort -k1,1g | head -1')
f = line.split("\t")
A1, A2, Z, NrSamp = f[4], f[5], float(f[6]), float(f[12])
esi = sh(f'zcat {MR}/eqtlgen_esi_freq.txt.gz | awk \'$2=="rs6486423"\'')
p2 = esi.split("\t")
freq_a1 = float(p2[6])
assert (p2[4], p2[5]) == (A1, A2)
denom = math.sqrt(2 * freq_a1 * (1 - freq_a1) * (NrSamp + Z ** 2))
beta_x, se_x = Z / denom, 1 / denom
print(f"eQTLGen raw: A1={A1} A2={A2} Z={Z} N={int(NrSamp)} freq_a1={freq_a1}")
print(f"eQTL beta={beta_x:.6f} se={se_x:.6f} F={Z ** 2:.1f}")

# FI 无 rsid 列 → 按 chr:pos 匹配
fi = sh(f'zcat {MR}/FI_GCST90002238.tsv.gz | awk -F"\\t" \'$1==11 && $2==18415057\'')
x = fi.split("\t")
FE, FO, FF, FB, FSE, FP = x[2], x[3], float(x[4]), float(x[5]), float(x[6]), x[7]
print(f"FI raw:      chr11:18415057 ea={FE} oa={FO} freq_ea={FF} beta={FB} se={FSE} p={FP}")
# 期望：effect_allele == a2(T) → flipped，freq_ea 应 ≈ 1-freq_a1
print(f"预期: ea==a2({A2}) → flipped; freq 核对 |{FF:.4f} - {1-freq_a1:.4f}| = {abs(FF-(1-freq_a1)):.4f}")
if FE == A2 and FO == A1:
    flip = -1.0
elif FE == A1 and FO == A2:
    flip = 1.0
else:
    raise AssertionError(f"unexpected FI alleles {FE}/{FO} vs {A1}/{A2}")
ratio, wse, wz, wp = wald_ratio(FB * flip, FSE, beta_x, se_x)
print(f"INDEP Wald:  beta={ratio:.6f} se={wse:.6f} z={wz:.4f} p={wp:.6g}")

r = res[(res["gene"] == "LDHA") & (res["outcome"] == "FI_GCST90002238")].iloc[0]
print(f"SCRIPT Wald: beta={r['wald_beta']:.6f} se={r['wald_se']:.6f} "
      f"z={r['wald_z']:.4f} p={r['wald_p']:.6g} allele={r['allele_status']}")
ok = abs(ratio - r["wald_beta"]) < 1e-6 and abs(wp - r["wald_p"]) < 1e-9
print("VERDICT:", "MATCH ✓" if ok else "MISMATCH ✗")

print()
print("=" * 72)
print("独立重算 3：BH-FDR 复核（statsmodels 独立调用）")
print("=" * 72)
from statsmodels.stats.multitest import multipletests
pvals = res["wald_p"].to_numpy()
q = multipletests(pvals, method="fdr_bh")[1]
print(f"max |Δq| = {abs(q - res['p_fdr'].to_numpy()).max():.2e}")
print(f"FDR<0.05 对数: {(res['p_fdr'] < 0.05).sum()} (脚本) / {(q < 0.05).sum()} (独立)")
