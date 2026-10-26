# -*- coding: utf-8 -*-
"""
docking_pipeline.py — Kpy 轴靶点分子对接流水线（任务 #10）
项目：丙酮酰化（Kpy）与糖尿病足（DFU）

设计（遵循项目计划书 P9）：
  受体：人源 PDB 实验结构，带共晶配体（阳性对照即共晶配体本身）
  候选：ChEMBL 生物活性排序化合物（IC50/Kd/Ki，pChEMBL>=5）+ 文献已核实化合物
        （PubChem/ChEMBL 按名检索核验身份；活性值经 ChEMBL API 逐一核实）
  对照：每体系 1 个阳性对照（共晶配体）+ 1 个负面对照（蔗糖，无已知结合）
  对接：AutoDock Vina 1.2.7；盒子以共晶配体单链质心为中心；exhaustiveness=8；seed=42
  验证：阳性对照重对接 RMSD（<2 Å 视为协议可靠）；结合能以阳性对照为基准
  措辞纪律：对接结果仅提示"需实验验证的候选"，不作为结合证据的唯一依据

  体系（8 个）：
    PKM2  C-pocket  4JPG（1OX 激活剂）         PKM2  A-pocket  4JPG（FBP）
    EP300 AcCoA位点 4PZS（ACO）                SIRT3 NAD+位点  4BN4（ADP-核糖）
    LDHA 活性位点   4JNK（ZHK + NADH 保留）    PDK2 ATP口袋    4MPC（PV2）
    PDK2 丙酮酸位点 2BU8（DCA + ADP 保留）     PDK2 R域N端位点 2BU5（TF1）
  位点归属纪律（2026-10-25 核实）：
    - DCA 结合 PDK2 R 域丙酮酸位点（Knoechel 2006），非 ATP 口袋 → 移至 2BU8 体系
    - ChEMBL 派啶候选来自 2021 专利（功能数据、位点未注释），与 R 域配体 TF1/TF2
      共享丙酮酸等排体 warhead → 同时对接三个口袋，如实报告位点不确定性
    - M77976 无 ChEMBL PDK2 活性记录 → 剔除
输出：03_结果/药物预测/docking/*.csv、docking_results.csv、figures/
"""
import json

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import math
import os
import re
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(_ROOT)
PDB = BASE / "01_数据/分子对接/PDB"
TOOLS = BASE / "01_数据/工具"
WORK = BASE / "01_数据/分子对接/work"
OUT = BASE / "03_结果/药物预测/docking"
FIG = OUT / "figures"
for d in [WORK, OUT, FIG, WORK / "receptors", WORK / "ligands", WORK / "poses"]:
    d.mkdir(parents=True, exist_ok=True)

VINA = TOOLS / "vina.exe"
t0 = time.time()


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)
    with open(OUT / "docking.log", "a", encoding="utf-8") as fh:
        fh.write(str(msg) + "\n")


def sh(cmd):
    r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    return r.returncode, r.stdout, r.stderr


def sani(s):
    """文件命名的标签消毒：去掉 / \\ : * ? " < > | + 等非法字符。"""
    return re.sub(r"[^A-Za-z0-9._-]", "", s)


# ================= 1. 体系配置（全部可溯源） =================
# 受体信息来自 RCSB（2026-10-25 检索）；盒子中心=共晶配体质心（单链）
SYSTEMS = [
    {
        "target": "PKM2", "gene": "PKM", "pdb": "4JPG",
        "site": "C-pocket (allosteric)", "pos_lig": ("1OX", "602", "B"),
        "box": {"cx": 9.0, "cy": 4.4, "cz": 11.9, "sx": 22, "sy": 22, "sz": 22},
        "site_note": "LYS311/ALA388/PHE26/HIS29 (A-B 亚基界面变构位点)",
        "ligands": [
            ("TEPP-46", "Cn1c2cc([S+](C)[O-])sc2c2cnn(Cc3cccc(N)c3)c(=O)c21", "literature",
             "PKM2 变构激活剂（A/C pocket 结合型）"),
            ("shikonin", "CC(C)=CC[C@@H](O)C1=CC(=O)c2c(O)ccc(O)c2C1=O", "literature",
             "PKM2 不可逆抑制剂（C-pocket）"),
            ("CHEMBL5590652", "Clc1ccc(-c2ccn3c(-c4cnc5ccccc5n4)cnc3n2)cc1", "chembl",
             "EC50 0.6 nM (ChEMBL)"),
            ("CHEMBL5569778", "COc1ccccc1Nc1nc(-c2cnc3ccccn23)cs1", "chembl",
             "Kd 0.9 nM / IC50 23 nM (ChEMBL)"),
            ("CHEMBL5573729", "Nc1nc(-c2cnc3ccccn23)cs1", "chembl", "IC50 1.7 nM (ChEMBL)"),
        ],
    },
    {
        "target": "PKM2", "gene": "PKM", "pdb": "4JPG",
        "site": "A-pocket (FBP site)", "pos_lig": ("FBP", "601", "A"),
        "box": {"cx": -20.8, "cy": -16.5, "cz": 13.9, "sx": 22, "sy": 22, "sz": 22},
        "site_note": "THR432/SER434/SER437/LYS433/ARG516/GLY518 (FBP 位点)",
        "ligands": [
            ("TEPP-46", "Cn1c2cc([S+](C)[O-])sc2c2cnn(Cc3cccc(N)c3)c(=O)c21", "literature",
             "PKM2 变构激活剂（FBP 位点）"),
            ("shikonin", "CC(C)=CC[C@@H](O)C1=CC(=O)c2c(O)ccc(O)c2C1=O", "literature",
             "PKM2 抑制剂"),
            ("CHEMBL5590652", "Clc1ccc(-c2ccn3c(-c4cnc5ccccc5n4)cnc3n2)cc1", "chembl",
             "EC50 0.6 nM (ChEMBL)"),
            ("CHEMBL5569778", "COc1ccccc1Nc1nc(-c2cnc3ccccn23)cs1", "chembl",
             "Kd 0.9 nM / IC50 23 nM (ChEMBL)"),
        ],
    },
    {
        "target": "EP300", "gene": "EP300", "pdb": "4PZS",
        "site": "acetyl-CoA substrate site", "pos_lig": ("ACO", "1701", "A"),
        "box": {"cx": 14.4, "cy": 43.3, "cz": 1.4, "sx": 24, "sy": 26, "sz": 22},
        "site_note": "THR1411/SER1400/ARG1410/ILE1457/TRP1466/ARG1462 (AcCoA 位点)",
        "ligands": [
            ("A-485", "CNC(=O)Nc1ccc2c(c1)CC[C@@]21OC(=O)N(CC(=O)N(Cc2ccc(F)cc2)"
                      "[C@@H](C)C(F)(F)F)C1=O", "literature",
             "p300 催化抑制剂（临床候选）"),
            ("C646", "CC1=CC(=C(C=C1C)[N+](=O)[O-])C2=CC=C(O2)/C=C/C3=C(NN(C3=O)"
                     "C4=CC=C(C=C4)C(=O)O)C", "literature", "p300 HAT 抑制剂"),
            ("curcumin", "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
             "literature", "p300 抑制剂（泛、非特异）"),
            ("CHEMBL5948549", "CCNC(=O)c1ccc2c(C(=O)[C@H](NC[C@H](C)c3ccc(C#N)cc3)"
                              "c3ccccc3)c[nH]c2c1", "chembl", "IC50 8.8 nM (ChEMBL)"),
            ("CHEMBL6021203", "C[C@@H](CN[C@@H](C(=O)c1cnn2cc(-c3cnn(C)c3)ccc12)"
                              "c1ccccc1)c1ccc(C#N)cc1", "chembl", "IC50 9.4 nM (ChEMBL)"),
            ("CHEMBL5770805", "N#Cc1ccc(CCN[C@@H](C(=O)c2cnn3cc(-c4cn[nH]c4)ccc23)"
                              "c2ccccc2)cc1", "chembl", "IC50 15.6 nM (ChEMBL)"),
        ],
    },
    {
        "target": "SIRT3", "gene": "SIRT3", "pdb": "4BN4",
        "site": "ADP-ribose/NAD+ C-pocket", "pos_lig": ("AR6", "1395", "A"),
        "box": {"cx": 9.8, "cy": 3.3, "cz": 4.6, "sx": 24, "sy": 24, "sz": 24},
        "site_note": "SER321/THR320/ARG158/TYR165/ALA146/PHE157/HIS248 (NAD+ 位点)",
        "ligands": [
            ("nicotinamide", "NC(=O)c1cccnc1", "literature", "烟酰胺（NAM，产物抑制）"),
            ("honokiol", "C=CCc1ccc(O)c(-c2ccc(O)c(CC=C)c2)c1", "literature",
             "文献报道 SIRT3 激活剂"),
            ("resveratrol", "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1", "literature",
             "sirtuin 调节剂（争议性）"),
            ("quercetin", "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12", "literature",
             "sirtuin 调节剂"),
            ("sirtinol", "CC(NC(=O)c1ccccc1/N=C/c1c(O)ccc2ccccc12)c1ccccc1",
             "literature", "sirtuin 抑制剂"),
            ("CHEMBL2332041", "CC(C)(C)OC(=O)c1ccc(C(=O)NCCC2CCN(c3ncnc4cc(C(N)=O)"
                              "sc34)CC2)s1", "chembl", "IC50 3.2 nM (ChEMBL)"),
            ("CHEMBL2332039", "CCNC(=O)c1ccc(C(=O)NCCC2CCN(c3ncnc4cc(C(N)=O)sc34)"
                              "CC2)s1", "chembl", "IC50 4.0 nM (ChEMBL)"),
            ("CHEMBL2338810", "CS(=O)(=O)NCCC1CCN(c2ncnc3cc(C(N)=O)sc23)CC1",
             "chembl", "IC50 7.2 nM (ChEMBL)"),
        ],
    },
    {
        "target": "LDHA", "gene": "LDHA", "pdb": "4JNK",
        "site": "active site", "pos_lig": ("ZHK", "803", "A"),
        "box": {"cx": 43.3, "cy": 21.3, "cz": 48.1, "sx": 20, "sy": 20, "sz": 20},
        "site_note": "ARG105/HIS192/ASP140/GLU191 (活性位点)",
        "ligands": [
            ("oxamic acid", "NC(=O)C(=O)O", "literature",
             "经典 LDHA 抑制剂（底物位点）"),
            ("FX11", "CCCc1c(O)c(O)c(C(=O)O)c2cc(Cc3ccccc3)c(C)cc12", "literature",
             "gossypol 衍生物 LDHA 抑制剂"),
            ("galloflavin", "O=C1Cc2cc(O)c(=O)oc2-c2c1cc(O)c(O)c2O", "literature",
             "NADH 结合位点 LDHA 抑制剂"),
            ("CHEMBL4859303", "NS(=O)(=O)c1ccc(Cc2c(-c3cccc(OCc4ccccc4)c3)nn(-c3nc("
                              "C(=O)O)cs3)c2CC2CC2)cc1F", "chembl", "Kd 0.016 nM / IC50 19 nM (ChEMBL)"),
            ("CHEMBL4782843", "NS(=O)(=O)c1ccc(Cc2c(-c3ccc(F)c(C#CC4CCCC4)c3)nn(-c3nc("
                              "C(=O)O)cs3)c2CC2CC2)cc1F", "chembl", "Kd 0.019 nM / IC50 22 nM (ChEMBL)"),
            ("CHEMBL4850573", "NS(=O)(=O)c1ccc(Cc2c(-c3ccc(F)c(-c4ccc(Oc5nn[nH]c5C(=O)O)"
                              "cc4)c3)nn(-c3nc(C(=O)O)cs3)c2CC2CC2)cc1F", "chembl",
             "IC50 0.07 nM (ChEMBL)"),
        ],
    },
    {
        # ATP 口袋：Tso et al. 2014 JBC 289:4432（4MPC, 1.70 A），PS2/PV2 共晶 ATP 竞争性抑制剂
        "target": "PDK2", "gene": "PDK2", "pdb": "4MPC",
        "site": "ATP pocket", "pos_lig": ("PV2", "502", "A"),
        "box": {"cx": 13.9, "cy": 2.1, "cz": -37.3, "sx": 20, "sy": 20, "sz": 20},
        "site_note": "LEU252/ASN255/ARG258/ALA259/GLU262/ASP290/GLY294/VAL295/THR354 (ATP 口袋)",
        "ligands": [
            ("CHEMBL5896437", "C[C@@H]1CN(C(=O)[C@@](C)(O)C(F)(F)F)CC[C@@H]1C1("
                              "c2nnc(-c3ccc(Cl)cc3)o2)CC1", "chembl",
             "Kd 2.4 nM / IC50 130 nM (ChEMBL, 专利；位点未注释)"),
            ("CHEMBL6011956", "COc1ccc(-c2nnc(C3([C@H]4CCN(C(=O)[C@@](C)(O)C(F)(F)F)"
                              "C[C@H]4C)CC3)[nH]2)cc1", "chembl",
             "Kd 2.5 nM / IC50 88 nM (ChEMBL, 专利；位点未注释)"),
            ("CHEMBL5777543", "Cc1ccc(-c2nnc(C3([C@H]4CCN(C(=O)[C@@](C)(O)C(F)(F)F)"
                              "C[C@H]4C)CC3)[nH]2)cc1", "chembl",
             "Kd 2.6 nM / IC50 90 nM (ChEMBL, 专利；位点未注释)"),
        ],
    },
    {
        # 丙酮酸位点（R 域中心）：Knoechel et al. 2006 Biochemistry 45:402（2BU8, 2.50 A），
        # TF4 = 二氯乙酸（DCA）共晶；ADP 同在 ATP 位点（质心距 DCA 27.2 A，盒外，保留以维持晶体构象）
        "target": "PDK2", "gene": "PDK2", "pdb": "2BU8",
        "site": "pyruvate site (R-domain)", "pos_lig": ("TF4", "1379", "A"),
        "box": {"cx": 56.35, "cy": 44.81, "cz": 80.86, "sx": 20, "sy": 20, "sz": 20},
        "site_note": "LEU53/TYR80/ILE111/HIS115/ARG154 (R 域丙酮酸位点)",
        "ligands": [
            ("pyruvic acid", "CC(=O)C(=O)O", "literature",
             "丙酮酸位点生理配体（PDK2 激动剂）"),
            ("CHEMBL5896437", "C[C@@H]1CN(C(=O)[C@@](C)(O)C(F)(F)F)CC[C@@H]1C1("
                              "c2nnc(-c3ccc(Cl)cc3)o2)CC1", "chembl",
             "Kd 2.4 nM / IC50 130 nM (ChEMBL, 专利；位点未注释)"),
            ("CHEMBL6011956", "COc1ccc(-c2nnc(C3([C@H]4CCN(C(=O)[C@@](C)(O)C(F)(F)F)"
                              "C[C@H]4C)CC3)[nH]2)cc1", "chembl",
             "Kd 2.5 nM / IC50 88 nM (ChEMBL, 专利；位点未注释)"),
            ("CHEMBL5777543", "Cc1ccc(-c2nnc(C3([C@H]4CCN(C(=O)[C@@](C)(O)C(F)(F)F)"
                              "C[C@H]4C)CC3)[nH]2)cc1", "chembl",
             "Kd 2.6 nM / IC50 90 nM (ChEMBL, 专利；位点未注释)"),
        ],
    },
    {
        # R 域 N 端位点：Knoechel et al. 2006（2BU5, 2.35 A），TF1 共晶；
        # TF1 与专利派啶候选共享 (2R)-3,3,3-三氟-2-羟基-2-甲基丙酰基（丙酮酸等排体）
        "target": "PDK2", "gene": "PDK2", "pdb": "2BU5",
        "site": "R-domain N-terminal site", "pos_lig": ("TF1", "1386", "A"),
        "box": {"cx": 56.5, "cy": 60.92, "cz": 78.74, "sx": 22, "sy": 22, "sz": 22},
        "site_note": "LEU23/GLN27/PHE28/PHE31/SER41/PHE44/LEU45/LEU160/GLN163 (R 域 N 端位点)",
        "ligands": [
            ("CHEMBL5896437", "C[C@@H]1CN(C(=O)[C@@](C)(O)C(F)(F)F)CC[C@@H]1C1("
                              "c2nnc(-c3ccc(Cl)cc3)o2)CC1", "chembl",
             "Kd 2.4 nM / IC50 130 nM (ChEMBL, 专利；位点未注释)"),
            ("CHEMBL6011956", "COc1ccc(-c2nnc(C3([C@H]4CCN(C(=O)[C@@](C)(O)C(F)(F)F)"
                              "C[C@H]4C)CC3)[nH]2)cc1", "chembl",
             "Kd 2.5 nM / IC50 88 nM (ChEMBL, 专利；位点未注释)"),
            ("CHEMBL5777543", "Cc1ccc(-c2nnc(C3([C@H]4CCN(C(=O)[C@@](C)(O)C(F)(F)F)"
                              "C[C@H]4C)CC3)[nH]2)cc1", "chembl",
             "Kd 2.6 nM / IC50 90 nM (ChEMBL, 专利；位点未注释)"),
        ],
    },
]

NEGATIVE = ("sucrose", "OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)"
                       "[C@H]2O)[C@@H](O)[C@@H]1O")

log("=" * 70)
log("docking_pipeline.py — Kpy 轴分子对接（Vina 1.2.7）")
log("=" * 70)

# ================= 2. 受体准备 =================
# 金属/辅因子保留策略（科学依据，写入稿件方法学）：
#   4JNK (LDHA)：保留 NAI —— NADH 辅因子，构成底物/抑制剂口袋壁（与 ZHK 质心距 11.1 A）
#   2BU8 (PDK2)：保留 ADP —— 与 DCA 共晶的核苷酸（质心距 DCA 27.2 A，在对接盒外；
#                保留以维持 Knoechel 2006 所述 DCA+ADP 构象状态）
#   4BN4 (SIRT3)：结构性 ZN 距 AR6 口袋质心 34.0 A，且 obabel 无法处理裸金属行
#               （kekulize 失败，见工作日志）→ 不保留；对对接盒无影响
#   其余 HETATM（HOH/DMS/EPE/SO4/LAC/TLA/OP2/NA/MG/共晶配体）一律移除
#   （共晶配体另行提取为阳性对照，不得残留在受体中，否则重对接验证失效）
KEEP_HETS = {
    "4JNK": ("NAI",),
    "2BU8": ("ADP",),
}


def prep_receptor(pdb_id, drop_hets=("HOH",)):
    """提取蛋白链 + 指定保留的金属/辅因子，去配体/水，加 H，转 PDBQT。"""
    src = PDB / f"{pdb_id}.pdb"
    keep = KEEP_HETS.get(pdb_id, ())
    prot = WORK / "receptors" / f"{pdb_id}_protein.pdb"
    lines = []
    for line in open(src, encoding="utf-8", errors="ignore"):
        if line.startswith("ATOM"):
            lines.append(line)
        elif line.startswith("HETATM"):
            rn = line[17:20].strip()
            if rn in keep:
                lines.append(line)
        elif line.startswith("END"):
            break
    prot.write_text("".join(lines), encoding="utf-8")
    out = WORK / "receptors" / f"{pdb_id}.pdbqt"
    if out.exists() and out.stat().st_size > 10000:
        n_at = sum(1 for l in open(out, encoding="utf-8", errors="ignore")
                   if l.startswith("ATOM") or l.startswith("HETATM"))
        log(f"  受体 {pdb_id}: 已存在，跳过重转换（{n_at} atoms, "
            f"{out.stat().st_size/1e6:.1f} MB）")
        return out
    cmd = (f'obabel -ipdb "{prot}" -h -O "{out}" --partialcharge gasteiger '
           f'-p 7.4 2>/dev/null')
    rc, so, se = sh(cmd)
    if rc != 0 or not out.exists():
        log(f"  !! 受体转换失败 {pdb_id}: {se[:200]}")
        return None
    # obabel 会在 PDBQT 写入 ROOT/ENDROOT/BRANCH/ENDBRANCH/TORSDOF 等柔性对接标签，
    # 而 Vina 1.2.7 的刚性受体解析器会拒绝 ROOT 标签。刚性对接只需 ATOM/HETATM 行，
    # 故剥离这些柔性标签，生成干净的刚性受体 PDBQT（保留原子电荷与元素列）。
    raw = out.read_text(encoding="utf-8", errors="ignore").splitlines()
    rigid = [l for l in raw if l.startswith(("ATOM", "HETATM"))]
    out.write_text("\n".join(rigid) + "\n", encoding="utf-8")
    n_at = len(rigid)
    kept = f"，保留 {'/'.join(keep)}" if keep else ""
    log(f"  受体 {pdb_id}: {n_at} atoms → 刚性 PDBQT{kept}")
    return out


log("\n[1] 受体准备")
receptors = {}
for sysd in SYSTEMS:
    pdb = sysd["pdb"]
    if pdb not in receptors:
        receptors[pdb] = prep_receptor(pdb)

# ================= 3. 配体准备 =================
def smiles_to_pdbqt(name, smiles, tag):
    """SMILES → 3D (RDKit+MMFF) → meeko PDBQT"""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "SMILES 解析失败"
    mol = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 42
    if AllChem.EmbedMolecule(mol, ps) != 0:
        return None, "3D 嵌入失败"
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        pass
    prep = MoleculePreparation()
    try:
        setups = prep(mol)
    except Exception as e:
        return None, f"meeko prepare 失败: {e}"
    if not setups:
        return None, "meeko 无输出"
    pdbqt_string, is_ok, err = PDBQTWriterLegacy().write_string(setups[0])
    if not is_ok:
        return None, f"meeko: {err}"
    pdbqt = WORK / "ligands" / f"{tag}.pdbqt"
    pdbqt.write_text(pdbqt_string, encoding="utf-8")
    return pdbqt, None


def crystal_ligand_pdbqt(pdb_id, resn, resi, ch, tag):
    """从 PDB 提取共晶配体 → 干净的刚性 PDBQT（保留晶体坐标，残基名统一 UNL）"""
    src = PDB / f"{pdb_id}.pdb"
    lig_pdb = WORK / "ligands" / f"{tag}.pdb"
    lines = ["REMARK   crystal ligand extracted from PDB %s %s%s%s\n" % (pdb_id, resn, resi, ch)]
    for line in open(src, encoding="utf-8", errors="ignore"):
        if line.startswith("HETATM"):
            rn, c, ri = line[17:20].strip(), line[21], line[22:27].strip()
            if rn == resn and ri == resi and c == ch:
                lines.append(line)
    lines.append("END\n")
    lig_pdb.write_text("".join(lines), encoding="utf-8")
    raw_pdbqt = WORK / "ligands" / f"{tag}.raw.pdbqt"
    cmd = f'obabel -ipdb "{lig_pdb}" -h -O "{raw_pdbqt}" --partialcharge gasteiger -p 7.4 2>/dev/null'
    rc, so, se = sh(cmd)
    if rc != 0 or not raw_pdbqt.exists():
        return None, f"obabel 失败: {se[:150]}"
    # obabel 产物含 ROOT/ENDROOT/BRANCH/ENDBRANCH/TORSDOF 柔性标签，且共晶配体残基名
    # 可能以数字开头（如 1OX），二者均使 Vina 1.2.7 解析失败。处理：
    # ① 仅保留 ATOM/HETATM（含加 H 后的 H 原子）；② 残基名统一改为 UNL；
    # ③ 重新包成干净的刚性配体 ROOT…ENDROOT + TORSDOF 0（Vina 所需格式）。
    raw = raw_pdbqt.read_text(encoding="utf-8", errors="ignore").splitlines()
    atoms = [l[:17] + "UNL" + l[20:] for l in raw if l.startswith(("ATOM", "HETATM"))]
    if not atoms:
        return None, "obabel 未输出原子"
    rigid = "ROOT\n" + "\n".join(atoms) + "\nENDROOT\nTORSDOF 0\n"
    pdbqt = WORK / "ligands" / f"{tag}.pdbqt"
    pdbqt.write_text(rigid, encoding="utf-8")
    raw_pdbqt.unlink(missing_ok=True)
    return pdbqt, None


log("\n[2] 配体准备")
lig_files = {}   # tag -> path
failures = []
# 阳性对照（共晶配体）
for sysd in SYSTEMS:
    resn, resi, ch = sysd["pos_lig"]
    tag = sani(f"POS_{sysd['target']}_{sysd['site'].split()[0]}_{resn}")
    p, err = crystal_ligand_pdbqt(sysd["pdb"], resn, resi, ch, tag)
    if p:
        lig_files[tag] = p
        log(f"  阳性对照 {tag}: {p.name}")
    else:
        failures.append((tag, err))
        log(f"  !! 阳性对照 {tag} 失败: {err}")
# 负面对照
p, err = smiles_to_pdbqt(NEGATIVE[0], NEGATIVE[1], "NEG_sucrose")
if p:
    lig_files["NEG_sucrose"] = p
    log("  负面对照 NEG_sucrose ✓")
else:
    failures.append(("NEG_sucrose", err))
# 候选化合物
for sysd in SYSTEMS:
    for name, smiles, src, note in sysd["ligands"]:
        tag = sani(f"LIG_{name}")
        if tag in lig_files:
            continue
        p, err = smiles_to_pdbqt(name, smiles, tag)
        if p:
            lig_files[tag] = p
        else:
            failures.append((tag, err))
            log(f"  !! {name} 配体准备失败: {err}")
log(f"  配体就绪: {len(lig_files)} 个；失败 {len(failures)} 个")

# ================= 4. Vina 对接 =================
log("\n[3] Vina 对接")
results = []
for sysd in SYSTEMS:
    pdb = sysd["pdb"]
    rec = receptors.get(pdb)
    if rec is None:
        continue
    b = sysd["box"]
    resn, resi, ch = sysd["pos_lig"]
    pos_tag = sani(f"POS_{sysd['target']}_{sysd['site'].split()[0]}_{resn}")
    site_key = sani(sysd["site"].split()[0])
    # 该体系的对接对象：阳性对照 + 候选 + 负面对照
    tags = [pos_tag, "NEG_sucrose"]
    for name, smiles, src, note in sysd["ligands"]:
        t = sani(f"LIG_{name}")
        if t in lig_files:
            tags.append(t)
    for tag in tags:
        lig = lig_files[tag]
        out_pose = WORK / "poses" / f"{pdb}_{site_key}_{tag}.pdbqt"
        logp = WORK / "poses" / f"{pdb}_{site_key}_{tag}.log"
        # Vina 1.2.7 已移除 --log 选项；把 stdout（含结果表）重定向到 log 文件
        cmd = (f'"{VINA}" --receptor "{rec}" --ligand "{lig}" '
               f'--center_x {b["cx"]} --center_y {b["cy"]} --center_z {b["cz"]} '
               f'--size_x {b["sx"]} --size_y {b["sy"]} --size_z {b["sz"]} '
               f'--exhaustiveness 8 --num_modes 9 --energy_range 3 --seed 42 '
               f'--out "{out_pose}" '
               f'> "{logp}" 2>&1')
        rc, so, se = sh(cmd)
        # 解析 log 中的最佳亲和力（Vina 结果表打印到 stdout → 已重定向到 logp）
        aff = None
        rmsd_lb = None
        if logp.exists() and rc == 0:
            txt = logp.read_text(errors="ignore")
            # RMSD 列可能为 0（整数）或 0.000，故用 (\d+(?:\.\d+)?) 兼容两种写法
            m = re.findall(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)", txt, re.M)
            if m:
                aff = float(m[0][0]); rmsd_lb = float(m[0][1])
        else:
            log(f"  !! Vina 失败 {tag}: rc={rc} {se[:120]}")
        results.append({
            "target": sysd["target"], "gene": sysd["gene"], "pdb": pdb,
            "site": sysd["site"], "ligand": tag, "role": (
                "positive_control" if tag == pos_tag else
                "negative_control" if tag == "NEG_sucrose" else "candidate"),
            "affinity_kcal_mol": aff, "rmsd_lower_bound": rmsd_lb,
            "box": f'{b["cx"]},{b["cy"]},{b["cz"]} {b["sx"]}x{b["sy"]}x{b["sz"]}',
        })
        log(f"  {pdb} {sysd['site'][:12]:12s} {tag:42s} aff={aff}")

df = pd.DataFrame(results)
df.to_csv(OUT / "docking_raw.csv", index=False)
log(f"\n对接完成: {len(df)} runs ({time.time()-t0:.0f}s)")

# ================= 5. 阳性对照 RMSD（Kabsch） =================
log("\n[4] 阳性对照重对接 RMSD")


def parse_pdbqt_atoms(path, first_model_only=False):
    """解析 PDBQT 原子坐标。

    first_model_only=True 时只读取 pose 文件的 MODEL 1（Vina 已按能量排序，
    为最优 pose），避免把 9 个 pose 的全部原子累加导致参考/姿态原子数不匹配
    （2026-10-25 修复：原实现漏算 1OX/TF1(2 MODEL)、PV2/TF4(9 MODEL) 的 RMSD）。
    参考晶体坐标文件为单块（无 MODEL），first_model_only 对其无效（恒读全块）。
    """
    atoms = []
    model_seen = 0
    for line in open(path, encoding="utf-8", errors="ignore"):
        if line.startswith("MODEL"):
            model_seen += 1
            if first_model_only and model_seen > 1:
                break
            continue
        if line.startswith("ENDMDL"):
            if first_model_only and model_seen >= 1:
                break
            continue
        if line.startswith(("ATOM", "HETATM")):
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                el = line[77:79].strip()
                nm = line[12:16].strip()
                atoms.append((x, y, z, el, nm))
            except ValueError:
                continue
    return atoms


def kabsch_rmsd(P, Q):
    P = np.asarray(P); Q = np.asarray(Q)
    assert P.shape == Q.shape, f"shape mismatch {P.shape} vs {Q.shape}"
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return float(np.sqrt(((Pc - Qc @ R) ** 2).sum() / len(P)))


rmsd_rows = []
for sysd in SYSTEMS:
    resn, resi, ch = sysd["pos_lig"]
    pos_tag = sani(f"POS_{sysd['target']}_{sysd['site'].split()[0]}_{resn}")
    site_key = sani(sysd["site"].split()[0])
    ref = WORK / "ligands" / f"{pos_tag}.pdbqt"          # 晶体坐标（经 obabel 转换）
    pose = WORK / "poses" / f"{sysd['pdb']}_{site_key}_{pos_tag}.pdbqt"
    if not (ref.exists() and pose.exists()):
        continue
    A = [a for a in parse_pdbqt_atoms(ref) if a[3] != "H"]
    B = [a for a in parse_pdbqt_atoms(pose, first_model_only=True) if a[3] != "H"]
    nA, nB = len(A), len(B)
    if nA != nB or nA < 5:
        rmsd_rows.append({"target": sysd["target"], "site": sysd["site"],
                          "pos_lig": resn, "n_atoms": min(nA, nB), "rmsd": None,
                          "note": f"atom count mismatch (ref={nA}, pose={nB})"})
        continue
    # 参考与姿态同源（同一输入配体，原子顺序一致）→ 按索引配对计算 RMSD；
    # 若 Vina 改写原子名则记录警告但不阻断（索引配对仍有效）
    name_ok = all(A[i][4] == B[i][4] for i in range(nA))
    if not name_ok:
        log(f"  ! {sysd['target']} {sysd['site'][:20]} {resn}: 原子名被 Vina 改写"
            f"（ref={A[0][4]}…），按索引配对计算 RMSD")
    rmsd = kabsch_rmsd([a[:3] for a in A], [a[:3] for a in B])
    rmsd_rows.append({"target": sysd["target"], "site": sysd["site"],
                      "pos_lig": resn, "n_atoms": nA, "rmsd": round(rmsd, 2),
                      "note": "heavy-atom RMSD, model 1, index-paired"})
    log(f"  {sysd['target']} {sysd['site'][:20]:20s} {resn}: RMSD={rmsd:.2f} Å (n={nA})")

rdf = pd.DataFrame(rmsd_rows)
rdf.to_csv(OUT / "positive_control_rmsd.csv", index=False)

# ================= 6. 汇总 =================
# 每个体系的基准：阳性对照亲和力
bench = (df[df["role"] == "positive_control"]
         .set_index(["target", "site"])["affinity_kcal_mol"].to_dict())
df["pos_control_affinity"] = df.apply(
    lambda r: bench.get((r["target"], r["site"])), axis=1)
df["delta_vs_pos"] = df["affinity_kcal_mol"] - df["pos_control_affinity"]
df.to_csv(OUT / "docking_results.csv", index=False)

log("\n[5] 汇总（按体系）")
for (tgt, site), sub in df.groupby(["target", "site"]):
    log(f"--- {tgt} / {site} ---")
    for _, r in sub.sort_values("affinity_kcal_mol").iterrows():
        log(f"   {r['role']:18s} {r['ligand']:44s} {r['affinity_kcal_mol']:8.2f} "
            f"Δ={r['delta_vs_pos']:+6.2f}")

if failures:
    pd.DataFrame(failures, columns=["tag", "error"]).to_csv(
        OUT / "ligand_failures.csv", index=False)

log(f"\nDONE — docking_pipeline.py ({time.time()-t0:.0f}s)")
