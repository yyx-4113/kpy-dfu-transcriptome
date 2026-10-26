# -*- coding: utf-8 -*-
"""
build_repro_package.py — 组装 07_复现包/（GitHub 发布材料）

设计原则：
  - **脚本零修改**：复制 02_脚本/ 到 07_复现包/scripts/，同时对硬编码的
    绝对路径做**参数化改写**（`BASE = "D:/..."` → 读环境变量 KPYRN_PROJECT_ROOT
    或相对路径），保证复现者无需手工改 44 个文件。
  - **小体积入仓**：GitHub 仓库只放代码、结果 CSV/JSON、图件、稿件与附表；
    scRNA npz / 表达矩阵等大数据不入仓，由 MANIFEST + Zenodo 归档（DOI）。
  - **零重算**：本脚本只做复制与路径改写，不重新计算任何统计量。

输出：07_复现包/
"""
import os
import re
import shutil
import hashlib
import json
import datetime

BASE = _ROOT
PKG = f"{BASE}/07_复现包"
SCRIPTS = f"{PKG}/scripts"

# ── 需要排除的文件（体积大 / 不应入仓 / 机器相关） ────────────────
EXCLUDE_DIR_NAMES = {"__pycache__", "provenance", "figures", ".git", ".ipynb_checkpoints"}
EXCLUDE_FILE_SUFFIX = {".pyc", ".pyo"}
# 大文件阈值（字节）：单文件超过此值不进 GitHub 仓库（转 Zenodo 归档）
# 0.5 MB 阈值保证 GitHub 仓库整体保持在 ~30 MB 以内
SIZE_LIMIT = 512 * 1024

manifest_rows = []
zenodo_extra = []


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def copy_tree(src, dst, skip_dirs=(), skip_suffix=(), size_limit=None, base_note=""):
    """递归复制并记录清单。size_limit 非空时跳过超限文件。
    同时清理本目录中本轮不再产生的陈旧文件（保持幂等）。"""
    n = 0
    skipped = []
    written = set()
    for root, dirs, files in os.walk(src):
        rel = os.path.relpath(root, src)
        # 剪枝目录
        dirs[:] = [d for d in dirs
                   if d not in skip_dirs and d not in EXCLUDE_DIR_NAMES]
        for fn in files:
            if os.path.splitext(fn)[1] in skip_suffix or fn in EXCLUDE_FILE_SUFFIX:
                continue
            sp = os.path.join(root, fn)
            rp = os.path.join(dst, rel) if rel != "." else dst
            os.makedirs(rp, exist_ok=True)
            dp = os.path.join(rp, fn)
            size = os.path.getsize(sp)
            if size_limit and size > size_limit:
                skipped.append((os.path.relpath(sp, BASE).replace("\\", "/"), size))
                continue
            shutil.copy2(sp, dp)
            written.add(os.path.abspath(dp))
            manifest_rows.append({
                "path": os.path.relpath(dp, PKG).replace("\\", "/"),
                "bytes": size,
                "sha256": sha256(dp),
                "origin": os.path.relpath(sp, BASE).replace("\\", "/"),
            })
            n += 1
    # 清理陈旧文件（上一轮复制、本轮因超限/来源删除而不再出现的）
    if os.path.isdir(dst):
        for _root, _dirs, _files in os.walk(dst):
            for _fn in _files:
                _p = os.path.abspath(os.path.join(_root, _fn))
                if _p not in written:
                    os.remove(_p)
    return n, skipped


def copy_file(src, dst_rel, note=""):
    dst = os.path.join(PKG, dst_rel)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst)
    manifest_rows.append({
        "path": dst_rel.replace("\\", "/"),
        "bytes": os.path.getsize(dst),
        "sha256": sha256(dst),
        "origin": os.path.relpath(src, BASE).replace("\\", "/"),
    })
    return dst


print("=== 组装 07_复现包 ===")
os.makedirs(PKG, exist_ok=True)

# ── 1. 脚本（含路径参数化改写） ──────────────────────────────────
# 注意：不使用 shutil.rmtree 删除整目录（长中文路径下回收站机制会失败）；
#       采用「就地覆盖 + 清理陈旧文件」的方式同步，保持幂等。
os.makedirs(SCRIPTS, exist_ok=True)
_expected_scripts = set()

# 路径参数化：把硬编码的项目根路径改为读环境变量 KPYRN_PROJECT_ROOT。
# 需要覆盖 4 种声明形态（实测项目中 46 个脚本全部命中）：
#   1) BASE = "D:/..."                    → Python 字符串
#   2) BASE = Path("D:/...")              → Python pathlib
#   3) OUT  = "D:/.../03_结果/药物预测"   → 子目录常量
#   4) BASE <- "D:/..."                   → R
#   5) file.path("D:/...", ...)           → R 复合
_PY_STR = (
    'import os as _os\n'
    '_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(\n'
    '    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))\n'
)
_PY_BASE = (
    'BASE = Path(_os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(\n'
    '    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..")))\n'
)
_R_BASE = (
    'BASE <- Sys.getenv("KPYRN_PROJECT_ROOT",\n'
    '                   file.path(dirname(dirname(normalizePath(sys.frame(1)$ofile)))))'
)

n_script = 0
path_fixed = []
for root, dirs, files in os.walk(f"{BASE}/02_脚本"):
    dirs[:] = [d for d in dirs if d not in EXCLUDE_DIR_NAMES]
    for fn in files:
        if not fn.endswith((".py", ".R")):
            continue
        sp = os.path.join(root, fn)
        rel = os.path.relpath(sp, f"{BASE}/02_脚本")
        dp = os.path.join(SCRIPTS, rel)
        os.makedirs(os.path.dirname(dp), exist_ok=True)
        with open(sp, encoding="utf-8") as fh:
            txt = fh.read()
        before = txt
        is_r = fn.endswith(".R")

        # 形态 1/3：Python 字符串常量（含子目录）
        # 形如 OUT = "D:/.../丙酮酰化与糖尿病足/03_结果/药物预测"
        #   → OUT = _ROOT + "/03_结果/药物预测"
        # 形如 BASE = "D:/.../丙酮酰化与糖尿病足"（无后缀）
        #   → BASE = _ROOT
        def _repl_str(m):
            suf = m.group(4) or ""
            return (f'{m.group(1)}{m.group(2)}{m.group(3)}'
                    + (f'_ROOT + "{suf}"' if suf else '_ROOT'))

        txt = re.sub(
            r'(?m)^(\s*)([A-Z_][A-Z0-9_]*)(\s*=\s*)"D:/2026\.10/丙酮酰化与神经病理性疼痛/丙酮酰化与糖尿病足(/[^"]*)?"',
            _repl_str, txt)
        # 形态 2：Python pathlib
        #   BASE = Path("D:/.../根/01_数据")  →  BASE = Path(_ROOT, "01_数据")
        #   BASE = Path("D:/.../根")          →  BASE = Path(_ROOT)
        def _repl_path(m):
            suf = (m.group(4) or "").lstrip("/")
            return (f'{m.group(1)}{m.group(2)}{m.group(3)}'
                    + (f'Path(_ROOT, "{suf}")' if suf else 'Path(_ROOT)'))

        txt = re.sub(
            r'(?m)^(\s*)([A-Z_][A-Z0-9_]*)(\s*=\s*)Path\("D:/2026\.10/丙酮酰化与神经病理性疼痛/丙酮酰化与糖尿病足([^"]*)?"\)',
            _repl_path, txt)
        # 形态 4：R 字符串常量
        #   BASE <- "D:/.../根"           → BASE <- Sys.getenv("KPYRN_PROJECT_ROOT")
        #   OUT  <- "D:/.../根/03_结果/x" → OUT  <- file.path(Sys.getenv(...), "03_结果/x")
        def _repl_r(m):
            suf = m.group(4) or ""
            if suf:
                return (f'{m.group(1)}{m.group(2)} {m.group(3)} '
                        f'file.path(Sys.getenv("KPYRN_PROJECT_ROOT"), '
                        f'"{suf.lstrip("/")}")')
            return (f'{m.group(1)}{m.group(2)} {m.group(3)} '
                    f'Sys.getenv("KPYRN_PROJECT_ROOT")')

        txt = re.sub(r'(?m)^(\s*)([A-Za-z_.][A-Za-z0-9_.]*)\s*(<-|=)\s*"D:/2026\.10/丙酮酰化与神经病理性疼痛/丙酮酰化与糖尿病足(/[^"]*)?"',
                      _repl_r, txt)
        # 形态 5：R file.path("D:/...", ...)
        def _repl_r_fp(m):
            suf = (m.group(1) or "").lstrip("/")
            return f'file.path(Sys.getenv("KPYRN_PROJECT_ROOT"), "{suf}")'

        txt = re.sub(r'file\.path\("D:/2026\.10/丙酮酰化与神经病理性疼痛/丙酮酰化与糖尿病足(/[^"]*)?"',
                      _repl_r_fp, txt)

        # 注入根路径定义（Python 需 _ROOT；R 用 Sys.getenv 内联，无需注入）
        if txt != before:
            if not is_r:
                # Python：在首个 import 之后插入 _ROOT 定义
                if "_ROOT = _os.environ.get" not in txt:
                    m_imp = re.search(r'(?m)^(import .+$|from .+ import .+$)', txt)
                    if m_imp:
                        txt = (txt[:m_imp.end()] + "\n\n" + _PY_STR
                               + txt[m_imp.end():])
                    else:
                        txt = _PY_STR + "\n" + txt
                if "Path(_ROOT" in txt and "from pathlib import Path" not in txt:
                    txt = "from pathlib import Path\n" + txt
            path_fixed.append(rel.replace("\\", "/"))
        with open(dp, "w", encoding="utf-8") as fh:
            fh.write(txt)
        _expected_scripts.add(os.path.abspath(dp))
        manifest_rows.append({
            "path": f"scripts/{rel.replace(chr(92), '/')}",
            "bytes": os.path.getsize(dp),
            "sha256": sha256(dp),
            "origin": f"02_脚本/{rel.replace(chr(92), '/')}",
        })
        n_script += 1

# 清理上一轮遗留、本轮不再产生的脚本文件（保持幂等）
for _root, _dirs, _files in os.walk(SCRIPTS):
    for _fn in _files:
        _p = os.path.abspath(os.path.join(_root, _fn))
        if _p not in _expected_scripts:
            os.remove(_p)

print(f"  scripts/       : {n_script} 个脚本，{len(path_fixed)} 个完成路径参数化")

# ── 2. 产物：小体积结果 / 图件 / 稿件 / 附表 ─────────────────────
n_res, skip_res = copy_tree(f"{BASE}/03_结果", f"{PKG}/results", size_limit=SIZE_LIMIT)
print(f"  results/       : {n_res} 个文件"
      + (f"，{len(skip_res)} 个大文件转入 Zenodo" if skip_res else ""))

n_fig, _ = copy_tree(f"{BASE}/04_图表", f"{PKG}/figures")
print(f"  figures/       : {n_fig} 个文件")

n_genes, _ = copy_tree(f"{BASE}/01_数据/基因集", f"{PKG}/data/gene_sets")
print(f"  data/gene_sets/: {n_genes} 个文件")

n_meta, skip_meta = copy_tree(f"{BASE}/01_数据/元数据", f"{PKG}/data/metadata",
                              size_limit=SIZE_LIMIT)
print(f"  data/metadata/ : {n_meta} 个文件"
      + (f"，{len(skip_meta)} 个大文件转入 Zenodo" if skip_meta else ""))
for _p, _s in skip_meta:
    zenodo_extra.append((_p, _s))

# 基因集里的大文件（msigdb gmt 6.7MB）不入 git
_msig = f"{PKG}/data/gene_sets/免疫/msigdb_c7_immunesigdb_v2023.2.Hs.symbols.gmt"
if os.path.exists(_msig) and os.path.getsize(_msig) > SIZE_LIMIT:
    os.remove(_msig)
    for _r in list(manifest_rows):
        if _r["path"].endswith("msigdb_c7_immunesigdb_v2023.2.Hs.symbols.gmt"):
            manifest_rows.remove(_r)
    zenodo_extra.append(("data/gene_sets/免疫/"
                         "msigdb_c7_immunesigdb_v2023.2.Hs.symbols.gmt", 7067060))

# 稿件与审计
for src, rel in [
    (f"{BASE}/05_稿件/manuscript_en.md", "manuscript/manuscript_en.md"),
    (f"{BASE}/05_稿件/manuscript_cn.md", "manuscript/manuscript_cn.md"),
    (f"{BASE}/05_稿件/audit_facts.json", "manuscript/audit_facts.json"),
    (f"{BASE}/05_稿件/audit_facts.csv", "manuscript/audit_facts.csv"),
    (f"{BASE}/05_稿件/audit_facts_meta.json", "manuscript/audit_facts_meta.json"),
    (f"{BASE}/05_稿件/audit_verification.json", "manuscript/audit_verification.json"),
    (f"{BASE}/05_稿件/审计追踪表.csv", "manuscript/audit_traceability.csv"),
    (f"{BASE}/05_稿件/审计报告.md", "manuscript/audit_report.md"),
    (f"{BASE}/06_补充材料/S4_RT-qPCR_验证方案.md",
     "manuscript/S4_RT-qPCR_validation_protocol.md"),
    # 稿件 §2.12/§4.7 引用 "Supplementary Material S4b"（最小可行集设计留档），
    # 复现包必须自洽包含该文件，否则包内引用悬空。
    (f"{BASE}/06_补充材料/S4b_湿实验最小集SOP.md",
     "manuscript/S4b_minimal_wetlab_set_SOP.md"),
]:
    if os.path.exists(src):
        copy_file(src, rel)
print(f"  manuscript/    : 已复制稿件、审计台账、S4 方案与 S4b 设计留档")

# 附表是稿件的正式组成部分（S1–S10），必须完整入 git，
# 不受通用 SIZE_LIMIT 约束（单文件最大 22 MB，远低于 GitHub 100 MB 限制）。
n_supp, _ = copy_tree(f"{BASE}/05_稿件/附表", f"{PKG}/supplementary_tables")
print(f"  supplementary_tables/: {n_supp} 个文件（附表为稿件组成，全量保留）")

# 数据溯源清单（原始数据 URL + SHA256）
copy_file(f"{BASE}/01_数据/MANIFEST.csv", "data/MANIFEST.csv")
if os.path.exists(f"{BASE}/01_数据/MR/MANIFEST_mr.csv"):
    copy_file(f"{BASE}/01_数据/MR/MANIFEST_mr.csv", "data/MANIFEST_mr.csv")

# ── 3. 写包内 MANIFEST ──────────────────────────────────────────
import csv as _csv
with open(f"{PKG}/MANIFEST.csv", "w", encoding="utf-8", newline="") as fh:
    w = _csv.DictWriter(fh, fieldnames=["path", "bytes", "sha256", "origin"])
    w.writeheader()
    for r in sorted(manifest_rows, key=lambda x: x["path"]):
        w.writerow(r)

# Zenodo 归档清单（未入仓的大文件）
zenodo_rows = [{"path": p, "bytes": s} for p, s in skip_res]
if os.path.exists(f"{BASE}/01_数据/处理后数据"):
    for root, dirs, files in os.walk(f"{BASE}/01_数据/处理后数据"):
        rel = os.path.relpath(root, BASE).replace("\\", "/")
        for fn in files:
            sp = os.path.join(root, fn)
            sz = os.path.getsize(sp)
            if sz > SIZE_LIMIT:
                zenodo_rows.append({"path": f"{rel}/{fn}", "bytes": sz})
for _p, _s in zenodo_extra:
    zenodo_rows.append({"path": _p, "bytes": _s})
with open(f"{PKG}/ZENODO_ARCHIVE_LIST.csv", "w", encoding="utf-8", newline="") as fh:
    w = _csv.DictWriter(fh, fieldnames=["path", "bytes"])
    w.writeheader()
    for r in sorted(zenodo_rows, key=lambda x: x["path"]):
        w.writerow(r)

print(f"  MANIFEST.csv   : {len(manifest_rows)} 条")
print(f"  ZENODO_ARCHIVE_LIST.csv : {len(zenodo_rows)} 条（转 Zenodo 的大文件）")
print("DONE — build_repro_package.py")
