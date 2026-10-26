#!/usr/bin/env bash
# run_all.sh — Kpy × DFU 全流程复现脚本
#
# 用法：
#   bash run_all.sh                # 全部阶段
#   bash run_all.sh --stage 4      # 只跑第 4 阶段
#   bash run_all.sh --from 4 --to 7
#
# 环境要求：
#   export KPYRN_PROJECT_ROOT=/abs/path/to/project   # 指向含 01_数据/ 的工作副本
#   python (3.13) / Rscript (4.4.3) / obabel / vina.exe
#
# 依赖顺序为硬约束，不可跳序：
#   meta_analysis.R → wgcna.R → ml_biomarker.py
#   export_estimate_sig.R → immune_estimate_correlate.py
#   immune_deconv.R → immune_estimate_correlate.py → immune_figures.py
#   scrna_load_qc.py → scrna_annotate_score.py → scrna_composition_stats.py
#   docking_pipeline.py → recompute_rmsd.py → docking_summary.py / docking_figures.py
#   build_audit_facts.py → verify_manuscript_audit.py

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="$SCRIPT_DIR/scripts"
export KPYRN_PROJECT_ROOT="${KPYRN_PROJECT_ROOT:-$SCRIPT_DIR}"

PY="${PYTHON:-python}"
RS="${RSCRIPT:-Rscript}"

FROM=1
TO=11
if [ "${1:-}" = "--stage" ]; then FROM="$2"; TO="$2"; fi
if [ "${1:-}" = "--from" ]; then FROM="$2"; TO="${3:-11}"; fi

stage() {
  local n="$1"; shift
  if [ "$n" -lt "$FROM" ] || [ "$n" -gt "$TO" ]; then
    echo "--- 跳过阶段 $n ---"; return
  fi
  echo ""
  echo "==================== 阶段 $n ===================="
  for cmd in "$@"; do
    echo ">>> $cmd"
    eval "$cmd"
  done
}

# ── 阶段 1：数据获取（联网） ──────────────────────────────────
stage 1 \
  "$PY $S/geo_enumerate.py" \
  "$PY $S/geo_sample_table.py" \
  "$PY $S/download_data.py" \
  "$PY $S/rcsb_receptor_search.py" \
  "$PY $S/chembl_candidates.py"

# ── 阶段 2：基因集构建 ───────────────────────────────────────
stage 2 "$PY $S/build_kpyrn.py"

# ── 阶段 3：预处理 ───────────────────────────────────────────
stage 3 \
  "$PY $S/preprocess_rnaseq.py" \
  "$PY $S/preprocess_gse230426.py" \
  "$PY $S/scrna_load_qc.py"

# ── 阶段 4：核心统计主链（严格串行） ─────────────────────────
stage 4 \
  "$RS $S/deg_limma.R" \
  "$RS $S/meta_analysis.R" \
  "$RS $S/redraw_meta_volcano.R" \
  "$RS $S/wgcna.R" \
  "$PY $S/ml_biomarker.py" \
  "$PY $S/make_figures_ml.py" \
  "$RS $S/deg_gse230426.R"

# ── 阶段 5：免疫反卷积 ───────────────────────────────────────
stage 5 \
  "$RS $S/export_estimate_sig.R" \
  "$RS $S/immune_deconv.R" \
  "$PY $S/immune_estimate_correlate.py" \
  "$PY $S/immune_figures.py"

# ── 阶段 6：单细胞与空间（内部串行，彼此独立） ───────────────
stage 6 \
  "$PY $S/scrna_annotate_score.py" \
  "$PY $S/scrna_composition_stats.py" \
  "$PY $S/scrna_gse223964.py" \
  "$PY $S/xcheck_gse223964.py" \
  "$PY $S/geommx_gse166120.py"

# ── 阶段 7：孟德尔随机化 ─────────────────────────────────────
stage 7 \
  "$PY $S/mr_eqtl_glycemic.py" \
  "$PY $S/mr_independent_check.py"

# ── 阶段 8：分子对接 ─────────────────────────────────────────
stage 8 \
  "$PY $S/docking_pipeline.py" \
  "$PY $S/recompute_rmsd.py" \
  "$PY $S/docking_summary.py" \
  "$PY $S/docking_figures.py"

# ── 阶段 9：RT-qPCR 方案设计 ─────────────────────────────────
stage 9 \
  "$PY $S/design_rtqpcr_primers.py" \
  "$PY $S/verify_primers.py" \
  "$PY $S/verify_primers_blast_v3.py" \
  "$PY $S/fix_refseq_provenance.py" \
  "$PY $S/compute_qpcr_sample_size.py"

# ── 阶段 10：图件汇总 ────────────────────────────────────────
stage 10 \
  "$PY $S/figures_main.py" \
  "$PY $S/figures_gse230426.py" \
  "$PY $S/figures_geommx_gse166120.py"

# ── 阶段 11：审计与附表收尾（严格最后） ───────────────────────
stage 11 \
  "$PY $S/build_audit_facts.py" \
  "$PY $S/build_supplementary_tables.py" \
  "$PY $S/verify_manuscript_audit.py"

echo ""
echo "==================== 全部完成 ===================="
echo "审计结果：05_稿件/审计报告.md"
