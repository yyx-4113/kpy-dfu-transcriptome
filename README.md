# Transcriptomic Footprint of the Lysine-Pyruvylation Regulatory Network in Diabetic Foot Ulcer Tissue

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23677954.svg)](https://doi.org/10.5281/zenodo.23677954)

**Reproducibility package** for the manuscript *"The lysine-pyruvylation regulatory network is transcriptionally dysregulated in diabetic foot ulcer tissue"*.

This repository contains the full analysis code, derived result tables, figures, manuscript text, and a machine-readable **audit ledger** that traces every numeric claim in the manuscript back to the file it was computed from.

> **Single-author notice.** This work is authored by a single investigator. The ORCID, affiliation, and citation metadata in `CITATION.cff` are authoritative for authorship and version-of-record purposes.

---

## 1. What this study did (and did not) claim

| | |
|---|---|
| **Question** | Does the lysine-pyruvylation (Kpy) regulatory network show a reproducible transcriptional signature in diabetic foot ulcer (DFU) tissue? |
| **Design** | Public-data re-analysis of **10 GEO cohorts** (5 bulk, 2 single-cell, 1 spatially resolved, 1 longitudinal infection, 1 independent outcome cohort), plus cis-eQTL Mendelian randomisation and molecular docking |
| **Key positive result** | 26 of 128 scored KpyRN genes meta-analysis significant (FDR < 0.05) with **26/26 cross-cohort direction concordance**; glycolytic arm up, IFN-I read-out down *in part* |
| **Key negative results** | No significant module-level signal at donor/subject level in GSE223964 (single-cell) or GSE166120 (GeoMx DSP); all whole-genome DEG contrasts in the infection cohort negative; no docked candidate outperformed its crystallographic reference cofactor |
| **Critical scope limit** | **Kpy modification itself was never measured.** No Kpy antibody and no mass spectrometry were available. Every finding concerns the *transcriptional footprint* of the Kpy regulatory network — upstream pyruvate supply and downstream IFN-I read-out — and is reported as such throughout. |
| **Not performed** | The pre-specified tissue RT-qPCR protocol was **designed but not executed**. It is published as a protocol with sample size derived from the observed pooled effect sizes. |

---

## 2. Repository map

```
.
├── README.md                     ← this file
├── CITATION.cff                  ← citation metadata (version of record)
├── LICENSE                       ← MIT
├── MANIFEST.csv                  ← SHA-256 + byte size of every file in this repo
├── ZENODO_ARCHIVE_LIST.csv       ← large files excluded from git; archived on Zenodo
├── requirements.txt              ← Python dependencies (pinned)
├── environment.yml               ← R environment (conda)
├── run_all.sh                    ← full pipeline, staged, with dependency ordering
├── scripts/                      ← 46 analysis scripts (all paths parameterised)
├── data/
│   ├── MANIFEST.csv              ← source-data provenance: URL + SHA-256 + download date
│   ├── MANIFEST_mr.csv           ← provenance for MR source files
│   ├── gene_sets/                ← KpyRN five-module gene set + immune signatures
│   └── metadata/                 ← GEO sample-level metadata (26 series)
├── results/                      ← all derived result tables (CSV/JSON/MD)
├── figures/                      ← manuscript main figures (PNG)
├── manuscript/
│   ├── manuscript_en.md          ← English manuscript (IMRaD)
│   ├── manuscript_cn.md          ← Chinese manuscript (equivalent structure)
│   ├── audit_facts.json          ← ★ the single source of truth for manuscript numbers
│   ├── audit_facts.csv           ← same, flat table
│   ├── audit_traceability.csv    ← fact id → manuscript location, bidirectional
│   ├── audit_report.md           ← coverage statistics and residual gaps
│   └── audit_verification.json   ← machine-readable audit result
└── supplementary_tables/         ← 28 supplementary tables S1–S10 (+ index & notes)
```

---

## 3. Quick start

### 3.1 Requirements

| Component | Version used | Notes |
|---|---|---|
| Python | 3.13.14 | see `requirements.txt` |
| R | 4.4.3 | see `environment.yml` (Bioconductor 3.20) |
| AutoDock Vina | 1.2.7 | `01_数据/工具/vina.exe` in the source project; any 1.2.x build works |
| Open Babel | 3.1.0 | must be on `PATH` (`obabel`) for receptor/ligand PDBQT preparation |
| Git Bash / bash + awk + zcat | — | required by the independent MR re-computation path |

Key pinned Python packages: pandas 2.3.3, numpy 2.4.6, scipy 1.18.1, statsmodels 0.15.0, scikit-learn 1.9.0, matplotlib 3.11.1, scanpy 1.12.4, anndata 0.13.3, primer3-py 2.4.0, rdkit 2026.3.6, meeko 0.8.0, openpyxl 3.1.5.

Key R packages: limma 3.62.2, metafor 5.2.1, ggplot2 4.0.3, WGCNA 1.74, edgeR 4.4.2, jsonlite 2.0.0, EPIC 1.1.7 (GitHub `GfellerLab/EPIC`, not on CRAN).

### 3.2 Path configuration

All 46 scripts in `scripts/` read the project root from the `KPYRN_PROJECT_ROOT` environment variable. If it is unset, they fall back to the parent directory of `scripts/`.

```bash
export KPYRN_PROJECT_ROOT=/absolute/path/to/this/repo   # POSIX
set KPYRN_PROJECT_ROOT=D:\path\to\this\repo              # Windows cmd
```

The upstream project working copy (which contains `01_数据/原始数据/` and the large processed matrices) must be visible at that root for the full pipeline to run. This repository ships the code and the derived results, not the multi-gigabyte raw downloads.

### 3.3 Running

```bash
bash run_all.sh              # full pipeline, staged
bash run_all.sh --stage 4    # single stage
```

`run_all.sh` enforces the dependency order; stages cannot be run out of sequence without breaking downstream inputs. The strict serial constraints are:

```
meta_analysis.R  →  wgcna.R  →  ml_biomarker.py
export_estimate_sig.R  →  immune_estimate_correlate.py
immune_deconv.R  →  immune_estimate_correlate.py  →  immune_figures.py
scrna_load_qc.py  →  scrna_annotate_score.py  →  scrna_composition_stats.py
docking_pipeline.py  →  recompute_rmsd.py  →  docking_summary.py / docking_figures.py
build_audit_facts.py  →  verify_manuscript_audit.py
```

### 3.4 Verify the shipped results before re-running anything

```bash
python scripts/verify_manuscript_audit.py
```

Expected output:

```
事实覆盖（英稿）：600/617 = 97.2%
事实覆盖（中稿）：588/617 = 95.3%
未溯源 token：0 处 / 0 个不同取值
```

Every number in the manuscript body is matched against `manuscript/audit_facts.json` in both directions: (a) coverage — does each fact appear in the text? (b) traceability — does each number in the text resolve to a fact? The expected values are **derived from the fact ledger at runtime**, never hard-coded, so the audit cannot silently protect a stale number.

To confirm file integrity against the ledger:

```bash
python - <<'PY'
import csv, hashlib, pathlib
bad = 0
for r in csv.DictReader(open('MANIFEST.csv', encoding='utf-8')):
    p = pathlib.Path(r['path'])
    h = hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None
    if h != r['sha256']:
        print('MISMATCH', r['path']); bad += 1
print(f'{bad} mismatches')
PY
```

---

## 4. Pipeline stages and expected outputs

| Stage | Scripts | Key outputs (under `results/`) |
|---|---|---|
| 1. Data acquisition | `download_data.py`, `geo_enumerate.py`, `geo_sample_table.py`, `rcsb_receptor_search.py`, `chembl_candidates.py` | `data/MANIFEST.csv`, `data/metadata/` |
| 2. Gene set | `build_kpyrn.py` | `data/gene_sets/KpyRN_modules.csv` (140 entries, 128 scored) |
| 3. Preprocessing | `preprocess_rnaseq.py`, `preprocess_gse230426.py`, `scrna_load_qc.py` | processed matrices, `results/QC/` |
| 4. Core statistics | `deg_limma.R`, `meta_analysis.R`, `wgcna.R`, `ml_biomarker.py`, `deg_gse230426.R` | `results/DEG/`, `results/meta/`, `results/WGCNA/`, `results/ML/` |
| 5. Immune infiltration | `export_estimate_sig.R`, `immune_deconv.R`, `immune_estimate_correlate.py`, `immune_figures.py` | `results/免疫浸润/` |
| 6. Single-cell & spatial | `scrna_annotate_score.py`, `scrna_composition_stats.py`, `scrna_gse223964.py`, `xcheck_gse223964.py`, `geommx_gse166120.py` | `results/单细胞/`, `results/单细胞GSE223964/`, `results/空间GSE166120/` |
| 7. Mendelian randomisation | `mr_eqtl_glycemic.py`, `mr_independent_check.py` | `results/MR/` |
| 8. Molecular docking | `docking_pipeline.py`, `recompute_rmsd.py`, `docking_summary.py`, `docking_figures.py` | `results/药物预测/docking/` |
| 9. RT-qPCR design | `design_rtqpcr_primers.py`, `verify_primers.py`, `verify_primers_blast_v3.py`, `fix_refseq_provenance.py`, `compute_qpcr_sample_size.py` | `results/实验验证设计/` |
| 10. Figures | `figures_main.py`, `make_figures_ml.py`, `figures_gse230426.py`, `figures_geommx_gse166120.py`, `docking_figures.py`, `immune_figures.py` | `figures/` |
| 11. Audit & tables | `build_audit_facts.py`, `build_supplementary_tables.py`, `verify_manuscript_audit.py` | `manuscript/audit_*`, `supplementary_tables/` |

---

## 5. Methodological decisions worth knowing before you re-run this

These are the choices that most affect reproducibility. Each is enforced in code, not just described.

**5.1 Nested designs are analysed at the level of the independent unit.** Single-cell ROIs are nested in donors and GeoMx ROIs are nested in subjects. Donor-level and subject-level inference is primary throughout; ROI/cell-level results are reported only as sensitivity analyses, and a nominal result that disappears when re-analysed at the correct level is classified as a false positive rather than a finding. See `results/空间GSE166120/summary.json` (`note` field).

**5.2 The audit ledger is generated from products, never re-derived.** `build_audit_facts.py` reads numbers out of result files; it does not recompute statistics. Values such as the docking candidate count (29), the infection-cohort sample sizes (20/18/15), the MCP-counter signature size (111), and the GeoMx negative-probe count (88) are each derived from their authoritative source file at build time and carry the extraction expression in the `expr` field.

**5.3 The positive-control docking RMSD gate.** Eight receptor–ligand systems were validated by re-docking the crystallographic ligand; heavy-atom RMSD was 0.00 Å in every system (acceptance threshold < 2 Å). Candidate affinities are always reported relative to the co-crystal reference in the same pocket (`delta_vs_pos`), never as absolute "better than" claims.

**5.4 BLAST specificity was established against `core_nt`, not `refseq_rna`.** The `refseq_rna` database does not contain a complete set of human transcripts — the LDHA forward primer returns zero human hits there while correctly hitting `NG_008185` in `core_nt` — which produces false negatives. `verify_primers_blast_v3.py` therefore queries `core_nt` with XML output (the Text format truncates hit descriptions, which silently hid the gene symbol) and classifies hits into true off-target (human RefSeq transcript), paralog, pseudogene, and genome-level tiers.

**5.5 `primer_verification.csv`'s `overall_pass` column is not a pass/fail verdict.** It inherits the `V1_specific` flag from the v1 BLAST run, which used the truncating Text format and therefore returns `False` for all 12 pairs. Specificity is adjudicated by `verify_primers_blast_v3.py`; the structural and GC/Tm criteria are recorded separately in the same file.

**5.6 The infection-cohort `n_pairs_inferred` is an inference, not metadata.** The source series lacks explicit donor identifiers; pairing was inferred from sample naming and is flagged as such in `data/metadata/` and in the manuscript limitations.

---

## 6. Data provenance and licensing

All primary data are public GEO / PDB / ChEMBL resources. `data/MANIFEST.csv` records, for every downloaded file: the GEO series accession, the file name, the full FTP URL, the byte size, the SHA-256 checksum, and the download timestamp. MR source files (eQTLGen, GWAS Catalog) are tracked separately in `data/MANIFEST_mr.csv`.

Large files (single-cell count matrices, bulk expression matrices) are not stored in git. Their paths, sizes, and checksums are listed in `ZENODO_ARCHIVE_LIST.csv` and are deposited with the Zenodo record referenced in the manuscript's data availability statement. The small derived tables, figures, and metadata that carry the actual findings **are** in this repository.

Third-party signature resources retain their original licenses; `data/gene_sets/免疫/MCPcounter_License.txt` is included verbatim.

## 7. Citation

See `CITATION.cff`. If you use these materials, please cite the manuscript and the Zenodo deposit for the version of the data you used, and cite the original GEO series for any cohort you reuse.

## 8. Limitations reproduced from the manuscript

(1) Kpy was not measured — all findings concern the transcriptional footprint of the network. (2) Two bulk cohorts (GSE80178, GSE68183/86) could not be processed at sample level and were not used for computation. (3) The discovery cohorts are small (n = 21 each) and the WGCNA soft-threshold fit is unstable at that n; module results were used only as a filter. (4) The infection cohort lacks explicit donor identifiers and healing-outcome metadata. (5) GeoMx coverage is partial (42.2 % of the scoring set) with no upstream normalisation metrics. (6) Single-cell annotation is marker-score based with 12.7 % ambiguous cells; sub-lineage claims are exploratory. (7) The RT-qPCR validation protocol was pre-specified but not executed.

## 9. License

Code and derived material: MIT (see `LICENSE`). Primary data remain under their original terms; the KpyRN gene set is released under MIT with per-gene provenance recorded in `data/gene_sets/KpyRN_modules.csv`.
