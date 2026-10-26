# Release notes — v1.0.0

**Repository:** `yyx-4113/kpy-dfu-transcriptome`
**Title:** Transcriptomic footprint of the lysine-pyruvylation regulatory network in diabetic foot ulcer tissue
**Date:** 2026-10-25 (package assembled) · tag created for *BMC Medical Genomics* submission
**License:** MIT
**Corresponding Zenodo DOI:** _to be assigned on first Zenodo archival; paste the version DOI back into the manuscript Data Availability statement and the SNAPP submission field._

## What this release contains

A fully source-traceable reproduction package for the manuscript
*"The lysine-pyruvylation regulatory network is transcriptionally dysregulated in diabetic foot ulcer tissue"*
(submitted to *BMC Medical Genomics* as a Research Article).

- **`scripts/`** — 60+ analysis scripts (all paths parameterised via `KPYRN_PROJECT_ROOT`), covering data acquisition, gene-set construction, preprocessing, differential-expression meta-analysis, WGCNA, consensus machine learning, immune deconvolution, single-cell & spatial transcriptomics, cis-eQTL Mendelian randomisation, molecular docking, RT-qPCR primer design (protocol, unexecuted), and figure/audit generation.
- **`data/`** — gene sets (KpyRN five-module set, immune signatures), GEO sample-level metadata (14 series), and provenance MANIFESTs (SHA-256 + download URL + date) for source and MR inputs.
- **`results/`** — all derived result tables (CSV/JSON/MD) for every analysis stage.
- **`figures/`** — manuscript main figures (PNG, 300 dpi).
- **`manuscript/`** — English and Chinese manuscript text plus the machine-readable audit ledger (`audit_facts.json`) and its verification report.
- **`supplementary_tables/`** — supplementary tables S1–S11 (+ index & notes).
- **`CITATION.cff`**, **`LICENSE`** (MIT), **`README.md`**, **`MANIFEST.csv`** (SHA-256 of every shipped file), **`ZENODO_ARCHIVE_LIST.csv`** (large matrices archived on Zenodo).

## Scope limits (reproduced from the manuscript)

1. **Kpy modification itself was never measured** — no Kpy antibody, no mass spectrometry. Every result concerns the *transcriptional footprint* of the Kpy regulatory network.
2. Two bulk series (**GSE80178**, **GSE68183/86**) are inventoried in S2 as sample metadata but could not be processed at sample level and were **not used for computation**; the analysed set is 10 cohorts.
3. The pre-specified tissue RT-qPCR protocol was designed but **not executed** (published as a protocol with sample size derived from observed effect sizes).
4. Discovery cohorts are small (n = 21 each); WGCNA soft-threshold fit is unstable at that n, so module results were used only as a filter.
5. Single-cell annotation is marker-score based (12.7 % ambiguous cells); sub-lineage claims are exploratory.
6. GeoMx spatial coverage is partial (42.2 % of the scoring set); infection cohort lacks explicit donor identifiers (pairing inferred from sample naming, flagged).

## How to verify

```bash
python scripts/verify_manuscript_audit.py     # expect: 未溯源 token = 0
python - <<'PY'                                # MANIFEST integrity
import csv, hashlib, pathlib
bad=0
for r in csv.DictReader(open('MANIFEST.csv', encoding='utf-8')):
    p=pathlib.Path(r['path'])
    h=hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None
    if h!=r['sha256']:
        print('MISMATCH', r['path']); bad+=1
print(f'{bad} mismatches')
PY
```

## Citation

Cite the manuscript and the Zenodo deposit for the version of the data you used; cite the original GEO series for any cohort you reuse. See `CITATION.cff`.
