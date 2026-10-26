# -*- coding: utf-8 -*-
"""
apply_reference_fixes.py — 把《参考文献核验报告》的结论落到稿件

原则：
  1. 每条 (old, new) 必须**恰好命中一次**，否则报错退出（防止静默漏改）。
  2. 改写值全部来自 Crossref / PubMed / GEO 的实证查询（见 03_结果/参考文献核验.md）。
  3. 参考文献表整体替换；正文引用标记同步重编号（旧 34 → 新 33）。
"""

import re
import sys

MS = "05_稿件/manuscript_en.md"
text = open(MS, encoding="utf-8").read()
orig_len = len(text)

# ─────────────────────────────────────────────────────────────
# 一、参考文献表整体替换（38 条 → 35 条，删除 3 条疑似不存在/无落点，
#      新增 Guyon 2002；全部补 DOI/卷期页；9 条元数据错误按实证改写）
# ─────────────────────────────────────────────────────────────
NEW_REFS = """1. Chinese Burn Association, the Yangtze River Delta Integrated Diabetic Foot Alliance, and the Editorial Committee of the Chinese Journal of Burns and Wound Repair. Practical guidelines for the prevention and management of diabetic foot disease in China. Burns Trauma. 2025;13:tkaf064. doi:10.1093/burnst/tkaf064
2. Zuo Y, Wang Q, Tian W, et al. Pyruvate is a natural suppressor of interferon signaling by inducing STAT1 protein pyruvylation. Cell. 2026;189(7):1975–1989.e19. doi:10.1016/j.cell.2026.01.023
3. Song X, Peng P, Zheng H, et al. Lysine pyruvylation couples glycolytic flux to epigenetic regulation. Nat Metab. 2026;8(7):1495–1507. doi:10.1038/s42255-026-01556-2
4. Park JS, Lee SJ. Pyruvylation expands the redox repertoire of lysine acylation. Trends Biochem Sci. 2026. doi:10.1016/j.tibs.2026.08.008
5. GEO GSE134431. Gene Expression Omnibus, NCBI. Accessed 2026-10-25.
6. Jiang N, Xu C, Xu Y, et al. Comprehensive transcriptomic analysis of immune-related genes in diabetic foot ulcers: new insights into mechanisms and therapeutic targets. Int Immunopharmacol. 2024;139:112638. doi:10.1016/j.intimp.2024.112638
7. GEO GSE143735. Gene Expression Omnibus, NCBI. Accessed 2026-10-25.
8. GEO GSE80178. Gene Expression Omnibus, NCBI. Accessed 2026-10-25.
9. GEO GSE68183. Gene Expression Omnibus, NCBI. Accessed 2026-10-25.
10. GEO GSE68186. Gene Expression Omnibus, NCBI. Accessed 2026-10-25.
11. GEO GSE165816. Gene Expression Omnibus, NCBI. Accessed 2026-10-25.
12. Theocharidis G, Thomas BE, Sarkar D, Mumme HL, et al. Single cell transcriptomic landscape of diabetic foot ulcers. Nat Commun. 2022;13(1):181. doi:10.1038/s41467-021-27801-8
13. GEO GSE166120. Gene Expression Omnibus, NCBI. Accessed 2026-10-25.
14. Agidigbi TS, Kwon HK, Knight JR, Zhao D, et al. Transcriptomic identification of genes expressed in invasive S. aureus diabetic foot ulcer infection. Front Cell Infect Microbiol. 2023;13:1198115. doi:10.3389/fcimb.2023.1198115
15. GEO GSE223964. Gene Expression Omnibus, NCBI. Accessed 2026-10-25.
16. GEO GSE230426. Gene Expression Omnibus, NCBI. Accessed 2026-10-25.
17. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B. 1995;57(1):289–300. doi:10.1111/j.2517-6161.1995.tb02031.x
18. Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43(7):e47. doi:10.1093/nar/gkv007
19. Viechtbauer W. Conducting meta-analyses in R with the metafor package. J Stat Softw. 2010;36(3):1–48. doi:10.18637/jss.v036.i03
20. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics. 2008;9:559. doi:10.1186/1471-2105-9-559
21. Breiman L. Random forests. Mach Learn. 2001;45(1):5–32. doi:10.1023/A:1010933404324
22. Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent. J Stat Softw. 2010;33(1):1–22. doi:10.18637/jss.v033.i01
23. Chen T, Guestrin C. XGBoost: a scalable tree boosting system. In: Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining. 2016:785–794. doi:10.1145/2939672.2939785
24. Guyon I, Weston J, Barnhill S, Vapnik V. Gene selection for cancer classification using support vector machines. Mach Learn. 2002;46(1–3):389–422. doi:10.1023/A:1012487302797
25. R Core Team. R: A Language and Environment for Statistical Computing. Version 4.4.3. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/
26. Racle J, de Jonge K, Baumgaertner P, Speiser DE, Gfeller D. Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. eLife. 2017;6:e26476. doi:10.7554/eLife.26476
27. Becht E, Giraldo NA, Lacroix L, et al. Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biol. 2016;17(1):218. doi:10.1186/s13059-016-1070-5
28. Hunter MV, Moncada R, Weiss JM, Yanai I. Spatially resolved transcriptomics reveals the architecture of the tumor-microenvironment interface. Nat Commun. 2021;12(1):6278. doi:10.1038/s41467-021-26614-z
29. Davey Smith G, Hemani G. Mendelian randomization: genetic anchors for causal inference in epidemiological studies. Hum Mol Genet. 2014;23(R1):R89–R98. doi:10.1093/hmg/ddu328
30. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, energy minimization, and multithreading. J Comput Chem. 2010;31(2):455–461. doi:10.1002/jcc.21334
31. Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0: new docking methods, expanded force field, and Python bindings. J Chem Inf Model. 2021;61(8):3891–3898. doi:10.1021/acs.jcim.1c00203
32. Open Babel Development Team. Open Babel: The Open Source Chemistry Toolbox. Version 3.1.0. https://openbabel.org/
33. Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen SG. Primer3—new capabilities and interfaces. Nucleic Acids Res. 2012;40(15):e115. doi:10.1093/nar/gks596
34. Bustin SA, Benes V, Garson JA, et al. The MIQE guidelines: minimum information for publication of quantitative real-time PCR experiments. Clin Chem. 2009;55(4):611–622. doi:10.1373/clinchem.2008.112797
35. Vandesompele J, De Preter K, Pattyn F, et al. Accurate normalization of real-time quantitative RT-PCR data by geometric averaging of multiple internal control genes. Genome Biol. 2002;3(7):research0034. doi:10.1186/gb-2002-3-7-research0034"""

lines = text.split("\n")
i = next(k for k, l in enumerate(lines) if re.match(r"^1\. Chinese Medical", l))
j = next(k for k, l in enumerate(lines) if re.match(r"^38\. Therneau", l))
assert j > i, "参考文献表定位失败"
lines[i:j + 1] = NEW_REFS.split("\n")
text = "\n".join(lines)
print(f"[1] 参考文献表已替换：38 条 → 35 条（行 {i+1}–{j+1}）")

# ─────────────────────────────────────────────────────────────
# 二、正文修正（每条须恰好命中一次）
# ─────────────────────────────────────────────────────────────
EDITS = [
    # ── (A) 引言：四个数字口径按 PMC 全文实证改写 ──
    ("In China, the reported prevalence of DFU is 8.1%, the lifetime risk of foot ulceration reaches 34%, "
     "the overall amputation rate is 19.03% and the annual case fatality is 14.4% [1]. Infection is the "
     "leading cause of prolonged non-healing and amputation, and current management — debridement combined "
     "with systemic therapy — yields a mean healing time of approximately 113 days with a 31.6% recurrence "
     "rate at one year [1].",
     "In China, the reported incidence of diabetic foot is 8.1%, the lifetime risk of developing foot "
     "lesions reaches 34%, and among patients with diabetic foot the overall amputation rate is 19.03% "
     "(2.14% major and 16.88% minor) and the annual mortality rate is 14.4% [1]. Infection is the leading "
     "cause of prolonged non-healing and amputation; healing rates at 12 weeks of treatment range from 24% "
     "to 82%, and the recurrence rate within one year after healing is 31.6% in individuals older than 50 [1]."),

    # ── (B) §2.1 数据集清单：修引用编号、修"Ten datasets"口径、补未计算使用的说明 ──
    ("Ten GEO datasets were used (Table 1). Five bulk cohorts (GSE134431 [5], GSE199939 [6], GSE143735 [7], "
     "GSE80178 [8], GSE68183/GSE68186 [9,10]), two single-cell cohorts (GSE165816 [11], GSE223964 [12]), "
     "one spatially resolved cohort (GSE166120 [13]), one longitudinal infection cohort (GSE230426 [14]) "
     "and one additional bulk cohort used for cross-validation.",
     "Ten GEO series were retrieved (Supplementary Table S2). Five bulk cohorts (GSE134431 [5], "
     "GSE199939 [6], GSE143735 [7], GSE80178 [8], GSE68183/GSE68186 [9,10]), two single-cell cohorts "
     "(GSE165816 [11], GSE223964 [15]), one spatially resolved cohort (GSE166120 [13]) and one "
     "longitudinal infection cohort (GSE230426 [14,16]); of these, GSE80178 and GSE68183/GSE68186 could "
     "not be processed at sample level and were not used computationally (Section 4.7, item 2)."),

    # ── (C) §2.1：Theocharidis 归属补引文号，并写明 GEO 系列标题（愈合设计） ──
    ("inspection of the GEO series matrix and the primary publication (Theocharidis et al., "
     "*Nature Communications*, 2022) established that it is a **NanoString GeoMx Digital Spatial "
     "Profiling (DSP) experiment**",
     "inspection of the GEO series matrix and the primary publication (Theocharidis et al., "
     "*Nature Communications*, 2022 [12]; GEO series title: *Spatial transcriptomics of healing and "
     "non-healing diabetic foot ulcers*) established that it is a **NanoString GeoMx Digital Spatial "
     "Profiling (DSP) experiment**"),

    # ── (D) §2.3：limma / metafor / BH 引文 ──
    ("Differential expression was performed per cohort with limma (R 4.4.3). Effect sizes were then "
     "combined across cohorts by random-effects meta-analysis with restricted maximum likelihood (REML) "
     "in `metafor`, retaining I² as a heterogeneity measure",
     "Differential expression was performed per cohort with limma (R 4.4.3) [18, 25]. Effect sizes were "
     "then combined across cohorts by random-effects meta-analysis with restricted maximum likelihood "
     "(REML) in `metafor` [19], retaining I² as a heterogeneity measure"),
    ("The three prespecified outcome contrasts yielded zero DEG at FDR < 0.05",
     "The three prespecified outcome contrasts yielded zero DEG at FDR < 0.05 (Benjamini–Hochberg [17])"),

    # ── (E) §2.4：WGCNA 引文 ──
    ("Weighted gene co-expression network analysis was run on GSE134431 (top 5,000 expressed genes).",
     "Weighted gene co-expression network analysis was run on GSE134431 (top 5,000 expressed genes) [20]."),

    # ── (F) §2.5：四个算法各自的引文（Breiman 替代 Geurts；新增 Guyon） ──
    ("Four algorithms (LASSO with 10-fold cross-validated penalty, random forest, SVM-RFE, XGBoost) "
     "were run in parallel",
     "Four algorithms (LASSO with 10-fold cross-validated penalty [22], random forest [21], SVM-RFE [24], "
     "XGBoost [23]) were run in parallel"),

    # ── (G) §2.6：EPIC/MCP-counter 引文 + 归档声明据实改写 ──
    ("**EPIC** (6 BRef non-malignant reference populations, NNLS with unit-sum constraint), "
     "**MCP-counter** (10 populations from the official 111-gene signature, implemented from the "
     "ebecht/MCPcounter source), and **ESTIMATE** (140 stromal + 140 immune genes). Group comparisons "
     "used Wilcoxon rank-sum tests with Benjamini–Hochberg FDR correction.",
     "**EPIC** (6 BRef non-malignant reference populations, NNLS with unit-sum constraint) [26], "
     "**MCP-counter** (10 populations from the official 111-gene signature, implemented from the "
     "ebecht/MCPcounter source) [27], and **ESTIMATE** (140 stromal + 140 immune genes). Group "
     "comparisons used Wilcoxon rank-sum tests with Benjamini–Hochberg FDR correction [17]. The "
     "MCP-counter gene signature and the ESTIMATE signature matrix are archived under "
     "`01_数据/基因集/免疫/`; the EPIC reference matrix is distributed inside the `GfellerLab/EPIC` R "
     "package (v1.1.7) and is therefore not duplicated as a flat file, while the calling script "
     "(`02_脚本/immune_deconv.R`) and its outputs (`03_结果/免疫浸润/*_epic_fractions.csv`) are archived."),

    # ── (H) §2.8：空间转录组方法学引文 ──
    ("### 2.8 Spatially resolved analysis (GSE166120, GeoMx DSP)\n\nGeoMx DSP DCC files were parsed",
     "### 2.8 Spatially resolved analysis (GSE166120, GeoMx DSP)\n\nSpatially resolved profiling was used "
     "to ask whether the KpyRN signal varies between anatomical regions of the same ulcer, in the sense "
     "of regional composition rather than whole-transcriptome spatial mapping [28]. GeoMx DSP DCC files "
     "were parsed"),

    # ── (I) §2.9：GSE230426 主要出版物引文 ──
    ("GSE230426 is a longitudinal PBMC cohort of patients with *S. aureus* diabetic foot infection "
     "(20 baseline and 18 post-treatment samples",
     "GSE230426 is a longitudinal PBMC cohort of patients with *S. aureus* diabetic foot infection [14] "
     "(20 baseline and 18 post-treatment samples"),

    # ── (J) §2.10：MR 引文 ──
    ("cis-eQTL instruments for hub genes were taken from eQTLGen lead SNPs and used in two-sample MR "
     "against type 2 diabetes",
     "cis-eQTL instruments for hub genes were taken from eQTLGen lead SNPs and used in two-sample MR [29] "
     "against type 2 diabetes"),

    # ── (K) §2.11：Vina / Vina 1.2.0 / Open Babel 引文 ──
    ("AutoDock Vina 1.2.7 was used with exhaustiveness 8 and seed 42. Receptors were human experimental "
     "structures (PKM2 4JPG, EP300 4PZS, SIRT3 4BN4, LDHA 4JNK, PDK2 4MPC/2BU8/2BU5), protonated and "
     "converted to PDBQT with Open Babel after stripping flexible-residue tags;",
     "AutoDock Vina 1.2.7 was used with exhaustiveness 8 and seed 42 [30, 31]. Receptors were human "
     "experimental structures (PKM2 4JPG, EP300 4PZS, SIRT3 4BN4, LDHA 4JNK, PDK2 4MPC/2BU8/2BU5), "
     "protonated and converted to PDBQT with Open Babel [32] after stripping flexible-residue tags;"),

    # ── (L) §2.12：primer3 编号随重编号改为 33，并补 MIQE / geNorm 落点 ──
    ("Primer sequences (designed with primer3 [34]), Tm and GC windows, secondary-structure checks, "
     "BLAST specificity tables and the LDHA 3'UTR anchoring required by the processed pseudogene "
     "LOC105375980 are documented in Supplementary Materials S4 and Table S10 and are not repeated here",
     "Primer sequences (designed with primer3 [33]), Tm and GC windows, secondary-structure checks, "
     "BLAST specificity tables and the LDHA 3'UTR anchoring required by the processed pseudogene "
     "LOC105375980 are documented in Supplementary Materials S4 and Table S10, reported against the "
     "MIQE guidelines [34] with reference-gene selection by geNorm [35], and are not repeated here"),

    # ── (M) §2.13：软件版本与统计纪律引文 ──
    ("Analyses used fixed seeds; all software versions and file checksums are recorded.",
     "Analyses used fixed seeds; all software versions (including R 4.4.3 [25]) and file checksums are recorded."),

    # ── (N) §3.5 / §3.6 / §3.7：各队列主要出版物引文 ──
    ("In GSE165816, 174,962 cells passed to 167,837 after QC (95.93%).",
     "In GSE165816 [11, 12], 174,962 cells passed to 167,837 after QC (95.93%)."),
    ("**GSE223964 (diabetic versus non-diabetic skin, 27,998 cells, 8 donors).** All 32 donor-level",
     "**GSE223964 (diabetic versus non-diabetic skin, 27,998 cells, 8 donors) [15].** All 32 donor-level"),
    ("In GSE230426 (PBMC, *S. aureus* DFU infection), the whole-transcriptome comparison yielded",
     "In GSE230426 (PBMC, *S. aureus* DFU infection) [14, 16], the whole-transcriptome comparison yielded"),

    # ── (O) §3.8：MR 方法引文 ──
    ("Thirteen of fifteen candidate genes had a significant cis-eQTL lead instrument",
     "Under the two-sample MR framework [29], thirteen of fifteen candidate genes had a significant "
     "cis-eQTL lead instrument"),

    # ── (P) §4.7(8)：删除与自有数据集矛盾的"无数据可答"断言 ──
    ("Relatedly, the healing-outcome arm of the question cannot be answered with currently available "
     "public data: among 36 GEO records matching DFU and healing terms, none provides a "
     "healed-versus-unhealed stratification of human DFU tissue transcriptomes (the two per-cohort "
     "healed-versus-unhealed comparisons yielded zero DEGs and the meta-analysis only 12 meta-DEGs). "
     "Closing that arm would require prospective tissue sampling with 12-week follow-up, which the "
     "design record",
     "Relatedly, the healing-outcome arm of the question remains **unanswered rather than unanswerable**: "
     "exactly one public dataset provides a healed-versus-unhealed stratification of human DFU tissue — "
     "GSE166120 (2 healers versus 2 non-healers, 12-week outcome) — and it was analysed here and was "
     "null at subject level (Section 3.6); the two bulk per-cohort healed-versus-unhealed comparisons "
     "yielded zero DEGs and their meta-analysis only 12 meta-DEGs. With n = 4 subjects and a panel "
     "covering 42.2% of the scoring set, that dataset cannot resolve the outcome question, and closing "
     "the arm would require prospective tissue sampling with 12-week follow-up, which the design record"),

    # ── (Q) 数据可用性声明：Table 1 不存在，改指真实存在的附表 ──
    ("All primary data are public and were obtained from the NCBI Gene Expression Omnibus under the "
     "accessions listed in Table 1;",
     "All primary data are public and were obtained from the NCBI Gene Expression Omnibus under the "
     "accessions listed in Supplementary Table S2;"),
]

fails = []
for k, (old, new) in enumerate(EDITS, 1):
    n = text.count(old)
    if n != 1:
        fails.append((k, n, old[:70]))
        continue
    text = text.replace(old, new, 1)
    print(f"  [{k:2d}] OK")

if fails:
    print("\n★ 以下条目未恰好命中一次，已中止（未写盘）：")
    for k, n, snip in fails:
        print(f"  [{k}] 命中 {n} 次 :: {snip}")
    sys.exit(1)

open(MS, "w", encoding="utf-8").write(text)
print(f"\n[写盘] {MS}（{orig_len} → {len(text)} 字符）")
