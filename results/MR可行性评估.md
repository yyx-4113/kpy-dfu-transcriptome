# 孟德尔随机化（MR）可行性评估与数据源决策记录

> 项目：丙酮酰化（Kpy）与糖尿病足（DFU）｜ 任务 #9（P8 阶段）
> 评估日期：2026-10-24 ｜ 执行人：杨永新（单作者）
> 依据：`00_计划书/项目计划书.md` P8 阶段预设规则——"若 DFU 专用 GWAS 不可得，则采用预设的替代方案 (a) 枢纽基因 cis-pQTL/eQTL → T2D/HbA1c/空腹胰岛素；(b) T2D → DFU 相关中间表型仅作叙述性证据，不作无据 MR；若确实无合格工具变量，如实报告'MR 不可行'，不制造伪 MR。"

---

## 1. DFU 作为结局的 GWAS 检索证据链

| # | 检索/核查 | 工具/来源 | 结果 | 证据状态 |
|---|-----------|-----------|------|----------|
| 1 | GWAS Catalog v2 REST API `?efo_trait=diabetic foot` | `https://www.ebi.ac.uk/gwas/rest/api/v2/studies?efo_trait=diabetic%20foot` | **仅 2 项研究**，均出自同一论文（Meng 2017，GoDARTS）：GCST004725（Foot ulcer and neuropathy in diabetes）、GCST004724（Foot ulcer in diabetes and neuropathy） | ✅ 已核实 |
| 2 | GCST004724/GCST004725 全文统计可得性 | Catalog API `full_summary_stats_available` 字段 | **两项均为 False**（无公开全套 summary statistics；Meng 2017 为 CC BY 论文但未向公共库存放松开统计量） | ✅ 已核实 |
| 3 | FinnGen DFU 端点 | r7/r8/r9 站点不可达（HTTP 000）；r10（2408 端点）、r11（2444）、r12（2470）全量端点表本地检索 | **现行 FinnGen 无 DFU 专用端点**；最接近者为 L12_ULCERLOWLIMB（下肢溃疡 NEC，5324 例，病因混合）与 DM_NEUROPATHY（糖尿病神经病变，3862 例） | ✅ 已核实（端点表已存档 `01_数据/元数据/finngen_r*_phenos.json`） |
| 4 | BMJ Open Diab Res Care 2023（doi:10.1136/bmjdrc-2023-003523）所用结局 | 全文（WebFetch 核实） | 其 "DFU" 结局实为 **Pan UKBB 2020，仅 130 例 vs 420,343 对照**（自我报告）；文中 "DPAD" 为**糖尿病外周动脉疾病**（111,970 例，FinnGen R7），并非 DFU | ✅ 已核实 |
| 5 | Pan UKBB DFU（130 例） | 文献（BMJ 2023 引用） | 病例数 130，**MR 效能严重不足**（SE 不可接受），不作为结局 | ✅ 已核实 |
| 6 | J Diabetes Res 2025（PMID 40190410，HbA1c→下肢溃疡） | 摘要核实 | 该研究为 UKB **单样本内 MR**（需个体级数据），非两样本 MR 可用源 | ✅ 已核实 |
| 7 | DD2/AfterEU/SDR/GoShare/DCS 四队列 DFU meta（980 例，COL4A2/ALK 提示位点） | Dundee/DD2 项目页 + CORE 摘要 | 论文自称 "ongoing study"，**未发表且统计量未公开** | ✅ 已核实 |
| 8 | 2025-2026 更新文献检索 | WebSearch（2026-10-24） | 无新的公开 DFU GWAS 全文统计 | ✅ 已核实 |

### 结论 1（结局侧）
**DFU 作为 MR 结局不可行**：现存唯一 DFU 专用 GWAS（Meng 2017，699 例）不公开全套统计量；次优公开源（Pan UKBB，130 例）效能不足；FinnGen 现行版本无 DFU 端点。**不制造伪 MR**——本稿件不报告任何"X → DFU"的因果估计。

---

## 2. 预设替代方案 (a) 的数据源决策（暴露：枢纽基因 cis-eQTL；结局：糖代谢性状）

计划书 P8 替代方案 (a)："枢纽基因 cis-pQTL/eQTL → T2D/HbA1c/空腹胰岛素（MR，评估'Kpy 机器→糖代谢紊乱'方向）"。

### 2.1 暴露侧：eQTLGen cis-eQTL（公开下载，无需 JWT）
- **来源**：eQTLGen 2018-10-17 significant cis-eQTLs（Vosa et al. 2021 Sci Adv 之前版本，31,684 受试者，欧洲血统）
- **URL**：`https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/cis-eQTL_significant_20181017.txt.gz`
- **许可/存档**：下载后 SHA256 记录于 `01_数据/MR/MANIFEST_mr.csv`
- **工具变量构建规则**（预注册）：
  1. 对每个枢纽基因，取其 1 Mb 窗口内 FDR<0.05 的 cis-eQTL；
  2. **取窗口内 P 值最小的 lead SNP 作为唯一工具变量**（单 IV Wald ratio 设计）——原因：公开 eQTLGen 显著集不携带 LD 信息，无法在无 LD 参考面板下做标准 clumping；此设计为文献常用保守做法，局限性在稿件中如实声明；
  3. 报告每个 IV 的 F 统计量；F<10 者剔除该基因分析。

### 2.2 结局侧：公共 GWAS Catalog 全文统计（无需 JWT）
| 结局 | Accession | PMID | 样本量（EUR） | FTP 路径 | 状态 |
|------|-----------|------|---------------|----------|------|
| 2 型糖尿病 | GCST90018926 | 34594039 | 38,841 例 / 451,248 对照 | `.../GCST90018001-GCST90019000/GCST90018926/GCST90018926_buildGRCh37.tsv.gz` | ✅ 已下载（523 MB） |
| 空腹胰岛素 | GCST90002238 | 34059833 | 151,013 个体（MAGIC/Lagou 2021） | `.../GCST90002001-GCST90003000/GCST90002238/GCST90002238_buildGRCh37.tsv.gz` | ✅ 已下载（480 MB） |
| 2 型糖尿病（第二来源） | GCST006867 | 30054458 | 61,714 例（Xue 2018 Nat Commun） | `.../GCST006001-GCST007000/GCST006867/Xue_et_al_T2D_META_Nat_Commun_2018.gz` | ✅ 已下载（空格分隔格式） |
| eQTLGen cis 显著集 | —（官方 release 2018-10-17） | 31,684 样本 | 19,250 基因 | `.../cis-eqtl/cis-eQTL_significant_20181017.txt.gz` | ✅ 已下载（260 MB） |
| eQTLGen SMR `.esi`（等位基因频率） | 官方 SMR_formatted 包 | — | — | `.../SMR_formatted/cis-eQTL-SMR_20191212.tar.gz` 内 `.esi.gz` | ✅ 已流式提取（110 MB） |

- **HbA1c 结局不可得**：GWAS Catalog 中以 "HbA1c/glycated hemoglobin/haemoglobin A1c" 检索**无任何带全文统计的研究**；MAGIC 官网（magicinvestigators.org）本机不可达（HTTP 000）。→ HbA1c 结局**放弃**并在稿件方法学中如实说明。

### 2.3 分析设计（预注册）
- 主分析：**Wald ratio**（每基因单 IV）；效应方向报告 OR（T2D）与 SD 单位变化（FI）
- 敏感性：F 统计量；Steiger 方向过滤（若暴露/结局均为同一方向可推断）；多效性无法用单 IV 的 MR-Egger 检验——如实声明为方法学局限
- 多重检验：16 枢纽基因 × 2-3 结局，BH-FDR 校正
- **预期局限性**（稿件中声明）：单 IV 设计无法做正式多效性/异质性检验；eQTL 组织来源为血液（非皮肤），外推至皮肤需谨慎；MR 检验的是"枢纽基因表达 → 糖代谢结局"方向，不能反向证明"糖代谢 → DFU"。

---

## 3. 与主假说的关系声明

替代方案 (a) 的 MR 回答的是 **H4 的前半**（"枢纽基因的表达受遗传变异影响并与糖代谢结局存在因果方向"）：若某枢纽基因的 cis-eQTL 对 T2D/FI 有因果效应，则该基因处于糖代谢紊乱的上游通路，支持其作为"Kpy 机器 ↔ 糖代谢"轴的候选节点；若全部为阴性，则如实报告"未发现枢纽基因表达对糖代谢结局的因果效应"，与反证条件一致（转向"微环境代谢状态读数"解释框架）。

**H4 的后半**（药物对接）由任务 #10 独立完成，不依赖本 MR。

---

## 5. 执行结果（2026-10-25，脚本 `02_脚本/mr_eqtl_glycemic.py`）

**工具变量**（13/15 枢纽基因；ADH1B 与 PKM 在 eQTLGen 无显著 cis-eQTL，按预注册规则如实 SKIPPED）：
每基因取 1Mb 窗口内 P 最小的 lead SNP，beta/se 按 eQTLGen 官方公式换算
（beta = Z/√(2f(1-f)(N+Z²))，se = 1/√(2f(1-f)(N+Z²))，f 取自 SMR `.esi` 文件，
N = NrSamples）。全部 F = Z² 在 51–2186 之间（>10 准入达标）。

**结局命中与协调**：T2D GCST90018926 命中 12/13（EP300 rs6002278 不在该 GWAS → SKIPPED）；
FI GCST90002238 命中 13/13（无 rsid 列，按 chr:pos 匹配）；T2D Xue 2018 命中 9/13。
等位基因协调含 strandflip 情形；频率一致性阈值 |Δfreq|<0.2 全部通过（最大 0.033）。

**结果**（34 个基因×结局对，BH-FDR 校正）：

| 基因 → 结局 | Wald β | SE | p | FDR | OR | 方向一致性 |
|---|---|---|---|---|---|---|
| **DHX58 → T2D (Xue 2018)** | −0.0712 | 0.0195 | 2.6e-4 | **0.0088** | 0.931 | 两 T2D 来源同向（−0.071/−0.045） |
| DHX58 → T2D (GCST90018926) | −0.0446 | 0.0168 | 0.0081 | 0.138 | 0.956 | ↑ |
| DHX58 → FI | −0.0104 | 0.0052 | 0.046 | 0.260 | — | 同向（名义） |
| PGAM1 → T2D | −0.0266 | 0.0130 | 0.040 | 0.260 | 0.974 | 两来源同向（−0.027/−0.025） |
| LDHA → FI | +0.0269 | 0.0119 | 0.023 | 0.260 | — | 名义 |
| SIRT3 → FI | −0.0261 | 0.0122 | 0.032 | 0.260 | — | 名义 |

**结论**：34 对中 **1 对通过 FDR<0.05**（DHX58 → T2D，OR=0.931，跨来源方向一致），
6 对名义显著但 FDR 校正后不显著。DHX58（RIG-I，IFN-I 通路传感器）表达越高、
T2D 风险越低——与"IFN-I 信号代谢保护作用"的既有文献方向一致，且方向在两个独立
T2D GWAS 中可复现。PGAM1（糖酵解酶）→ T2D 亦为两来源一致的负向名义趋势。
按预注册反证条件如实报告：**未发现糖酵解枢纽基因表达对糖代谢结局的 FDR 显著因果效应**
（除 DHX58），MR 证据强度有限，稿件中定位为"支持性/探索性"证据而非核心结论。

**验证纪律**：`02_脚本/mr_independent_check.py` 以完全独立代码路径（subprocess+awk
直读原始 gz + 手写 Wald 公式）重算 3 项：① DHX58→T2D Xue（p=0.000259782 逐位一致）；
② LDHA→FI flipped 等位基因场景（p=0.0231395 一致）；③ BH-FDR（max|Δq|=1.1e-15）。
产物：`03_结果/MR/mr_results.csv`、`mr_instruments.csv`、`mr_harmonization_audit.csv`、
`MANIFEST_mr.csv`（5 文件 SHA256）、`figures/fig_mr_forest.png`。

**如实记录的方法学事件**（开发中发现并修复，不影响最终数值）：
1. 频率对齐优先级 bug——回文 SNP（G/C）上 `COMPLEMENT[a1]==a2`，曾把 a2 误判为链翻转
   而使用未翻转频率，致 DHX58/STAT1 在 FI 中被伪剔除；修复后频率差仅 0.018/0.036。
2. Xue 2018 文件为空格分隔（非 TAB），分块读取需自动探测分隔符。
3. eQTLGen significant 文件无 beta/se 列——按官方 README 公式换算并经独立重算验证。

---

## 4. 备选增强（需用户授权）

若用户提供 OpenGWAS JWT（既往慢性疼痛 MR 项目持有），可扩展：LD clumping 参考（`POST /ld/clump`）、pQTL 暴露源、更多糖代谢结局（HbA1c/FG/2hGlu）、以及正式 MR-Egger/加权中位数多 IV 敏感性分析。**未获得 JWT 前不使用任何 JWT 端点**，也不从记忆臆造任何 GWAS ID。
