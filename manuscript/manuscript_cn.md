# 丙酮酰化调控网络在糖尿病足溃疡中的转录组改变：多队列、空间与因果推断研究

**作者**：杨永新

**单位**：福建中医药大学附属第二医院麻醉科，福建 福州 350003

**通信作者**：杨永新（960856791@qq.com）；ORCID 0009-0004-9698-6552

---

## 摘要

### 背景

赖氨酸丙酮酰化（lysine pyruvylation, Kpy）是由 2026 年两篇研究定义的新型蛋白翻译后修饰：其底物为糖酵解终产物丙酮酸，由 HAT1/p300 写入、SIRT3 擦除，作用机制为通过 STAT1 丙酮酰化抑制 I 型干扰素（IFN-I）信号并偶联糖酵解通量与表观遗传转录。糖尿病足溃疡（DFU）创面兼具高糖、缺氧、高糖酵解通量与反复感染，理论上构成 Kpy 高负荷的最适病理场景，但 Kpy 与 DFU 的关系此前完全空白。

### 方法

我们由两篇 Kpy 原始文献加经典通路资源构建 128 个计分基因、五模块、逐基因可溯源的 KpyRN 基因集，在 10 个人类 DFU 相关公共数据集（5 个 bulk、2 个单细胞、1 个空间、1 个纵向感染队列、1 个独立结局队列；其中 7 个进入计算）中评估：逐队列 limma + REML 随机效应 meta 分析；WGCNA 模块–性状关联；四算法共识机器学习（LASSO、随机森林、SVM-RFE、XGBoost，10 个预设随机种子，嵌套交叉验证）；三种独立免疫解卷积方法（EPIC、MCP-counter、ESTIMATE）；cis-eQTL 孟德尔随机化；AutoDock Vina 1.2.7 分子对接（8 个受体–配体体系，含阳性与负面对照）；以及基于实测合并效应量推导样本量、**仅作设计层面预设而未实施**的组织 RT-qPCR 验证方案。

### 结果

两个 DFU 队列 meta 分析中，进入 META1 的 128 个计分 KpyRN 基因有 97 个被覆盖，其中 26 个达到 FDR<0.05（26.8%，非 KpyRN 基因为 24.6%；Fisher 精确 OR 1.12（95% CI 0.71–1.77），p=0.64，即无模块级富集），且 **26/26 在两队列方向一致**：13 个糖酵解酶上调（PGAM1 +1.25 log2FC，padj=1.5×10⁻³²；TPI1 +1.29；GAPDH +1.32；LDHA +1.00；PFKP +0.90；SLC16A1 +0.83），线粒体丙酮酸转运体 MPC2 下调（−0.45）；而 5 个 IFN-I 读出基因下调（IFI44L −1.09；IFI27 −1.02；DHX58 −0.74；DDX58 −0.47；STAT2 −0.41），同时 ISG20、GBP1 上调——提示 ISG 内部异质性，而非 IFN-I 通路整体受抑。三层交集（meta-DEG ∩ WGCNA 模块 ∩ KpyRN）完全由 6 个糖酵解酶构成（BPGM、LDHA、PFKP、PGAM1、PGM2、TPI1），按预设独立性检验全部归于丙酮酸供应模块；IFN-I 模块在该池中贫化（21 个背景基因 0 个进入，期望 3.3；单侧 p=0.026，双侧 p=0.052，BH q=0.18），按方向性观察而非显著结果报告。六基因共识逻辑回归在 GSE199939 外部验证 DFU vs 对照达到 AUC 1.000（完美分离，零假设下精确概率 1/C(21,10)=2.8×10⁻⁶），但前臂皮肤两个预设对比分别为 AUC 0.05（n=9，显著低于随机，精确单侧 p=0.016）与 0.69（n=8），愈合结局任务无交叉验证信号。免疫解卷积在两队列方向并不一致：GSE134431 显示细胞毒性 T 细胞/NK/髓样树突细胞减少，而 GSE199939 表现为单核系、CD4 T、中性粒细胞与成纤维细胞升高、髓样树突细胞降低，且没有任何细胞毒性或 NK 指标显著（FDR 0.30–0.79）。GSE165816 单细胞 pseudobulk 在角质细胞中证实糖酵解协同上调；另两个公共队列（GSE223964 单细胞、GSE166120 GeoMx DSP）在供体/受试体层面**均无显著模块级信号**；纵向金黄色葡萄球菌感染队列（GSE230426）PBMC 中 M4 IFN-I 与 M3 Kpy 靶基因模块评分随感染缓解显著下降，而 6 个枢纽基因全部阴性。孟德尔随机化发现唯一 FDR 显著因果效应：DHX58→2 型糖尿病（Wald β=−0.0712，SE 0.0195，p=2.6×10⁻⁴，FDR=0.0088，OR=0.931）。对接协议验证 8/8 体系通过（阳性对照重对接 RMSD<0.001 Å；刚性对照仅验证盒子与受体准备），29 个候选无一在原始亲和力上超过共晶对照，但按配体效率（每重原子亲和力）归一化后 8 个体系中有 5 个反转排序。

### 结论

KpyRN 的转录组足迹在糖尿病足溃疡组织中确实改变，但方向具有异质性：可重现的信号是组织驻留的糖酵解优势与**选择性**（而非整体性）的 IFN-I 读出衰减，且细胞毒性浸润减少仅见于两个 bulk 队列之一。KpyRN 作为一个基因集合在 meta-DEG 中未富集（26/97 vs 24.6% 背景率；Fisher OR 1.12，p=0.64），可重现的是特定基因（糖酵解核心与部分 ISG）的行为，而非整个网络。阳性证据集中于 bulk 转录组 meta 分析与独立 PBMC 感染队列的模块层；在两个独立单细胞与空间队列的供体/受试体层面不可评估，且不延伸至循环血液。我们识别出一个通过校正的因果估计（DHX58→2 型糖尿病，OR 0.931，95%CI 0.896–0.968），对接覆盖 5 个靶蛋白 8 个口袋，原始亲和力无候选超过共晶，但按分子大小归一化后 8 个中有 5 个反转。Kpy 本身未被测量，组织 RT-qPCR 方案已预设但未执行且不计划执行；本研究将 KpyRN 定位为 DFU 微环境候选**代谢状态读数**，并明确界定公共数据推断的边界，而非已确立的机制。

**关键词**：糖尿病足溃疡；赖氨酸丙酮酰化；丙酮酸；糖酵解；I 型干扰素；STAT1；生物信息学；孟德尔随机化；分子对接

---

## 1. 引言

糖尿病足溃疡（diabetic foot ulcer, DFU）在糖尿病患者中患病率高、负担沉重。依据《中国糖尿病足防治实践指南（Ⅰ）》（2025），我国 DFU 发病率为 8.1%，一生足溃疡风险达 34%，总截肢率 19.03%（2.14% 大截、16.88% 小截），年病死率 14.4%；感染是创面迁延不愈与截肢的首要原因，12 周愈合率 24%–82%，愈合后 1 年复发率（50 岁以上）31.6%[1]。其创面微环境特征为高糖、缺氧、高糖酵解通量与反复细菌感染，并伴局部免疫应答受损。

赖氨酸丙酮酰化（Kpy）是赖氨酸残基的一种翻译后修饰，其独特之处在于底物为糖酵解终产物丙酮酸。2026 年两项研究奠定了其分子基础：Zuo 等（*Cell*）证明糖酵解通量升高经丙酮酸供应驱动 STAT1 在 K201 位点发生丙酮酰化，阻碍 STAT1–STAT2 二聚化，削弱 IFN-I 信号与抗菌抗病毒免疫；高血糖患者外周血单个核细胞 STAT1 丙酮酰化升高、干扰素刺激基因表达受抑（r≈−0.75~−0.82）[2]。Song 等（*Nature Metabolism*）证明 Kpy 是酶促可逆修饰——HAT1 与 p300 以丙酮酰辅酶 A 写入、SIRT3 擦除——CUT&Tag 显示 Kpy 富集于启动子并激活 FGFR4、LGALS8、RSPH4A 等靶基因；缺氧、呼吸链抑制与丙酮酸补充均升高 Kpy[3]。后续综述将"糖尿病与心血管代谢疾病"列为 Kpy 研究的未开发方向[4]。

迄今尚无研究考察 Kpy 调控网络（KpyRN）在糖尿病并发症中的状态，亦无 DFU 转录组研究从 Kpy 或翻译后修饰网络角度切入。这一缺口并非仅是增量性的，原因有二：其一，Kpy 的分子机器明确（"写入–擦除"对 HAT1/p300 与 SIRT3，作用于代谢定义的底物），使网络层面的转录组分析可具机制解释力，而非一份无预设的基因清单；其二，Kpy 的功能读出（IFN-I 介导的抗感染免疫抑制）与 DFU 的临床失效模式高度契合——感染驱动的迁延不愈与截肢。DFU 生信研究中最接近的竞争性翻译后修饰轴是赖氨酸乳酸化，其同样由糖酵解中间体驱动，但缺乏"写入–擦除"架构与直接的抗病毒免疫读出。

据此我们提出四个问题：（Q1）KpyRN 在 DFU 溃疡组织中是否失调，失调方向是否与丙酮酸可利用性升高一致？（Q2）是否存在跨队列、跨算法稳定重现的 KpyRN 核心基因子集，可区分 DFU 与对照并分层愈合结局？（Q3）KpyRN 异常是否与免疫浸润重塑及特定细胞类型耦联？（Q4）KpyRN 枢纽基因是否有影响糖代谢结局的因果证据，其蛋白产物是否可被小分子靶向？

我们预设了反证条件：若 KpyRN 的 meta 信号不高于全基因组背景率、若枢纽基因无法在外部队列重现、或若验证队列呈阴性，我们将如实报告阴性结果并将 KpyRN 重新定位为 DFU 代谢状态的读数，而非机制。我们报告所有此类结局，包括三个验证队列的阴性结果，不做事后重构。

---

## 2. 材料与方法

### 2.1 研究设计与数据来源

本研究为公开、去标识化人类转录组数据的二次分析，未招募新受试者，未接触可识别个人隐私信息，故本节分析无需伦理审查。**预设**湿实验验证方案（2.12 节）仅为设计文件、并未实施，其执行所需的伦理路径见补充材料 S4。本报告遵循 STROBE（观察性队列分析）、STROBE-MR（孟德尔随机化部分）、TRIPOD（预测模型）与 MIQE（预设 RT-qPCR 方案）等社区报告规范。

共使用 10 个 GEO 数据集（附表 S2）：5 个 bulk 队列（GSE134431[5]、GSE199939[6,36]、GSE143735[7]、GSE80178[8]、GSE68183/GSE68186[9,10]）、2 个单细胞队列（GSE165816[11]、GSE223964[15]）、1 个空间队列（GSE166120[13]）、1 个纵向感染队列（GSE230426[14,16]）及 1 个交叉验证用 bulk 队列；其中 GSE80178 与 GSE68183/GSE68186 因 pdInfo 包不可用而无法在样本级处理、未纳入计算。14 个已编目 GEO 系列（283 个样本）的样本级元数据经 NCBI E-utilities 获取，并在项目 MANIFEST 中留存 SHA-256 校验和。

执行过程中对原始数据清单做了两处技术纠正，并在全文中标注：其一，GSE166120 原被登记为 Visium 空间转录组，经 GEO series matrix 与原始文献（Theocharidis 等，*Nature Communications*，2022）核实，实为 **NanoString GeoMx Digital Spatial Profiling（DSP）**实验，采用 Cancer Transcriptome Atlas（CTA）v2.0 靶向板（1,825 个 target），每个 GSM 对应一个 ROI 而非空间点。其二，GSE165816 的基因 ID 经文件级核实即基因符号（ENSG 前缀标识符为 0），故基因标识符映射为经验证的空操作，而非假设。

### 2.2 KpyRN 基因集构建

构建五模块基因集，每个条目均标注具体来源（附表 S1；完整表见 `01_数据/基因集/KpyRN_modules.csv`）：

- **M1 Kpy 机器（写入/擦除/底物；4 基因）**：HAT1、EP300、SIRT3、STAT1，来源 *Nature Metabolism* 2026[3] 与 *Cell* 2026[2]。
- **M2 丙酮酸供应与分流（90 基因）**：糖酵解核心酶系（HK1/HK2/GCK、GPI、PFKM/PFKL/PFKP、ALDOA、TPI1、GAPDH、PGK1、PGAM1、ENO1/ENO2、PKM、LDHA/LDHB、PGAM2、PGAM5、PGM2、BPGM 等）、线粒体丙酮酸输入与氧化（SLC16A1、MPC1、MPC2、PDHA1、PDK1–PDK4、PDHB、PC、ALDH4A1），来源 KEGG hsa00010、hsa00620 及文献[3]的丙酮酸供应分支。
- **M3 已实验鉴定 Kpy 靶基因（3 基因）**：FGFR4、LGALS8、RSPH4A，来自文献[3]的 CUT&Tag 与 RNA-seq。
- **M4 IFN-I 读出（31 基因）**：IFNAR1/2、JAK1、TYK2、STAT2、IRF9、IRF1 及经典 ISG（IFIT1/2/3、MX1/2、OAS1/2/3/L、ISG15、ISG20、IFI6/16/27/44/44L、IFITM1/3、RSAD2、DDX58、DHX58、EIF2AK2、CXCL10、GBP1），来源文献[2]及 KEGG hsa04064、Reactome JAK-STAT。
- **M5 组蛋白载体（12 基因）**：仅保留溯源，**不纳入任何评分**，因组蛋白修饰具细胞周期依赖性，会在增殖旺盛的创缘组织中造成混杂。

故计分集合含 **128 个基因**；M5 仅出现在溯源表中。

### 2.3 预处理、差异表达与随机效应 meta 分析

各队列独立处理，**不做跨平台矩阵拼接**以避免平台构成伪影。GSE134431 为 RPKM、GSE199939 为 TPM、GSE143735 为原始 counts，均做 log2 转换，基因符号按 HGNC 统一并按预设冲突规则（保留平均表达最大值）处理，队列内过滤低表达基因。GSE165816 逐样本处理，QC 条件为 200≤n_genes≤8,000 且 pct_mt≤20；174,962 个原始细胞经 QC 保留 167,837 个（95.93%）。

差异表达用 limma（R 4.4.3），效应量再经 `metafor` 以 REML 随机效应跨队列合并，报告 I² 异质性并保留每基因的队列特异 log2FC 与标准误。事前指定两个 meta 分析：**META1**，GSE134431 与 GSE199939 的 DFU vs 非溃疡/对照皮肤（2 队列）；**META2**，GSE134431 与 GSE143735 的愈合 vs 不愈合。逐队列 DEG 计数：GSE134431 DFU vs 糖尿病足皮肤 2,814/12,450（817 上调、1,997 下调）；GSE199939 DFU vs 对照 9,773/17,239（3,958 上调、5,815 下调），与该队列原始文献报道的 8,800 个 DEG 一致[6]。三个预设结局对比均为零 FDR<0.05（GSE134431 愈合 vs 不愈合；GSE143735 愈合 vs 不愈合；GSE143735 DFU vs 糖尿病无溃疡）。

### 2.4 WGCNA 与三层交集

在 GSE134431（表达量前 5,000 基因）上运行加权基因共表达网络分析。因 n=21 时软阈值拟合不稳，软阈值 power 取 6 作为兜底，此点已记录为局限；模块结果仅用作交集层，不作独立生物学论断。共检出 10 个模块（9 个非 grey），其中 4 个与 DFU 状态显著相关。机器学习候选池定义为两层交集（meta-DEG ∩ KpyRN）与三层交集（meta-DEG ∩ DFU 相关 WGCNA 模块 ∩ KpyRN）的并集，共 26 个基因。

### 2.5 共识机器学习

在发现集 GSE134431 上，以 10 个预设随机种子并行运行四种算法（LASSO 以 10 折交叉验证选惩罚参数、随机森林、SVM-RFE、XGBoost），设两个任务：Task A，DFU vs 糖尿病足皮肤（诊断）；Task B，愈合 vs 不愈合（结局）。被至少 3/4 算法选中者定为共识枢纽基因。性能以训练集重复 5 折交叉验证与外部验证（GSE199939 承担 Task A，GSE143735 承担 Task A 与 Task B）评估，报告 AUC 与 DeLong 置信区间。所有六个共识基因在两个外部队列均存在，故无评估被跳过（缺失基因本将记为 SKIPPED 而非静默填补）。

### 2.6 免疫浸润解卷积

因 `immunedeconv` 已从 Bioconductor 下线、`xCell2` 需 R≥5.0（本机 R 4.4.3 / Bioc 3.20 不提供），我们直接依据已发表签名与源码自行实现三种解卷积方法，签名与源码全部归档于 `01_数据/基因集/免疫/`：**EPIC**（6 类 BRef 非恶性参考群，NNLS 且比例和为 1）、**MCP-counter**（官方 111 基因签名，10 群，按 ebecht/MCPcounter 源码实现）、**ESTIMATE**（140 基质 + 140 免疫基因）。组间比较用 Wilcoxon 秩和检验 + BH 校正；枢纽基因与 KpyRN 模块评分与免疫估计量的关联用 Spearman 相关 + FDR 校正。

### 2.7 单细胞分析

GSE165816（54 样本）以标记评分法注释（CP10K + log1p 标度下 23 个谱系标记集均值的 argmax），top1−top2<0.25 判为 Ambiguous（12.7%）。组织学 sanity check 通过（PBMC 84.7% 血液谱系、足皮肤 73.4%、前臂皮肤 79.8% 实质谱系）。KpyRN 模块与枢纽基因评分由 donor × 谱系的 pseudobulk CPM（log2）计算，仅纳入细胞数 ≥100 的 donor×谱系组合，在供体层面以 Wilcoxon 秩和检验 + BH 校正比较。此为标记评分注释而非参考转录组投影（未使用 SingleR/Azimuth），亚谱系结论仅为探索性。

GSE223964 为已处理 h5ad（27,998 细胞 × 18,287 基因；X 已 log1p，无 raw counts），pseudobulk 取每位供体每个细胞类型 log1p 均值（不做 counts 求和），以 Mann–Whitney U（exact）比较糖尿病供体（n=5）与非糖尿病供体（n=3）皮肤的模块评分。另以独立脚本实现第二条供体级路径（中位 pseudobulk + 10,000 次标签置换），以防计算缺陷伪装成阴性结果。

### 2.8 空间/区域分析（GSE166120，GeoMx DSP）

GeoMx DSP DCC 文件按 probe 级解析，经 CTA v2.0 PKC 的 `SystematicName` 字段映射到基因：取每个 probe 自身 `SystematicName` 的首个非空值（弃用 `DisplayName`，因其使 24 个 target 混编多基因 probe；亦不能只取 `SystematicName` 首元素，因部分 probe 该字段为空列表会致索引错位），最终 8,670/8,670 probe、1,825/1,825 target 全部命中、零遗漏。88 个 ERCC 负对照 probe 予以剔除。因 DCC 仅提供 probe 级读段计数而无上游 DSQ2 质控，未做 GC/面积/核数归一化，仅跨 ROI 比较相对丰度（构成比），不作任何绝对定量声明。模块评分先**逐基因跨 ROI** 标准化（而非逐 ROI 内跨基因——后者数学上恒为 0），再跨基因取均值。推断主单位为**受试体**（n=4；ROI 嵌套于受试体）；ROI 级检验（n=23）仅作敏感性分析。

### 2.9 纵向感染队列分析（GSE230426）

GSE230426 为金黄色葡萄球菌糖尿病足感染患者的 PBMC 纵向队列（基线 20 例、治疗后 18 例，其中 8 周 15 例、4 周 3 例）。公开元数据未提供愈合结局，故不构建愈合/不愈合对比。公开数据无显式 donor 列，配对依"连续样本编号 + 同年龄 + 互补时点"推断（17 对、4 单例），故**主模型采用非配对 limma-voom**，推断配对与仅 8 周后模型作敏感性分析。模块层用 `camera()` 对 128 个计分 KpyRN 基因作竞争性基因集检验。

### 2.10 孟德尔随机化

按事前登记的可行性评估，取 eQTLGen lead SNP 作为枢纽基因的 cis-eQTL 工具变量，以 Wald 比值法对 2 型糖尿病（两个独立来源）与空腹胰岛素作两样本 MR；等位基因协调过程全部记录。工具变量强度以 F 统计评估，要求 F>10。

### 2.11 分子对接

使用 AutoDock Vina 1.2.7，exhaustiveness=8、seed=42。受体为人源实验结构（PKM2 4JPG、EP300 4PZS、SIRT3 4BN4、LDHA 4JNK、PDK2 4MPC/2BU8/2BU5），经 Open Babel 加氢并剥离柔性标签转 PDBQT，对接盒以共晶配体质心为中心。配体由 meeko/RDKit 制备，来源为 pChEMBL≥5 的 ChEMBL 记录与文献核实化合物。**协议验证**：在全部 8 个受体–配体体系中重对接共晶配体并计算重原子 RMSD（Kabsch，仅取 MODEL 1，按索引配对），验收标准为 RMSD<2 Å。结合能以同一口袋的共晶配体为基准报告（Δ = 候选 − 共晶），并在每个体系加入蔗糖作负面对照。全程将对接定位为产生假说的筛选手段，绝不作为结合证据。

### 2.12 预设组织 RT-qPCR 验证方案（设计留档）

该方案为设计层面预设、**未实施且不计划实施**，仅作为设计留档，不向本稿提供任何数据、对任何论断均无证据权重。原设计覆盖 6 个共识枢纽基因与 6 个内参候选（ACTB、GAPDH、TBP、HPRT1、RPLP0、SDHA），样本量基于 **META1 实测**合并效应（Cohen's d = pooled logFC / pooled SE），α=0.05（双侧）、power=0.80，在配对与独立两种设计及"效应量减半"保守情景下分别计算，目标 n=10/组（按 20% 脱落招募 12/组）；后因设计层面收缩，最小可行集改为两个轴代表靶基因——LDHA（糖酵解轴）与 DHX58（IFN-I 轴）——加两个内参，由 HPRT1、RPLP0、TBP 三个通过 BLAST 的候选中事后选 2 个，并改写为补充材料 S4b 的 H1–H4 四条定向假设。内参设 6 个候选，待事后以 geNorm 与 NormFinder 排序而非预先指定——因 GBP1 与 DHX58 本身为干扰素刺激基因、不可作内参，且 GAPDH/ACTB 在 DFU 组织中受炎症细胞浸润干扰。引物用 primer3-py 2.4.0 在 RefSeq 规范 NM_ 记录上设计（accession 由 NCBI E-utilities 程序化解析，非凭记忆书写），目标产物 80–180 bp、Tm 58–62 °C、GC 40–60%、Poly-X≤4、预测二级结构低于 47 °C。原设计另含一个预设阳性判定标准（多数目标基因方向与 meta 分析一致且 BH-FDR<0.05），随靶基因集收缩至 3 个而失去区分度，留档不再执行。

**设计层面的收缩（其结果本身就是发现）**：独立性检验提示 KpyRN 为两条部分独立的轴而非一个共调控模块（此解读仅作假设，见 §4.1）；BLAST 最终重跑（`core_nt` + 人源限定 + 分层判据，12 对全量复核）判定 **6 对可用、6 对 OFFTARGET_P0**——三个代谢靶基因（FAHD2B、ADH1B、ALDH3A2）与三个内参候选（ACTB、GAPDH、SDHA）因 3' 端（近）完美脱靶匹配被淘汰（ADH1B→ADH1A/ADH1C 旁系、ACTB→POTEM/ACTBL2、FAHD2B/GAPDH/SDHA→未定性 LOC 位点、ALDH3A2 反向引物 210 条 95% 命中），任何重设计均无法消除；合格靶基因只剩 3 个，原"≥3/6"标准随之失效，留档不再执行。经典内参 ACTB/GAPDH 的落选本身即说明：内参不能凭文献沿用，必须逐引物实测。详见补充材料 S4/S4b 与附表 S10。

### 2.13 统计原则、可复现性与审计

本稿每一个数字均登记于机器可读审计台账（`05_稿件/audit_facts.json`），其中记录每条数值的来源文件与提取表达式；审计脚本据此从台账**动态派生**期望值并反向核对稿件，不硬编码任何期望值。本次研究在执行中依此纪律发现并修复了若干上游缺陷（见 2.14 节），如实报告。分析使用固定随机种子，软件版本与文件校验和全部留存。

---

### 2.14 分析过程中发现并修复的缺陷

本研究在执行中依审计纪律发现并修复四类上游缺陷，因其影响上述数字的解读方式而在此如实报告。（i）PDBQT 解析器累加了 Vina 输出中所有 MODEL 块的原子，导致含 2 或 9 个 MODEL 的受体（1OX、TF1、PV2、TF4）索引配对失败，4 个体系（含全部 3 个 PDK2 位点）的阳性对照 RMSD 缺失，而 4 个单 MODEL 体系侥幸通过掩盖了该缺陷；改为仅解析 MODEL 1 后，8 个体系全部重算为 <0.001 Å。（ii）某图件脚本对两个原子数不同的分子施加 Kabsch 叠加，属维度非法；改为在受体坐标系中直接作 2D 投影并在图注中明示。（iii）GSE223964 中分组变量为分类变量，pandas 默认 `observed=False` 将其展开为 16 组笛卡尔积（其中 8 组为空），产生 NaN 均值、虚报分组计数（8 vs 8）与完全空白的 p 值列；改为 `observed=True` 后恢复真实供体数（5 vs 3）与有效统计量。同一分析另曾以配对 Wilcoxon 检验比较两组独立供体，纠正为 Mann–Whitney U（exact），并在同一次提交中同步修正下游列引用。（iv）在每个 ROI 内跨基因做 z-score 在数学上恒为 0（z 行均值必为 0），约 10⁻¹⁷ 的量级曾被误读为"运行成功"；标准化方向改为逐基因跨 ROI，并以手工复算一个 ROI–基因对（1.011253）验证。GeoMx 系列另有三处：重复出现且行标签相同的 `!Sample_characteristics_ch1` 行被按标签存为字典（静默丢字段）、GSM 标识被从标题行而非 accession 行读取、`Non-ulcer`/`Non-Ulcer` 拼写变体把同一解剖分层拆成两组。此外数据清单本身有一处错误：GSE166120 被登记为 Visium 而实为 GeoMx DSP，这使该队列的证据权重由空间转录组降为部分覆盖的靶向分析。

---

## 3. 结果

### 3.1 KpyRN 在 DFU 组织中的表达特征：方向一致但无模块级富集

META1（DFU vs 非溃疡皮肤，2 队列，REML 随机效应）检验 11,615 个基因，FDR<0.05 得 2,863 个 meta-DEG。128 个计分 KpyRN 基因中 97 个进入该 meta 分析，其中 26 个显著，且 **26 个全部在两队列方向一致**（26/97=26.8%，非 KpyRN 基因为 24.6%；Fisher 精确 OR 1.12（95% CI 0.71–1.77），p=0.64，即无模块级富集）。信号由两个方向相反的模块主导：13 个糖酵解酶上调（M2），5 个 IFN-I 读出基因下调（M4）。代表性糖酵解效应为 PGAM1 +1.25 log2FC（padj=1.5×10⁻³²）、TPI1 +1.29（1.4×10⁻⁶）、GAPDH +1.32（1.0×10⁻²）、PGK1 +1.07（1.4×10⁻¹⁵）、LDHA +1.00（1.1×10⁻⁹）、PFKP +0.90（3.3×10⁻⁷）、SLC16A1 +0.83（2.4×10⁻⁴）、LDHB +0.75、PGM2 +0.67（2.4×10⁻³）、ALDOA +0.66、ENO1 +1.24（3.2×10⁻²）、BPGM +0.83、PDHB +0.26。IFN-I 读出中，IFI44L −1.09（1.3×10⁻⁶）、IFI27 −1.02（1.8×10⁻⁶）、DHX58 −0.74（8.6×10⁻⁵）、DDX58 −0.47（9.1×10⁻³）、STAT2 −0.41（3.5×10⁻³），而 ISG20（+0.88，3.9×10⁻³）与 GBP1（+0.84，1.3×10⁻³）上调；线粒体丙酮酸转运体 MPC2 下调（−0.45，3.0×10⁻³），ADH1B（−4.42）、ALDH3A2（−1.63）、FAHD2B（−1.12）亦下调。三者附明确异质性警告：ADH1B 队列特异估计 −2.79（SE 0.44）vs −6.29（SE 1.00）（I²=90.2%）、ALDH3A2 −1.00（SE 0.28）vs −2.24（SE 0.23）（I²=91.2%）、FAHD2B −0.87（SE 0.12）vs −1.38（SE 0.14）（I²=87.1%），其合并值平均了两个实质不同的队列效应，不应读作单一共同效应。全部 11,615 个受检基因中位 I² 为 61.2%。

该结果有两点约束其解释。第一，计分 KpyRN 基因的显著比例（26/97=26.8%）**并未**超过同一 meta 分析的全基因组背景率（2,837/11,518=24.6%），故不存在模块级富集；KpyRN 应被理解为"机制上有预设、部分存活的假说"，而非统计上富集的基因集。第二，ISG 异质性具生物学意义，不支持"IFN-I 整体受抑"的简化叙事：GBP1 与 ISG20 升高而 IFI44L、IFI27、DHX58、DDX58、STAT2 下降，更符合干扰素刺激转录本**组成**的改变（例如转向 NF-κB 驱动的 GBP1）而非 IFN-I 信号通路的整体阻断。META2（愈合 vs 不愈合）在 10,300 个受检基因中仅 12 个 meta-DEG，提示结局对比的信号远弱于发生对比，与逐队列阴性结果一致。

### 3.2 三层交集由丙酮酸供应模块驱动，而非整个网络

与 DFU 显著相关的模块为 magenta（r=−0.79，p=2.1×10⁻⁵）、yellow（−0.74）、red（−0.73）、blue（+0.70）。meta-DEG ∩ DFU 相关 WGCNA 模块 ∩ KpyRN 三层交集为 **6 个基因，全部为糖酵解酶**：BPGM（+0.83，padj 0.026）、LDHA（+1.00，1.1×10⁻⁹）、PFKP（+0.90，3.3×10⁻⁷）、PGAM1（+1.25，1.5×10⁻³²）、PGM2（+0.67，2.4×10⁻³）、TPI1（+1.29，1.4×10⁻⁶）。与两层交集取并集得 26 基因候选池。

六个基因均标注于 KEGG hsa00010，丙酮酸供应模块本身即由该通路构建，故收敛本身权重有限，我们预设独立性检验而非把交集当 corroboration（附表 S11；`03_结果/WGCNA/intersection_independence.json`）。剔除丙酮酸供应模块的 90 个基因后三层交集为空，Kpy 机器（M1）、Kpy 靶基因（M3）与 IFN-I 读出（M4）模块均无成员进入。KpyRN 基因进入 DFU 相关池的数量不超随机期望——WGCNA 限于 5,000 高变异基因，128 个计分基因中 56 个位于该背景，其中 6 个进入 vs 期望 8.9（超几何 p=0.19）。IFN-I 模块表现相反：它在池中贫化（21 个背景基因 0 个进入，期望 3.3，单侧 p=0.026、双侧 p=0.052、跨七子集 BH q=0.18），其中 14 个落入 turquoise 模块（与 DFU 状态无关，r=−0.22，p=0.33）；该方向未预设，故按方向性观察而非显著结果报告。六基因交集说明糖酵解臂与 DFU 转录组共变，而 IFN-I 臂不共变——这与"一个共调控 KpyRN"相反。因 n=21 已是 WGCNA 实际下限且软阈值拟合不稳，这些模块仅用作交集过滤。

### 3.3 六基因共识分类器可区分足皮肤 DFU，但不可迁移

Task A（诊断）发现集交叉验证 AUC：LASSO 0.983±0.022、随机森林 0.969±0.021、SVM-RFE 0.955±0.029、XGBoost 0.911±0.067。在外部 GSE199939 中，LASSO、随机森林、SVM-RFE 均达 AUC 0.991（95%CI 0.967–1.000），而 XGBoost 为 0.500：其特征重要性被单一基因 FAHD2B 支配（importance 0.97），而该基因在 GSE199939 中近乎恒定。这实证了单算法选基因的不可靠性，也说明多算法共识步骤是必要的而非装饰性的。六基因共识逻辑回归模型（ADH1B、ALDH3A2、DHX58、FAHD2B、GBP1、LDHA）在 GSE199939 达 **AUC 1.000**，模型将这 21 例样本完全分离，DeLong 方差为零、区间塌缩至 1.000–1.000；有信息量的量是零假设下完美分离的精确概率 1/C(21,10)=2.8×10⁻⁶。共识基因集与共识模型均在完整训练集上选定与拟合，故该数字背后无内部验证；应读作"本样本近乎完美可分"这种由构造决定的大效应，而非经校准的性能估计。所有六个共识基因在两个外部队列均存在，故无评估被跳过或替换。

三项限定仍须说明。其一，前臂皮肤迁移在 GSE143735 以两个预设对比评估：DFU 愈合者 vs 非溃疡糖尿病皮肤（n=9；LASSO 0.10、随机森林 0.35、SVM-RFE 0.20、共识 0.05）与 DFU 未愈合者 vs 非溃疡糖尿病皮肤（n=8；0.50、0.31、0.81、0.69）。第一个对比中共识 AUC 0.05 显著低于随机（精确单侧 p=0.016），即签名被系统性反转而非仅无信息；第二个对比则否。迁移失败因而是依赖对比的，模型无论何种异质性均须限定于足皮肤材料。其二，我们只报告判别力，不报校准统计量（Brier 分数、校准斜率）：一个在 21 样本上拟合的六预测因子逻辑模型按构造即过拟合，其校准估计无意义。其三，共识模型不应与临床可用检测混淆；它只是一个用来探问 KpyRN 基因是否携带任何组间信息的工具。

Task B（愈合结局）无交叉验证信号（四算法 AUC 0.39–0.55）。该任务外部队列仅 9 例，单个模型达 AUC 1.000、另一个 0.850，在无内部信号锚定时这些值不可解释——它们只是 9 个观测上统计量的散布，仅列于附表 S4 备查。我们**不主张** KpyRN 可预测愈合结局。

### 3.4 免疫解卷积显示 DFU 组织细胞毒性与单核浸润减少

GSE134431（DFS vs DFU，8 vs 13）中 6 项免疫指标达 FDR<0.05 且方向一致：EPIC CD8 T 细胞（FDR=0.0063）、MCP-counter CD8 T 细胞（0.0063）、MCP-counter 细胞毒性淋巴细胞（0.0063）、MCP-counter T 细胞（0.0092）、MCP-counter NK 细胞（0.0101）、MCP-counter 髓系树突细胞（0.0151）。GSE199939（DFU vs 对照，10 vs 11）中同样有 6 项显著：MCP-counter 单核细胞系（0.0025）、EPIC CD4 T 细胞（0.0115）、MCP-counter 成纤维细胞（0.0115）、EPIC otherCells（0.0268）、MCP-counter 髓系树突细胞（0.0268）、EPIC 中性粒细胞（0.0414）。GSE134431 愈合 vs 不愈合无 FDR 显著差异。

两个 bulk 队列因此并未趋同于单一免疫表型，细胞毒性浸润减少只是单队列发现。GSE134431 中细胞毒性浸润减少与 3.1 节的 IFN-I 读出下调方向一致，但这只是 GSE134431 内部的自洽，并非复制；GSE199939 中没有任何细胞毒性或 NK 指标显著（FDR 0.30–0.79）。两队列也不构成跨方法互证——平台不同，比较二者检验的是可重复性而非方法；真正的跨方法比较是同队列内 EPIC 与 MCP-counter，而该一致性随指标而异。枢纽基因–免疫相关广泛（GSE134431 60/300 对、GSE199939 62/300 对 FDR<0.05），最强为 GSE134431 中 DHX58 与 MCP-counter CD8 T 细胞（ρ=+0.840）。EPIC 与 MCP-counter 一致性：B 细胞高（Spearman 0.87/0.90）、CD8 T 细胞高（0.91/0.95）、中性粒细胞较低（0.54/0.93）、NK 弱（−0.40/0.29）；EPIC CD4 与 MCP T 细胞的负相关是 EPIC 比例和约束的预期后果，已如实报告。本文任何免疫方向论断都不依赖单一队列。

### 3.5 单细胞分析将糖酵解偏移定位到角质细胞

GSE165816 中 174,962 个细胞经 QC 保留 167,837 个（95.93%）。足皮肤供体层面（14 DFU vs 19 非 DFU）4 个谱系达 FDR<0.05：Mono_CD14（中位占比 0.0099 vs 0.0010，q=6.0×10⁻³）、Mono_CD16（0.0071 vs 0.0002，q=6.0×10⁻³）、Schwann/Neural（0.0005 vs 0.0026，q=0.017）、NK（0.0047 vs 0.0017，q=0.037）。前臂皮肤（5 vs 6）、PBMC（5 vs 5）、足皮肤 4 组 Kruskal–Wallis 及愈合 vs 不愈合对比均无 FDR 显著差异。

供体 × 谱系层面的 KpyRN pseudobulk 呈现一致图景：**角质细胞糖酵解酶在 DFU 中协同上调**——PGAM1（log2 CPM 8.37 vs 7.29）、TPI1（9.26 vs 8.52）、LDHA（9.99 vs 9.30）、ENO1（9.85 vs 9.08）、PKM（9.86 vs 8.87）、GAPDH（11.97 vs 11.27）、HK2（4.12 vs 2.58）、PGM2（5.51 vs 4.86），均 q<0.05——在单细胞分辨率与特定细胞类型中重现了 bulk meta 方向，构成跨平台互证。成纤维细胞与周细胞/SMC 中 ADH1B 降低（分别为 2.37 vs 6.57、0.67 vs 3.40）而 ISG 转录本升高（IFITM3、ISG20、IFI16）。ISG 模式同样异质（IFITM3、IFITM1、ISG20、IFI16、IFIT2 升高；IRF9、JAK1、IFNAR1、EIF2AK2 降低），与 bulk 观察一致，进一步支持"改变具选择性而非 IFN-I 整体关闭"。

两项方法学限定已记录。全部细胞中位数因大量零表达细胞而被压低（PGAM1 仅 53% 细胞可检出、ADH1B 18%、IFI44L 5%），故表达模式必须以 pseudobulk 聚合值或"表达细胞比例"解读，绝不可与中位数混谈。此外 bulk 解卷积（细胞毒/单核浸润减少）与单细胞直接计数（单核/NK 比例增加）方向不同：因对照组成不同（bulk 对照为糖尿病非溃疡足皮肤；单细胞对照为非糖尿病与非溃疡糖尿病皮肤）且方法不同（解卷积 vs 直接计数），我们并列呈现而不作单向裁决。

### 3.6 两个独立验证队列未提供阳性模块级支持

**GSE223964**（糖尿病 vs 非糖尿病皮肤，27,998 细胞，8 位供体）：全部 32 个供体级"细胞类型×模块"对比均为阴性（Mann–Whitney U exact，无 FDR<0.10）。独立的第二条供体级路径（中位 pseudobulk，10,000 次标签置换，统计量不同）同样无 FDR<0.10（3 个名义 perm p<0.10），确认阴性为真实结果而非计算缺陷。个别对比中可见方向收敛（Blood Vessels 的 M4 IFN-I 读出糖尿病皮肤更高，均值路径 +0.0423、置换路径 +0.0351；Keratinocytes M4 +0.0348/+0.0338）。最有信息量的细胞层面观察是：**在糖尿病皮肤中 LDHA 可在免疫细胞中检出**（淋巴细胞中位 1.061、巨噬细胞 0.914、肥大细胞 1.243），而相应非糖尿病细胞全部为 0.000；成纤维细胞（1.925 vs 1.937）、角质形成细胞（2.228 vs 2.457）、周细胞（1.696 vs 2.058）中糖尿病组反而略低——即 LDHA 在本队列呈**免疫细胞限制型**升高，而非全皮肤升高。细胞层面（我们标注为描述性，因细胞嵌套于供体、该分析存在伪重复）显示 ALDH3A2 阳性细胞比例在糖尿病皮肤中下降（19.0% vs 26.1%），FAHD2B（3.95% vs 3.24%）与 DHX58（1.51% vs 1.10%）略升。本数据集的限度须强调：在 5 vs 3 位供体下，32 个"细胞类型×模块"对比经 BH 校正后最小可达 q 为 1.0，故全阴性是设计属性而非实证发现；6 个枢纽基因中 5 个在多数细胞类型中位数为 0；矩阵为 log1p 且无 raw counts，pseudobulk 无法与 bulk RNA-seq 直接比较。

**GSE166120**（GeoMx DSP，4 位受试体 23 个 ROI，2 愈合、2 不愈合）：CTA v2.0 板仅覆盖 128 个计分 KpyRN 基因中的 54 个（42.2%），模块覆盖为 M4 IFN-I 读出 26/31（83.9%）、M1 2/4（50.0%）、M3 1/3（33.3%）、M2 丙酮酸供应 25/90（27.8%）；6 个枢纽基因中仅 LDHA、ADH1B、GBP1 在板内，FAHD2B、DHX58、ALDH3A2 确实缺失（已逐 probe 核实）。受试体层面（n=4）全部对比为阴性（Healer vs Non-Healer MWU exact p=1.0；无 FDR<0.10），且因 epidermis 仅 1 位受试体、非溃疡仅 2 位，ROI 位置分层不可检验。ROI 层面 M4 IFN-I 读出的名义 p=0.0405（ulcer −0.29 vs non-ulcer +0.45）未通过 FDR 校正（0.1731），且在受试体层面重算时完全消失；驱动力为 Subject 3 的单个极端 ROI（−1.45）。我们判定其为伪重复驱动的假阳性，不作为阳性结论报告。由于板覆盖不足评分集一半，板内子集模块分数不可与全 KpyRN 分数比较，本实验只能提供探索性的区域定位证据。

### 3.7 纵向感染队列显示模块层变化但无枢纽基因信号

GSE230426（PBMC，金黄色葡萄球菌 DFU 感染）中，全转录组比较在**全部三个模型**（非配对、推断配对、仅 8 周后）均为**零个 FDR<0.05 的 DEG**（受检 57,059 个基因）。预设的模块层竞争性检验则有信息：主非配对模型中 M4 IFN-I 读出模块（31 基因）下调（FDR=0.0094）、M3 Kpy 靶基因模块（3 基因）下调（FDR=0.0094），M2 丙酮酸供应名义下调（FDR=0.061）；推断配对敏感性模型方向一致（M3 FDR=0.021；M4 FDR=0.052）；仅 8 周后模型证据更强（M4 FDR=0.0016；M2 FDR=0.0152；M3 FDR=0.0159，均下调）。主模型与配对模型在 128 个 KpyRN 基因上的 logFC 一致性高（Spearman r=0.930；方向一致率 88.3%）。

其解释受两项事实约束。第一，方向与组织中观察**相反**：IFN-I 读出在 DFU 溃疡组织相对非溃疡皮肤为下调（3.1 节），而在 PBMC 活动性感染期高于缓解期。这是两个不同生物学过程——组织局部 IFN-I 抑制与循环系统感染相关免疫激活——不可混同表述。第二，**6 个共识枢纽基因在 PBMC 中全部不显著**（|logFC| 均低于 0.52，最小 adj.p 为 0.41，三模型一致）。这是一项实质性的阴性结果：KpyRN 糖酵解枢纽信号呈**组织驻留型**而非循环型，不支持其作为血循环标志物，反向支持 RT-qPCR 验证取**组织**样本而非外周血。

### 3.8 孟德尔随机化：一个 FDR 显著估计，位于 DHX58

15 个候选基因中 13 个有显著 cis-eQTL lead 工具变量（ADH1B 与 PKM 无，记为 SKIPPED），F 统计 51 至 2,186，全部 >10。34 个基因–结局配对经 BH-FDR 校正后恰有 1 对显著：**DHX58 → 2 型糖尿病**，Wald β=−0.0712（SE 0.0195），p=2.6×10⁻⁴，FDR=0.0088，OR=0.931，工具变量 F=1,189。第二个独立 T2D 数据集方向一致（β=−0.0446，名义 p=0.0081，FDR=0.138）。另 4 个基因–结局配对名义显著（p<0.05）但未通过校正：PGAM1→T2D 仅两个 T2D 来源之一显著（β=−0.0266，p=0.040；另一来源 p=0.11），LDHA→空腹胰岛素（β=+0.0269，p=0.023）、SIRT3→空腹胰岛素（β=−0.0261，p=0.032）与 DHX58→空腹胰岛素（β=−0.0104，p=0.046）。故我们将 MR 定位为**支持性但较弱的**证据：它与"组织中 M4 IFN-I 读出减弱"这一发现方向自洽（DHX58 表达越高，T2D 风险越低，而 DHX58 表达降低伴随 IFN-I 受抑、易感感染的织状态），但 34 次检验中仅 1 对 FDR 显著并不构成因果机制。

### 3.9 对接筛选五个靶蛋白，但无候选优于共晶对照

协议验证在全部 8 个体系通过：共晶配体重对接在每个体系均以 <0.001 Å 重原子 RMSD 重现实验构象（PKM2 C-pocket、PKM2 A-pocket、EP300 乙酰辅酶 A 位点、SIRT3 ADP-核糖/NAD⁺ 口袋、LDHA 活性位点、PDK2 ATP 口袋、PDK2 丙酮酸位点、PDK2 R 域 N 端位点），远优于预设的 <2 Å 验收标准；但刚性共晶配体由其自身晶体构象在以其质心为中心的盒子中对接，构象恢复近乎同义反复，故该验证仅确认盒子与受体准备，并不验证对柔性候选施加的扭转搜索（候选具 11–15 个可旋转键，而刚性参照不带分子内应变项，使 Δ 偏向共晶参照）。55 次对接（29 个独特候选、8 个阳性对照、每体系 1 个负面对照）中，各体系最优候选为：LDHA CHEMBL4782843（−10.35 kcal/mol；共晶 −13.16，Δ=+2.81）与 CHEMBL4859303（−10.29，Δ=+2.87）；PDK2 ATP 口袋 CHEMBL6011956（−8.01；共晶 −9.58，Δ=+1.57）；PDK2 R 域 CHEMBL5896437（−9.06；共晶 −11.39，Δ=+2.33）；PKM2 A-pocket CHEMBL5590652（−8.92；共晶 −12.22，Δ=+3.30）；EP300 CHEMBL6021203（−10.79；共晶 −17.10，Δ=+6.31）；SIRT3 CHEMBL2332041（−11.50；共晶 −21.26，Δ=+9.76）；PKM2 C-pocket TEPP-46（−10.25；共晶 1OX −20.30，Δ=+10.05）。文献参照配体表现与已知定位一致：LDHA 的 FX11（−7.98）与 galloflavin（−7.88）、EP300 的 A-485（−10.58）与 C646（−10.17）、SIRT3 的 quercetin（−9.61）与 sirtinol（−8.70），以及弱效的 oxamic acid（−4.07，与其已知弱抑制一致）。

负面对照（蔗糖）在每个体系均属最弱结合者之一（8 个体系亲和力自 −4.01 至 −7.45 kcal/mol），提示无非特异性假阳性信号。有一个候选名义上优于其共晶配体（CHEMBL5777543 于 PDK2 丙酮酸位点，−4.82 vs −4.08，Δ=−0.73）；但在同一浅口袋中蔗糖（−4.01）与共晶不可区分，口袋深度低（Δ≈−0.73 kcal/mol，对应 22×22×22 Å 盒），该协议的分辨力不足。故我们**降级**该观察，不主张"优于共晶"。同理，PKM2 C-pocket 体系由 42 重原子的共晶配体主导，其绝对分值异常强，该体系以 A-pocket 交叉印证而非绝对值解读。3 个 ChEMBL 候选虽有功能活性但位点未注释，已在三个 PDK2 口袋同时对接并如实报告位点不确定性。所有候选均表述为需实验验证（SPR/ITC/酶活/细胞实验）；对接是筛选手段，不是结合证据。

另一项约束支配 Δ 比较并改变结论。AutoDock Vina 分数随分子大小单调上升，故 Δ = 候选 − 共晶既衡量结合也衡量大小。按每重原子计数归一化（配体效率，附表 S9c）后，8 个体系中有 5 个反转排序——PKM2 C-pocket、EP300、SIRT3、LDHA 与 PDK2 丙酮酸位点。两项指标各有已知偏向（原始亲和力偏袒大分子、配体效率偏袒小分子），且因每体系仅以 exhaustiveness 8、单一种子运行一次，无不确定度估计（Vina 典型打分误差 2–3 kcal/mol，足以涵盖若干上述 Δ 值）。因此"无候选优于其晶体参照"只是本协议、本精度下的陈述，不应读作关于分子结合优劣的证据。

---

## 4. 讨论

### 4.1 数据支持什么

本研究核心且可重现的观察是：赖氨酸丙酮酰化调控网络的转录组足迹在糖尿病足溃疡组织中发生偏移，且偏移是有方向的、并非均匀——糖酵解酶丰度升高，而 IFN-I 读出仅部分下降。两条跨平台证据支持糖酵解分支：两队列 bulk meta 分析（26/26 方向一致）与角质细胞单细胞 pseudobulk。第三条是同源性观察——线粒体丙酮酸转运体 MPC2 的下降，把胞质丙酮酸蓄积与 Kpy 写入底物池连接起来，但它反映的是同一糖酵解轴的信号，而非独立的证据线。免疫信号仅在一个 bulk 队列指向同一方向：GSE134431 把细胞毒性浸润减少定位在与 IFN-I 读出下降相同的方向，而后者恰是文献[2]中 Kpy 所影响的生物学轴；GSE199939 未重复这一结果，故该对齐仅建立在单队列之上。

我们的数据与文献[2]提出的模型一致——糖酵解通量升高提高丙酮酸可利用性，进而驱动 Kpy 机器、削弱 IFN-I 信号与抗菌免疫——但同样与一个更简单的模型一致：创面缺氧与高糖各自独立提升糖酵解，并各自独立抑制 IFN-I 应答，而完全无 Kpy 参与。区分二者需要直接测量丙酮酰化，而我们没有做：既无抗体法也无质谱，本分析涉及的完全是网络的**转录**足迹，而非修饰本身。读者应将"KpyRN 失调"理解为"与丙酮酸可利用性改变及选择性 IFN-I 缺陷相一致的转录组证据"。

### 4.2 三个验证队列为何阴性，及其含义

四个验证队列中两个在模块层无阳性信号（GSE223964 供体层、GSE166120 受试体层），第三个（GSE230426）在相反区室给出模块层而非基因层信号。这些尝试并非"未能复制"的失败，而是各自收窄了可主张的范围。

**生物学定位。** KpyRN 糖酵解枢纽信号见于**溃疡组织**与**溃疡组织的角质细胞**，在感染队列与单细胞队列的 PBMC 中均不存在。这一组织限制性是实质且可操作的结论：它排除了循环替代指标，也是预设验证方案采用组织而非血液的原因。

**统计效力层级。** 信号在何层级可检出本身即含信息：bulk meta 分析（每队列 12,000–17,000 基因量级、21 与 21 例样本）可解析 KpyRN 模块；8 位供体的单细胞队列与 4 位受试体的空间队列不能，因由 3–5 位供体计算的模块评分对 bulk 所观测到的效应量几乎无统计效力。因此这两个队列的结果应读作"相对于 bulk 观察到的效应量而言效力不足"，而这些问题的细胞类型特异与区域特异版本所需样本量，远超当前任何公开数据集。

**失败的结构性原因。** 每个阴性结果都有可辨认的成因。GeoMx 实验仅能覆盖评分集的 42.2% 与 6 个枢纽基因中的 3 个，故未被覆盖的三分之二网络中即便存在很大真实效应也将不可见。单细胞数据集缺 raw counts，限制了 pseudobulk 保真度。感染队列缺 donor 标识与愈合结局元数据，被迫采用非配对主模型并放弃了结局对比。脱离其结构性原因而把这些结果仅称为"阴性"会失真；不加限定地报告为阴性则会低估证据。

### 4.3 IFN-I 异质性是一项发现而非麻烦

bulk meta 分析与单细胞 pseudobulk 均显示 ISG 图景混杂：IFI44L、IFI27、DHX58、DDX58、STAT2 下降，而 ISG20、GBP1 升高；单细胞层面 IFITM3、IFITM1、ISG20、IFI16、IFIT2 升高而 IRF9、JAK1、IFNAR1、EIF2AK2 下降。统一的 IFN-I 阻断无法产生这种模式。最经济的解释是干扰素刺激转录组**组成**的改变——例如从依赖 IRF9/STAT2 复合物的 ISG 转向由炎症性 NF-κB 或 IFN-γ 驱动者——而非干扰素信号总强度的改变。这有实际含义：单个 ISG 将是衡量 Kpy 轴的糟糕读出；最有信息量的组织水平测量是"经典 I 型 ISG vs 炎症性 ISG"的配对对比，而非任选一个 ISG。

免疫表型并非两个队列都讲述同一故事。GSE134431 显示溃疡组织中细胞毒性 T 细胞与 NK 浸润减少；GSE199939 中任何细胞毒性指标均无显著变化，单核系评分反而升高。GSE134431 中组织单核浸润减少、单细胞供体层面单核**比例**增加，二者最可能的调和是浸润物**组成**改变（单位单核细胞所含淋巴细胞减少，或向髓系主导的浸润谱偏移），而非单核细胞绝对丢失。我们如实并列，不作单向裁决。

### 4.4 机器学习结果须保守解读

六基因共识模型在外部队列的 AUC 1.000 表面可观，但不应被读作临床性能估计：它来自 21 例样本（10 DFU vs 11 对照），因模型对本样本完全可分，置信区间塌缩至边界，且同一模型在前臂皮肤上得分 0.05。XGBoost 的失败更具启发性而非尴尬：被单基因主导（FAHD2B，importance 0.97）且外部恰为随机水平（AUC 0.500）的模型，展示了看似自信的高内样本交叉验证 AUC（0.911±0.067）如何与零外部信号并存。若我们报告最优单算法，就会报告一个带自信置信区间的随机性能。这正是四算法共识、10 种子重复、外部队列与前臂皮肤失败的诚实报告必须并列呈现的原因。

签名由 ADH1B（组织中 −4.42 log2FC）、ALDH3A2、FAHD2B、DHX58 与 GBP1 连同 LDHA 组成一事亦值得注意：仅 LDHA 与 GBP1 是 M2/M4 的典型成员，ADH1B、ALDH3A2、FAHD2B 三个成员为代谢基因，其组织效应量大但跨队列高度异质（三者 I² 均 >87%），意味着驱动原始 RT-qPCR 样本量计算的合并效应量对这些基因并不稳定——且这三者的引物已被 BLAST 最终重跑淘汰（见 §2.12）。它们在组织中复现失败不会推翻核心糖酵解结论，而复现成功则属强支持性证据。

### 4.5 因果推断：一个信号，恰当的保留

DHX58→2 型糖尿病结果（cis-eQTL 驱动的 DHX58 表达每升高一个单位，OR 0.931）方向耐人寻味：KpyRN 在 DFU 组织中衰减的 IFN-I 读出包含 DHX58 降低，而 DHX58 表达与 T2D 风险呈负相关，即表达越高越具保护性。一个"DHX58 降低标志易感感染的组织状态、其表达与代谢病存在因果关系"的模型，与整体框架自洽。但 34 次检验中仅 1 对 FDR 显著、3 个名义显著的代谢配对全部未通过校正，证据偏弱。我们亦注意到 DHX58 与 DDX58 系 RIG-I 的两个别名：在感染队列的 ENSG→symbol 映射中，DDX58 曾"缺失"仅因映射记录的是别名 RIGI，按别名规范化后回收——这是提醒：公共数据中看似基因层面的缺失常常只是命名学伪影。

### 4.6 对接：它确立了什么，未确立什么

全部 8 个体系以 <0.001 Å 通过阳性对照重对接，确立了协议能重现晶体构象；并未确立任何候选化合物能够结合。在每一体系中，全部候选的得分都弱于其共晶参照，这是针对深而明确的口袋筛选药样小分子的预期且诚实的结局；实际结论不是"无化合物可结合"，而是"这些候选是需生物物理验证的先导化合物"。唯一表面上更优的结果（PDK2 丙酮酸位点，Δ=−0.73）因另一原因更具启发性：在该浅口袋中蔗糖负面对照与共晶不可区分（−4.01 vs −4.08），显示分辨力不佳的体系中的数值"胜出"并不携带信息。若不附蔗糖对照而报告该数字，将是本研究最可能误导读者的做法。

### 4.7 局限

除各队列已述局限外，以下约束适用于全部分析。

**测量内容。** Kpy 本身未被测量，全部发现均属网络的转录足迹层面。组织 RT-qPCR 验证方案已预设但**未执行且不计划执行**，仅以设计文件呈现，故本文主张相应地局限于公共数据推断。相关地，愈合结局这一分支仍是**未答而非不可答**：唯一提供 DFU 组织愈合/未愈合分层的公开数据集是 GSE166120（2 名愈合者 vs 2 名未愈合者，12 周结局），本稿对其作过分析且在受试体层面为阴性；两个 bulk 队列的愈合/未愈合对比零 DEG，其 meta 分析仅 12 个 meta-DEG。在 n=4 受试体、面板覆盖评分集 42.2% 的条件下，该数据集无法了结结局问题；关闭此分支需要具 12 周随访的前瞻性组织采样，这正是设计文件（补充材料 S4b）所规定而本工作未尝试的。

**数据与样本量。** 两个 bulk 队列（GSE80178、GSE68183/86）因 pdInfo 包不可用而无法在样本级处理、未纳入计算，削弱了 meta 分析的实际广度。发现队列样本量小（各 21 例），n=21 且软阈值拟合不稳的 WGCNA 处于有效性下限，模块结果仅用作过滤。感染队列缺显式 donor 标识与愈合结局元数据。GeoMx 覆盖不完整（评分集的 42.2%）且无上游归一化指标。单细胞注释基于标记评分，12.7% 为 Ambiguous 细胞，若干亚谱系仅 1–2 个标记，故亚谱系结论为探索性。

**推断层级。** GeoMx 比较为嵌套结构（4 位受试体下 23 个 ROI），单细胞比较亦为嵌套（8 位供体下若干细胞）；我们全程强制以供体/受试体层面为推断主口径，并将那些在正确层级重算后消失的 ROI 级名义结果归类为假阳性而非发现。

### 4.8 何种实验能推进本研究

四项实验将实质改变本工作的地位。（i）以质谱或经验证抗体直接定量 DFU 与非溃疡组织的 Kpy，最好配以方案中已指定的配对丙酮酸测定，使修饰本身（而非仅其转录）变得可检验。（ii）执行预设的组织 RT-qPCR，并纳入"经典 I 型 ISG vs 代谢基因"的配对对比，以分别检验失调的两个分支而非合并检验。（iii）更大规模、供体层面且面板全覆盖的 DFU 单细胞或空间队列，以检验我们在 n=5 vs 3 条件下观察到的细胞类型限制性 LDHA 模式。（iv）针对 Kpy 靶基因启动子（FGFR4、LGALS8、RSPH4A）的亚硫酸氢盐或 CUT&Tag 分析，直接检验机制的表观遗传分支。

---

## 5. 结论

赖氨酸丙酮酰化调控网络在糖尿病足溃疡组织的转录组中发生改变：糖酵解酶表达可重现地升高，I 型干扰素读出**选择性而非均匀**地减弱，且细胞毒性浸润减少仅见于两个 bulk 队列之一。KpyRN 作为一个基因集合在 meta-DEG 中未富集（26/97 vs 24.6% 背景率；Fisher OR 1.12，p=0.64），可重现的是特定基因——糖酵解核心与部分 ISG——的行为，而非整个网络。证据在 bulk 组织 meta 分析与独立 PBMC 感染队列的模块层最强；另两个公共队列在供体与受试体层面不可评估，且不延伸至循环血液。一个因果估计通过校正（DHX58→2 型糖尿病，OR 0.931，95%CI 0.896–0.968），但其工具变量为单一 cis-eQTL，假定无法检验。对接覆盖 5 个靶蛋白 8 个口袋，协议阳性对照为刚性、原始亲和力无候选超过共晶——而按分子大小归一化后 8 个中有 5 个反转，故该排序本身随协议而定。Kpy 本身未被测量，组织 RT-qPCR 方案已预设但未执行且不计划执行。我们将 KpyRN 定位为 DFU 微环境候选代谢状态读数，并界定其证据边界，而非已确立的机制。

---

## 声明

**伦理审批与知情同意**：不适用。本研究全部分析使用公开、去标识化的二次数据，未招募新受试者，未接触可识别个人隐私信息。预设湿实验验证方案（补充材料 S4）未实施，其执行所需的伦理审批与知情同意要求由该文件说明；本工作不执行该方案。

**发表同意**：不适用。

**数据与材料可用性**：所有原始数据均来自 NCBI Gene Expression Omnibus，编号见附表 S2；所用补充文件的精确路径、大小与 SHA-256 校验和记录于项目 MANIFEST，并以版本标签归档于下述公共仓库。所有衍生结果文件（meta 分析表、机器学习性能表、免疫解卷积输出、单细胞与空间输出、MR 表、对接表）连同其生成脚本一并存放于同一仓库；本稿每一个数字均登记于机器可读审计台账，其中记录源文件与提取表达式。Zenodo DOI 将在归档存入时分配。

**利益冲突**：作者声明无利益冲突。

**基金资助**：本研究未获得任何资助机构（公共、商业或非营利部门）的专项经费支持。

**作者贡献**：YY 构思本研究、设计分析方案、执行全部分计算分析、对输出进行独立重算验证并撰写稿件。YY 对本文内容负全部责任。

**致谢**：作者感谢本文所用公共 GEO 数据集的研究者与受试者。

**AI 辅助技术声明**：撰写过程中使用了基于大语言模型的助手（WorkBuddy/CodeBuddy）提供有限且有人监督的辅助，包括（i）低风险的语体与结构润色，以及（ii）较高风险的辅助——依据分析台账核对数值一致性、降 AI 痕迹处理，以及覆盖统计、报告规范、编辑与语言等维度的结构化多轮评审。未使用任何 AI 工具生成、编造或修改数据，亦未在未经独立核对的情况下采纳任何 AI 输出。全部分析、统计决策、解释及科学内容均由作者完成并验证，作者保留全部责任；未将任何作者身份或责任归于 AI 系统。未向任何 AI 系统共享机密或可识别数据——所有输入均为公开、去标识化的二次数据。

---

## 图表说明

**图 1** 研究流程与分析管线。检索到 10 个公共人类 DFU 相关数据集；其中 2 个 bulk 系列（GSE80178、GSE68183/86）不能在样本层面处理，故 8 个进入五模块、128 基因的 KpyRN 框架，并沿五条并行路径分析：(i) 逐队列 limma + 随机效应 meta 分析；(ii) WGCNA 模块–性状关联；(iii) 四算法共识机器学习 + 外部验证；(iv) 三算法免疫解卷积 + 单细胞与 GeoMx DSP 验证；(v) cis-eQTL 孟德尔随机化与 AutoDock Vina 对接（含阳性与负面对照）。实线箭头示数据流向，虚线箭头示独立验证路径。

**图 2** KpyRN 证据层的概要。左：三层证据韦恩图——meta 分析显著 KpyRN 基因（26）、四算法 ML 共识枢纽基因（6）、WGCNA 5,000 基因背景中出现的候选池基因（18），及两两交集（meta∩ML=6、ML∩WGCNA=3、meta∩WGCNA=18）与三者交集（3）。右：共识枢纽基因的合并 META1 效应量，FDR<0.05 以星号标注，Task A/Task B 归属以颜色区分。

**图 3** 发现队列的 WGCNA 模块–性状相关。GSE134431（n=21，软阈值 power 6，9 个非 grey 模块）中模块特征基因与 DFU 状态、愈合结局的相关热图，与 DFU 状态负相关者示在溃疡组织中下调；4 个对 DFU 显著的模块（magenta、yellow、red、blue）标注其相关系数与 p 值。

**图 4** 计分 KpyRN 基因的 meta 效应量。META1（DFU vs 非溃疡皮肤，2 队列，REML 随机效应）中 26 个显著 KpyRN 基因的合并 log2FC 森林式图，按效应量排序并标注模块，附队列特异估计以评估方向一致性。

**图 5** 共识机器学习性能与外部验证。(A) 发现集重复 5 折交叉验证 AUC（10 种子，均值±SD）；(B) GSE199939 中各算法与六基因共识逻辑回归模型的外部 ROC；(C) 两个 bulk 队列中共识基因在 DFU 与对照的表达箱线图；(D) 外部队列 AUC 与 DeLong 置信区间，含前臂皮肤 GSE143735 的迁移失败。

**图 6** 免疫浸润与枢纽基因–免疫关联。(A) GSE134431 的 DFS vs DFU 与 GSE199939 的 DFU vs 对照中显著免疫估计量的箱线图；(B) 枢纽基因与免疫指标的 Spearman 相关热图，FDR 显著对以标记；(C) EPIC 与 MCP-counter 在 B 细胞、CD8 T 细胞、单核细胞、中性粒细胞与 NK 细胞上的跨方法一致性散点。

**图 7** GSE165816 中 KpyRN 信号的单细胞定位。(A) 足皮肤供体层面 DFU vs 非 DFU 的细胞类型组成，4 个 FDR 显著谱系高亮；(B) 按谱系与条件分层的糖酵解与 IFN-I 读出基因 pseudobulk log2 CPM 热图；(C) 角质细胞糖酵解酶按条件的表达；(D) 伪重复对照示意：全部单细胞检验在供体层面执行。

**图 8** 三个验证队列的阴性与边界发现。(A) GSE223964：供体层面按细胞类型与条件的模块评分，并列置换检验 p 值与 Mann–Whitney U 结果（全部阴性）；(B) GSE166120：CTA 板对 KpyRN 计分集的按模块覆盖，标注枢纽基因覆盖；(C) GSE166120：ROI 级与受试体级 p 值对照，演示伪重复驱动的显著性消失；(D) GSE230426：三模型的 PBMC camera 模块评分，显示 M3/M4 一致下降与枢纽基因阴性结果。

**图 9** 孟德尔随机化森林图。按基因给出对 2 型糖尿病与空腹胰岛素的 Wald 比值估计及 95% 置信区间；DHX58→T2D 作为唯一 FDR 显著配对高亮，并与 XGBoost 主导的单基因模型对照。

**图 10** 分子对接结果。(A) 8 个受体–配体体系的阳性对照重对接 RMSD，标注 2 Å 验收阈值；(B) 各体系亲和力总览：候选亲和力、晶体共晶参照与蔗糖负面对照；(C) 四个最优候选体系的受体坐标系内代表性姿态投影。

---

## 补充材料

**附表 S1** KpyRN 基因集（模块、亚模块、角色、来源、证据说明、计分标志；140 条目，128 计分）。
**附表 S2** 全部 26 个 GEO 系列的样本级元数据与分组变量。
**附表 S3** 随机效应 meta 分析结果（META1 与 META2），含合并 logFC、SE、z、p、FDR、I²、队列特异 logFC 与 SE 及两队列方向一致性。
**附表 S4** 机器学习性能：逐算法交叉验证 AUC（10 种子均值±SD）、外部队列 AUC 与 DeLong 置信区间、基因选择频次与 SKIPPED 记录日志。
**附表 S5** 免疫浸润统计（三种解卷积方法、全部对比、FDR 校正）、枢纽基因×免疫相关矩阵与跨方法一致性表。
**附表 S6** 单细胞（GSE165816、GSE223964）与 GeoMx DSP（GSE166120）结果表：供体/受试体层面组成统计、pseudobulk 模块评分、ROI 元数据、板覆盖与伪重复对照。
**附表 S7** 感染队列（GSE230426）结果表：三模型全基因组 DEG 计数、camera 模块结果、枢纽基因结果与模型一致性统计。
**附表 S8** 孟德尔随机化表：工具变量与 F 统计、等位基因协调审计（含全部排除记录）、逐配对结果与 FDR 校正。
**附表 S9** 分子对接表：受体与位点定义、对接盒几何、参数、全部 55 次对接的亲和力与相对共晶 Δ、阳性对照 RMSD 验证表。
**附表 S10** RT-qPCR 设计与验证：引物序列、Tm、GC 与产物长度；二级结构与异源二聚体检查；BLAST 特异性结果；样本量计算（配对、独立与保守情景）。
**补充材料 S4** 预设组织 RT-qPCR 验证方案：科学目的、伦理与知情同意路径、样本量推导、组织采集与丙酮酸酶法、引物表、统计分析计划与预设阳性判定标准、与已完成公共数据的关系、溯源清单。**该方案未实施**；撰稿时不存在组织标本、伦理批准、qPCR 或丙酮酸数据。

---

## 参考文献

1. Chinese Burn Association, the Yangtze River Delta Integrated Diabetic Foot Alliance, and the Editorial Committee of the Chinese Journal of Burns and Wound Repair. Practical guidelines for the prevention and management of diabetic foot disease in China. Burns Trauma. 2025;13:tkaf064. doi:10.1093/burnst/tkaf064
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
35. Vandesompele J, De Preter K, Pattyn F, et al. Accurate normalization of real-time quantitative RT-PCR data by geometric averaging of multiple internal control genes. Genome Biol. 2002;3(7):research0034. doi:10.1186/gb-2002-3-7-research0034
36. GEO GSE199939. Gene Expression Omnibus, NCBI. Accessed 2026-10-25.
37. GEO GSE114248. Gene Expression Omnibus, NCBI. Accessed 2026-10-25.
38. GEO GSE146028. Gene Expression Omnibus, NCBI. Accessed 2026-10-25.
39. GEO GSE231643. Gene Expression Omnibus, NCBI. Accessed 2026-10-25.
40. GEO GSE245703. Gene Expression Omnibus, NCBI. Accessed 2026-10-25.
