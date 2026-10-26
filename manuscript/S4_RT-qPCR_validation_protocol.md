# 补充材料 S4：DFU 组织 RT-qPCR 验证实验方案（预设方案 / pre-specified protocol）

> **文件性质**：本文件是**预设实验方案（pre-specified protocol）**，在实验实施前完成撰写与冻结。
> 截至本稿撰写时，本研究**尚未采集任何湿实验样本**（伦理审批未提交/未获批，无组织标本，
> 无 qPCR 数据，无丙酮酸酶法数据）。因此本节内容为**方法学设计**，**不含任何实验结果**，
> 稿件正文中对应位置的表述必须明确区分"已完成的公共数据分析"与"预设的湿实验验证方案"。
> 方案中所有数值（引物序列、样本量、效应量、Tm、GC、产物长度）均可回溯到本项目的产物文件。

> **★ 状态更新（2026-10-25）：本方案不执行。** 项目定位为**纯公共数据再分析**，不开展
> 湿实验：无伦理申请、无组织采集、无 qPCR、无随访，稿件中亦不作任何实验验证声明。
> 本 S4 与 `S4b_湿实验最小集SOP.md` 均为**未执行的设计留档**。
>
> 留档价值有二：① 因 §3.2 的独立性检验把叙事改为"两条部分独立轴"，且 BLAST v6 重跑
> （2026-10-25）判定 12 对引物中 **6 对 OFFTARGET_P0**（FAHD2B、ADH1B、ALDH3A2、
> ACTB、GAPDH、SDHA——败于 class I ADH / 肌动蛋白旁系家族、未定性 LOC 位点的
> （近）完美匹配、及 ALDH3A2 反向引物的 210 条 95% 命中），最小可行集在 S4b 中收缩为
> **LDHA（糖酵解轴）+ DHX58（IFN-I 轴）+ 2 个事后选定内参 = 4 对引物**；
> ② 明确公共数据无法闭合的缺口——S4b §3.2 指明"愈合者 vs 未愈合者"分层在现有 GEO 中不存在，
> 只能前瞻性采集，而本文不采集。
>
> 两份文件不冲突：S4 记录完整设计（6 基因 / 12 对引物及其 BLAST 判定，其中 6 对已判
> OFFTARGET_P0 不得进湿实验），S4b 记录最小可行子集与预注册判据。**任何"将执行/已排期"的表述均属错误。**

---

## S4.1 科学目的

在组织水平验证本研究公共数据分析所提出的核心结论：**丙酮酰化相关网络（KpyRN）在
糖尿病足溃疡（DFU）创面组织中失调，且该失调与丙酮酸可利用性升高方向一致**。

具体检验三个方向性预设（依据 `03_结果/meta/META1_DFU_vs_control.csv` 的 pooled 效应量）：

| 目标 | 方向 | META1 pooled_logFC | I² | 队列方向一致 |
|---|---|---|---|---|
| LDHA（糖酵解/丙酮酸供应） | ↑ DFU | +0.9987 | 0.0 | 2/2 一致 |
| GBP1（IFN-I 读出） | ↑ DFU | +0.8406 | 0.0 | 2/2 一致 |
| DHX58（RIG-I，IFN-I 读出） | ↓ DFU | −0.7372 | 57.2 | 2/2 一致 |
| FAHD2B（丙酮酸/糖代谢） | ↓ DFU | −1.1177 | 87.1 | 2/2 一致 |
| ALDH3A2（乙醇/乙醛代谢） | ↓ DFU | −1.6264 | 91.2 | 2/2 一致 |
| ADH1B（乙醇代谢） | ↓ DFU | −4.4223 | 90.2 | 2/2 一致 |

**预设说明（诚实边界）**：
1. 上表 6 个基因中 **LDHA 与 GBP1 的 I²=0（两队列完全一致）**，是**最可靠**的验证靶点；
   ADH1B / ALDH3A2 / FAHD2B 的 **I²>87%（队列间高度异质）**，其 pooled 效应量不稳定，
   湿实验若复现失败**不视为本研究核心结论被推翻**，但**若复现成功则是强支持**。
2. META1 效应量来自**两个不同平台/数据批次**的 bulk RNA-seq（GSE134431 RPKM、
   GSE199939 TPM），**与 RT-qPCR 绝对表达量不可直接类比**；本节只检验**方向**，
   不预设倍数变化与 RNA-seq 一致。
3. 本研究**未检测任何 Kpy 修饰本身**（无 Kpy 抗体、无质谱），因此**不能**声称验证了
   "Kpy 修饰水平改变"。可验证的是**其上游底物供应（丙酮酸）与下游转录读出（IFN-I）**。
   这一点必须在正文与本节中明确。

---

## S4.2 伦理与知情同意（先行，未获批不得实施任何样本采集）

| 项目 | 要求 |
|---|---|
| 伦理审查 | 本单位（福建中医药大学附属第二医院）医学伦理委员会；或具备资质的合作单位 IRB |
| 提交时点 | **先于任何样本采集与实验操作**；方案冻结后提交，随机对照/单臂观察性设计 |
| 研究类型 | 观察性研究（组织/血液样本采集 + 分子检测），不涉及治疗干预 |
| 知情同意 | 书面知情同意；由具备 DFU 诊疗资质的医师在床旁获取，避免影响常规诊疗 |
| 样本来源 | DFU 创面组织（术中/换药时获得的**废弃**组织，不额外增加创伤）+ 同期外周血 + 同一患者配对非溃疡皮肤（内对照） |
| 匿名化 | 样本编码由独立第三方（不接触数据者）完成，实验室仅接收编码，身份对应表加密保存 |
| 退出与不良事件 | 受试者可随时退出，不影响常规医疗；组织采集不额外增加创口 |
| 数据留存 | 编码数据保存 ≥5 年；身份对应表在研究结束后 5 年销毁 |

**特别声明（与既往项目一致的红线）**：本研究若**无法**获得伦理批准或组织样本，
则**不实施**湿实验，稿件改为完全依赖公共数据的表述，并在 Limitations 中如实声明
"缺乏组织水平湿实验验证"。**不得**以任何形式把预设方案表述为已完成的验证。

---

## S4.3 样本量论证（基于本项目 meta 分析的真实效应量，非假设值）

效应量来源：`03_结果/meta/META1_DFU_vs_control.csv` 的 `pooled_logFC` 与 `pooled_se`，
Cohen's d = pooled_logFC / pooled_se；配对检验按配对标准化效应量近似。

**配对设计（DFU 溃疡组织 vs 同一患者配对非溃疡皮肤），α=0.05 双侧，power=0.80：**

| 基因 | Cohen's d | 所需配对数 | 保守情景（效应量×0.5） |
|---|---|---|---|
| LDHA | 6.70 | 5 | 5 |
| GBP1 | 3.80 | 5 | 7 |
| DHX58 | 4.51 | 5 | 6 |
| FAHD2B | 4.39 | 5 | 6 |
| ALDH3A2 | 2.63 | 6 | 9 |
| ADH1B | 2.53 | 6 | 9 |

**独立组设计（DFU vs 健康对照），α=0.05 双侧，power=0.80：**

| 基因 | Cohen's d | 所需 n/组 | 保守情景 |
|---|---|---|---|
| LDHA | 6.70 | 1 | 2 |
| GBP1 | 3.80 | 2 | 5 |
| DHX58 | 4.51 | 1 | 4 |
| FAHD2B | 4.39 | 1 | 4 |
| ALDH3A2 | 2.63 | 3 | 10 |
| ADH1B | 2.53 | 3 | 10 |

**最终目标样本量（综合计划书 n≥8/组 要求与功效计算、并预留脱落/技术失败）**：

- **DFU 创面组织组 n = 10**（≥8/组要求；对全部 6 基因的保守情景均充足——原最严者为
  ALDH3A2/ADH1B（配对 n=9、独立 n/组=10），该两基因已因 BLAST v6 OFFTARGET_P0 退出；
  仅计合格靶基因 LDHA/GBP1/DHX58 时保守情景配对 n ≤ 7，n=10 依然充足，无需下调）
- **配对非溃疡皮肤组 n = 10**（与 DFU 组同患者配对；配对设计是主分析）
- **健康对照皮肤组 n = 8**（≥8/组要求；独立组分析为次要分析）
- **外周血（配对）n = 10**（与组织组同患者；用于丙酮酸负荷与基因表达的循环层面参照）
- **招募量 = 分析目标 ×1.2（预留 20% 取材失败/RNA 质量不合格/失访）**：
  DFU 组织 **12**、配对非溃疡皮肤 **12**、健康对照皮肤 **10**、配对外周血 **12**
- **组织样本合计 34 例**（12 + 12 + 10），血样 12 例

> 数值来源：`03_结果/实验验证设计/sample_size_calculation.csv` 与 `.json`
> （由 `02_脚本/compute_qpcr_sample_size.py` 依 `META1_DFU_vs_control.csv` 的真实
> pooled 效应量计算生成，非人工设定）

**统计功效的诚实说明**：以上功效计算基于 **RNA-seq pooled 效应量**，而 RT-qPCR 的
测量误差通常**小于** bulk RNA-seq，效应量可能更大；但 RNA-seq 效应量本身受
批次/平台/细胞组成混杂影响。**因此 n=10 是在"复现方向"这一目标下的合理规模，
若目标升级为"精确估计效应量"则需更大样本**——本方案不声称后者。

---

## S4.4 样本采集与处理

| 环节 | 标准操作 |
|---|---|
| 组织 | 溃疡中心（避开坏死与感染灶外缘，取 2 份 3 mm 皮活检）与**同一患者**足部非溃疡皮肤（距溃疡 ≥5 cm，取 1 份）；离体后 30 min 内液氮速冻，−80 °C 保存 |
| 血液 | EDTA 抗凝管，采血后 2 h 内 300×g 离心 10 min 分离血浆，−80 °C 保存（用于丙酮酸测定） |
| RNA 提取 | TRIzol 提取 → RQI ≥ 7.0（Bioanalyzer）；DNase I 去除基因组 DNA；**残余基因组 DNA 用无 RT 对照验证** |
| cDNA 合成 | 500 ng 总 RNA，RevertAid 第一链合成；**每批样本同时做 −RT 对照**（无逆转录酶） |
| 丙酮酸测定 | 组织匀浆液与血浆各取，用**丙酮酸酶法（pyruvate kinase/lactate dehydrogenase 偶联，340 nm 吸光度）**测定；组织值以总蛋白归一（µmol/g protein），血样以 mmol/L 报告 |

**丙酮酸作为替代指标的科学定位（诚实边界）**：丙酮酸是 Kpy 的**直接前体**，
但丙酮酸水平反映的是**瞬时代谢状态**而非 Kpy 修饰稳态；且组织匀浆中的丙酮酸
极易在取材后被内源酶消耗。**因此丙酮酸测定仅作为"Kpy 前体可利用性"的辅助指标，
不能替代 Kpy 修饰的直接检测**。若条件允许，建议增加 **SIRT3/PDK2 蛋白水平
（Western blot）** 作为 Kpy 机器组分的功能性旁证——这仍不等同于 Kpy 修饰的直接测量。

---

## S4.5 RT-qPCR 引物（Primer3 primer3-py 2.4.0 设计，模板为 NCBI RefSeq 规范 NM_ 记录）

引物设计参数：产物 80–180 bp；Tm 58–62 °C（最适 60 °C）；GC 40–60%；
最大 Poly-X = 4；拒绝任何含 N 的引物；自互补/发卡/异源二聚体 Tm < 47 °C。

### 目标基因（6 个 ML 共识枢纽基因）

**BLAST v6 重跑后（2026-10-25）：6 个靶基因中 3 个通过（LDHA / DHX58 / GBP1），3 个判 OFFTARGET_P0 不得进湿实验。**

| 基因 | 模块 | RefSeq 模板 | 正向引物 (5'→3') | 反向引物 (5'→3') | Tm F/R (°C) | GC F/R (%) | 产物 (bp) | BLAST v6 判定 |
|---|---|---|---|---|---|---|---|---|
| **LDHA** ‡ | M2 丙酮酸供应 | NM_005566.4 | ACCACTTCCATTGTAAGTCCCA | CCAGCTTCCAGAGGACAAGAT | 59.3/59.4 | 45.5/52.4 | 174 | **USABLE**（双端 P0=0，LOC105375980 命中=0） |
| **DHX58** | M4 IFN-I 读出 | NM_024119.3 | TGACCTGCTCATCTGCACAG | TGTAGACGGTGTCCTTGTGC | 60.0/60.0 | 55/55 | 137 | **USABLE**（双端 P0=0，tgtTx 8/8） |
| **GBP1** | M4 IFN-I 读出 | NM_002053.3 | AAGGGAACAGCCTGGACATG | ATTGCCACCACCACCATAGG | 60.0/60.0 | 55/55 | 145 | **USABLE**（双端 P0=0，tgtTx 1/1） |
| ~~FAHD2B~~ | M2 丙酮酸供应 | NM_001320848.2 | CAATGGGAAACAGTGGCTGC | CATCTGGTTGGTGTTGCTGC | 60.0/60.0 | 55/55 | 157 | OFFTARGET_P0（F→LOC124900512，100% identity、0 错配） |
| ~~ADH1B~~ | M2 丙酮酸供应 | NM_000668.5 | TGGCCAAAATTGATGCAGCC | AACACAGCACAGGTAGAGCC | 60.0/60.0 | 50/55 | 124 | OFFTARGET_P0（F→ADH1C、R→ADH1A，均 100%——class I ADH 旁系家族） |
| ~~ALDH3A2~~ | M2 丙酮酸供应 | NM_000382.3 | TGTGTTCCCAGCTGCTTAGG | TTTCGGAGACAGGGTCTTGC | 60.0/60.0 | 55/55 | 129 | OFFTARGET_P0（R→DYM/A1CF/SUN2/TMEM67 等 **210 条** 95%、3' 端完美） |

> ‡ LDHA 改为 **3'UTR 锚定**（偏离其余基因的 CDS 锚定范式），理由与偏差声明见下方第 5 条迭代。
> 该行 Tm 为 primer3-py 2.4.0 在 NM_005566.4 上的**实测复算值**（59.29/59.44 °C，ΔTm 0.15 °C），
> 其余各行保留设计目标值 60.0 °C。该对引物落在 LDHA 3'UTR 的多转录本共有区段（与 6 条
> LDHA RefSeq 转录本全长匹配），因此报告的是 **LDHA 转录本总体丰度**而非单一异构体；
> 其相对优势在于**不与 LOC105375980 共扩增**（BLAST 命中 = 0）。

### 内参基因候选（最终选择以实测稳定性为准）

**BLAST v6 重跑后（2026-10-25），6 个内参候选中仅 3 个通过（TBP / HPRT1 / RPLP0）；ACTB、GAPDH、SDHA 均判 OFFTARGET_P0，退出候选。**

| 基因 | RefSeq 模板 | 正向引物 (5'→3') | 反向引物 (5'→3') | Tm F/R (°C) | GC F/R (%) | 产物 (bp) | BLAST v6 判定 |
|---|---|---|---|---|---|---|---|
| TBP | NM_003194.5 | TGCTGAGAAGAGTGTGCTGG | TGGCAGGCTGTTGTTCTGAT | 60.0/59.9 | 55/50 | 124 | **USABLE** |
| HPRT1 | NM_000194.3 | TTGCTTTCCTTGGTCAGGCA | ATCCAACACTTCGTGGGGTC | 60.1/60.0 | 50/55 | 85 | **USABLE** |
| RPLP0 | NM_001002.4 | AATCTCCAGGGGCACCATTG | TCGAACACCTGCTGGATGAC | 60.0/60.0 | 55/55 | 144 | **USABLE** |
| ~~ACTB~~ | NM_001101.3 | ATTGCCGACAGGATGCAGAA | CGGACTCGTCATACTCCTGC | 60.0/60.0 | 50/60 | 172 | OFFTARGET_P0（F→POTEM、ACTBL2，各 95%、3' 端完美） |
| ~~GAPDH~~ | NM_002046.6 | GATTTGGTCGTATTGGGCGC | TTCCCGTTCTCAGCCTTGAC | 60.0/60.0 | 55/55 | 169 | OFFTARGET_P0（R→LOC124901337，95%、3' 端完美） |
| ~~SDHA~~ | NM_004168.4 | GCTGTGGCCCTGAGAAAGAT | TCTCTGAAATGCCAGGCAGG | 60.0/60.0 | 55/55 | 90 | OFFTARGET_P0（F→LOC124904260，95%、3' 端完美） |

**内参选择规则（不预先指定唯一内参）**：对通过 BLAST 的 **3 个候选（TBP / HPRT1 / RPLP0）全部检测**，再用 **geNorm** 与 **NormFinder** 对本研究的组织与血样分别做稳定性排序，取**几何平均数（geometric mean）最靠前的 2 个**作为最终内参。理由：DFU 组织中炎症细胞浸润会显著改变白细胞相关基因（如 RPLP0）的表达，而 GBP1、DHX58 本身即为 ISG（会受 IFN-I 状态直接影响，**故不可作为内参**）——这一判断直接来自本研究的免疫浸润分析结果，是方案设计的实质性约束。**经典内参 ACTB/GAPDH 的落选（POTEM/ACTBL2/LOC124901337 的 3' 端完美匹配）本身即为例证：通用内参不能凭文献沿用，必须逐引物实测。**

**引物验证状态（提交前必须完成）**：
- 二级结构与自互补/异源二聚体：已由 primer3-py 独立复算，全部 < 47 °C（12/12 通过）
- GC 锁温与 Tm 匹配：全部合规（正反向 GC 差 ≤10%、Tm 差 ≤0.2 °C；LDHA 实测 GC 差 6.9%、ΔTm 0.15 °C）
- **引物在设计模板上的落位**：已由 `02_脚本/fix_refseq_provenance.py` 逐对 efetch
  设计 accession 并实测，**12/12 对正反向引物均精确命中模板**（产物长度 85–174 bp，
  坐标见 `provenance/primer_template_mapping.csv`；LDHA 因改为 3'UTR 锚定，产物 174 bp，
  为区间内最大值）
- **单靶特异性：BLAST（blastn vs `core_nt` + 人源物种限定 `txid9606[Organism:exp]`）
  验证结果见 `primer_blast_v6.csv`**（由 `02_脚本/regen_ldha_v5_blast.py`
  以 **v6 判据 + 重试/INCONCLUSIVE 修复**统一重跑全部 12 对生成，2026-10-25；
  分层明细 `primer_blast_v5_ldha_all.json`）——**重跑终判：6 对 USABLE、
  6 对 OFFTARGET_P0**（见上方两张表的判定列；`primer_blast_specificity_v3.csv`
  为 S10d 附表的同源视图）
  —— 判据：人源 RefSeq 转录本层面**无全长匹配的真脱靶**（全长 = `align_len` 等于引物
  长度；BLAST 的局部比对不计），且目标基因存在**全长匹配的人源 RefSeq 转录本**
  且 identity ≥ 90%；**LDHA 额外要求 LOC105375980 命中 = 0**（见下方第 5 条迭代）
  —— 若某对引物存在真脱靶，**须更换该对引物**（`primer_design_report.json` 中每个基因
  备选了 5 对候选引物，可直接替换）
- **熔解曲线（Tm）验证**：每对引物做单一熔解峰确认；出现双峰即判定为非特异
- **标准曲线**：5 点系列稀释，效率 90–110%、R² ≥ 0.99 方可纳入

> **方法学说明（迭代记录，提交时应保留以示透明）**：本项特异性验证经历多轮
> 版本迭代，多数修正的都是**会导致错误结论的工具层缺陷**，值得记录以备复核：
>
> 1. **v1（Text 输出）→ v2（XML 输出）**：v1 使用 `FORMAT_TYPE=Text`，实测 NCBI
>    返回的命中列表中 description 列被**截断**（形如 `... uncharacterized... 39.2 1.2
>    100%`），基因符号常落在截断点之后，导致正则提取不到符号而**被误判为"非特异"**
>    （FAHD2B 前向人源 29 条命中但 `genes=` 为空）。v2 改用 `FORMAT_TYPE=XML`
>    的完整 `<Hit_def>`，并修正 `Hsp_identity` 为**匹配碱基数**而非百分比
>    （identity% = identity / align-len）。
> 2. **v2（refseq_rna）→ v3（core_nt）**：v2 发现 LDHA 前向引物"人源 0 命中"，
>    初判为引物不可用；**独立复核发现是数据库覆盖缺陷**——`refseq_rna`
>    并不完整收录人源转录本（该库对人源仅收录部分记录，且优先呈现 PREDICTED:
>    的 XM_ 记录），而 `core_nt`（NCBI Core Nucleotide，完整收录人源 RefSeq/GenBank）
>    正确命中 `NG_008185`（LDHA RefSeqGene）。该引物在 NM_001135239.1 第 136–155 位
>    精确匹配、产物 120 bp，**完全可用**。v3 改用 core_nt，并把判据细化为
>    "人源 RefSeq 转录本层面无真脱靶"，同时将同源旁系同源（如 ACTB↔ACTG1）、
>    假基因与基因组非转录区命中**分层单列**而非笼统计为脱靶。
> 3. **v3（core_nt 无物种限定）→ v4（core_nt + 人源限定 + 全长匹配筛选）**：
>    v3 改用 core_nt 后，ACTB 反向引物仍报"人源 0 命中"。**独立复核发现是
>    命中列表截断**——`HITLIST_SIZE=500` 且排序按 bit-score 优先，**进化上更近的
>    非人源物种会挤占前 500 名**。v3 留档 XML 证据：第 3 个查询起命中数**全部
>    等于 500**（恰为上限，确证截断）；ACTB 反向查询中 `Pseudopipra pipra` 24 条、
>    `Neopelma chrysocephalum` 22 条、`Corapipo altera` 22 条占主导，人源仅剩
>    7 条；另一查询中 `Vibrio cholerae` 独占 250/500 条。修复：加
>    `ENTREZ_QUERY=txid9606[Organism:exp]` 物种限定。**加限定后同一引物立即得到
>    752 条人源命中、202 条 RefSeq 转录本、100% identity**（含 `NG_007992`
>    ACTB RefSeqGene），确证为查询层缺陷。
>    v4 进一步发现第二个判读缺陷：**BLAST 允许局部比对**，`align_len < 引物长度`
>    的部分匹配不构成可扩增脱靶却会被计入脱靶数。实测 LDHA 前向引物 2,261 条
>    人源命中中 **1,081 条仅 13 bp 匹配**、342 条 16 bp（多为基因组克隆记录）；
>    **全长 20 bp 匹配仅 27 条**，其中 100% identity 者**全部**为 LDHA 位点及其
>    邻近增强子（`NG_008185`/`NG_122038`/`NG_122039`/`AC084117`），**人源 RefSeq
>    转录本层面无 100% 全长非目标命中**；ZNF518B/PDS5B 的命中均为 **18/20 = 90%**
>    （含 2 个错配、反向链），不构成有效扩增。故 v4 判定层加**全长匹配**筛选，
>    并新增"目标基因须存在全长匹配的人源 RefSeq 转录本"这一必要条件
>    （否则该引物不可用）。
> 4. **溯源文件缺陷（与引物设计本身无关）**：复核时发现
>    `provenance/refseq_mrna_sequences.json` 对 LDHA 存的是 NM_001165414.2
>    （variant 3）、对 ACTB 存的是 **XM_054343706.1（POTEF）**——即 NCBI
>    Gene→nuccore elink 返回列表中被误取的记录。设计脚本实际使用的模板
>    （LDHA NM_001135239.1、ACTB NM_001101.3）自始正确，缺陷仅在**溯源文件**，
>    但会导致"引物是否落在模板上"的复核出现假阴性。已由
>    `fix_refseq_provenance.py` 按设计 accession 重取并落位实测修复。
>
> 5. **LDHA 改为 3'UTR 锚定（偏离 CDS 范式，诚实偏差）**：LOC105375980
>    （NR_202040.1）携带 LDHA CDS 的 96.16% 拷贝且带 polyA 尾，进入 cDNA 池后
>    任何 CDS 锚定候选都会出现 P0 完美匹配脱靶。故放弃 CDS 锚定，改为靶向规范
>    转录本 **NM_005566.4** 的 3'UTR 区段 [1667,2047]（GC 31.76%）。
>    最终选中候选 8：F=`ACCACTTCCATTGTAAGTCCCA`（复用已验证 P0=0 的正向）/
>    R=`CCAGCTTCCAGAGGACAAGAT`（21-mer，3' 端加 T 逃逸原 20-mer 与
>    PIGS/NM_033198 的 3' 端 5 nt 全匹配脱靶），产物 **174 bp**（落 80–180
>    窗口）、GC 差 **6.9%**（≤10）。特异性由 `02_脚本/verify_ldha_c8_finalize.py`
>    （仅对新反向做 BLAST，正向复用 prior 实测 P0=0/tgt_tx=6/LOC=0）与
>    `regen_ldha_v5_blast.py`（全 12 对统一复核）双重确认。**偏差范围（须精确表述）**：
>    该区段并非 NM_005566.4 独有——正反向引物各自与 **6 条 LDHA RefSeq 转录本**
>    全长匹配（`target_full_match_refseq_tx = 6/6`），独立复核也显示
>    NM_001135239.1（variant 2）上落位与产物长度同样为 F@1832 / 174 bp。
>    因此该对定量的是 **LDHA 转录本总体丰度**，而非单一异构体；
>    真正的偏离点是**不与 LOC105375980 共扩增**（BLAST 对该位点命中 = 0），
>    而非"只测单转录本"。已在 `qpcr_primers.csv` 的 `scope_note` 与稿件
>    §4.7 Limitations 如实标注。
>
> 6. **v6 判据脚本的假阴性缺陷与修复**：`verify_primers_blast_v6.py` 在
>    BLAST 状态≠READY 时返回 `hits=[]` 并被误读为零命中，导致 FAHD2B-F、
>    ADH1B-F、ALDH3A2-R、ACTB-R 等多对被**误报**为 NO_TARGET/tgtTx=0
>    （明显 implausible）。新写 `regen_ldha_v5_blast.py`：状态≠READY 重试至多
>    3 轮，仍失败显式标记 **INCONCLUSIVE**（绝不把"取不到结果"冒充"零命中"），
>    一次生成同源的 `primer_blast_v5.csv` / `primer_blast_v6.csv` /
>    `primer_blast_specificity_v3.csv` / `primer_blast_v5_ldha_all.json`，
>    供 `verify_primers.py`、`build_supplementary_tables.py`（S10d）与审计复用。
>
> 7. **2026-10-25 全量重跑与文件间矛盾修复**：12 对引物经修复后的
>    `regen_ldha_v5_blast.py` 完整重跑（约 2.5 h，逐基因 checkpoint 防半途丢失），
>    终判 **6 对 USABLE（LDHA、DHX58、GBP1、TBP、HPRT1、RPLP0）/ 6 对
>    OFFTARGET_P0（FAHD2B、ADH1B、ALDH3A2、ACTB、GAPDH、SDHA）**。
>    重跑中发现 `primer_blast_v6.json` 停留在 18:25 修复前判定（旧 LDHA 引物判
>    OFFTARGET_P0、GAPDH 误判 USABLE），与 22:51 的新 CSV **互相矛盾**——
>    已由 `sync_blast_v6_json.py` 从 `primer_blast_v5_ldha_all.json` 重建
>    （12 基因 CSV↔JSON 判定一致性自检通过），`regen_ldha_v5_blast.py` 亦改为
>    同时写两个 JSON 防止复发；引物判定列由 `join_blast_verdict_into_primers.py`
>    写回 `qpcr_primers.csv`（`blast_v6_verdict` / `blast_v6_P0_genes`）。
>
> **★ 字段污染与重建记录（2026-10-25）**：`primer_verification.csv` 的旧版
> `V1_specific` / `overall_pass` 曾依赖有缺陷的 v1 BLAST（description 截断 →
> 12/12 的 `V1_specific` 均为 `False` → `overall_pass` 恒为 `False`），且旧文件
> 停留在 LDHA 重设计之前的引物。已由 `verify_primers.py`（v2：读
> `primer_blast_v5.csv` 的 v6 同源重跑结果 + primer3-py 2.4.0 独立复算二级结构）
> **全量重建**：现 V1 单靶 6/12、V2 无强二级结构 12/12、V3 GC/Tm 合规 12/12、
> 总通过 6/12，与 BLAST v6 的 USABLE 集合一致。该 CSV 现为可信产物。
> 结构与 GC/Tm 合规数见附表 S10b 及 `audit_facts.json` 的
> `qpcr.v2v3_pass_n` / `qpcr.v3_gc_tm_pass_n`。

---

## S4.6 统计方案（预先指定，避免事后择优）

| 比较 | 设计 | 检验 | 多重校正 |
|---|---|---|---|
| **主分析**：DFU 溃疡 vs 配对非溃疡皮肤 | 配对（n=10 配对） | 配对 t 检验（Shapiro-Wilk 正态性；不服从则 Wilcoxon 符号秩） | BH-FDR，6 个目标基因家族内 |
| **次要分析**：DFU vs 健康对照皮肤 | 独立（10 vs 8） | 两样本 t 检验 / Mann-Whitney U | BH-FDR |
| 血液 vs 组织 | 配对 | 配对 t / Wilcoxon | BH-FDR |
| 基因表达 vs 丙酮酸负荷 | 相关 | Spearman（n 太小不做 Pearson 假设） | BH-FDR |
| 结局（healed vs non-healed） | 如样本量允许（预设为探索性，**不做检验**，仅描述） | — | — |

**预设的阳性判定标准（避免事后择优）**：
1. **主要终点**：6 个目标基因中，**≥3 个**在配对比较中方向与 META1 预设一致
   **且** BH-FDR < 0.05，视为"方向性复现成功"；
   ⚠️ **2026-10-25 注**：该标准以"6 个靶基因"为前提，而 BLAST v6 重跑已淘汰
   FAHD2B/ADH1B/ALDH3A2 三对引物（OFFTARGET_P0），合格靶基因仅剩 LDHA/GBP1/DHX58
   三个，"≥3/6"在合格子集上等价于"3/3 全部复现"，标准失去原有区分度——
   已由 **S4b 的 H1–H4 四个定向假设**取代（S4b 为最小可行集的预注册判据）。
   本条款留档以记录标准演化，不再作为执行判据；
2. **次要终点**：KpyRN **模块评分**（由板内可检测的 M2/M4 基因构成，用本项目
   `KpyRN_modules.csv` 的评分规则计算）在配对比较中差异 FDR < 0.05；
3. **丙酮酸**：组织与血浆丙酮酸在 DFU vs 配对非溃疡间升高，且与 LDHA 表达正相关；
4. **阴性结果的解读**：若未复现，**如实报告为"部分/未复现"**，
   并结合本项目已完成的公共数据结论讨论**批次效应、细胞组成混杂与 RT-qPCR
   检测窗口**等可能原因，**不得**改写统计阈值或更换检验方法以追求显著。

**技术质量门槛（不达标样本剔除，须在结果中报告剔除数与原因）**：
RQI ≥ 7.0；−RT 对照无扩增；内参 Cq 变异系数 < 0.5；技术复孔差值 < 0.5 Cq。

---

## S4.7 本方案与已完成公共数据验证的关系（诚实定位）

本研究的**公共数据验证已全部完成**（见正文 Methods/Results），构成独立的证据层：

| 数据集 | 层次 | 已完成结论 |
|---|---|---|
| GSE134431 / GSE199939 / GSE143735 | bulk meta | KpyRN 26/128 基因 FDR<0.05，26/26 两队列方向一致（DFU vs 对照） |
| GSE165816 | 单细胞（DFU vs 对照） | 细胞类型定位结果见正文 5.6 |
| GSE223964 | 单细胞（糖尿病 vs 非糖尿病皮肤） | 供体级 32 个模块对比全阴性（两条独立路径一致） |
| GSE230426 | PBMC（感染期 vs 缓解期） | M4 IFN-I 读出与 M3 Kpy 靶基因显著下调（三模型稳健） |
| GSE166120 | GeoMx DSP 空间/区域 | 受试体级全部阴性（ROI 级名义 p 经 FDR 与受试体级重算后不成立） |
| MR | 因果推断 | DHX58→T2D 1 对 FDR 显著，三重独立验证通过 |
| 分子对接 | 药物重定位 | 8/8 阳性对照 RMSD 0.00 Å；55 次对接，Δ 基准完整 |

**本节（湿实验方案）不改变上述任何已完成的结论，也不为其背书**；
其作用是提供**可执行、可复现、可被同行审查的组织水平验证路径**。
若后续实施，其结果将作为**独立证据层**加入，届时稿件的证据等级需重新评估。

---

## S4.8 溯源清单（本节所有数值的来源文件）

| 内容 | 来源文件 |
|---|---|
| 6 个目标基因的 META1 pooled logFC / SE / I² | `03_结果/meta/META1_DFU_vs_control.csv` |
| 引物序列、Tm、GC、产物长度、模板 accession | `03_结果/实验验证设计/qpcr_primers.csv` |
| 引物设计参数与每基因 5 对备选 | `03_结果/实验验证设计/primer_design_report.json` |
| 二级结构/二聚体独立复算、GC/Tm 合规 | `03_结果/实验验证设计/primer_verification.csv` |
| BLAST 单靶特异性（v6 重跑终判，6 USABLE / 6 OFFTARGET_P0） | `03_结果/实验验证设计/primer_blast_v6.csv`（分层明细 `primer_blast_v5_ldha_all.json`） |
| 引物 BLAST 判定列（写回） | `03_结果/实验验证设计/qpcr_primers.csv`（`blast_v6_verdict` / `blast_v6_P0_genes`） |
| KpyRN 模块定义与评分规则 | `01_数据/基因集/KpyRN_modules.csv` |
| 公共数据验证各数据集结论 | `03_结果/分析摘要.md`（5.1–5.11 节） |
| 引物设计/验证脚本 | `02_脚本/design_rtqpcr_primers.py`、`02_脚本/verify_primers.py`、`02_脚本/verify_primers_blast.py` |

---

## S4.9 补充材料编号与稿件位置对照

| 补充材料 | 文件 | 稿件位置 | 状态 |
|---|---|---|---|
| S1 | 样本分组总表 | Supplementary Table S1 | 已完成 |
| S2 | KpyRN 基因集溯源表 | Supplementary Table S2 | 已完成 |
| S3 | 对接参数与结合能表 | Supplementary Table S3 | 已完成 |
| **S4** | **本文件：DFU 组织 RT-qPCR 验证实验方案** | Supplementary Methods / 补充材料 | **预设方案（未实施）** |
| S5 | GSE166120 GeoMx 板覆盖清单 | Supplementary Table S5 | 已完成 |
| S6 | 全部脚本与产物清单（可复现性） | Supplementary Table S6 | 已完成 |

**状态声明（本节最重要的诚实性要求）**：S4 是**方案**，不是结果。稿件正文若引用
本节，只能表述为"we pre-specified a tissue-level RT-qPCR validation protocol
(Supplementary Methods S4)"，**不得**使用"we validated / we confirmed /
RT-qPCR showed"等表述。相应的 Limitations 必须包含："No wet-lab tissue validation
was performed; the tissue-level RT-qPCR protocol is pre-specified but not executed,
and the Kpy modification itself was not directly measured in any dataset."
