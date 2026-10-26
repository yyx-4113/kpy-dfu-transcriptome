# 稿件数字审计追踪报告

- 稿件：`05_稿件/manuscript_en.md`（366 行）
- 事实库：`05_稿件/audit_facts.json`（660 条）
- **事实覆盖率**：635/660 = **96.2%**
- **未溯源数字 token**：0 处，0 个不同取值

## 方法说明

期望值全部由 `audit_facts.json` 动态派生，本审计脚本内**不硬编码任何期望数字**。
事实库每条记录包含：值、类型、源文件路径、提取表达式。正文数字的匹配采用多形式
指纹（整数、1–6 位小数、科学计数法、百分数、负号变体），避免格式差异造成的假阴性。

未溯源 token 主要来自：统计方法学常数（如 0.05、0.80、10000 次置换）、
软件版本、阈值设定与生物学常数（呼吸链复合体编号、氨基酸残基位置）。
这些是**方法参数**而非结果数字，故不进入事实库；逐条列于下表供人工判读。

## 未溯源数字 token 明细（去重，按出现次数排序）

（无）

## 逐事实覆盖明细

完整表见 `05_稿件/审计追踪表.csv`。未在稿件中直接出现的 fact 分为两类：
（i）支撑性事实（中间量、敏感性分析、清单项），在稿件中以聚合或表格形式引用；
（ii）待补充事实——需在终稿中显式报告或明确标注为未使用。

未被正文直接引用的事实：25 条

| fact_id | 值 | 源文件 |
|---|---|---|
| `wgcna.module.black.num` | 7 | 03_结果/WGCNA/module_trait_cor.csv |
| `wgcna.module.green.num` | 5 | 03_结果/WGCNA/module_trait_cor.csv |
| `scrna165816.n_cells_ambiguous` | 21235 | 01_数据/处理后数据/scRNA/scRNA_GSE165816_cell_scores.csv |
| `scrna165816.pseudobulk.Fibroblast.PGAM1.nondfu_median` | 7.54 | 03_结果/单细胞/pseudobulk_Foot_skin_DFU_vs_nonDFU_by_lineage.csv |
| `scrna165816.pseudobulk.Fibroblast.ENO1.nondfu_median` | 8.17 | 03_结果/单细胞/pseudobulk_Foot_skin_DFU_vs_nonDFU_by_lineage.csv |
| `scrna165816.pseudobulk.Pericyte_SMC.PGAM1.dfu_median` | 7.74 | 03_结果/单细胞/pseudobulk_Foot_skin_DFU_vs_nonDFU_by_lineage.csv |
| `scrna165816.pseudobulk.Pericyte_SMC.TPI1.nondfu_median` | 8.41 | 03_结果/单细胞/pseudobulk_Foot_skin_DFU_vs_nonDFU_by_lineage.csv |
| `scrna165816.pseudobulk.Pericyte_SMC.PKM.nondfu_median` | 8.18 | 03_结果/单细胞/pseudobulk_Foot_skin_DFU_vs_nonDFU_by_lineage.csv |
| `mr.n_nominal_only` | 5 | 03_结果/MR/mr_results.csv |
| `docking.n_pdb` | 7 | 03_结果/药物预测/docking/docking_results.csv |
| `docking.PKM2__C_pocket_allosteric.negative_affinity` | -6.96 | 03_结果/药物预测/docking/docking_results.csv |
| `docking.PKM2__C_pocket_allosteric.n_candidates` | 5 | 03_结果/药物预测/docking/docking_results.csv |
| `docking.PKM2__A_pocket_FBP_site.negative_affinity` | -5.38 | 03_结果/药物预测/docking/docking_results.csv |
| `docking.SIRT3__ADP_ribose_NADplus_C_pocket.negative_affinity` | -6.92 | 03_结果/药物预测/docking/docking_results.csv |
| `docking.LDHA__active_site.negative_affinity` | -6.63 | 03_结果/药物预测/docking/docking_results.csv |
| `scrna223964.n_donors_diabetic` | 5 | 03_结果/单细胞GSE223964/pseudobulk_condition_stats.csv |
| `geommx.roi.ulcer_edge` | 5 | 03_结果/空间GSE166120/summary.json |
| `qpcr.status` | pre-specified, not performed | 06_补充材料/S4_RT-qPCR_验证方案.md |
| `qpcr.amplicon_max_bp` | 180 | 02_脚本/design_rtqpcr_primers.py |
| `qpcr.v4.n_pairs_with_3p_anchor_mismatch` | 5 | 03_结果/实验验证设计/final_primers_v4.csv |
| `qpcr.v6.n_P1_hits_total` | 141 | 03_结果/实验验证设计/primer_blast_v6.csv |
| `qpcr.cds_map.total_repeat_masked_nt` | 3880 | 03_结果/实验验证设计/cds_duplication_map.json |
| `method.estimate.total_genes` | 280 | 01_数据/基因集/免疫/estimate_signatures.csv |
| `method.geommx.pk_c_total_probes` | 8684 | 01_数据/原始数据/GSE166120/GSE166120_Alpha_CTA_v2.0.pkc.txt.gz |
| `dataset.n_bulk` | 5 | 01_数据/元数据/样本表/ |