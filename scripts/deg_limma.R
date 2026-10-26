# deg_limma.R — 逐队列差异表达分析（limma, trend=TRUE）
# 项目：丙酮酰化（Kpy）与糖尿病足（DFU）
# 输入：01_数据/处理后数据/<GSE>_log2_expr.csv + _annotation.csv
# 输出：03_结果/DEG/<GSE>_<contrast>_limma.csv（全基因排序统计量）
# 过滤规则：基因在较小组的 >=50% 样本中 log2 表达 > 1（预注册，防止低表达噪声）

# Windows R 中文路径支持：切换到 UTF-8 代码页 locale（须在任何中文路径使用前执行）
Sys.setlocale("LC_ALL", "Chinese_China.65001")

library(limma)

BASE <- Sys.getenv("KPYRN_PROJECT_ROOT")
PROC <- file.path(BASE, "01_数据/处理后数据")
OUT  <- file.path(BASE, "03_结果/DEG")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

run_deg <- function(gse, contrast_label, group_col, level1, level2,
                    subset_col = NULL, subset_level = NULL) {
  expr <- read.csv(file.path(PROC, paste0(gse, "_log2_expr.csv")),
                   row.names = 1, check.names = FALSE)
  ann <- read.csv(file.path(PROC, paste0(gse, "_annotation.csv")),
                  row.names = 1, check.names = FALSE)
  # GSE143735 派生溃疡状态列
  if (gse == "GSE143735" && "ulcer_status" %in% names(ann)) {
    ann$ulcer_status <- ifelse(ann$group2 %in% c("DFU_healer", "DFU_nonhealer"),
                               "DFU", "DM_no_ulcer")
  }
  # GSE143735 派生溃疡状态列（DFU 各愈合亚组合并为 DFU）
  if (gse == "GSE143735") {
    ann$ulcer_status <- ifelse(ann$group2 %in% c("DFU_healer", "DFU_nonhealer"),
                               "DFU", "DM_no_ulcer")
  }
  if (!is.null(subset_col)) {
    keep <- ann[[subset_col]] == subset_level
    ann <- ann[keep, , drop = FALSE]
  }
  # 仅保留两对比组的样本，并剔除分组列为空/NA 者
  ann <- ann[ann[[group_col]] %in% c(level1, level2), , drop = FALSE]
  expr <- expr[, rownames(ann), drop = FALSE]
  grp <- factor(ann[[group_col]], levels = c(level2, level1))
  stopifnot(length(levels(grp)) == 2)
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- levels(grp)
  # 低表达过滤：较小组 >=50% 样本 log2>1
  keep <- rowSums(expr[, grp == level1, drop = FALSE] > 1) >= 0.5 * sum(grp == level1) |
          rowSums(expr[, grp == level2, drop = FALSE] > 1) >= 0.5 * sum(grp == level2)
  n_before <- nrow(expr)
  expr <- expr[keep, , drop = FALSE]
  fit <- lmFit(expr, design)
  cont <- makeContrasts(contrasts = paste0(level1, "-", level2), levels = design)
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2, trend = TRUE)
  tt <- topTable(fit2, number = Inf, sort.by = "p")
  tt$gene <- rownames(tt)
  out_file <- file.path(OUT, paste0(gse, "_", contrast_label, "_limma.csv"))
  write.csv(tt, out_file, row.names = FALSE)
  n_up <- sum(tt$adj.P.Val < 0.05 & tt$logFC > 0)
  n_dn <- sum(tt$adj.P.Val < 0.05 & tt$logFC < 0)
  cat(sprintf("[DEG] %s %s: %d->%d genes; up=%d down=%d (FDR<0.05)\n",
              gse, contrast_label, n_before, nrow(expr), n_up, n_dn))
}

# ---- C1: DFU/溃疡 vs 非溃疡皮肤（发生） ----
run_deg("GSE134431", "DFU_vs_DFS",   "group2", "DFU", "DFS")
run_deg("GSE199939", "DFU_vs_control", "group2", "DFU", "control")

# ---- C2: 愈合 vs 不愈合（结局） ----
run_deg("GSE134431", "Healer_vs_Nonhealer", "healing", "Healer", "Nonhealer")
run_deg("GSE143735", "Healer_vs_Nonhealer", "group2", "DFU_healer", "DFU_nonhealer")

# ---- C3: 溃疡 vs 糖尿病非溃疡皮肤（GSE143735，糖尿病内部对比） ----
run_deg("GSE143735", "DFU_vs_DM_no_ulcer", "ulcer_status", "DFU", "DM_no_ulcer")

cat("DEG DONE\n")
