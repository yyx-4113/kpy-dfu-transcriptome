# meta_analysis.R — 多队列效应量 meta 分析 + KpyRN 交集（假说 H1 检验）
# 项目：丙酮酰化（Kpy）与糖尿病足（DFU）
# 输入：03_结果/DEG/*.csv；01_数据/基因集/KpyRN_modules.csv
# 输出：03_结果/meta/*.csv + KpyRN_meta_overview.csv + 火山图

Sys.setlocale("LC_ALL", "Chinese_China.65001")

library(metafor)
library(ggplot2)

BASE <- Sys.getenv("KPYRN_PROJECT_ROOT")
DEG  <- file.path(BASE, "03_结果/DEG")
OUT  <- file.path(BASE, "03_结果/meta")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

meta_deg <- function(files, labels, out_name) {
  # files: 各队列 limma 结果路径；labels: 队列名
  mats <- lapply(files, function(f) {
    d <- read.csv(f)
    d <- d[!duplicated(d$gene), c("gene", "logFC", "P.Value")]
    names(d) <- c("gene", "logFC", "p")
    # 由 limma 的 t 与自由度反推 SE：改用 adj.P.Val 不可得 SE；用 p 与 logFC 近似不可靠
    # -> 重新读取含 t 统计量的列（topTable 输出含 t 与 P.Value，B）
    full <- read.csv(f)
    full <- full[!duplicated(full$gene), ]
    data.frame(gene = full$gene, logFC = full$logFC,
               se = abs(full$logFC / full$t), p = full$P.Value)
  })
  common <- Reduce(intersect, lapply(mats, function(m) m$gene))
  res <- data.frame(gene = common)
  for (i in seq_along(mats)) {
    m <- mats[[i]][match(common, mats[[i]]$gene), ]
    res[[paste0("logFC_", labels[i])]] <- m$logFC
    res[[paste0("se_", labels[i])]] <- m$se
  }
  # 随机效应 meta（REML）
  pooled <- mapply(function(fc, se) {
    if (all(is.na(fc))) return(c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_))
    r <- tryCatch(rma.uni(yi = fc, vi = se^2, method = "REML", control = list(maxiter = 200)),
                  error = function(e) NULL)
    if (is.null(r)) return(c(NA, NA, NA, NA, NA))
    c(r$b[[1]], r$se, r$zval, r$pval, r$I2)
  }, res[[paste0("logFC_", labels[1])]], res[[paste0("se_", labels[1])]],
     SIMPLIFY = FALSE)
  # 逐基因 meta
  out <- do.call(rbind, lapply(seq_len(nrow(res)), function(i) {
    y <- as.numeric(res[i, paste0("logFC_", labels)])
    v <- as.numeric(res[i, paste0("se_", labels)])^2
    ok <- !is.na(y) & !is.na(v)
    if (sum(ok) < 1) return(data.frame(pooled_logFC = NA, pooled_se = NA, z = NA, p = NA, I2 = NA))
    r <- tryCatch(rma.uni(yi = y[ok], vi = v[ok], method = "REML"), error = function(e) NULL)
    if (is.null(r)) return(data.frame(pooled_logFC = NA, pooled_se = NA, z = NA, p = NA, I2 = NA))
    data.frame(pooled_logFC = r$b[[1]], pooled_se = r$se, z = r$zval, p = r$pval,
               I2 = ifelse(is.finite(r$I2), r$I2, NA))
  }))
  out$gene <- res$gene
  out$padj <- p.adjust(out$p, method = "BH")
  # 方向一致性
  fcs <- res[, paste0("logFC_", labels)]
  out$n_cohorts_pos <- rowSums(fcs > 0, na.rm = TRUE)
  out$n_cohorts_neg <- rowSums(fcs < 0, na.rm = TRUE)
  # 审计：保留各队列原始 logFC / SE
  for (lb in labels) {
    out[[paste0("logFC_", lb)]] <- res[[paste0("logFC_", lb)]]
    out[[paste0("se_", lb)]] <- res[[paste0("se_", lb)]]
  }
  out <- out[, c("gene", "pooled_logFC", "pooled_se", "z", "p", "padj", "I2",
                 "n_cohorts_pos", "n_cohorts_neg",
                 paste0("logFC_", labels), paste0("se_", labels))]
  write.csv(out, file.path(OUT, paste0(out_name, ".csv")), row.names = FALSE)
  n_deg <- sum(out$padj < 0.05, na.rm = TRUE)
  cat(sprintf("[META] %s: %d genes tested, %d meta-DEG (FDR<0.05)\n",
              out_name, nrow(out), n_deg))
  invisible(out)
}

# ---- Meta-C1: DFU vs 非溃疡皮肤 ----
m1 <- meta_deg(
  c(file.path(DEG, "GSE134431_DFU_vs_DFS_limma.csv"),
    file.path(DEG, "GSE199939_DFU_vs_control_limma.csv")),
  c("GSE134431", "GSE199939"), "META1_DFU_vs_control")

# ---- Meta-C2: 愈合 vs 不愈合 ----
m2 <- meta_deg(
  c(file.path(DEG, "GSE134431_Healer_vs_Nonhealer_limma.csv"),
    file.path(DEG, "GSE143735_Healer_vs_Nonhealer_limma.csv")),
  c("GSE134431", "GSE143735"), "META2_Healer_vs_Nonhealer")

# ---- KpyRN 交集总览（假说 H1 检验） ----
kpyrn <- read.csv(file.path(BASE, "01_数据/基因集/KpyRN_modules.csv"))
kpy_genes <- unique(kpyrn$gene_symbol[kpyrn$scoring_include == "Yes"])

overview <- kpyrn[kpyrn$scoring_include == "Yes", ]
overview$META1_logFC <- m1$pooled_logFC[match(overview$gene_symbol, m1$gene)]
overview$META1_padj  <- m1$padj[match(overview$gene_symbol, m1$gene)]
overview$META2_logFC <- m2$pooled_logFC[match(overview$gene_symbol, m2$gene)]
overview$META2_padj  <- m2$padj[match(overview$gene_symbol, m2$gene)]
write.csv(overview, file.path(OUT, "KpyRN_meta_overview.csv"), row.names = FALSE)

cat("\n===== KpyRN × DFU 假说检验（H1）=====\n")
for (mod in unique(overview$module)) {
  sub <- overview[overview$module == mod, ]
  up <- sum(sub$META1_logFC > 0, na.rm = TRUE)
  dn <- sum(sub$META1_logFC < 0, na.rm = TRUE)
  sig <- sum(sub$META1_padj < 0.05, na.rm = TRUE)
  cat(sprintf("  %-22s n=%3d  META1: up=%3d down=%3d  meta-DEG=%2d\n", mod, nrow(sub), up, dn, sig))
}
cat(sprintf("\nKpyRN 总体：META1 显著 %d/%d；META2 显著 %d/%d\n",
            sum(overview$META1_padj < 0.05, na.rm = TRUE), nrow(overview),
            sum(overview$META2_padj < 0.05, na.rm = TRUE), nrow(overview)))

# 方向检验（M2 丙酮酸供应上调 / M4 IFN-I 读出下调）
m2mod <- overview[overview$module == "M2_Pyruvate_supply", ]
m4mod <- overview[overview$module == "M4_IFNI_readout", ]
cat(sprintf("  M2 logFC>0 比例: %d/%d (%.0f%%)\n", sum(m2mod$META1_logFC > 0, na.rm = TRUE),
            sum(!is.na(m2mod$META1_logFC)), 100 * mean(m2mod$META1_logFC > 0, na.rm = TRUE)))
cat(sprintf("  M4 logFC<0 比例: %d/%d (%.0f%%)\n", sum(m4mod$META1_logFC < 0, na.rm = TRUE),
            sum(!is.na(m4mod$META1_logFC)), 100 * mean(m4mod$META1_logFC < 0, na.rm = TRUE)))

# ---- 火山图（META1 + META2，标注 KpyRN 基因） ----
plot_volcano <- function(m, title, out_file) {
  # 注：metafor 的 rma.uni 输出 p 值列名为 `p`（不是 `p_raw`）。
  #     旧版此处写 `m$p_raw`，对 NULL 列做 is.na() 过滤会使 subset 全部丢弃，
  #     产出完全空白的 PNG（META1_volcano.png 非白像素=0）。已修正为 `p`。
  m <- m[!is.na(m$pooled_logFC) & !is.na(m$p), ]
  m$sig <- !is.na(m$padj) & m$padj < 0.05
  m$kpy <- as.factor(m$gene %in% kpy_genes)
  p <- ggplot(m, aes(x = pooled_logFC, y = -log10(pmax(p, 1e-300)))) +
    geom_point(aes(color = sig, size = kpy), alpha = 0.5) +
    scale_color_manual(values = c("grey50", "#C0392B")) +
    scale_size_manual(values = c(0.8, 2.2)) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey30") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
    labs(title = title, x = "pooled log2FC", y = "-log10(p)") +
    theme_classic(base_size = 12)
  ggssave <- function(...) ggplot2::ggsave(...)
  ggplot2::ggsave(out_file, p, width = 7, height = 6, dpi = 150)
}
plot_volcano(m1, "META1: DFU vs non-ulcer skin (red=meta-DEG, large=KpyRN)",
             file.path(OUT, "META1_volcano.png"))
plot_volcano(m2, "META2: healer vs nonhealer (red=meta-DEG, large=KpyRN)",
             file.path(OUT, "META2_volcano.png"))
cat("\nMETA DONE\n")
