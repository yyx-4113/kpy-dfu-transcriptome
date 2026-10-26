# redraw_meta_volcano.R — 仅重绘 META1/META2 火山图（零统计重算）
#
# 背景：
#   meta_analysis.R 的 plot_volcano() 误用 `m$p_raw`。metafor::rma.uni 的输出
#   列名为 `p`，`m$p_raw` 为 NULL；`!is.na(NULL)` 返回 logical(0)，使
#   `m[...]` 子集全部丢弃 → 产出**完全空白**的 PNG（实测 META1_volcano.png
#   非白像素 = 0，而 META2 有 49,909，源于两者生成时间不同）。
#
# 本脚本：
#   - **不重算任何统计量**：直接读取已审计的 META1/META2 结果 CSV
#     （03_结果/meta/META1_DFU_vs_control.csv、META2_Healer_vs_Nonhealer.csv）
#     以及 KpyRN 基因集，重绘火山图。
#   - 修正后的 plot_volcano 逻辑与 meta_analysis.R 中已修复的版本一致。
#   - 幂等：可重复运行，输出确定。
#
# 输入：03_结果/meta/META{1,2}_*.csv；01_数据/基因集/KpyRN_modules.csv
# 输出：03_结果/meta/META1_volcano.png；03_结果/meta/META2_volcano.png

Sys.setlocale("LC_ALL", "Chinese_China.65001")
suppressPackageStartupMessages({
  library(ggplot2)
})

BASE <- Sys.getenv("KPYRN_PROJECT_ROOT")
OUT  <- file.path(BASE, "03_结果/meta")

kpyrn <- read.csv(file.path(BASE, "01_数据/基因集/KpyRN_modules.csv"))
kpy_genes <- unique(kpyrn$gene_symbol[kpyrn$scoring_include == "Yes"])

plot_volcano <- function(m, title, out_file) {
  # 修正点：p 值列名为 `p`（metafor::rma.uni 输出），非 `p_raw`
  m <- m[!is.na(m$pooled_logFC) & !is.na(m$p), ]
  m$sig <- !is.na(m$padj) & m$padj < 0.05
  m$kpy <- as.factor(m$gene %in% kpy_genes)
  cat(sprintf("  [%s] 参与绘图的基因数 = %d；meta-DEG(padj<0.05) = %d；KpyRN 计分基因 = %d\n",
              basename(out_file), nrow(m), sum(m$sig), sum(m$kpy == "TRUE")))
  p <- ggplot(m, aes(x = pooled_logFC, y = -log10(pmax(p, 1e-300)))) +
    geom_point(aes(color = sig, size = kpy), alpha = 0.5) +
    scale_color_manual(values = c("grey50", "#C0392B")) +
    scale_size_manual(values = c(0.8, 2.2)) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey30") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
    labs(title = title, x = "pooled log2FC", y = "-log10(p)") +
    theme_classic(base_size = 12)
  ggplot2::ggsave(out_file, p, width = 7, height = 6, dpi = 150)
}

cat("重绘 META 火山图（仅绘图，不重算统计量）\n")

m1 <- read.csv(file.path(OUT, "META1_DFU_vs_control.csv"))
plot_volcano(m1, "META1: DFU vs non-ulcer skin (red=meta-DEG, large=KpyRN)",
             file.path(OUT, "META1_volcano.png"))

m2 <- read.csv(file.path(OUT, "META2_Healer_vs_Nonhealer.csv"))
plot_volcano(m2, "META2: healer vs nonhealer (red=meta-DEG, large=KpyRN)",
             file.path(OUT, "META2_volcano.png"))

cat("REDRAW DONE\n")
