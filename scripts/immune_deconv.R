# immune_deconv.R — 免疫浸润反卷积（EPIC + MCPcounter 双方法交叉）
# 项目：丙酮酰化（Kpy）与糖尿病足（DFU）
#
# 方法说明（2026-10-24 更新）：
#   immunedeconv 已从 Bioconductor 移除（官方 oldstats 页确认"不属于当前 release/devel"），
#   且 xCell2 需 R>=4.5（本机 R-4.4.3 对应 Bioc 3.20 无此包）。
#   改用两个独立的底层方法直接实现：
#   1) EPIC 1.1.7（GfellerLab/EPIC, GitHub 安装）——非负最小二乘解卷积，参考谱 BRef，
#      输出 7 类免疫/间质细胞比例 + total immune/stromal/other。
#   2) MCPcounter（按 ebecht/MCPcounter 官方 Source/R/MCPcounter.R 算法直接实现，
#      签名取自官方 Signatures/genes.txt，HUGO symbols 版 111 个标记基因）。
#   3) ESTIMATE（stromal/immune/ESTIMATE score）在 Python 端实现
#      （算法复刻 omnideconv/immunedeconv R/estimate.R，签名取自其 extdata）。
#
# 输入：01_数据/处理后数据/<GSE>_log2_expr.csv（行=基因符号，列=样本，log2 标度）
# 输出：03_结果/免疫浸润/<GSE>_<method>_fractions.csv + session_info.txt

Sys.setlocale("LC_ALL", "Chinese_China.65001")

BASE <- Sys.getenv("KPYRN_PROJECT_ROOT")
PROC <- file.path(BASE, "01_数据/处理后数据")
OUT  <- file.path(BASE, "03_结果/免疫浸润")
SIG  <- file.path(BASE, "01_数据/基因集/免疫")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

library(EPIC)

GSES <- c("GSE134431", "GSE199939")

# ---- MCPcounter 官方算法直实现（行=基因符号，列=样本；输入须为 log2 标度） ----
mcp_counter_estimate <- function(expression, markers_file) {
  markers <- read.table(markers_file, sep = "\t", stringsAsFactors = FALSE,
                        header = TRUE, colClasses = "character", check.names = FALSE)
  markersG <- markers
  features <- subset(markersG, get("HUGO symbols") %in% rownames(expression))
  markers.names <- unique(features[, "Cell population"])
  features <- split(features[, "HUGO symbols"], features[, "Cell population"])
  features <- features[intersect(markers.names, names(features))]
  res <- as.data.frame(do.call(cbind, lapply(features, function(x) {
    apply(expression[intersect(rownames(expression), x), , drop = FALSE], 2,
          mean, na.rm = TRUE)
  })))
  res
}

sink(file.path(OUT, "immune_deconv.log"))
cat("=== immune_deconv.R ===\n")
cat("R version:", R.version.string, "\n")
cat("EPIC version:", as.character(packageVersion("EPIC")), "\n")
cat("MCPcounter: direct implementation of ebecht/MCPcounter Source/R/MCPcounter.R\n")
cat("ESTIMATE: implemented in Python (immune_estimate_correlate.py)\n\n")

for (gse in GSES) {
  expr <- read.csv(file.path(PROC, paste0(gse, "_log2_expr.csv")),
                   row.names = 1, check.names = FALSE)
  mat <- as.matrix(expr)
  cat(sprintf("[%s] matrix: %d genes x %d samples\n", gse, nrow(mat), ncol(mat)))

  # ---- 方法 1: EPIC ----
  epic_res <- tryCatch({
    r <- EPIC(mat, reference = "BRef")  # scaleExprs=TRUE 为默认；BRef=非恶性参考谱
    frac <- t(r$cellFractions)
    write.csv(cbind(sample = rownames(frac), frac),
              file.path(OUT, paste0(gse, "_epic_fractions.csv")),
              row.names = FALSE)
    cat(sprintf("  [%s] EPIC: OK (%d cell types + %s)\n",
                gse, ncol(frac) - 3, paste(colnames(frac)[(ncol(frac) - 2):ncol(frac)], collapse = "/")))
  }, error = function(e) {
    cat(sprintf("  [%s] EPIC: FAILED %s\n", gse, conditionMessage(e)))
  })

  # ---- 方法 2: MCPcounter ----
  mcp_res <- tryCatch({
    r <- mcp_counter_estimate(mat, file.path(SIG, "MCPcounter_genes.txt"))
    write.csv(cbind(sample = rownames(r), r),
              file.path(OUT, paste0(gse, "_mcpcounter_fractions.csv")),
              row.names = FALSE)
    cat(sprintf("  [%s] MCPcounter: OK (%d populations: %s)\n",
                gse, ncol(r), paste(colnames(r), collapse = ", ")))
  }, error = function(e) {
    cat(sprintf("  [%s] MCPcounter: FAILED %s\n", gse, conditionMessage(e)))
  })
  cat("\n")
}
cat("IMMUNE DECONV DONE\n")
sink()
