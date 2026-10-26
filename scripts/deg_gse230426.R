# deg_gse230426.R — GSE230426 感染队列次级分析（任务 #11a）
# 项目：丙酮酰化（Kpy）与糖尿病足（DFU）
#
# 设计（series matrix 2026-10-25 核实）：
#   S. aureus 感染 DFU 患者 PBMC；BEFORE THERAPY(0wk, n=20) vs AFTER THERAPY(n=18: 15×8wk+3×4wk)
#   主分析：非配对 limma-voom（保守；配对结构为推断性，不作为主分析依据）
#   敏感性 A：17 推断配对（连续ID+同年龄+互补时点），~donor + stage
#   敏感性 B：仅 AFTER 8wk（剔除 3 个 4wk）vs BEFORE，非配对
#   模块检验：camera（robust，inter.roast 基准），M1–M4（评分模块）
#   结局对比（healed vs non-healed）：公开元数据无该字段 → 不构建（如实记录）
# 输出：03_结果/感染队列GSE230426/

Sys.setlocale("LC_ALL", "Chinese_China.65001")
library(limma); library(edgeR); library(jsonlite)

BASE <- Sys.getenv("KPYRN_PROJECT_ROOT")
PROC <- file.path(BASE, "01_数据/处理后数据")
OUT  <- file.path(BASE, "03_结果/感染队列GSE230426")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cnt <- read.csv(file.path(PROC, "GSE230426_counts.csv"), row.names = 1, check.names = FALSE)
ann <- read.csv(file.path(PROC, "GSE230426_annotation.csv"), row.names = 1,
                check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(identical(colnames(cnt), rownames(ann)))
cat(sprintf("counts: %d genes × %d samples\n", nrow(cnt), ncol(cnt)))

K <- read.csv(file.path(BASE, "01_数据/基因集/KpyRN_modules.csv"),
              stringsAsFactors = FALSE)
K_scored <- K[K$scoring_include == "Yes" & K$gene_symbol %in% rownames(cnt), ]
idx <- split(K_scored$gene_symbol, K_scored$module)
cat(sprintf("评分模块: %s\n",
    paste(sprintf("%s=%d", names(idx), lengths(idx)), collapse = ", ")))

# ML 共识枢纽基因（TaskA，n_algorithms>=3）
hub <- c("LDHA", "FAHD2B", "ADH1B", "DHX58", "ALDH3A2", "GBP1")
stopifnot(all(hub %in% rownames(cnt)))

# ---------- voom：TMM 归一化 + 精度权重 ----------
dge <- DGEList(counts = cnt, samples = ann)
dge <- calcNormFactors(dge, method = "TMM")

fit_hub_table <- function(fit, cont, label) {
  f <- eBayes(contrasts.fit(fit, cont), robust = TRUE)
  tb <- topTable(f, number = Inf, sort.by = "none")
  tb$gene <- rownames(tb)
  tb$model <- label
  tb
}
summarize_model <- function(tb, label) {
  list(
    label = label,
    n_genes = nrow(tb),
    n_fdr05 = sum(tb$adj.P.Val < 0.05),
    dir_up = sum(tb$adj.P.Val < 0.05 & tb$logFC > 0),
    dir_down = sum(tb$adj.P.Val < 0.05 & tb$logFC < 0)
  )
}

results <- list(); deg_tables <- list()

# ---------- 主分析：非配对 AFTER vs BEFORE ----------
ann$stage_f <- factor(ann$stage, levels = c("BEFORE", "AFTER"))
des_main <- model.matrix(~ 0 + stage_f, data = ann)
colnames(des_main) <- c("BEFORE", "AFTER")
v_main <- voom(dge, design = des_main, plot = FALSE)
fit_main <- lmFit(v_main, des_main)
cont_main <- makeContrasts(AFTER - BEFORE, levels = des_main)
tb_main <- fit_hub_table(fit_main, cont_main, "main_unpaired")
deg_tables$main <- tb_main
results$main <- summarize_model(tb_main, "main_unpaired (AFTER vs BEFORE, n=20/18)")

# ---------- 敏感性 A：推断配对（~donor + stage）----------
ann_p <- ann[ann$donor_inferred != "unpaired", ]
cnt_p <- cnt[, rownames(ann_p)]
ann_p$donor_f <- factor(ann_p$donor_inferred)
ann_p$stage_f <- factor(ann_p$stage, levels = c("BEFORE", "AFTER"))
des_p <- model.matrix(~ 0 + donor_f + stage_f, data = ann_p)
# 17 对完美配对（每 donor 恰 1×BEFORE+1×AFTER）→ stage_fBEFORE 与 donor 哑元
# 完全共线被自动丢弃（秩亏去共线，18 列）。保留的 AFTER 列系数即"配对内
# AFTER−BEFORE 差"（等价配对 t 检验），故对比取单项 AFTER。
cn <- gsub("donor_f|stage_f", "", colnames(des_p))
colnames(des_p) <- cn
stopifnot("AFTER" %in% cn, !"BEFORE" %in% cn)
dge_p <- DGEList(counts = cnt_p, samples = ann_p)
dge_p <- calcNormFactors(dge_p, method = "TMM")
v_p <- voom(dge_p, design = des_p, plot = FALSE)
fit_p <- lmFit(v_p, des_p)
cont_p <- makeContrasts(AFTER, levels = des_p)
tb_p <- fit_hub_table(fit_p, cont_p, "paired_inferred")
deg_tables$paired <- tb_p
results$paired <- summarize_model(tb_p, "paired_inferred (17 对, ~donor+stage)")

# ---------- 敏感性 B：AFTER 8wk only ----------
ann8 <- ann[!(ann$stage == "AFTER" & ann$week == 4), ]
cnt8 <- cnt[, rownames(ann8)]
ann8$stage_f <- factor(ann8$stage, levels = c("BEFORE", "AFTER"))
des8 <- model.matrix(~ 0 + stage_f, data = ann8)
colnames(des8) <- c("BEFORE", "AFTER")
dge8 <- DGEList(counts = cnt8, samples = ann8)
dge8 <- calcNormFactors(dge8, method = "TMM")
v8 <- voom(dge8, design = des8, plot = FALSE)
fit8 <- lmFit(v8, des8)
cont8 <- makeContrasts(AFTER - BEFORE, levels = des8)
tb8 <- fit_hub_table(fit8, cont8, "after8wk_only")
deg_tables$after8 <- tb8
results$after8 <- summarize_model(tb8, "after8wk_only (BEFORE n=20 vs AFTER-8wk n=15)")

# ---------- 模块级 camera（主模型 + 配对敏感性 + 8wk-only 敏感性）----------
run_camera <- function(v, des, cont, label) {
  cm <- camera(v, idx, design = des, contrast = cont)
  data.frame(module = rownames(cm), ngenes = cm$NGenes,
             direction = cm$Direction, p = cm$PValue, fdr = cm$FDR,
             model = label, row.names = NULL)
}
cam_main <- run_camera(v_main, des_main, cont_main, "main_unpaired")
cam_p    <- run_camera(v_p, des_p, cont_p, "paired_inferred")
cam_8    <- run_camera(v8, des8, cont8, "after8wk_only")
cam <- rbind(cam_main, cam_p, cam_8)
write.csv(cam, file.path(OUT, "camera_modules.csv"), row.names = FALSE)
cat("\n== camera 模块检验 ==\n"); print(cam[, c("model","module","ngenes","direction","fdr")])

# ---------- 枢纽基因表（三个模型）----------
hub_tab <- do.call(rbind, lapply(deg_tables, function(tb)
  tb[tb$gene %in% hub, c("gene", "logFC", "t", "P.Value", "adj.P.Val", "model")]))
hub_tab$hub <- hub_tab$gene %in% c("LDHA","FAHD2B","ADH1B","DHX58","ALDH3A2","GBP1")
write.csv(hub_tab[order(hub_tab$gene), ], file.path(OUT, "hub_genes.csv"), row.names = FALSE)
cat("\n== 枢纽基因（6 个 ML 共识基因）==\n")
print(hub_tab[order(hub_tab$gene, hub_tab$model), ], digits = 3)

# ---------- KpyRN 全基因方向一致性（主 vs 配对敏感性）----------
m1 <- tb_main[tb_main$gene %in% K_scored$gene_symbol, ]
m2 <- tb_p[tb_p$gene %in% K_scored$gene_symbol, ]
m1 <- m1[order(m1$gene), ]; m2 <- m2[order(m2$gene), ]
stopifnot(identical(m1$gene, m2$gene))
conc <- cor(m1$logFC, m2$logFC, method = "spearman")
dir_agree <- mean(sign(m1$logFC) == sign(m2$logFC))
results$concordance_main_vs_paired <- list(
  spearman_logFC = round(conc, 4),
  direction_agreement = round(dir_agree, 4),
  n_genes = nrow(m1))
cat(sprintf("\n主 vs 配对敏感性（KpyRN %d 基因）：Spearman logFC r=%.3f，方向一致率 %.1f%%\n",
            nrow(m1), conc, 100 * dir_agree))

# ---------- 落盘 ----------
for (nm in names(deg_tables))
  write.csv(deg_tables[[nm]][order(deg_tables[[nm]]$adj.P.Val), ],
            file.path(OUT, sprintf("deg_%s.csv", nm)), row.names = FALSE)
write_json(results, file.path(OUT, "summary.json"), auto_unbox = TRUE,
           pretty = TRUE, digits = 6)
cat("\nDONE — deg_gse230426.R\n")
