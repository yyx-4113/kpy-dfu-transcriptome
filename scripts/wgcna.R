# wgcna.R — WGCNA 模块挖掘（GSE134431，n=21）
# 项目：丙酮酰化（Kpy）与糖尿病足（DFU）
# 输入：01_数据/处理后数据/GSE134431_*；03_结果/meta/META1/META2
# 输出：03_结果/WGCNA/（模块归属、模块-性状相关、三层交集候选池）
# 说明：n=21 为 WGCNA 可行下限；软阈值拟合不佳时如实报告并采用文献常用 power=6 兜底。

Sys.setlocale("LC_ALL", "Chinese_China.65001")

library(WGCNA)
library(ggplot2)

BASE <- Sys.getenv("KPYRN_PROJECT_ROOT")
PROC <- file.path(BASE, "01_数据/处理后数据")
OUT  <- file.path(BASE, "03_结果/WGCNA")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

enableWGCNAThreads()
options(stringsAsFactors = FALSE)

expr <- read.csv(file.path(PROC, "GSE134431_log2_expr.csv"), row.names = 1, check.names = FALSE)
ann  <- read.csv(file.path(PROC, "GSE134431_annotation.csv"), row.names = 1, check.names = FALSE)
expr <- expr[, rownames(ann)]

# 性状
datTraits <- data.frame(
  DFU = ifelse(ann$group2 == "DFU", 1, 0),
  Healer = ifelse(ann$healing == "Healer", 1, ifelse(ann$healing == "Nonhealer", 0, NA))
)
rownames(datTraits) <- rownames(ann)

# 取表达量最高的前 5000 基因（按 MAD），降低 n=21 的维度灾难
mads <- apply(expr, 1, mad)
top_genes <- names(sort(mads, decreasing = TRUE))[1:5000]
datExpr <- t(expr[top_genes, ])
cat(sprintf("WGCNA input: %d genes x %d samples\n", ncol(datExpr), nrow(datExpr)))

# 软阈值
powers <- c(1:10, 12, 14, 16, 18, 20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)
fit <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
best_power <- sft$fitIndices$power[which.max(fit)]
best_r2 <- max(fit, na.rm = TRUE)
cat(sprintf("Soft threshold: best power=%d (signed R^2=%.3f); n=21 小样本，若 R^2<0.8 采用 power=6 兜底并如实报告\n",
            best_power, best_r2))
use_power <- ifelse(best_r2 >= 0.8, best_power, 6)

# 网络构建
net <- blockwiseModules(datExpr, power = use_power, networkType = "signed",
                        TOMType = "signed", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = FALSE, verbose = 0)
moduleColors <- labels2colors(net$colors)
cat(sprintf("Modules detected: %d\n", length(table(net$colors))))

# 模块-性状相关
MEs <- net$MEs
MEs0 <- MEs[, -match("ME0", colnames(MEs))]
moduleTraitCor <- cor(MEs0, datTraits, use = "p")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# 输出模块-性状相关热图数据
mt <- data.frame(module = sub("ME", "", colnames(MEs0)),
                 color = gsub("ME", "", colnames(MEs0)),
                 cor_DFU = moduleTraitCor[, "DFU"], p_DFU = moduleTraitP[, "DFU"],
                 cor_Healer = moduleTraitCor[, "Healer"], p_Healer = moduleTraitP[, "Healer"])
mt <- mt[order(-abs(mt$cor_DFU)), ]
write.csv(mt, file.path(OUT, "module_trait_cor.csv"), row.names = FALSE)
print(mt)

# 基因-模块归属 + kME
geneModules <- data.frame(gene = colnames(datExpr), module = moduleColors)
kme <- signedKME(datExpr, MEs0)
for (m in colnames(kme)) geneModules[[m]] <- kme[, m]
write.csv(geneModules, file.path(OUT, "gene_module_membership.csv"), row.names = FALSE)

# 与 DFU 显著相关的模块（数值编号 -> 颜色名）
sig_mod_nums <- as.numeric(mt$module[mt$p_DFU < 0.05])
sig_mods <- as.character(labels2colors(sig_mod_nums))
cat(sprintf("\n与 DFU 显著相关的模块(编号->颜色): %s\n",
            paste(sprintf("%s->%s", sig_mod_nums, sig_mods), collapse = ", ")))

# ---- 三层交集：META-DEG ∩ DFU 相关模块 ∩ KpyRN ----
kpyrn <- read.csv(file.path(BASE, "01_数据/基因集/KpyRN_modules.csv"))
kpy_genes <- unique(kpyrn$gene_symbol[kpyrn$scoring_include == "Yes"])
m1 <- read.csv(file.path(BASE, "03_结果/meta/META1_DFU_vs_control.csv"))
deg1 <- m1$gene[m1$padj < 0.05]
m2 <- read.csv(file.path(BASE, "03_结果/meta/META2_Healer_vs_Nonhealer.csv"))
deg2 <- m2$gene[m2$padj < 0.05]

mod_genes <- geneModules$gene[geneModules$module %in% sig_mods]
pool <- list(
  L1_metaDEG = deg1,
  L2_WGCNA_module = mod_genes,
  L3_KpyRN = kpy_genes
)
three_layer <- Reduce(intersect, pool)
cat(sprintf("\n三层交集（META1-DEG ∩ DFU模块 ∩ KpyRN）: %d 基因\n", length(three_layer)))
cat(paste(sort(three_layer), collapse = ", "), "\n")

# 两层交集（备用候选池扩大）：DEG ∩ KpyRN
deg_kpy <- intersect(deg1, kpy_genes)
cat(sprintf("两层交集（META1-DEG ∩ KpyRN）: %d 基因\n", length(deg_kpy)))

# 候选池定义（预注册）：三层交集 ∪ (两层交集中 kME 最高的基因)
cand <- unique(c(three_layer, deg_kpy))
cand_df <- data.frame(gene = cand)
cand_df$module <- geneModules$module[match(cand, geneModules$gene)]
cand_df$META1_logFC <- m1$pooled_logFC[match(cand, m1$gene)]
cand_df$META1_padj <- m1$padj[match(cand, m1$gene)]
cand_df$META2_logFC <- m2$pooled_logFC[match(cand, m2$gene)]
cand_df$META2_padj <- m2$padj[match(cand, m2$gene)]
cand_df$kpy_module <- kpyrn$module[match(cand, kpyrn$gene_symbol)]
write.csv(cand_df, file.path(OUT, "ML_candidate_pool.csv"), row.names = FALSE)
cat(sprintf("\n候选池（ML 输入）: %d 基因 -> ML_candidate_pool.csv\n", nrow(cand_df)))

# ---- 图：模块-性状热图 ----
pdf(file.path(OUT, "module_trait_heatmap.pdf"), width = 8, height = 6)
par(mfrow = c(1, 1))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(datTraits),
               yLabels = colnames(MEs0), colorLabels = TRUE,
               main = "Module-trait relationships (GSE134431)")
dev.off()

# ---- 图：候选池 KpyRN 基因热图（GSE134431） ----
if (length(cand) >= 2 && requireNamespace("pheatmap", quietly = TRUE)) {
  sub <- expr[cand, ]
  sub <- sub - rowMeans(sub)
  annotation <- data.frame(Group = ann$group2,
                           Healing = ifelse(is.na(ann$healing), "NA", ann$healing),
                           row.names = rownames(ann))
  pheatmap::pheatmap(sub, annotation_col = annotation, show_colnames = FALSE,
                     filename = file.path(OUT, "candidate_heatmap.png"),
                     width = 10, height = 8, dpi = 150,
                     color = colorRampPalette(c("#2166AC", "white", "#C0392B"))(100),
                     main = "KpyRN candidate genes (GSE134431)")
} else {
  cat("[note] pheatmap 不可用或候选基因不足，跳过热图\n")
}
cat("WGCNA DONE\n")
