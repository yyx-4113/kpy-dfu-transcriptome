# 导出 ESTIMATE 签名到 CSV（供 Python 使用）
Sys.setlocale("LC_ALL", "Chinese_China.65001")
base <- Sys.getenv("KPYRN_PROJECT_ROOT")
x <- readRDS(file.path(base, "01_数据/基因集/免疫/estimate_genes.rds"))
sig <- x$signature_geneset
out <- data.frame(
  gene = as.character(unlist(sig[1, 3:ncol(sig)])),
  signature = "StromalSignature"
)
out2 <- data.frame(
  gene = as.character(unlist(sig[2, 3:ncol(sig)])),
  signature = "ImmuneSignature"
)
out <- rbind(out, out2)
write.csv(out, file.path(base, "01_数据/基因集/免疫/estimate_signatures.csv"),
          row.names = FALSE, quote = FALSE)
cat("exported:", nrow(out), "genes\n")
print(table(out$signature))
cat("OK-DONE\n")
