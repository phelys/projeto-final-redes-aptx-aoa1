# ============================================================
# 07_classification.R — Classificador APTX-KO x WT + MÉTRICAS
# Gera AUC, precisão, recall, F1 e matriz de confusão.
# Compara DOIS conjuntos de atributos:
#   (A) BASELINE: top DEGs (efeito do genótipo)  -> referência
#   (B) REDE:     eigengenes dos módulos WGCNA    -> proposto
# Validação cruzada leave-one-out (LOO) dado o n pequeno (n=12).
# ============================================================
source("R/config.R")
suppressMessages({ library(glmnet); library(caret); library(pROC); library(matrixStats) })
set.seed(SEED)

vsd     <- readRDS(file.path(DIR_PROCESSED, "vsd.rds"))
coldata <- readRDS(file.path(DIR_PROCESSED, "coldata.rds"))
y <- factor(coldata$genotype, levels = c("WT", "KO"))

# Conjunto A: top DEGs do efeito do genótipo (baseline)
deg <- read.csv(file.path(DIR_TABLES, "de_genotype_KO_vs_WT_degs.csv"))
top_deg <- intersect(head(deg$gene[order(deg$padj)], 50), rownames(vsd))
Xa <- t(vsd[top_deg, rownames(coldata), drop = FALSE])

# Conjunto B: eigengenes dos módulos (rede)
wg <- readRDS(file.path(DIR_PROCESSED, "wgcna.rds"))
Xb <- as.matrix(wg$MEs[rownames(coldata), , drop = FALSE])

loo_metrics <- function(X, y) {
  n <- nrow(X); probs <- numeric(n)
  for (i in seq_len(n)) {
    fit <- cv.glmnet(X[-i, , drop = FALSE], y[-i], family = "binomial",
                     alpha = 0.5, nfolds = min(5, n - 1))
    probs[i] <- predict(fit, X[i, , drop = FALSE], s = "lambda.min", type = "response")
  }
  pred <- factor(ifelse(probs >= 0.5, "KO", "WT"), levels = c("WT", "KO"))
  cm  <- caret::confusionMatrix(pred, y, positive = "KO")
  roc <- pROC::roc(response = y, predictor = probs, levels = c("WT", "KO"), quiet = TRUE)
  list(auc = as.numeric(pROC::auc(roc)),
       precision = cm$byClass["Precision"], recall = cm$byClass["Recall"],
       f1 = cm$byClass["F1"], confusion = cm$table, roc = roc)
}

res_A <- loo_metrics(Xa, y)
res_B <- loo_metrics(Xb, y)

metrics <- data.frame(
  modelo   = c("Baseline (DEGs)", "Rede (eigengenes WGCNA)"),
  AUC      = c(res_A$auc, res_B$auc),
  Precisao = c(res_A$precision, res_B$precision),
  Recall   = c(res_A$recall, res_B$recall),
  F1       = c(res_A$f1, res_B$f1)
)
write.csv(metrics, file.path(DIR_TABLES, "classificacao_metricas.csv"), row.names = FALSE)
print(metrics)

capture.output(
  cat("== Baseline (DEGs) ==\n"); print(res_A$confusion),
  cat("\n== Rede (eigengenes) ==\n"); print(res_B$confusion),
  file = file.path(DIR_TABLES, "matriz_confusao.txt")
)

png(file.path(DIR_FIGURES, "roc_comparativo.png"), width = 1100, height = 1000, res = 150)
plot(res_A$roc, col = "grey40", main = "ROC — APTX-KO x WT (LOO-CV)")
plot(res_B$roc, col = "firebrick", add = TRUE)
legend("bottomright",
       legend = c(sprintf("Baseline DEGs (AUC=%.2f)", res_A$auc),
                  sprintf("Rede eigengenes (AUC=%.2f)", res_B$auc)),
       col = c("grey40", "firebrick"), lwd = 2)
dev.off()

message("Classificação concluída. Métricas em results/tables/classificacao_metricas.csv")
