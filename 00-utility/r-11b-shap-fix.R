# Diagnose and fix the SHAP computation for the XGBoost model

suppressPackageStartupMessages({ library(xgboost); library(glmnet) })

core <- read.csv("codex_output/results/ml-consensus-genes.csv", stringsAsFactors = FALSE)$gene
e <- readRDS("codex_output/data/GSE96804-expr.rds")
g <- read.csv("codex_output/data/GSE96804-group.csv", stringsAsFactors = FALSE)
keep <- g$group %in% c("Control", "DN")
e <- e[, g$sample[keep]]; g <- g[keep, ]
xr <- t(scale(t(e[core, ])))
xr[is.na(xr)] <- 0
x <- t(xr)
y <- as.numeric(factor(g$group, levels = c("Control", "DN"))) - 1
cat("input:", dim(x), "| cases:", sum(y), "\n")

dm <- xgb.DMatrix(x, label = y, feature_names = core)
set.seed(2026)
params <- list(max_depth = 3, learning_rate = 0.05, objective = "binary:logistic", nthread = 1)
bst <- xgb.train(params, dm, nrounds = 300)

cat("\n--- importance (gain) ---\n")
imp <- xgb.importance(model = bst)
print(head(imp, 15))

cat("\n--- predcontrib via DMatrix ---\n")
sh <- predict(bst, dm, predcontrib = TRUE)
cat("dim:", paste(dim(sh), collapse = " x "), "\n")
cat("colnames head:", paste(head(colnames(sh), 5), collapse = ", "), "\n")

sh <- as.matrix(sh)
drop <- which(colnames(sh) %in% c("BIAS", "(Intercept)"))
if (length(drop)) sh <- sh[, -drop, drop = FALSE]
ma <- colMeans(abs(sh))
df <- data.frame(gene = names(ma), mean_abs_shap = round(as.numeric(ma), 4))
df <- df[order(-df$mean_abs_shap), ]
print(head(df, 20), row.names = FALSE)
write.csv(df, "codex_output/results/ml-shap-importance.csv", row.names = FALSE)
cat("\nnonzero features:", sum(df$mean_abs_shap > 0), "of", nrow(df), "\n")
