# Compare the LASSO risk score with a ridge version: LASSO coefficients are extreme because
# the training cohort separates perfectly, which makes the probability cut-off untransferable.

suppressPackageStartupMessages({ library(glmnet); library(pROC) })

log_con <- file("codex_output/logs/r-14b-ridge.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 12b: LASSO vs ridge risk score ===\n\n")
set.seed(2026)

sel <- read.csv("codex_output/results/specificity-filtered-genes.csv", stringsAsFactors = FALSE)$gene

load_cohort <- function(expr_path, grp_path, pos) {
  e <- readRDS(expr_path); g <- read.csv(grp_path, stringsAsFactors = FALSE)
  keep <- g$group %in% c("Control", pos)
  avail <- intersect(sel, rownames(e))
  m <- t(scale(t(e[avail, g$sample[keep], drop = FALSE]))); m[is.na(m)] <- 0
  nm <- make.names(sel, unique = TRUE)
  full <- matrix(0, nrow = ncol(m), ncol = length(sel), dimnames = list(colnames(m), nm))
  full[, make.names(avail, unique = TRUE)] <- t(m)
  list(x = full, y = as.numeric(g$group[keep] == pos))
}
tr  <- load_cohort("codex_output/data/GSE96804-expr.rds", "codex_output/data/GSE96804-group.csv", "DN")
va1 <- load_cohort("codex_output/data/GSE30528-expr.rds", "codex_output/data/GSE30528-group.csv", "DKD")
va2 <- load_cohort("codex_output/data/GSE99325-expr.rds", "codex_output/data/GSE99325-group.csv", "DN")

eval_fit <- function(alpha, sname) {
  fit <- cv.glmnet(tr$x, tr$y, family = "binomial", alpha = alpha, nfolds = 10)
  cf <- as.matrix(coef(fit, s = "lambda.1se"))
  nz <- sum(cf[-1, 1] != 0)
  p_tr  <- as.numeric(predict(fit, tr$x,  s = "lambda.1se", type = "response"))
  p_va1 <- as.numeric(predict(fit, va1$x, s = "lambda.1se", type = "response"))
  p_va2 <- as.numeric(predict(fit, va2$x, s = "lambda.1se", type = "response"))
  data.frame(penalty = sname, n_nonzero = nz,
             max_abs_coef = round(max(abs(cf[-1, 1])), 2),
             AUC_training = round(as.numeric(pROC::auc(pROC::roc(tr$y,  p_tr,  quiet = TRUE))), 3),
             AUC_GSE30528 = round(as.numeric(pROC::auc(pROC::roc(va1$y, p_va1, quiet = TRUE))), 3),
             AUC_GSE99325 = round(as.numeric(pROC::auc(pROC::roc(va2$y, p_va2, quiet = TRUE))), 3))
}

res <- rbind(eval_fit(1, "LASSO (alpha=1)"), eval_fit(0, "Ridge (alpha=0)"))
print(res, row.names = FALSE)
write.csv(res, "codex_output/results/risk-score-lasso-vs-ridge.csv", row.names = FALSE)

cat("\n=== ridge coefficients (lambda.1se), sorted ===\n")
fit0 <- cv.glmnet(tr$x, tr$y, family = "binomial", alpha = 0, nfolds = 10)
cf0 <- as.matrix(coef(fit0, s = "lambda.1se"))
cf0 <- cf0[cf0[, 1] != 0, , drop = FALSE]
print(round(sort(cf0[-1, 1], decreasing = TRUE), 4))
write.csv(data.frame(gene = rownames(cf0), coefficient = round(cf0[, 1], 5)),
          "codex_output/results/risk-score-ridge-coefficients.csv", row.names = FALSE)

cat("\n=== ridge score: probability distribution per cohort ===\n")
p <- function(x) sprintf("%.3f", x)
for (nm in c("training","GSE30528","GSE99325")) {
  d <- switch(nm, training = list(x = tr$x, y = tr$y), GSE30528 = va1, GSE99325 = va2)
  pr <- as.numeric(predict(fit0, d$x, s = "lambda.1se", type = "response"))
  cat(sprintf("%-10s cases median %s | controls median %s\n",
              nm, p(median(pr[d$y == 1])), p(median(pr[d$y == 0]))))
}

cat("\n=== classification at 0.5 for the ridge model ===\n")
for (nm in c("training","GSE30528","GSE99325")) {
  d <- switch(nm, training = list(x = tr$x, y = tr$y), GSE30528 = va1, GSE99325 = va2)
  pr <- as.numeric(predict(fit0, d$x, s = "lambda.1se", type = "response"))
  pred <- as.integer(pr >= 0.5)
  tp <- sum(pred == 1 & d$y == 1); tn <- sum(pred == 0 & d$y == 0)
  fp <- sum(pred == 1 & d$y == 0); fn <- sum(pred == 0 & d$y == 1)
  cat(sprintf("%-10s sens %.2f | spec %.2f | acc %.2f\n",
              nm, tp/(tp+fn), tn/(tn+fp), (tp+tn)/length(d$y)))
}

cat("\ndone\n")
sink()
close(log_con)
