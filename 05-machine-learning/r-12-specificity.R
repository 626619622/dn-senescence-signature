# Step 10: two robustness checks
#   (a) can the signature tell DKD apart from OTHER kidney diseases, not just healthy tissue?
#   (b) how much of the performance comes from a single gene (FOS)?

suppressPackageStartupMessages({
  library(glmnet); library(randomForest); library(xgboost); library(pROC)
})

log_con <- file("codex_output/logs/r-12-specificity.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 10: specificity and single-gene checks ===\n")
cat("date:", format(Sys.time()), "\n\n")
set.seed(2026)

core <- read.csv("codex_output/results/ml-consensus-genes.csv", stringsAsFactors = FALSE)$gene

prep <- function(expr_path, grp_path, case_label) {
  e <- readRDS(expr_path); g <- read.csv(grp_path, stringsAsFactors = FALSE)
  m <- t(scale(t(e[core, , drop = FALSE])))
  m[is.na(m)] <- 0
  list(x = as.data.frame(t(m)),
       y = factor(ifelse(g$group == case_label, "Case", "Control"),
                  levels = c("Control", "Case")),
       group = g$group)
}

tr  <- prep("codex_output/data/GSE96804-expr.rds", "codex_output/data/GSE96804-group.csv", "DN")
va1 <- prep("codex_output/data/GSE30528-expr.rds", "codex_output/data/GSE30528-group.csv", "DKD")
va2 <- prep("codex_output/data/GSE99325-expr.rds", "codex_output/data/GSE99325-group.csv", "DN")

# ---------- (a) DKD vs other nephropathies in GSE99325 --------------------
e <- readRDS("codex_output/data/GSE99325-expr.rds")
g <- read.csv("codex_output/data/GSE99325-group.csv", stringsAsFactors = FALSE)
sel <- g$group %in% c("DN", "Other")
m <- t(scale(t(e[core, g$sample[sel], drop = FALSE])))
m[is.na(m)] <- 0
dko <- list(x = as.data.frame(t(m)),
            y = factor(ifelse(g$group[sel] == "DN", "Case", "Control"),
                       levels = c("Control", "Case")))
cat("specificity test: DN", sum(dko$y == "Case"), "vs other nephropathies",
    sum(dko$y == "Control"), "\n\n")

# ---------- fit the models on the training cohort -------------------------
dat <- data.frame(y = tr$y, tr$x)
x_tr <- as.matrix(tr$x); y_tr <- tr$y

logit_fit <- cv.glmnet(x_tr, y_tr, family = "binomial", alpha = 1, nfolds = 5)
rf_fit <- randomForest(y ~ ., data = dat, ntree = 1000)
xgb_fit <- xgb.train(list(max_depth = 3, learning_rate = 0.05, objective = "binary:logistic"),
                     xgb.DMatrix(x_tr, label = as.numeric(y_tr) - 1), nrounds = 300)

predict_all <- function(d) {
  list(
    Logistic = as.numeric(predict(logit_fit, as.matrix(d$x), s = "lambda.min", type = "response")),
    RandomForest = as.numeric(predict(rf_fit, d$x, type = "prob")[, "Case"]),
    XGBoost = as.numeric(predict(xgb_fit, xgb.DMatrix(as.matrix(d$x)))),
    FOS_only = as.numeric(scale(d$x$FOS))
  )
}

res <- list()
for (nm in c("GSE30528","GSE99325","GSE99325_DKD_vs_other")) {
  d <- switch(nm, GSE30528 = va1, GSE99325 = va2, GSE99325_DKD_vs_other = dko)
  ps <- predict_all(d)
  for (mdl in names(ps)) {
    r <- pROC::roc(d$y, ps[[mdl]], quiet = TRUE, direction = "<")
    ci <- as.numeric(pROC::ci.auc(r))
    res[[length(res) + 1]] <- data.frame(
      cohort = nm, model = mdl, n_control = sum(d$y == "Control"), n_case = sum(d$y == "Case"),
      AUC = round(as.numeric(r$auc), 3),
      CI_low = round(ci[1], 3), CI_high = round(ci[3], 3), stringsAsFactors = FALSE)
  }
}
rdf <- do.call(rbind, res)
write.csv(rdf, "codex_output/results/ml-specificity.csv", row.names = FALSE)
cat("=== AUC by cohort and model ===\n")
print(rdf, row.names = FALSE)

cat("\n=== note ===\n")
cat("If the signature also separates DKD from other nephropathies it carries disease-specific\n")
cat("information rather than merely detecting abnormal kidney tissue. The FOS-only row shows\n")
cat("how much of the performance is attributable to a single gene.\n")

r_fos <- pROC::roc(tr$y, tr$x$FOS, quiet = TRUE, direction = "<")
cat(sprintf("\nFOS alone in the training cohort: AUC = %.3f\n", as.numeric(r_fos$auc)))

sink()
close(log_con)
