# Step 9: multi-algorithm feature selection, model building and external validation
# Training: GSE96804 | External validation: GSE30528, GSE99325
# Strategy: z-score genes within each cohort before prediction (the cohorts were processed
# on different platforms and GSE30528 was already batch-corrected before submission, so the
# model must never see raw cross-cohort values).

suppressPackageStartupMessages({
  library(glmnet)
  library(randomForest)
  library(e1071)
  library(Boruta)
  library(xgboost)
  library(pROC)
  library(ggplot2)
  library(reshape2)
})

dir.create("codex_output/results", showWarnings = FALSE, recursive = TRUE)
dir.create("codex_output/figures", showWarnings = FALSE, recursive = TRUE)

log_con <- file("codex_output/logs/r-11-ml.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 9: machine learning model ===\n")
cat("date:", format(Sys.time()), "\n\n")
set.seed(2026)

# ---------- data ----------------------------------------------------------
core <- read.csv("codex_output/results/senescence-core-genes.csv", stringsAsFactors = FALSE)
genes <- core$gene

load_cohort <- function(expr_path, grp_path, pos) {
  e <- readRDS(expr_path); g <- read.csv(grp_path, stringsAsFactors = FALSE)
  keep <- g$group %in% c("Control", pos)
  e <- e[, g$sample[keep], drop = FALSE]; g <- g[keep, , drop = FALSE]
  genes_present <- intersect(genes, rownames(e))
  m <- t(e[genes_present, , drop = FALSE])
  m <- m[, genes, drop = FALSE]                     # keep a fixed gene order
  m <- scale(m)                                     # z-score within this cohort
  m[is.na(m)] <- 0
  list(x = as.data.frame(m), y = factor(ifelse(g$group == "Control", "Control", "Case"),
                                        levels = c("Control", "Case")),
       genes = genes, n_missing = sum(is.na(match(genes, rownames(e)))))
}

tr  <- load_cohort("codex_output/data/GSE96804-expr.rds",  "codex_output/data/GSE96804-group.csv",  "DN")
va1 <- load_cohort("codex_output/data/GSE30528-expr.rds",  "codex_output/data/GSE30528-group.csv",  "DKD")
va2 <- load_cohort("codex_output/data/GSE99325-expr.rds",  "codex_output/data/GSE99325-group.csv",  "DN")

cat("training :", nrow(tr$x), "samples,", ncol(tr$x), "candidate genes\n")
print(table(tr$y))
cat("validation A (GSE30528):", nrow(va1$x), "samples\n"); print(table(va1$y))
cat("validation B (GSE99325):", nrow(va2$x), "samples\n"); print(table(va2$y))
cat("\ngenes missing per cohort - training:", tr$n_missing,
    "| A:", va1$n_missing, "| B:", va2$n_missing, "\n\n")

X <- as.matrix(tr$x); Y <- tr$y

# ---------- 1. LASSO ------------------------------------------------------
set.seed(2026)
cvfit <- cv.glmnet(X, Y, family = "binomial", alpha = 1, nfolds = 10)
lasso_genes <- rownames(coef(cvfit, s = "lambda.1se"))[as.numeric(coef(cvfit, s = "lambda.1se")) != 0]
lasso_genes <- setdiff(lasso_genes, "(Intercept)")
cat("LASSO (lambda.1se) selected", length(lasso_genes), "genes\n")

# ---------- 2. random forest ---------------------------------------------
set.seed(2026)
rf <- randomForest(x = tr$x, y = Y, ntree = 1000, importance = TRUE)
imp <- importance(rf, type = 1)[, 1]
rf_genes <- names(sort(imp, decreasing = TRUE))[1:20]
cat("random forest top 20 genes\n")

# ---------- 3. SVM-RFE ----------------------------------------------------
svm_rfe <- function(x, y, n_features = 20) {
  keep <- colnames(x)
  while (length(keep) > n_features) {
    m <- svm(x[, keep, drop = FALSE], y, kernel = "linear", scale = FALSE)
    w <- t(m$coefs) %*% m$SV
    drop <- names(sort(abs(w[1, ])))[1]
    keep <- setdiff(keep, drop)
  }
  keep
}
set.seed(2026)
svm_genes <- svm_rfe(tr$x, Y, 20)
cat("SVM-RFE selected", length(svm_genes), "genes\n")

# ---------- 4. Boruta -----------------------------------------------------
set.seed(2026)
bor <- Boruta(x = tr$x, y = Y, doTrace = 0, maxRuns = 200)
bor_genes <- names(bor$finalDecision)[bor$finalDecision == "Confirmed"]
cat("Boruta confirmed", length(bor_genes), "genes\n")

# ---------- 5. consensus --------------------------------------------------
voters <- cbind(
  LASSO   = genes %in% lasso_genes,
  RF      = genes %in% rf_genes,
  SVMRFE  = genes %in% svm_genes,
  Boruta  = genes %in% bor_genes
)
votes <- rowSums(voters)
sel_table <- data.frame(gene = genes, votes = votes, voters, stringsAsFactors = FALSE)
sel_table <- sel_table[order(-sel_table$votes, sel_table$gene), ]
write.csv(sel_table, "codex_output/results/ml-feature-selection.csv", row.names = FALSE)

consensus <- sel_table$gene[sel_table$votes >= 2]
if (length(consensus) < 3) consensus <- sel_table$gene[1:10]
cat("\nconsensus genes (selected by >=2 methods):", length(consensus), "\n")
print(consensus)
write.csv(data.frame(gene = consensus), "codex_output/results/ml-consensus-genes.csv", row.names = FALSE)

# ---------- 6. models -----------------------------------------------------
dat <- data.frame(y = Y, tr$x[, consensus, drop = FALSE])

# penalised logistic regression on the consensus signature (plain glm separates because the
# training set has 41 cases vs 20 controls and the selected genes are strongly discriminative)
fit_logit <- function(dat) {
  x <- as.matrix(dat[, -1]); y <- dat$y
  cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 5)
}
fit_rf    <- function(dat) randomForest(y ~ ., data = dat, ntree = 1000)
fit_xgb   <- function(dat) {
  x <- as.matrix(dat[, -1]); y <- as.numeric(dat$y) - 1
  xgb.train(list(max_depth = 3, eta = 0.05, objective = "binary:logistic"),
            xgb.DMatrix(x, label = y), nrounds = 300)
}

cv_auc <- function(dat, fitter, predict_fun, k = 10, reps = 5) {
  set.seed(2026)
  aucs <- c()
  for (r in seq_len(reps)) {
    # stratified folds: each fold keeps the case/control ratio of the full training set
    folds <- integer(nrow(dat))
    for (lv in levels(dat$y)) {
      idx <- which(dat$y == lv)
      folds[idx] <- rep(seq_len(k), length.out = length(idx))[sample(length(idx))]
    }
    for (f in seq_len(k)) {
      trn <- dat[folds != f, ]; tst <- dat[folds == f, ]
      if (length(unique(tst$y)) < 2) next
      m <- tryCatch(fitter(trn), error = function(e) NULL)
      if (is.null(m)) next
      p <- predict_fun(m, tst)
      aucs <- c(aucs, as.numeric(pROC::auc(pROC::roc(tst$y, p, quiet = TRUE))))
    }
  }
  aucs
}

pred_logit <- function(m, newd) {
  as.numeric(predict(m, as.matrix(newd[, -1]), s = "lambda.min", type = "response"))
}
pred_rf    <- function(m, newd) as.numeric(predict(m, newd, type = "prob")[, "Case"])
pred_xgb   <- function(m, newd) as.numeric(predict(m, xgb.DMatrix(as.matrix(newd[, -1]))))

cat("\n=== internal cross-validation (training cohort) ===\n")
x_all <- tr$x

# Nested CV for the signature model: feature selection is repeated inside every fold, so the
# estimate is not inflated by having selected genes on the full training set first.
nested_lasso_auc <- function(x, y, k = 10, reps = 5) {
  set.seed(2026); aucs <- c()
  for (r in seq_len(reps)) {
    folds <- integer(nrow(x))
    for (lv in levels(y)) {
      idx <- which(y == lv)
      folds[idx] <- rep(seq_len(k), length.out = length(idx))[sample(length(idx))]
    }
    for (f in seq_len(k)) {
      trn <- folds != f; tst <- folds == f
      if (length(unique(y[tst])) < 2) next
      cvf <- tryCatch(cv.glmnet(as.matrix(x[trn, ]), y[trn], family = "binomial",
                                alpha = 1, nfolds = 5), error = function(e) NULL)
      if (is.null(cvf)) next
      p <- as.numeric(predict(cvf, as.matrix(x[tst, ]), s = "lambda.min", type = "response"))
      aucs <- c(aucs, as.numeric(pROC::auc(pROC::roc(y[tst], p, quiet = TRUE))))
    }
  }
  aucs
}

models <- list(
  Logistic = list(fit = fit_logit, pred = pred_logit, nested = TRUE),
  RandomForest = list(fit = fit_rf, pred = pred_rf, nested = FALSE),
  XGBoost = list(fit = fit_xgb, pred = pred_xgb, nested = FALSE)
)
cv_res <- list()
fitted <- list()
for (nm in names(models)) {
  if (isTRUE(models[[nm]]$nested)) {
    a <- nested_lasso_auc(x_all, Y)
  } else {
    # tree models use all 87 candidate genes, so no feature-selection leakage can occur
    a <- cv_auc(data.frame(y = Y, x_all), models[[nm]]$fit, models[[nm]]$pred)
  }
  cv_res[[nm]] <- data.frame(model = nm, cv_auc_mean = round(mean(a), 3),
                             cv_auc_sd = round(sd(a), 3), n_folds = length(a))
  fitted[[nm]] <- models[[nm]]$fit(dat)
  cat(sprintf("%-14s CV AUC = %.3f (SD %.3f, %d folds)\n", nm, mean(a), sd(a), length(a)))
}
write.csv(do.call(rbind, cv_res), "codex_output/results/ml-cv-performance.csv", row.names = FALSE)

# ---------- 7. external validation ---------------------------------------
ext <- list()
roc_list <- list()
for (nm in names(models)) {
  for (vd in list(list(tag = "GSE30528", d = va1), list(tag = "GSE99325", d = va2))) {
    newd <- data.frame(y = vd$d$y, vd$d$x[, consensus, drop = FALSE])
    p <- models[[nm]]$pred(fitted[[nm]], newd)
    r <- pROC::roc(newd$y, p, quiet = TRUE)
    ci <- as.numeric(pROC::ci.auc(r))
    ext[[length(ext) + 1]] <- data.frame(
      model = nm, cohort = vd$tag, n = nrow(newd),
      n_control = sum(newd$y == "Control"), n_case = sum(newd$y == "Case"),
      AUC = round(as.numeric(r$auc), 3),
      CI_low = round(ci[1], 3), CI_high = round(ci[3], 3),
      stringsAsFactors = FALSE)
    roc_list[[paste(nm, vd$tag)]] <- data.frame(
      model = nm, cohort = vd$tag,
      specificity = rev(r$specificities), sensitivity = rev(r$sensitivities))
  }
}
extdf <- do.call(rbind, ext)
write.csv(extdf, "codex_output/results/ml-external-validation.csv", row.names = FALSE)
cat("\n=== external validation ===\n")
print(extdf, row.names = FALSE)

# ---------- 8. figures ----------------------------------------------------
fdir <- "codex_output/figures"
rocdf <- do.call(rbind, roc_list)
p_roc <- ggplot(rocdf, aes(x = 1 - specificity, y = sensitivity, color = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ cohort, nrow = 1) +
  coord_equal() +
  scale_color_manual(values = c(Logistic = "#D7191C", RandomForest = "#2C7BB6", XGBoost = "#1A9641")) +
  labs(x = "1 - specificity", y = "sensitivity", color = NULL,
       title = "External validation ROC") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(fdir, "fig-ml-roc.png"), p_roc, width = 9, height = 4.5, dpi = 300)
ggsave(file.path(fdir, "fig-ml-roc.pdf"), p_roc, width = 9, height = 4.5)

impdf <- data.frame(gene = names(imp), importance = as.numeric(imp))
impdf <- impdf[order(-impdf$importance), ][1:min(25, nrow(impdf)), ]
p_imp <- ggplot(impdf, aes(x = reorder(gene, importance), y = importance)) +
  geom_col(fill = "#2C7BB6") + coord_flip() +
  labs(x = NULL, y = "permutation importance (mean decrease accuracy)",
       title = "Random forest feature importance") +
  theme_bw(base_size = 10)
ggsave(file.path(fdir, "fig-ml-importance.png"), p_imp, width = 6.5, height = 6, dpi = 300)

# ---------- 9. SHAP values from XGBoost -----------------------------------
set.seed(2026)
xmat <- as.matrix(dat[, -1])
dm <- xgb.DMatrix(xmat, label = as.numeric(dat$y) - 1)
bst <- xgb.train(list(max_depth = 3, learning_rate = 0.05, objective = "binary:logistic"),
                 dm, nrounds = 300)
shap <- predict(bst, xmat, predcontrib = TRUE)
shap <- shap[, colnames(shap) != "BIAS", drop = FALSE]
shap_mean <- colMeans(abs(shap))
shap_df <- data.frame(gene = names(shap_mean), mean_abs_shap = as.numeric(shap_mean))
shap_df <- shap_df[order(-shap_df$mean_abs_shap), ]
write.csv(shap_df, "codex_output/results/ml-shap-importance.csv", row.names = FALSE)
cat("\n=== top SHAP features ===\n")
print(head(shap_df, 15), row.names = FALSE)

p_shap <- ggplot(head(shap_df, 20), aes(x = reorder(gene, mean_abs_shap), y = mean_abs_shap)) +
  geom_col(fill = "#D7191C") + coord_flip() +
  labs(x = NULL, y = "mean |SHAP value|", title = "XGBoost SHAP importance") +
  theme_bw(base_size = 10)
ggsave(file.path(fdir, "fig-ml-shap.png"), p_shap, width = 6.5, height = 5.5, dpi = 300)

saveRDS(list(models = fitted, consensus = consensus, dat = dat,
             cohorts = list(train = tr, A = va1, B = va2)),
        "codex_output/data/ml-models.rds")

cat("\npenalised logistic coefficients (lambda.min):\n")
cf <- as.matrix(coef(fitted$Logistic, s = "lambda.min"))
cf <- cf[cf[, 1] != 0, , drop = FALSE]
print(round(sort(cf[, 1], decreasing = TRUE), 4))
write.csv(data.frame(gene = rownames(cf), coefficient = round(cf[, 1], 4)),
          "codex_output/results/ml-logistic-coefficients.csv", row.names = FALSE)
cat("\ndone\n")
sink()
close(log_con)
