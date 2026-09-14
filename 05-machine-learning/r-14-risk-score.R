# Step 12: the final risk-score model
# Produces the paper-ready formula, the score distribution, the cut-off, calibration and
# decision-curve analysis for the 49-gene specificity-filtered signature.

suppressPackageStartupMessages({
  library(glmnet); library(pROC); library(ggplot2); library(reshape2)
})

dir.create("codex_output/results", showWarnings = FALSE, recursive = TRUE)
dir.create("codex_output/figures", showWarnings = FALSE, recursive = TRUE)

log_con <- file("codex_output/logs/r-14-risk-score.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 12: risk score model ===\n")
cat("date:", format(Sys.time()), "\n\n")
set.seed(2026)

sel <- read.csv("codex_output/results/specificity-filtered-genes.csv",
                stringsAsFactors = FALSE)$gene
cat("signature size:", length(sel), "genes\n\n")

load_cohort <- function(expr_path, grp_path, pos, label) {
  e <- readRDS(expr_path); g <- read.csv(grp_path, stringsAsFactors = FALSE)
  keep <- g$group %in% c("Control", pos)
  avail <- intersect(sel, rownames(e))
  m <- t(scale(t(e[avail, g$sample[keep], drop = FALSE])))
  m[is.na(m)] <- 0
  nm <- make.names(sel, unique = TRUE)
  full <- matrix(0, nrow = ncol(m), ncol = length(sel), dimnames = list(colnames(m), nm))
  full[, make.names(avail, unique = TRUE)] <- t(m)
  list(x = full, y = as.numeric(g$group[keep] == pos),
       label = label, samples = colnames(m))
}

tr  <- load_cohort("codex_output/data/GSE96804-expr.rds", "codex_output/data/GSE96804-group.csv", "DN",  "GSE96804 (training)")
va1 <- load_cohort("codex_output/data/GSE30528-expr.rds", "codex_output/data/GSE30528-group.csv", "DKD", "GSE30528 (validation)")
va2 <- load_cohort("codex_output/data/GSE99325-expr.rds", "codex_output/data/GSE99325-group.csv", "DN",  "GSE99325 (validation)")

# ---------- fit the penalised logistic model -------------------------------
cvfit <- cv.glmnet(tr$x, tr$y, family = "binomial", alpha = 1, nfolds = 10)
cf <- as.matrix(coef(cvfit, s = "lambda.1se"))
coefs <- cf[cf[, 1] != 0, , drop = FALSE]
cat("non-zero coefficients at lambda.1se:", nrow(coefs) - 1, "genes\n\n")
print(round(sort(coefs[-1, 1], decreasing = TRUE), 4))
write.csv(data.frame(gene = rownames(coefs), coefficient = round(coefs[, 1], 5)),
          "codex_output/results/risk-score-coefficients.csv", row.names = FALSE)

# ---------- paper-ready formula --------------------------------------------
beta <- coefs[-1, 1]
pos <- sort(beta[beta > 0], decreasing = TRUE)
neg <- sort(beta[beta < 0])
fmt <- function(v) sprintf("%+.4f", v)
formula_txt <- c(
  "RiskScore = ",
  paste(sprintf("%s x z(%s)", fmt(pos), names(pos)), collapse = "\n          + "),
  if (length(neg)) paste0("\n          ", paste(sprintf("%s x z(%s)", fmt(neg), names(neg)), collapse = "\n          ")) else NULL
)
formula_txt <- paste(formula_txt, collapse = "")
cat("\n=== risk score formula ===\n")
cat(formula_txt, "\n")
cat(sprintf("\nintercept = %.4f\n", coefs[1, 1]))
cat("z() denotes the gene expression value standardised within each cohort (z-score).\n")
cat("Probability = 1 / (1 + exp(-(intercept + RiskScore)))\n")

writeLines(c(
  "Specificity-filtered senescence risk score",
  paste0("Genes in the signature: ", length(sel)),
  paste0("Genes with non-zero coefficients: ", nrow(coefs) - 1),
  "",
  formula_txt,
  "",
  sprintf("Intercept = %.4f", coefs[1, 1]),
  "",
  "Notes:",
  "- z() is the expression value standardised across samples within a cohort (mean 0, SD 1).",
  "- The score was fitted on GSE96804 with LASSO-penalised logistic regression (lambda.1se).",
  "- Probability of diabetic nephropathy = 1 / (1 + exp(-(intercept + RiskScore)))."
), "codex_output/results/risk-score-formula.txt")

# ---------- apply the score to every cohort --------------------------------
score_of <- function(d) {
  z <- d$x[, names(beta), drop = FALSE]
  linpred <- as.numeric(z %*% beta)
  p <- 1 / (1 + exp(-(coefs[1, 1] + linpred)))
  list(score = linpred, prob = p, y = d$y)
}

sc <- list(train = score_of(tr), A = score_of(va1), B = score_of(va2))

# a second independent set: DKD versus other nephropathies (used for orientation only)
e <- readRDS("codex_output/data/GSE99325-expr.rds")
g <- read.csv("codex_output/data/GSE99325-group.csv", stringsAsFactors = FALSE)
sel2 <- g$group %in% c("DN", "Other")
avail2 <- intersect(sel, rownames(e))
m2 <- t(scale(t(e[avail2, g$sample[sel2], drop = FALSE]))); m2[is.na(m2)] <- 0
nm <- make.names(sel, unique = TRUE)
full2 <- matrix(0, nrow = ncol(m2), ncol = length(sel), dimnames = list(colnames(m2), nm))
full2[, make.names(avail2, unique = TRUE)] <- t(m2)
dko <- score_of(list(x = full2, y = as.numeric(g$group[sel2] == "DN")))

# ---------- cut-off from the training cohort (Youden index) ----------------
roc_tr <- pROC::roc(tr$y, sc$train$prob, quiet = TRUE)
best <- pROC::coords(roc_tr, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"))
cutoff <- as.numeric(best$threshold[1])
cat(sprintf("\nYouden cut-off on the training cohort: probability >= %.3f\n", cutoff))
cat(sprintf("training sensitivity = %.3f, specificity = %.3f\n", best$sensitivity[1], best$specificity[1]))

classify <- function(s, cutoff) {
  pred <- as.integer(s$prob >= cutoff)
  tp <- sum(pred == 1 & s$y == 1); tn <- sum(pred == 0 & s$y == 0)
  fp <- sum(pred == 1 & s$y == 0); fn <- sum(pred == 0 & s$y == 1)
  data.frame(n = length(s$y), tp = tp, tn = tn, fp = fp, fn = fn,
             accuracy = round((tp + tn) / length(s$y), 3),
             sensitivity = round(tp / max(1, tp + fn), 3),
             specificity = round(tn / max(1, tn + fp), 3),
             PPV = round(tp / max(1, tp + fp), 3),
             NPV = round(tn / max(1, tn + fn), 3))
}
perf <- rbind(
  data.frame(cohort = "GSE96804 (training)", classify(sc$train, cutoff)),
  data.frame(cohort = "GSE30528 (validation)", classify(sc$A, cutoff)),
  data.frame(cohort = "GSE99325 (validation)", classify(sc$B, cutoff)),
  data.frame(cohort = "GSE99325 (vs other nephropathy, circular)", classify(dko, cutoff))
)
cat("\n=== classification at the training cut-off ===\n")
print(perf, row.names = FALSE)
write.csv(perf, "codex_output/results/risk-score-performance.csv", row.names = FALSE)

aucs <- rbind(
  data.frame(cohort = "GSE96804 (training)", AUC = round(as.numeric(pROC::auc(roc_tr)), 3)),
  data.frame(cohort = "GSE30528 (validation)", AUC = round(as.numeric(pROC::auc(pROC::roc(va1$y, sc$A$prob, quiet = TRUE))), 3)),
  data.frame(cohort = "GSE99325 (validation)", AUC = round(as.numeric(pROC::auc(pROC::roc(va2$y, sc$B$prob, quiet = TRUE))), 3)),
  data.frame(cohort = "GSE99325 (vs other nephropathy, circular)", AUC = round(as.numeric(pROC::auc(pROC::roc(dko$y, dko$prob, quiet = TRUE))), 3))
)
cat("\nAUC:\n"); print(aucs, row.names = FALSE)

# ---------- calibration -----------------------------------------------------
calib <- function(s, label) {
  d <- data.frame(p = s$prob, y = s$y)
  d$bin <- cut(d$p, breaks = unique(quantile(d$p, probs = seq(0, 1, 0.2))),
               include.lowest = TRUE)
  agg <- aggregate(cbind(p, y) ~ bin, data = d, FUN = mean)
  agg$cohort <- label
  agg
}
cal <- do.call(rbind, list(calib(sc$train, "GSE96804 (training)"),
                           calib(sc$A, "GSE30528 (validation)"),
                           calib(sc$B, "GSE99325 (validation)")))
p_cal <- ggplot(cal, aes(x = p, y = y, color = cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(size = 2) + geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("#D7191C","#2C7BB6","#1A9641")) +
  labs(x = "predicted probability", y = "observed proportion of cases", color = NULL,
       title = "Calibration of the risk score") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave("codex_output/figures/fig-risk-calibration.png", p_cal, width = 6, height = 5, dpi = 300)

# ---------- decision curve analysis ----------------------------------------
dca <- function(s, label) {
  pts <- seq(0.05, 0.8, by = 0.01)
  n <- length(s$y); prev <- mean(s$y)
  nb_model <- sapply(pts, function(pt) {
    pred <- as.integer(s$prob >= pt)
    sum(pred == 1 & s$y == 1) / n - sum(pred == 1 & s$y == 0) / n * (pt / (1 - pt))
  })
  nb_all <- prev - (1 - prev) * (pts / (1 - pts))
  data.frame(threshold = rep(pts, 3),
             net_benefit = c(nb_model, nb_all, rep(0, length(pts))),
             strategy = rep(c("risk score", "treat all", "treat none"), each = length(pts)),
             cohort = label)
}
dca_df <- rbind(dca(sc$train, "GSE96804 (training)"),
                dca(sc$A, "GSE30528 (validation)"),
                dca(sc$B, "GSE99325 (validation)"))
p_dca <- ggplot(dca_df, aes(x = threshold, y = net_benefit, color = strategy)) +
  geom_line(linewidth = 0.7) + facet_wrap(~ cohort, nrow = 1) +
  scale_color_manual(values = c("risk score" = "#D7191C", "treat all" = "#2C7BB6", "treat none" = "grey40")) +
  labs(x = "threshold probability", y = "net benefit", color = NULL,
       title = "Decision curve analysis") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave("codex_output/figures/fig-risk-dca.png", p_dca, width = 11, height = 4.2, dpi = 300)

# ---------- score distribution ---------------------------------------------
dist <- rbind(
  data.frame(cohort = "GSE96804 (training)", prob = sc$train$prob, group = ifelse(sc$train$y == 1, "Case", "Control")),
  data.frame(cohort = "GSE30528 (validation)", prob = sc$A$prob, group = ifelse(sc$A$y == 1, "Case", "Control")),
  data.frame(cohort = "GSE99325 (validation)", prob = sc$B$prob, group = ifelse(sc$B$y == 1, "Case", "Control"))
)
dist$cohort <- factor(dist$cohort, levels = c("GSE96804 (training)", "GSE30528 (validation)", "GSE99325 (validation)"))
p_dist <- ggplot(dist, aes(x = group, y = prob, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
  geom_hline(yintercept = cutoff, linetype = "dashed", color = "grey30", linewidth = 0.4) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_fill_manual(values = c("Control" = "#9ecae1", "Case" = "#de2d26")) +
  labs(x = NULL, y = "predicted probability", fill = NULL,
       title = "Risk score distribution (dashed line = training cut-off)") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave("codex_output/figures/fig-risk-distribution.png", p_dist, width = 9, height = 4.2, dpi = 300)

saveRDS(list(coefs = coefs, cutoff = cutoff, scores = sc, perf = perf, aucs = aucs),
        "codex_output/data/risk-score-model.rds")
cat("\ndone\n")
sink()
close(log_con)
