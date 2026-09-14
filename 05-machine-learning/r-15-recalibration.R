# Step 13: probability recalibration
# The ridge score discriminates well across cohorts but its absolute probabilities are shifted:
# the training cohort (tumour-nephrectomy controls) is easy to separate, so probabilities in
# other cohorts are inflated. Here we quantify the shift with the standard calibration
# framework (calibration-in-the-large and calibration slope) and show how a simple intercept
# correction restores usable probabilities.

suppressPackageStartupMessages({
  library(glmnet); library(pROC); library(ggplot2)
})

log_con <- file("codex_output/logs/r-15-recalibration.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 13: probability recalibration ===\n")
cat("date:", format(Sys.time()), "\n\n")
set.seed(2026)

sel <- read.csv("codex_output/results/specificity-filtered-genes.csv", stringsAsFactors = FALSE)$gene
nm <- make.names(sel, unique = TRUE)
model <- readRDS("codex_output/data/final-risk-model.rds")
beta <- model$beta; intercept <- model$intercept

build <- function(expr_path, grp_path, groups_used, case_label, tag) {
  e <- readRDS(expr_path); g <- read.csv(grp_path, stringsAsFactors = FALSE)
  keep <- g$group %in% groups_used
  avail <- intersect(sel, rownames(e))
  m <- t(scale(t(e[avail, g$sample[keep], drop = FALSE]))); m[is.na(m)] <- 0
  full <- matrix(0, nrow = ncol(m), ncol = length(sel), dimnames = list(colnames(m), nm))
  full[, make.names(avail, unique = TRUE)] <- t(m)
  lp <- as.numeric(full[, names(beta), drop = FALSE] %*% beta) + intercept
  list(lp = lp, p = 1 / (1 + exp(-lp)),
       y = as.numeric(g$group[keep] == case_label), tag = tag)
}

coh <- list(
  build("codex_output/data/GSE96804-expr.rds", "codex_output/data/GSE96804-group.csv", c("Control","DN"),  "DN",  "GSE96804 (training)"),
  build("codex_output/data/GSE30528-expr.rds", "codex_output/data/GSE30528-group.csv", c("Control","DKD"), "DKD", "GSE30528 (validation)"),
  build("codex_output/data/GSE99325-expr.rds", "codex_output/data/GSE99325-group.csv", c("Control","DN"),  "DN",  "GSE99325 (validation)")
)
names(coh) <- vapply(coh, function(z) z$tag, character(1))

# ---------- 1. how far off is the calibration? ------------------------------
cat("=== calibration-in-the-large and calibration slope (logit(y) ~ a + b*logit(p)) ===\n")
cal_rows <- list()
for (k in names(coh)) {
  d <- coh[[k]]
  lp <- qlogis(pmin(pmax(d$p, 1e-6), 1 - 1e-6))
  fit <- glm(d$y ~ lp, family = binomial)
  a <- coef(fit)[1]; b <- coef(fit)[2]
  brier <- mean((d$p - d$y)^2)
  cal_rows[[k]] <- data.frame(cohort = k, n = length(d$y),
                              calibration_intercept = round(a, 3),
                              calibration_slope = round(b, 3),
                              mean_predicted = round(mean(d$p), 3),
                              observed_rate = round(mean(d$y), 3),
                              Brier = round(brier, 3),
                              AUC = round(as.numeric(pROC::auc(pROC::roc(d$y, d$p, quiet = TRUE))), 3))
}
cal <- do.call(rbind, cal_rows)
print(cal, row.names = FALSE)
write.csv(cal, "codex_output/results/recalibration-diagnostics.csv", row.names = FALSE)

cat("\nInterpretation:\n")
cat("- calibration intercept < 0 means the model over-predicts risk in that cohort\n")
cat("- calibration slope close to 1 means the spread of predictions is about right\n")
cat("- recalibration can fix the intercept; it cannot create discrimination\n\n")

# ---------- 2. cross-validated recalibration within each cohort -------------
# Recalibration parameters must never be estimated on the same data used to evaluate them,
# so they are fitted out-of-fold (5-fold) inside each cohort.
recalibrate_cv <- function(d, k = 5, slope_free = TRUE) {
  folds <- integer(length(d$y))
  for (lv in c(0, 1)) {
    idx <- which(d$y == lv)
    folds[idx] <- rep(seq_len(k), length.out = length(idx))[sample(length(idx))]
  }
  p_new <- rep(NA_real_, length(d$y))
  for (f in seq_len(k)) {
    tr <- folds != f; te <- folds == f
    lp_tr <- qlogis(pmin(pmax(d$p[tr], 1e-6), 1 - 1e-6))
    if (slope_free) {
      fit <- glm(d$y[tr] ~ lp_tr, family = binomial)
      p_new[te] <- plogis(coef(fit)[1] + coef(fit)[2] *
                            qlogis(pmin(pmax(d$p[te], 1e-6), 1 - 1e-6)))
    } else {
      # recalibration-in-the-large only: keep the slope at 1, shift the intercept
      off <- qlogis(mean(d$y[tr])) - mean(qlogis(pmin(pmax(d$p[tr], 1e-6), 1 - 1e-6)))
      p_new[te] <- plogis(qlogis(pmin(pmax(d$p[te], 1e-6), 1 - 1e-6)) + off)
    }
  }
  p_new
}

cat("=== before vs after cross-validated recalibration ===\n")
rows <- list()
for (k in names(coh)) {
  d <- coh[[k]]
  # the training cohort separates perfectly, so out-of-fold recalibration is undefined there
  if (grepl("training", k)) {
    p_large <- d$p
  } else {
    p_large <- recalibrate_cv(d, slope_free = FALSE)
  }
  p_free  <- recalibrate_cv(d, slope_free = TRUE)
  cls <- function(p) {
    pred <- as.integer(p >= 0.5)
    tp <- sum(pred == 1 & d$y == 1); tn <- sum(pred == 0 & d$y == 0)
    fp <- sum(pred == 1 & d$y == 0); fn <- sum(pred == 0 & d$y == 1)
    c(sens = tp / max(1, tp + fn), spec = tn / max(1, tn + fp), acc = (tp + tn) / length(d$y))
  }
  cm <- rbind(original = cls(d$p), intercept_only = cls(p_large), slope_free = cls(p_free))
  rows[[k]] <- data.frame(
    cohort = k,
    Brier_original = round(mean((d$p - d$y)^2), 3),
    Brier_intercept_only = round(mean((p_large - d$y)^2), 3),
    Brier_slope_free = round(mean((p_free - d$y)^2), 3),
    mean_pred_original = round(mean(d$p), 3),
    mean_pred_intercept_only = round(mean(p_large), 3),
    mean_pred_slope_free = round(mean(p_free), 3),
    sens_original = round(cm["original", "sens"], 2),
    spec_original = round(cm["original", "spec"], 2),
    acc_original = round(cm["original", "acc"], 2),
    sens_intercept_only = round(cm["intercept_only", "sens"], 2),
    spec_intercept_only = round(cm["intercept_only", "spec"], 2),
    acc_intercept_only = round(cm["intercept_only", "acc"], 2),
    sens_slope_free = round(cm["slope_free", "sens"], 2),
    spec_slope_free = round(cm["slope_free", "spec"], 2),
    acc_slope_free = round(cm["slope_free", "acc"], 2),
    stringsAsFactors = FALSE)
  d$p_recal <- p_large; coh[[k]] <- d
}
res <- do.call(rbind, rows)
print(res, row.names = FALSE)
write.csv(res, "codex_output/results/recalibration-performance.csv", row.names = FALSE)

# ---------- 3. calibration curves before and after --------------------------
calib_curve <- function(d, which_p, label) {
  p <- d[[which_p]]
  bins <- unique(quantile(p, probs = seq(0, 1, 0.25)))
  if (length(bins) < 3) return(NULL)
  g <- cut(p, breaks = bins, include.lowest = TRUE)
  agg <- aggregate(cbind(p, y) ~ g, data = data.frame(p = p, y = d$y), FUN = mean)
  data.frame(mean_predicted = agg$p, observed = agg$y, cohort = label,
             version = ifelse(which_p == "p", "before", "after"))
}
cur <- do.call(rbind, unlist(lapply(names(coh), function(k) {
  d <- coh[[k]]
  list(calib_curve(d, "p", k), calib_curve(d, "p_recal", k))
}), recursive = FALSE))

p_cal <- ggplot(cur, aes(x = mean_predicted, y = observed,
                         color = version, group = version)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(size = 1.8) + geom_line(linewidth = 0.6) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_color_manual(values = c(before = "#D7191C", after = "#2C7BB6")) +
  labs(x = "mean predicted probability", y = "observed proportion",
       color = NULL, title = "Calibration before and after intercept recalibration",
       subtitle = "recalibration fitted out-of-fold (5-fold CV within each cohort)") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave("codex_output/figures/fig-recalibration-curve.png", p_cal, width = 11, height = 4.4, dpi = 300)

# ---------- 4. distribution before and after --------------------------------
dist <- do.call(rbind, lapply(names(coh), function(k) {
  d <- coh[[k]]
  rbind(
    data.frame(cohort = k, version = "before", p = d$p,  group = ifelse(d$y == 1, "Case", "Control")),
    data.frame(cohort = k, version = "after",  p = d$p_recal, group = ifelse(d$y == 1, "Case", "Control")))
}))
dist$cohort <- factor(dist$cohort, levels = c("GSE96804 (training)","GSE30528 (validation)","GSE99325 (validation)"))
p_dist <- ggplot(dist, aes(x = group, y = p, fill = group)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30", linewidth = 0.4) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  facet_grid(version ~ cohort) +
  scale_fill_manual(values = c("Control" = "#9ecae1", "Case" = "#de2d26")) +
  labs(x = NULL, y = "probability", fill = NULL,
       title = "Probability distribution before and after recalibration") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave("codex_output/figures/fig-recalibration-distribution.png", p_dist, width = 10, height = 6, dpi = 300)

saveRDS(coh, "codex_output/data/recalibrated-cohorts.rds")
cat("\ndone\n")
sink()
close(log_con)
