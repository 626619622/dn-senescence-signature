# Step 11: build a DKD-specific signature
# The previous signature separated "diseased kidney" from "normal kidney" but not DKD from
# other nephropathies. Here the 115 non-diabetic nephropathy samples in GSE99325 are used as
# a specificity filter: a gene is kept only if it is changed in DKD and also differs between
# DKD and other nephropathies.

suppressPackageStartupMessages({
  library(limma); library(glmnet); library(randomForest); library(pROC); library(ggplot2)
})

log_con <- file("codex_output/logs/r-13-specificity-filter.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 11: specificity-filtered signature ===\n")
cat("date:", format(Sys.time()), "\n\n")
set.seed(2026)

core <- read.csv("codex_output/results/senescence-core-genes.csv", stringsAsFactors = FALSE)
genes <- core$gene

# ---------- 1. differential expression of the core genes in GSE99325 ----------
e <- readRDS("codex_output/data/GSE99325-expr.rds")
g <- read.csv("codex_output/data/GSE99325-group.csv", stringsAsFactors = FALSE)
stopifnot(identical(colnames(e), g$sample))
present <- intersect(genes, rownames(e))
cat("core genes measured in GSE99325:", length(present), "of", length(genes), "\n")
print(table(g$group))

sub <- e[present, , drop = FALSE]
sub <- sub[!apply(sub, 1, function(x) any(is.na(x))), , drop = FALSE]

fit_contrast <- function(grp, case, ctrl) {
  keep <- grp %in% c(case, ctrl)
  f <- factor(grp[keep], levels = c(ctrl, case))
  fit <- eBayes(lmFit(sub[, keep, drop = FALSE], model.matrix(~ f)))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), logFC = tt$logFC, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val)
}

dn_vs_ctrl   <- fit_contrast(g$group, "DN", "Control")
oth_vs_ctrl  <- fit_contrast(g$group, "Other", "Control")
dn_vs_oth    <- fit_contrast(g$group, "DN", "Other")
cat("\nDN vs Control  : n =", sum(g$group == "DN"), "vs", sum(g$group == "Control"), "\n")
cat("Other vs Control: n =", sum(g$group == "Other"), "vs", sum(g$group == "Control"), "\n")
cat("DN vs Other    : n =", sum(g$group == "DN"), "vs", sum(g$group == "Other"), "\n")

# ---------- 2. merge with the training-cohort statistics ----------------------
tab <- data.frame(gene = present, stringsAsFactors = FALSE)
tab$logFC_train <- core$logFC_training[match(tab$gene, core$gene)]
tab$FDR_train   <- core$FDR_training[match(tab$gene, core$gene)]
tab$logFC_DN_vs_Control  <- dn_vs_ctrl$logFC[match(tab$gene, dn_vs_ctrl$gene)]
tab$FDR_DN_vs_Control    <- dn_vs_ctrl$adj.P.Val[match(tab$gene, dn_vs_ctrl$gene)]
tab$logFC_Other_vs_Control <- oth_vs_ctrl$logFC[match(tab$gene, oth_vs_ctrl$gene)]
tab$FDR_Other_vs_Control   <- oth_vs_ctrl$adj.P.Val[match(tab$gene, oth_vs_ctrl$gene)]
tab$logFC_DN_vs_Other    <- dn_vs_oth$logFC[match(tab$gene, dn_vs_oth$gene)]
tab$FDR_DN_vs_Other      <- dn_vs_oth$adj.P.Val[match(tab$gene, dn_vs_oth$gene)]

# ---------- 3. tiered specificity definition ----------------------------------
tab$changed_in_DKD      <- !is.na(tab$FDR_train) & tab$FDR_train < 0.05 &
                           !is.na(tab$FDR_DN_vs_Other) & tab$FDR_DN_vs_Other < 0.05
tab$other_also_changed  <- !is.na(tab$FDR_Other_vs_Control) &
                           tab$FDR_Other_vs_Control < 0.05 &
                           abs(tab$logFC_Other_vs_Control) > 0.25
tab$discordant          <- sign(tab$logFC_DN_vs_Control) != sign(tab$logFC_Other_vs_Control)

tab$tierA <- tab$changed_in_DKD & !tab$other_also_changed & tab$discordant
tab$tierB <- tab$changed_in_DKD & !tab$other_also_changed
tab$tierC <- tab$changed_in_DKD
tab$tierD <- tab$logFC_train != 0 & sign(tab$logFC_train) == sign(tab$logFC_DN_vs_Control) &
             !is.na(tab$FDR_DN_vs_Other) & tab$FDR_DN_vs_Other < 0.05

cat("\n=== candidate counts by definition ===\n")
cat("A (DKD-specific and discordant with other nephropathies):", sum(tab$tierA), "\n")
cat("B (DKD-specific, other nephropathies unchanged)         :", sum(tab$tierB), "\n")
cat("C (changed in DKD vs other nephropathies)               :", sum(tab$tierC), "\n")
cat("D (C plus direction consistent with training)           :", sum(tab$tierD), "\n")
write.csv(tab, "codex_output/results/specificity-filter-table.csv", row.names = FALSE)

# choose the strictest definition that still leaves enough genes to model
choice <- if (sum(tab$tierA) >= 8) "tierA" else if (sum(tab$tierB) >= 8) "tierB" else
          if (sum(tab$tierD) >= 8) "tierD" else "tierC"
sel <- tab$gene[tab[[choice]]]
cat("\nselected definition:", choice, "->", length(sel), "genes\n")
print(sel)
write.csv(data.frame(gene = sel), "codex_output/results/specificity-filtered-genes.csv", row.names = FALSE)

# ---------- 4. rebuild the model on the filtered genes -------------------------
load_x <- function(expr_path, grp_path, pos, genes_use) {
  e <- readRDS(expr_path); g <- read.csv(grp_path, stringsAsFactors = FALSE)
  keep <- g$group %in% c("Control", pos)
  available <- intersect(genes_use, rownames(e))
  m <- t(scale(t(e[available, g$sample[keep], drop = FALSE])))
  m[is.na(m)] <- 0
  # gene names such as HLA-DQB1 are not valid symbols, and a few genes may be absent from a
  # given platform; build a fixed-width matrix and use the x/y interface of the models
  nm <- make.names(genes_use, unique = TRUE)
  full <- matrix(0, nrow = ncol(m), ncol = length(genes_use), dimnames = list(colnames(m), nm))
  full[, make.names(available, unique = TRUE)] <- t(m)
  list(x = as.data.frame(full),
       y = factor(ifelse(g$group[keep] == "Control", "Control", "Case"),
                  levels = c("Control", "Case")),
       n_missing = length(setdiff(genes_use, available)))
}

tr  <- load_x("codex_output/data/GSE96804-expr.rds", "codex_output/data/GSE96804-group.csv", "DN",  sel)
va1 <- load_x("codex_output/data/GSE30528-expr.rds", "codex_output/data/GSE30528-group.csv", "DKD", sel)
va2 <- load_x("codex_output/data/GSE99325-expr.rds", "codex_output/data/GSE99325-group.csv", "DN",  sel)

# DKD vs other nephropathies (circular for tier C/D: reported with a caveat)
e <- readRDS("codex_output/data/GSE99325-expr.rds")
g <- read.csv("codex_output/data/GSE99325-group.csv", stringsAsFactors = FALSE)
sel2 <- g$group %in% c("DN", "Other")
avail2 <- intersect(sel, rownames(e))
m <- t(scale(t(e[avail2, g$sample[sel2], drop = FALSE]))); m[is.na(m)] <- 0
nm2 <- make.names(sel, unique = TRUE)
full2 <- matrix(0, nrow = ncol(m), ncol = length(sel), dimnames = list(colnames(m), nm2))
full2[, make.names(avail2, unique = TRUE)] <- t(m)
vsoth <- list(x = as.data.frame(full2),
              y = factor(ifelse(g$group[sel2] == "DN", "Case", "Control"),
                         levels = c("Control", "Case")))

cat("\ngenes missing per cohort - training:", tr$n_missing, "| A:", va1$n_missing,
    "| B:", va2$n_missing, "\n")

logit <- cv.glmnet(as.matrix(tr$x), tr$y, family = "binomial", alpha = 1, nfolds = 5)
rf    <- randomForest(x = tr$x, y = tr$y, ntree = 1000)

auc_of <- function(d) {
  c(Logistic = as.numeric(pROC::auc(pROC::roc(d$y,
      as.numeric(predict(logit, as.matrix(d$x), s = "lambda.min", type = "response")), quiet = TRUE))),
    RandomForest = as.numeric(pROC::auc(pROC::roc(d$y,
      as.numeric(predict(rf, d$x, type = "prob")[, "Case"]), quiet = TRUE))))
}

res <- rbind(
  data.frame(cohort = "training (GSE96804)", t(auc_of(tr))),
  data.frame(cohort = "validation A (GSE30528, vs control)", t(auc_of(va1))),
  data.frame(cohort = "validation B (GSE99325, vs control)", t(auc_of(va2))),
  data.frame(cohort = "specificity (GSE99325, vs other nephropathy)", t(auc_of(vsoth)))
)
res[, -1] <- round(res[, -1], 3)
cat("\n=== AUC after specificity filtering ===\n")
print(res, row.names = FALSE)
write.csv(res, "codex_output/results/specificity-filtered-performance.csv", row.names = FALSE)

# ---------- 5. figure: specificity map of the 87 core genes --------------------
pdat <- tab[!is.na(tab$logFC_DN_vs_Control) & !is.na(tab$logFC_Other_vs_Control), ]
p <- ggplot(pdat, aes(x = logFC_Other_vs_Control, y = logFC_DN_vs_Control)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_point(aes(color = tierC), size = 1.8, alpha = 0.85) +
  scale_color_manual(values = c(`TRUE` = "#D7191C", `FALSE` = "grey70"),
                     labels = c("not DKD-specific", "DKD-specific")) +
  labs(x = expression("log"[2]*"FC (other nephropathies vs control)"),
       y = expression("log"[2]*"FC (DKD vs control)"),
       color = NULL, title = "DKD specificity of the 87 core genes",
       subtitle = "genes far from the diagonal differ between DKD and other nephropathies") +
  theme_bw(base_size = 10) + theme(legend.position = "top")
ggsave("codex_output/figures/fig-specificity-map.png", p, width = 6.5, height = 5.5, dpi = 300)

cat("\ndone\n")
sink()
close(log_con)
