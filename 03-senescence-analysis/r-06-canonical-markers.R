# Step 6e: do the canonical senescence markers actually move? (critical for the framing)

suppressPackageStartupMessages(library(limma))

log_con <- file("codex_output/logs/r-06-canonical-markers.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 6e: canonical senescence markers across cohorts ===\n\n")

markers <- c("CDKN1A", "CDKN2A", "CDKN2B", "TP53", "LMNB1", "GLB1", "MKI67",
             "SERPINE1", "IL6", "CXCL8", "HMGB1", "GDF15", "FN1", "SIRT1",
             "SIRT6", "ATM", "TERF2", "RB1", "E2F1", "CCNA2", "CCND2", "PTEN")

cohorts <- list(
  list(tag = "GSE96804", e = "codex_output/data/GSE96804-expr.rds",
       g = "codex_output/data/GSE96804-group.csv", pos = "DN"),
  list(tag = "GSE30528", e = "codex_output/data/GSE30528-expr.rds",
       g = "codex_output/data/GSE30528-group.csv", pos = "DKD"),
  list(tag = "GSE99325", e = "codex_output/data/GSE99325-expr.rds",
       g = "codex_output/data/GSE99325-group.csv", pos = "DN")
)

out <- list()
for (co in cohorts) {
  expr <- readRDS(co$e)
  grp <- read.csv(co$g, stringsAsFactors = FALSE)
  keep <- grp$group %in% c("Control", co$pos)
  grp <- grp[keep, , drop = FALSE]
  expr <- expr[, grp$sample, drop = FALSE]

  # mean expression per group for every marker that is measured (before any filtering)
  present <- intersect(markers, rownames(expr))
  cat("---", co$tag, "| markers measured:", length(present), "of", length(markers), "\n")
  cat("    not measured:", paste(setdiff(markers, present), collapse = ", "), "\n")

  grp_f <- factor(grp$group, levels = c("Control", co$pos))
  sub <- expr[present, , drop = FALSE]
  sub <- sub[!apply(sub, 1, function(x) any(is.na(x))), , drop = FALSE]
  fit <- eBayes(lmFit(sub, model.matrix(~ grp_f)))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)

  m <- aggregate(t(sub), by = list(group = grp_f), FUN = mean)
  ctrl_mean <- as.numeric(m[m$group == "Control", -1])
  case_mean <- as.numeric(m[m$group == co$pos, -1])
  names(ctrl_mean) <- names(case_mean) <- colnames(m)[-1]

  for (g in tt$gene) {
    out[[length(out) + 1]] <- data.frame(
      dataset = co$tag, gene = g,
      mean_control = round(ctrl_mean[g], 2), mean_case = round(case_mean[g], 2),
      logFC = round(tt$logFC[tt$gene == g], 3),
      p_value = signif(tt$P.Value[tt$gene == g], 3),
      stringsAsFactors = FALSE
    )
  }
  print(tt[order(tt$P.Value), c("gene", "logFC", "P.Value", "adj.P.Val")], row.names = FALSE)
  cat("\n")
}

res <- do.call(rbind, out)
write.csv(res, "codex_output/results/canonical-senescence-markers.csv", row.names = FALSE)

cat("=== direction of the classical senescence hallmarks ===\n")
cat("expected in senescence: CDKN1A up, CDKN2A up, LMNB1 down, GLB1 up, MKI67 down, IL6 up\n\n")
wide <- reshape(res[, c("dataset", "gene", "logFC")], idvar = "gene",
                timevar = "dataset", direction = "wide")
wide$consistent_up <- apply(wide[, -1], 1, function(x) all(x > 0, na.rm = TRUE))
wide$consistent_down <- apply(wide[, -1], 1, function(x) all(x < 0, na.rm = TRUE))
print(wide, row.names = FALSE)
sink()
close(log_con)
