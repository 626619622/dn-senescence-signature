# Step 6c: senescence pathway scoring (ssGSEA) across the training and validation cohorts

suppressPackageStartupMessages({
  library(GSVA)
  library(msigdbr)
  library(ggplot2)
  library(ggpubr)
  library(reshape2)
})

dir.create("codex_output/results", showWarnings = FALSE, recursive = TRUE)
dir.create("codex_output/figures", showWarnings = FALSE, recursive = TRUE)

log_con <- file("codex_output/logs/r-04-senescence-score.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 6c: senescence pathway scoring ===\n")
cat("date:", format(Sys.time()), "\n\n")

sets_wanted <- c("SAUL_SEN_MAYO", "REACTOME_CELLULAR_SENESCENCE",
                 "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
                 "GOBP_CELLULAR_SENESCENCE", "FRIDMAN_SENESCENCE_UP",
                 "RODWELL_AGING_KIDNEY_UP")
msig <- msigdbr(species = "Homo sapiens")
sen <- unique(msig[msig$gs_name %in% sets_wanted, c("gs_name", "gene_symbol")])
gene_sets <- split(sen$gene_symbol, sen$gs_name)

cohorts <- list(
  list(tag = "GSE96804",            expr = "codex_output/data/GSE96804-expr.rds",
       grp = "codex_output/data/GSE96804-group.csv",        pos = "DN"),
  list(tag = "GSE30528",            expr = "codex_output/data/GSE30528-expr.rds",
       grp = "codex_output/data/GSE30528-group.csv",        pos = "DKD"),
  list(tag = "GSE99325",            expr = "codex_output/data/GSE99325-expr.rds",
       grp = "codex_output/data/GSE99325-group.csv",        pos = "DN")
)

all_scores <- list()
summary_rows <- list()

for (co in cohorts) {
  cat("\n---", co$tag, "\n")
  expr <- readRDS(co$expr)
  grp <- read.csv(co$grp, stringsAsFactors = FALSE)
  # keep only samples that carry a case/control label
  keep <- grp$group %in% c("Control", co$pos)
  grp <- grp[keep, , drop = FALSE]
  expr <- expr[, grp$sample, drop = FALSE]

  # drop genes that are missing in every sample, then fill any remaining gaps with the
  # gene's own mean across the cohort (ssGSEA cannot handle NA values)
  if (anyNA(expr)) {
    cat("NA values found:", sum(is.na(expr)), "\n")
    all_na <- rowSums(is.na(expr)) == ncol(expr)
    expr <- expr[!all_na, , drop = FALSE]
    nas <- is.na(expr)
    if (any(nas)) {
      rm_vals <- rowMeans(expr, na.rm = TRUE)
      expr[nas] <- rm_vals[row(expr)[nas]]
    }
    cat("after handling:", sum(is.na(expr)), "NA remaining |", nrow(expr), "genes\n")
  }

  group <- factor(grp$group, levels = c("Control", co$pos))
  cat("samples used:", ncol(expr), "| genes:", nrow(expr), "\n")
  print(table(group))
  if (sum(group == "Control") < 3 | sum(group == co$pos) < 3) {
    cat("skipped: fewer than 3 samples in one group\n"); next
  }

  gs <- gene_sets
  gs <- lapply(gs, function(g) intersect(g, rownames(expr)))
  gs <- gs[lengths(gs) >= 10]
  cat("gene sets usable:", length(gs), "| sizes:",
      paste(names(gs), lengths(gs), sep = "=", collapse = ", "), "\n")

  # GSVA >= 2.0 uses parameter objects instead of the old method = "ssGSEA" argument
  if (exists("ssgseaParam", where = asNamespace("GSVA"))) {
    scores <- GSVA::gsva(GSVA::ssgseaParam(as.matrix(expr), gs, normalize = TRUE),
                         verbose = FALSE)
  } else {
    scores <- GSVA::gsva(as.matrix(expr), gs, method = "ssGSEA",
                         kcdf = "Gaussian", verbose = FALSE)
  }
  scores <- t(scores)
  scores <- as.data.frame(scores)
  scores$sample <- rownames(scores)
  scores$group <- group[match(scores$sample, grp$sample)]
  all_scores[[co$tag]] <- scores
  write.csv(scores, sprintf("codex_output/results/senescence-scores-%s.csv", co$tag),
            row.names = FALSE)

  for (s in setdiff(colnames(scores), c("sample", "group"))) {
    x <- scores[[s]][scores$group == "Control"]
    y <- scores[[s]][scores$group == co$pos]
    wt <- suppressWarnings(wilcox.test(y, x))
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      dataset = co$tag, gene_set = s,
      n_control = length(x), n_case = length(y),
      median_control = round(median(x), 4), median_case = round(median(y), 4),
      delta = round(median(y) - median(x), 4),
      p_value = signif(wt$p.value, 4), stringsAsFactors = FALSE
    )
  }
}

scoresum <- do.call(rbind, summary_rows)
scoresum$FDR <- ave(scoresum$p_value, scoresum$dataset,
                    FUN = function(p) p.adjust(p, method = "BH"))
scoresum <- scoresum[order(scoresum$dataset, scoresum$p_value), ]
write.csv(scoresum, "codex_output/results/senescence-score-comparison.csv", row.names = FALSE)
cat("\n=== case vs control comparison ===\n")
print(scoresum, row.names = FALSE)

# combined boxplot across datasets
long <- do.call(rbind, lapply(names(all_scores), function(tag) {
  d <- all_scores[[tag]]
  m <- reshape2::melt(d, id.vars = c("sample", "group"),
                      variable.name = "gene_set", value.name = "score")
  m$dataset <- tag
  m
}))
long$gene_set <- gsub("_", " ", long$gene_set)
long$gene_set <- factor(long$gene_set, levels = sort(unique(long$gene_set)))
long$dataset <- factor(long$dataset, levels = c("GSE96804", "GSE30528", "GSE99325", "GSE99325-GPL19109"))

p <- ggplot(long, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.size = 0.4, width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.5) +
  facet_grid(gene_set ~ dataset, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#9ecae1", "DN" = "#de2d26", "DKD" = "#de2d26")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 3) +
  labs(x = NULL, y = "ssGSEA score", fill = NULL,
       title = "Senescence pathway activity across cohorts") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0),
        strip.text.y = element_text(size = 7), panel.grid.minor = element_blank())

ggsave("codex_output/figures/fig-senescence-scores-all.pdf", p, width = 13, height = 12)
ggsave("codex_output/figures/fig-senescence-scores-all.png", p, width = 13, height = 12, dpi = 250)
cat("\nfigure saved\n")

saveRDS(all_scores, "codex_output/data/senescence-scores.rds")
cat("\ndone\n")
sink()
close(log_con)
