# Step 8a: sample-level pseudo-bulk analysis to rule out cell-composition artefacts
# Cell-level comparisons can be driven by how many cells each sample contributed;
# aggregating to sample x cell-type and testing at the sample level removes that bias.

suppressPackageStartupMessages({
  library(Seurat)
  library(limma)
  library(ggplot2)
  library(reshape2)
})

log_con <- file("codex_output/logs/r-10-pseudobulk.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 8a: pseudo-bulk analysis ===\n")
cat("date:", format(Sys.time()), "\n\n")

obj <- readRDS("codex_output/data/GSE195460-sc-annotated.rds")
score_cols <- c("SenMayo_score","CellularSen_score","SASP_score","KidneyAging_score","Core87_score")

# --- 1. cell composition per sample ---------------------------------------
comp <- as.data.frame.matrix(table(obj$sample, obj$celltype))
comp$sample <- rownames(comp)
comp$group <- obj$group[match(comp$sample, obj$sample)]
comp$total <- rowSums(comp[, levels(obj$celltype)])
for (ct in levels(obj$celltype)) comp[[paste0("pct_", ct)]] <- round(100 * comp[[ct]] / comp$total, 2)
write.csv(comp, "codex_output/results/sc-composition-per-sample.csv", row.names = FALSE)
cat("=== cell composition per sample (%) ===\n")
print(comp[, c("sample","group","total", paste0("pct_", levels(obj$celltype)))], row.names = FALSE)

# --- 2. pseudo-bulk scores: mean score per sample x cell type -------------
md <- obj@meta.data[, c("sample","group","celltype", score_cols)]
pb <- aggregate(md[, score_cols],
                by = list(sample = md$sample, group = md$group, celltype = md$celltype),
                FUN = mean)
n_cells <- aggregate(md$Core87_score, by = list(sample = md$sample, celltype = md$celltype),
                     FUN = length)
pb$n_cells <- n_cells$x[match(paste(pb$sample, pb$celltype), paste(n_cells$sample, n_cells$celltype))]
pb <- pb[pb$n_cells >= 10, ]
write.csv(pb, "codex_output/results/sc-pseudobulk-scores.csv", row.names = FALSE)

cat("\n=== sample-level test of senescence scores (pseudo-bulk, n>=10 cells per sample) ===\n")
res <- list()
for (ct in unique(pb$celltype)) {
  sub <- pb[pb$celltype == ct, ]
  if (length(unique(sub$group)) < 2) next
  if (sum(sub$group == "DKD") < 3 | sum(sub$group == "Control") < 3) next
  for (sc in score_cols) {
    a <- sub[[sc]][sub$group == "DKD"]; b <- sub[[sc]][sub$group == "Control"]
    w <- suppressWarnings(wilcox.test(a, b))
    res[[length(res) + 1]] <- data.frame(
      celltype = ct, score = sc,
      n_DKD = length(a), n_Control = length(b),
      mean_DKD = round(mean(a), 4), mean_Control = round(mean(b), 4),
      delta = round(mean(a) - mean(b), 4), p_value = signif(w$p.value, 3),
      stringsAsFactors = FALSE)
  }
}
rdf <- do.call(rbind, res)
rdf$FDR <- ave(rdf$p_value, rdf$score, FUN = function(p) p.adjust(p, "BH"))
rdf <- rdf[order(rdf$score, rdf$p_value), ]
write.csv(rdf, "codex_output/results/sc-pseudobulk-tests.csv", row.names = FALSE)
print(rdf, row.names = FALSE)

cat("\n=== direction agreement: cell-level vs sample-level ===\n")
cell_lv <- read.csv("codex_output/results/sc-senescence-score-tests.csv", stringsAsFactors = FALSE)
m <- merge(cell_lv[, c("celltype","score","delta","p_value","FDR")],
           rdf[, c("celltype","score","delta","p_value","FDR")],
           by = c("celltype","score"), suffixes = c("_cell","_sample"))
m$same_direction <- sign(m$delta_cell) == sign(m$delta_sample)
m$sample_sig <- m$p_value_sample < 0.05
write.csv(m, "codex_output/results/sc-cell-vs-sample-level.csv", row.names = FALSE)
cat("comparisons:", nrow(m),
    "| same direction:", sum(m$same_direction),
    "| significant at sample level:", sum(m$sample_sig), "\n")
print(m[order(m$p_value_sample), c("celltype","score","delta_cell","delta_sample",
                                   "p_value_sample","same_direction")], row.names = FALSE)

# --- 3. pseudo-bulk expression for the core genes --------------------------
counts <- LayerData(obj, assay = "RNA", layer = "counts")
core <- read.csv("codex_output/results/senescence-core-genes.csv", stringsAsFactors = FALSE)
core_genes <- intersect(core$gene, rownames(counts))

keys <- paste(obj$sample, obj$celltype, sep = "|")
ukeys <- unique(keys)
groups <- obj$group[match(ukeys, keys)]
pb_expr <- matrix(0, nrow = length(core_genes), ncol = length(ukeys),
                  dimnames = list(core_genes, ukeys))
for (i in seq_along(ukeys)) {
  cells <- which(keys == ukeys[i])
  pb_expr[, i] <- Matrix::rowSums(counts[core_genes, cells, drop = FALSE])
}
# library size per sample (not per cell) - each pseudo-bulk column inherits its sample's depth
lib_per_sample <- tapply(Matrix::colSums(counts), obj$sample, sum)
sample_of <- sub("\\|.*$", "", ukeys)
lib_vec <- as.numeric(lib_per_sample[sample_of])
stopifnot(all(!is.na(lib_vec)), all(lib_vec > 0))
pb_cpm <- t(t(pb_expr) / lib_vec * 1e6)
pb_log <- log2(pb_cpm + 1)

cell_n <- as.integer(table(keys)[ukeys])
ct_of <- sub("^.*\\|", "", ukeys)
keep <- cell_n >= 10 & ct_of %in% c("Proximal tubule", "Fibroblast", "Endothelium",
                                    "Podocyte", "Mesangium", "Tubule", "Cycling")
sel <- which(cell_n >= 10)
pb_log <- pb_log[, sel, drop = FALSE]
grp <- factor(groups[sel])
ctv <- ct_of[sel]

cat("\n=== pseudo-bulk differential expression of core genes (per cell type) ===\n")
lim_res <- list()
for (ct in unique(ctv)) {
  idx <- which(ctv == ct)
  if (length(unique(grp[idx])) < 2) next
  if (sum(grp[idx] == "DKD") < 3 | sum(grp[idx] == "Control") < 3) next
  if (length(idx) < length(levels(grp)) + 3) next
  fit <- eBayes(lmFit(pb_log[, idx, drop = FALSE], model.matrix(~ grp[idx])))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  lim_res[[ct]] <- data.frame(celltype = ct, gene = rownames(tt), logFC = tt$logFC,
                              p_value = tt$P.Value, adj_p = tt$adj.P.Val)
  cat(sprintf("%-22s genes tested: %3d | p<0.05: %3d | FDR<0.1: %3d\n",
              ct, nrow(tt), sum(tt$P.Value < 0.05), sum(tt$adj.P.Val < 0.1)))
}
lrd <- do.call(rbind, lim_res)
write.csv(lrd, "codex_output/results/sc-pseudobulk-core-gene-de.csv", row.names = FALSE)

if (!is.null(lrd)) {
  cat("\n=== core genes with FDR<0.1 in at least one cell type, training direction kept ===\n")
  sig <- lrd[lrd$adj_p < 0.1, ]
  train <- read.csv("codex_output/results/GSE96804-senescence-DEGs.csv", stringsAsFactors = FALSE)
  sig$direction_matches_bulk <- sign(sig$logFC) == sign(train$logFC[match(sig$gene, train$gene)])
  print(head(sig[order(sig$adj_p), ], 25), row.names = FALSE)
  write.csv(sig, "codex_output/results/sc-pseudobulk-core-significant.csv", row.names = FALSE)
}

# --- 4. figure ------------------------------------------------------------
fdir <- "codex_output/figures"
pb_long <- melt(pb, id.vars = c("sample","group","celltype","n_cells"),
                variable.name = "score", value.name = "value")
pb_long <- pb_long[pb_long$celltype %in% c("Proximal tubule","Podocyte","Endothelium",
                                           "Mesangium","Fibroblast","Collecting duct",
                                           "Distal tubule","Thick ascending limb"), ]
p <- ggplot(pb_long[pb_long$score == "Core87_score", ],
            aes(x = celltype, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.25) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1.1, alpha = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
  scale_fill_manual(values = c("Control" = "#9ecae1", "DKD" = "#de2d26")) +
  labs(x = NULL, y = "pseudo-bulk core-87 score", fill = NULL,
       title = "Sample-level senescence score by cell type",
       subtitle = "each point = one sample (3 DKD vs 5 control)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")
ggsave(file.path(fdir, "fig-sc-pseudobulk-core87.png"), p, width = 9, height = 5, dpi = 300)

cat("\ndone\n")
sink()
close(log_con)
