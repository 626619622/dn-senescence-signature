# Step 17b: senescence scoring and core-gene localisation in the GSE131882 replication cohort
suppressPackageStartupMessages({
  library(Seurat); library(msigdbr); library(ggplot2); library(reshape2)
})
dir.create("codex_output/figures", showWarnings = FALSE, recursive = TRUE)
log_con <- file("codex_output/logs/r-18c-sc-scores.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 17b: senescence scoring in GSE131882 ===\n")
cat("date:", format(Sys.time()), "\n\n")
obj <- readRDS("codex_output/data/GSE131882-sc-clustered.rds")
cat("nuclei:", ncol(obj), "| cell types:", length(levels(obj$celltype)), "\n")
sets <- list(SenMayo = "SAUL_SEN_MAYO",
             CellularSen = "REACTOME_CELLULAR_SENESCENCE",
             SASP = "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
             KidneyAging = "RODWELL_AGING_KIDNEY_UP")
msig <- msigdbr(species = "Homo sapiens")
for (nm in names(sets)) {
  g <- intersect(unique(msig$gene_symbol[msig$gs_name == sets[[nm]]]), rownames(obj))
  obj <- AddModuleScore(obj, features = list(g), name = nm, seed = 2026)
  colnames(obj@meta.data)[ncol(obj@meta.data)] <- paste0(nm, "_score")
  cat(sprintf("%-12s genes used: %d\n", nm, length(g)))
}
core <- read.csv("codex_output/results/senescence-core-genes.csv", stringsAsFactors = FALSE)$gene
core_in <- intersect(core, rownames(obj))
obj <- AddModuleScore(obj, features = list(core_in), name = "Core87", seed = 2026)
colnames(obj@meta.data)[ncol(obj@meta.data)] <- "Core87_score"
cat("core genes used:", length(core_in), "of", length(core), "\n\n")
score_cols <- c("SenMayo_score","CellularSen_score","SASP_score","KidneyAging_score","Core87_score")
md <- obj@meta.data[, c("celltype","group","sample", score_cols)]

# cell-level comparison
res <- list()
for (ct in levels(obj$celltype)) {
  sub <- md[md$celltype == ct, ]
  if (length(unique(sub$group)) < 2) next
  for (sc in score_cols) {
    a <- sub[[sc]][sub$group == "DKD"]; b <- sub[[sc]][sub$group == "Control"]
    if (length(a) < 10 || length(b) < 10) next
    w <- suppressWarnings(wilcox.test(a, b))
    res[[length(res) + 1]] <- data.frame(celltype = ct, score = sc,
      n_DKD = length(a), n_Control = length(b),
      mean_DKD = round(mean(a), 3), mean_Control = round(mean(b), 3),
      delta = round(mean(a) - mean(b), 3), p = signif(w$p.value, 3))
  }
}
rd <- do.call(rbind, res)
rd$FDR <- ave(rd$p, rd$score, FUN = function(p) p.adjust(p, "BH"))
rd <- rd[order(rd$score, rd$p), ]
cat("=== cell-level: DKD vs control (top 20) ===\n")
print(head(rd, 20), row.names = FALSE)
write.csv(rd, "codex_output/results/sc2-senescence-tests.csv", row.names = FALSE)

# sample-level (pseudo-bulk) comparison
pb <- aggregate(md[, score_cols], by = list(sample = md$sample, group = md$group,
                                            celltype = md$celltype), FUN = mean)
nc <- aggregate(md$Core87_score, by = list(sample = md$sample, celltype = md$celltype), FUN = length)
pb$n <- nc$x[match(paste(pb$sample, pb$celltype), paste(nc$sample, nc$celltype))]
pb <- pb[pb$n >= 10, ]
res2 <- list()
for (ct in unique(pb$celltype)) {
  sub <- pb[pb$celltype == ct, ]
  if (sum(sub$group == "DKD") < 3 | sum(sub$group == "Control") < 3) next
  for (sc in score_cols) {
    w <- suppressWarnings(wilcox.test(sub[[sc]][sub$group == "DKD"], sub[[sc]][sub$group == "Control"]))
    res2[[length(res2) + 1]] <- data.frame(celltype = ct, score = sc,
      delta_sample = round(mean(sub[[sc]][sub$group == "DKD"]) - mean(sub[[sc]][sub$group == "Control"]), 4),
      p_sample = signif(w$p.value, 3))
  }
}
rd2 <- do.call(rbind, res2)
rd2$FDR_sample <- ave(rd2$p_sample, rd2$score, FUN = function(p) p.adjust(p, "BH"))
rd2 <- rd2[order(rd2$p_sample), ]
cat("\n=== sample-level (pseudo-bulk) comparison ===\n")
print(head(rd2, 20), row.names = FALSE)
write.csv(rd2, "codex_output/results/sc2-pseudobulk-tests.csv", row.names = FALSE)

# direction agreement with the discovery cohort
disc <- read.csv("codex_output/results/sc-senescence-score-tests.csv", stringsAsFactors = FALSE)
m <- merge(disc[, c("celltype","score","delta")], rd[, c("celltype","score","delta")],
           by = c("celltype","score"), suffixes = c("_discovery","_replication"))
m$same_direction <- sign(m$delta_discovery) == sign(m$delta_replication)
cat("\n=== direction agreement with GSE195460 ===\n")
cat("comparable cell type x score pairs:", nrow(m),
    "| same direction:", sum(m$same_direction, na.rm = TRUE), "\n")
write.csv(m, "codex_output/results/sc2-direction-agreement.csv", row.names = FALSE)

# core gene localisation
avg_core <- as.data.frame(AverageExpression(obj, features = core_in, group.by = "celltype",
                                            assays = "RNA", layer = "data")$RNA)
zc <- t(scale(t(as.matrix(avg_core))))
top_ct <- colnames(zc)[apply(zc, 1, which.max)]
cat("\n=== core gene localisation in GSE131882 ===\n")
print(sort(table(top_ct), decreasing = TRUE))
write.csv(data.frame(gene = rownames(zc), top_celltype = top_ct, z = round(apply(zc, 1, max), 2)),
          "codex_output/results/sc2-core-gene-localisation.csv", row.names = FALSE)

# figures
png("codex_output/figures/fig-sc2-umap-celltypes.png", width = 7, height = 6, units = "in", res = 300)
print(DimPlot(obj, reduction = "umap", group.by = "celltype", label = TRUE, raster = FALSE) +
        ggtitle("GSE131882 (independent replication)"))
dev.off()
long <- melt(md, id.vars = c("celltype","group","sample"))
png("codex_output/figures/fig-sc2-senescence-scores.png", width = 9, height = 12, units = "in", res = 250)
print(ggplot(long, aes(x = celltype, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, linewidth = 0.2) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#9ecae1", "DKD" = "#de2d26")) +
  labs(x = NULL, y = "module score", fill = NULL,
       title = "Senescence scores in the replication cohort") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top"))
dev.off()
saveRDS(obj, "codex_output/data/GSE131882-sc-annotated.rds")
cat("\ndone\n")
sink()
close(log_con)
