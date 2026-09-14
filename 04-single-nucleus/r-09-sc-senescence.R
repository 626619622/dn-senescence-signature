# Step 7c: senescence markers and core-gene localisation at single-nucleus resolution

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(reshape2)
})

log_con <- file("codex_output/logs/r-09-sc-senescence.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 7c: senescence at single-nucleus resolution ===\n")
cat("date:", format(Sys.time()), "\n\n")

obj <- readRDS("codex_output/data/GSE195460-sc-clustered.rds")
cat("cells:", ncol(obj), "\n")

# --- cluster -> cell type mapping (from the marker panel scores) ----------
map <- c(
  "0"  = "Proximal tubule", "1"  = "Distal tubule",   "2"  = "Collecting duct",
  "3"  = "Thick ascending limb", "4"  = "Proximal tubule", "5"  = "Thick ascending limb",
  "6"  = "Collecting duct", "7"  = "Distal tubule",   "8"  = "Proximal tubule",
  "9"  = "Thick ascending limb", "10" = "Endothelium", "11" = "Cycling",
  "12" = "Podocyte",        "13" = "Collecting duct", "14" = "Endothelium",
  "15" = "Mesangium",       "16" = "Fibroblast",      "17" = "B/Plasma",
  "18" = "Proximal tubule", "19" = "T/NK"
)
obj$celltype <- factor(unname(map[as.character(obj$seurat_clusters)]),
                       levels = unique(unname(map)))
cat("cell types:\n")
print(table(obj$celltype, obj$group))

# --- senescence gene sets -------------------------------------------------
suppressPackageStartupMessages(library(msigdbr))
sets <- list(
  SenMayo        = "SAUL_SEN_MAYO",
  CellularSen    = "REACTOME_CELLULAR_SENESCENCE",
  SASP           = "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
  KidneyAging    = "RODWELL_AGING_KIDNEY_UP"
)
msig <- msigdbr(species = "Homo sapiens")
for (nm in names(sets)) {
  g <- unique(msig$gene_symbol[msig$gs_name == sets[[nm]]])
  g <- intersect(g, rownames(obj))
  obj <- AddModuleScore(obj, features = list(g), name = nm, seed = 2026)
  colnames(obj@meta.data)[ncol(obj@meta.data)] <- paste0(nm, "_score")
  cat(sprintf("%-12s genes used: %d\n", nm, length(g)))
}

core <- read.csv("codex_output/results/senescence-core-genes.csv", stringsAsFactors = FALSE)
core_genes <- intersect(core$gene, rownames(obj))
obj <- AddModuleScore(obj, features = list(core_genes), name = "Core87", seed = 2026)
colnames(obj@meta.data)[ncol(obj@meta.data)] <- "Core87_score"
cat("core genes used:", length(core_genes), "of", nrow(core), "\n\n")

# --- canonical senescence markers per cell type ---------------------------
markers <- c("CDKN1A","CDKN2A","CDKN2B","TP53","LMNB1","GLB1","SERPINE1",
             "IL6","CXCL8","HMGB1","GDF15","FN1","MKI67","COL1A1")
markers <- intersect(markers, rownames(obj))
counts <- LayerData(obj, assay = "RNA", layer = "counts")
dat    <- LayerData(obj, assay = "RNA", layer = "data")

det <- do.call(rbind, lapply(levels(obj$celltype), function(ct) {
  cells <- which(obj$celltype == ct)
  data.frame(celltype = ct, n_cells = length(cells),
             t(vapply(markers, function(g) {
               c(pct_expr = round(100 * mean(counts[g, cells] > 0), 1),
                 mean_expr = round(mean(dat[g, cells]), 3))
             }, numeric(2))), check.names = FALSE, stringsAsFactors = FALSE)
}))
write.csv(det, "codex_output/results/sc-marker-detection-by-celltype.csv", row.names = FALSE)

cat("=== detection rate (%) of canonical senescence markers by cell type ===\n")
pct <- do.call(rbind, lapply(levels(obj$celltype), function(ct) {
  cells <- which(obj$celltype == ct)
  round(100 * rowMeans(counts[markers, cells, drop = FALSE] > 0), 1)
}))
rownames(pct) <- levels(obj$celltype)
print(pct)

# --- senescence scores per cell type and condition -------------------------
score_cols <- c("SenMayo_score","CellularSen_score","SASP_score","KidneyAging_score","Core87_score")
md <- obj@meta.data[, c("celltype","group","sample", score_cols)]

agg <- aggregate(md[, score_cols], by = list(celltype = md$celltype, group = md$group), FUN = mean)
write.csv(agg, "codex_output/results/sc-senescence-scores-by-celltype.csv", row.names = FALSE)
cat("\n=== mean senescence score by cell type and condition ===\n")
print(agg, row.names = FALSE, digits = 3)

# Wilcoxon test per cell type: DKD vs Control
res <- list()
for (ct in levels(obj$celltype)) {
  sub <- md[md$celltype == ct, ]
  if (length(unique(sub$group)) < 2) next
  for (sc in score_cols) {
    a <- sub[[sc]][sub$group == "DKD"]; b <- sub[[sc]][sub$group == "Control"]
    if (length(a) < 10 || length(b) < 10) next
    w <- suppressWarnings(wilcox.test(a, b))
    res[[length(res) + 1]] <- data.frame(
      celltype = ct, score = sc, n_DKD = length(a), n_Control = length(b),
      mean_DKD = round(mean(a), 3), mean_Control = round(mean(b), 3),
      delta = round(mean(a) - mean(b), 3), p_value = signif(w$p.value, 3),
      stringsAsFactors = FALSE)
  }
}
rdf <- do.call(rbind, res)
rdf$FDR <- ave(rdf$p_value, rdf$score, FUN = function(p) p.adjust(p, "BH"))
rdf <- rdf[order(rdf$score, rdf$p_value), ]
write.csv(rdf, "codex_output/results/sc-senescence-score-tests.csv", row.names = FALSE)
cat("\n=== DKD vs Control per cell type (top 25 by p) ===\n")
print(head(rdf, 25), row.names = FALSE)

# --- core gene expression across cell types --------------------------------
avg_core <- AverageExpression(obj, features = core_genes, group.by = "celltype",
                              assays = "RNA", layer = "data")$RNA
avg_core <- as.data.frame(avg_core)
avg_core$gene <- rownames(avg_core)
write.csv(avg_core, "codex_output/results/sc-core87-by-celltype.csv", row.names = FALSE)

zcore <- t(scale(t(as.matrix(avg_core[, setdiff(colnames(avg_core), "gene")]))))
top_ct <- colnames(zcore)[apply(zcore, 1, which.max)]
cat("\n=== which cell type expresses each core gene most highly (top 25 genes) ===\n")
print(head(data.frame(gene = rownames(zcore), max_celltype = top_ct,
                      z = round(apply(zcore, 1, max), 2)), 25), row.names = FALSE)
cat("\ncore gene counts by top cell type:\n")
print(sort(table(top_ct), decreasing = TRUE))

# --- figures --------------------------------------------------------------
fdir <- "codex_output/figures"
long <- melt(md, id.vars = c("celltype","group","sample"),
             variable.name = "score", value.name = "value")
p1 <- ggplot(long, aes(x = celltype, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, linewidth = 0.2) +
  facet_wrap(~ score, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#9ecae1", "DKD" = "#de2d26")) +
  labs(x = NULL, y = "module score", fill = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
ggsave(file.path(fdir, "fig-sc-senescence-scores.png"), p1, width = 9, height = 12, dpi = 250)

p2 <- DotPlot(obj, features = markers, group.by = "celltype") + RotatedAxis() +
  ggtitle("Canonical senescence markers")
ggsave(file.path(fdir, "fig-sc-canonical-markers.png"), p2, width = 9, height = 5, dpi = 250)

p3 <- DotPlot(obj, features = core_genes[1:min(40, length(core_genes))],
              group.by = "celltype") + RotatedAxis() +
  theme(axis.text.x = element_text(size = 6)) + ggtitle("Core senescence-related genes")
ggsave(file.path(fdir, "fig-sc-core-genes.png"), p3, width = 13, height = 5, dpi = 250)

saveRDS(obj, "codex_output/data/GSE195460-sc-annotated.rds")
cat("\nfigures and results saved\ndone\n")
sink()
close(log_con)
