# Step 7a: load the GSE195460 single-nucleus matrices, apply QC and merge

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

dir.create("codex_output/logs", showWarnings = FALSE, recursive = TRUE)
dir.create("codex_output/results", showWarnings = FALSE, recursive = TRUE)

log_con <- file("codex_output/logs/r-07-sc-load.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 7a: load GSE195460 snRNA ===\n")
cat("date:", format(Sys.time()), "\n\n")

files <- list.files("codex_output/geo/GSE195460", pattern = "filtered_feature_bc_matrix\\.h5$",
                    full.names = TRUE)
files <- sort(files)
cat("files found:", length(files), "\n")

objs <- list()
qc_rows <- list()
for (f in files) {
  sample <- sub("^GSE195460_", "", sub("_filtered_feature_bc_matrix\\.h5$", "", basename(f)))
  group <- ifelse(grepl("^DN", sample), "DKD", "Control")
  cat("--- reading", sample, "(", group, ")\n")
  counts <- Read10X_h5(f)
  if (is.list(counts)) counts <- counts[["Gene Expression"]]
  cat("    genes:", nrow(counts), "| barcodes:", ncol(counts), "\n")

  o <- CreateSeuratObject(counts = counts, project = sample,
                          min.cells = 3, min.features = 200)
  o$sample <- sample
  o$group <- group
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
  o[["percent.hb"]] <- PercentageFeatureSet(o, pattern = "^HB[ABDEG]")

  qc_rows[[length(qc_rows) + 1]] <- data.frame(
    sample = sample, group = group, cells = ncol(o),
    median_features = round(median(o$nFeature_RNA), 0),
    median_counts = round(median(o$nCount_RNA), 0),
    median_mt = round(median(o$percent.mt), 2), stringsAsFactors = FALSE)

  objs[[sample]] <- o
  gc(verbose = FALSE)
}

qc <- do.call(rbind, qc_rows)
cat("\n=== per-sample QC (before filtering) ===\n")
print(qc, row.names = FALSE)
write.csv(qc, "codex_output/results/sc-qc-per-sample.csv", row.names = FALSE)

merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))
# Seurat v5 keeps one counts layer per sample after merge; join them into a single layer
# so that downstream normalisation and averaging work on the full matrix.
merged <- JoinLayers(merged)
cat("\nmerged object:", ncol(merged), "cells x", nrow(merged), "genes\n")
print(table(merged$group))

# QC filtering: nuclei with very few or very many features, high mitochondrial content
keep <- merged$nFeature_RNA > 300 & merged$nFeature_RNA < 6000 &
        merged$percent.mt < 5 & merged$percent.hb < 5
cat("\ncells passing QC:", sum(keep), "of", ncol(merged), "\n")
print(table(merged$group[keep]))
merged <- subset(merged, cells = colnames(merged)[keep])
cat("after filtering:", ncol(merged), "cells\n")

saveRDS(merged, "codex_output/data/GSE195460-sc-merged.rds")
cat("\nsaved: codex_output/data/GSE195460-sc-merged.rds\n")
cat("done\n")
sink()
close(log_con)
