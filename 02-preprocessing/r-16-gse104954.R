# Step 14: process GSE104954 (independent cohort for the specificity test)
# 195 CEL files across two platforms (GPL24120 HGU133A, GPL22945 HGU133Plus2).
# RMA is run per platform, probes are mapped to symbols, each gene is z-scored within its
# platform, and the two platforms are then pooled (both contain DKD and non-DKD samples).

suppressPackageStartupMessages({
  library(affy); library(AnnotationDbi); library(hgu133a.db); library(hgu133plus2.db)
})

dir.create("codex_output/logs", showWarnings = FALSE, recursive = TRUE)
dir.create("codex_output/data", showWarnings = FALSE, recursive = TRUE)

log_con <- file("codex_output/logs/r-16-gse104954.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 14: process GSE104954 ===\n")
cat("date:", format(Sys.time()), "\n\n")

meta <- read.csv("codex_output/geo/GSE104954/GSE104954-samples.csv", stringsAsFactors = FALSE)
cel_dir <- "codex_output/geo/GSE104954/cel"
cels <- list.files(cel_dir, pattern = "\\.CEL\\.gz$")
meta$file <- cels[match(meta$gsm, sub("_.*$", "", cels))]
stopifnot(!any(is.na(meta$file)))

meta$diagnosis <- sub("^.*diagnosis:\\s*", "", meta$characteristics)
meta$diagnosis[!grepl("diagnosis:", meta$characteristics)] <- NA
meta$group <- ifelse(grepl("Diabetic nephropathy", meta$diagnosis), "DN",
              ifelse(grepl("Tumor nephrectomy", meta$diagnosis), "Control", "Other"))
cat("group by platform:\n"); print(table(meta$group, meta$platform))

########################################################################
# 1. RMA per platform
########################################################################
rma_one <- function(platform, cdfname) {
  sub <- meta[meta$platform == platform, ]
  cat("\n--- RMA:", platform, "|", nrow(sub), "arrays | cdf:", cdfname, "\n")
  eset <- justRMA(filenames = file.path(cel_dir, sub$file), cdfname = cdfname, verbose = FALSE)
  cat("    probe sets:", nrow(eset), "| samples:", ncol(eset), "\n")
  cat("    value range:", paste(round(range(exprs(eset)), 2), collapse = " - "), "\n")
  e <- exprs(eset)
  colnames(e) <- sub$gsm          # justRMA names columns after the CEL files, not the GSMs
  list(expr = e, samples = sub$gsm, pkg = cdfname)
}
a <- rma_one("GPL24120", "hgu133a")
b <- rma_one("GPL22945", "hgu133plus2")

########################################################################
# 2. map probes to gene symbols and collapse
########################################################################
collapse <- function(expr, sym) {
  keep <- !is.na(sym) & sym != ""
  expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
  ord <- order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE)
  expr <- expr[ord, , drop = FALSE]; sym <- sym[ord]
  dedup <- !duplicated(sym)
  out <- expr[dedup, , drop = FALSE]; rownames(out) <- sym[dedup]
  out
}
map_a <- AnnotationDbi::mapIds(hgu133a.db, keys = rownames(a$expr), column = "SYMBOL",
                               keytype = "PROBEID", multiVals = "first")
map_b <- AnnotationDbi::mapIds(hgu133plus2.db, keys = rownames(b$expr), column = "SYMBOL",
                               keytype = "PROBEID", multiVals = "first")
ea <- collapse(a$expr, unname(map_a))
eb <- collapse(b$expr, unname(map_b))
cat("\ngene-level matrices:", nrow(ea), "(GPL24120) and", nrow(eb), "(GPL22945)\n")

########################################################################
# 3. z-score within platform, then pool
########################################################################
genes_common <- intersect(rownames(ea), rownames(eb))
cat("common genes:", length(genes_common), "\n")
za <- t(scale(t(ea[genes_common, , drop = FALSE])))
zb <- t(scale(t(eb[genes_common, , drop = FALSE])))
za[is.na(za)] <- 0; zb[is.na(zb)] <- 0
expr <- cbind(za, zb)

g <- meta[match(colnames(expr), meta$gsm), ]
stopifnot(all(!is.na(g$group)))
cat("\npooled matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")
print(table(g$group))

saveRDS(expr, "codex_output/data/GSE104954-expr.rds")
write.csv(data.frame(sample = g$gsm, group = g$group, platform = g$platform,
                     diagnosis = g$diagnosis),
          "codex_output/data/GSE104954-group.csv", row.names = FALSE)
cat("saved normalised matrix\n")
cat("done\n")
sink()
close(log_con)
