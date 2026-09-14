# Step 5a: load series matrices, map probes to genes, build per-cohort expression matrices
# Outputs: codex_output/data/<GSE>-expr.rds  +  <GSE>-group.csv  +  a QC summary

suppressPackageStartupMessages({
  library(GEOquery)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

dir.create("codex_output/data", showWarnings = FALSE, recursive = TRUE)
dir.create("codex_output/logs", showWarnings = FALSE, recursive = TRUE)

log_con <- file("codex_output/logs/r-01-prepare.log", open = "wt")
sink(log_con, split = TRUE)

cat("=== step 5a: prepare expression matrices ===\n")
cat("date:", format(Sys.time()), "\n\n")

qc <- list()

read_matrix <- function(path) {
  eset <- GEOquery::getGEO(filename = path, getGPL = FALSE)
  expr <- Biobase::exprs(eset)
  pd <- Biobase::pData(eset)
  list(expr = expr, pd = pd)
}

report <- function(tag, expr, groups, note = "") {
  q <- stats::quantile(expr, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  cat(sprintf("%-12s genes=%6d samples=%4d | quantiles: %s\n", tag, nrow(expr), ncol(expr),
              paste(sprintf("%.2f", q), collapse = " / ")))
  qc[[length(qc) + 1]] <<- data.frame(
    dataset = tag, genes = nrow(expr), samples = ncol(expr),
    min = q[1], q25 = q[2], median = q[3], q75 = q[4], max = q[5],
    groups = paste(names(table(groups)), table(groups), sep = "=", collapse = "; "),
    note = note, stringsAsFactors = FALSE
  )
}

collapse_to_symbol <- function(expr, symbols, method = "maxmean") {
  keep <- !is.na(symbols) & symbols != ""
  expr <- expr[keep, , drop = FALSE]
  symbols <- symbols[keep]
  # drop symbols mapping to many probes is handled by aggregation
  ord <- order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE)
  expr <- expr[ord, , drop = FALSE]
  symbols <- symbols[ord]
  dedup <- !duplicated(symbols)
  out <- expr[dedup, , drop = FALSE]
  rownames(out) <- symbols[dedup]
  out
}

########################################################################
# 1. GSE96804 (training, GPL17586 HTA 2.0, DN vs tumour-nephrectomy control)
########################################################################
cat("--- GSE96804\n")
g <- read_matrix("codex_output/geo/GSE96804/GSE96804_series_matrix.txt.gz")
pd <- g$pd
cat("pData columns:", paste(colnames(pd), collapse = ", "), "\n")
src <- as.character(pd$source_name_ch1)
grp96804 <- ifelse(grepl("diabetic", src, ignore.case = TRUE), "DN",
            ifelse(grepl("tumor|tumour|unaffected", src, ignore.case = TRUE), "Control", NA))
names(grp96804) <- rownames(pd)
cat("group counts:\n"); print(table(grp96804, useNA = "ifany"))

suppressPackageStartupMessages(library(hta20transcriptcluster.db))
sym <- AnnotationDbi::mapIds(hta20transcriptcluster.db, keys = rownames(g$expr),
                             column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
cat("probe->symbol mapped:", sum(!is.na(sym)), "of", length(sym), "\n")
expr96804 <- collapse_to_symbol(g$expr, unname(sym))
report("GSE96804", expr96804, grp96804)
saveRDS(expr96804, "codex_output/data/GSE96804-expr.rds")
write.csv(data.frame(sample = names(grp96804), group = unname(grp96804)),
          "codex_output/data/GSE96804-group.csv", row.names = FALSE)

########################################################################
# 2. GSE30528 (validation A, GPL571 U133A 2.0, DKD glomeruli vs control)
########################################################################
cat("\n--- GSE30528\n")
g <- read_matrix("codex_output/geo/GSE30528/GSE30528_series_matrix.txt.gz")
pd <- g$pd
ch <- apply(pd[, grep("characteristics", colnames(pd)), drop = FALSE], 1,
            function(x) paste(x, collapse = " | "))
grp30528 <- ifelse(grepl("diabetic kidney disease", ch, ignore.case = TRUE), "DKD",
            ifelse(grepl("control", ch, ignore.case = TRUE), "Control", NA))
names(grp30528) <- rownames(pd)
cat("group counts:\n"); print(table(grp30528, useNA = "ifany"))

suppressPackageStartupMessages(library(hgu133a2.db))
sym <- AnnotationDbi::mapIds(hgu133a2.db, keys = rownames(g$expr),
                             column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
cat("probe->symbol mapped:", sum(!is.na(sym)), "of", length(sym), "\n")
expr30528 <- collapse_to_symbol(g$expr, unname(sym))
report("GSE30528", expr30528, grp30528, note = "values are batch-corrected, centred near 0")
saveRDS(expr30528, "codex_output/data/GSE30528-expr.rds")
write.csv(data.frame(sample = names(grp30528), group = unname(grp30528)),
          "codex_output/data/GSE30528-group.csv", row.names = FALSE)

########################################################################
# 3. GSE99325 (validation B, GPL19184 custom CDF ENTREZG, tubulointerstitium)
########################################################################
cat("\n--- GSE99325 (spans two array platforms)\n")
meta <- read.csv("codex_output/geo/GSE99325/GSE99325-samples.csv", stringsAsFactors = FALSE)
target_gsm <- meta$gsm
src_map <- setNames(meta$source, meta$gsm)

# NOTE: all 18 DN samples sit on GPL19184, so merging platforms with ComBat would confound
# the batch term with the outcome. Keep the two platforms as separate validation sets.
for (pf in c("GPL19184", "GPL19109")) {
  path <- sprintf("codex_output/geo/GSE99325/GSE99340-%s_series_matrix.txt.gz", pf)
  g <- read_matrix(path)
  present <- intersect(target_gsm, colnames(g$expr))
  cat(sprintf("  %s: %d of %d target samples present\n", pf, length(present), length(target_gsm)))
  if (length(present) == 0) next

  src <- src_map[present]
  grp <- ifelse(grepl("Diabetic Nephropathy", src, ignore.case = TRUE), "DN",
         ifelse(grepl("Cadaveric Donor|Tumor Nephrectomy", src, ignore.case = TRUE), "Control", "Other"))
  names(grp) <- present

  # IDs are Entrez Gene IDs with an "_at" suffix (custom CDF v18 ENTREZG)
  entrez <- sub("_at$", "", rownames(g$expr))
  sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(entrez),
                               column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  e <- collapse_to_symbol(g$expr[, present, drop = FALSE],
                          unname(sym[match(entrez, names(sym))]))
  tag <- if (pf == "GPL19184") "GSE99325" else "GSE99325-GPL19109"
  report(tag, e, grp, note = paste("platform", pf, "; custom CDF ENTREZG v18"))
  saveRDS(e, sprintf("codex_output/data/%s-expr.rds", tag))
  write.csv(data.frame(sample = colnames(e), group = unname(grp), platform = pf),
            sprintf("codex_output/data/%s-group.csv", tag), row.names = FALSE)
}

########################################################################
# summary
########################################################################
cat("\n=== QC summary ===\n")
qcdf <- do.call(rbind, qc)
print(qcdf, row.names = FALSE)
write.csv(qcdf, "codex_output/data/qc-summary.csv", row.names = FALSE)

cat("\ndone\n")
sink()
close(log_con)
