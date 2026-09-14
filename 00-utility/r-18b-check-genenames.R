# Check the gene identifier format used by the GSE131882 Drop-seq matrices
x <- readRDS("codex_output/geo/GSE131882/GSM3823944_diabetes.s3.dgecounts.rds")
m <- x$umicount$inex$all
cat("dimensions:", dim(m), "\n")
cat("\nrowname examples:\n"); print(head(rownames(m), 15))
cat("\ncolname examples:\n"); print(head(colnames(m), 5))
cat("\nany gene symbols present (TP53, COL1A1, NPHS1):",
    sum(c("TP53","COL1A1","NPHS1") %in% rownames(m)), "of 3\n")
cat("any Ensembl IDs (start with ENSG):", sum(grepl("^ENSG", rownames(m))), "\n")
