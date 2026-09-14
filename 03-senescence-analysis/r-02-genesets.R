# Step 6a: inventory of senescence-related gene sets available in MSigDB

suppressPackageStartupMessages(library(msigdbr))

g <- msigdbr(species = "Homo sapiens")
cat("total gene sets:", length(unique(g$gs_name)), "\n")
cat("total gene-set memberships:", nrow(g), "\n\n")

hits <- grep("senescen|sen_?mayo|sasp|aging|ageing", unique(g$gs_name),
             ignore.case = TRUE, value = TRUE)
cat("senescence-related sets:", length(hits), "\n\n")
print(sort(hits))

cat("\n=== size of each candidate set ===\n")
sz <- table(g$gs_name[g$gs_name %in% hits])
print(sort(sz))

cat("\n=== collection breakdown ===\n")
print(table(g$gs_collection[g$gs_name %in% hits]))
