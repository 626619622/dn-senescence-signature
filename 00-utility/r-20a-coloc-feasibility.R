# Feasibility check: how many SNPs are available at the ITPR3 locus in the data we already have?
suppressPackageStartupMessages({ library(data.table) })

cat("=== eQTLGen significant cis-eQTLs for ITPR3 ===\n")
eq <- fread("codex_output/mr/eqtlgen-cis-eQTL-FDR0.05.txt.gz",
            select = c("SNP","SNPChr","SNPPos","GeneSymbol","Zscore","FDR"),
            showProgress = FALSE)
it <- eq[GeneSymbol == "ITPR3"]
cat("significant SNPs:", nrow(it), "\n")
cat("chromosome:", unique(it$SNPChr), "\n")
cat("position range:", min(it$SNPPos), "-", max(it$SNPPos),
    "(span", round((max(it$SNPPos) - min(it$SNPPos)) / 1000), "kb)\n")
cat("FDR range:", signif(min(it$FDR), 3), "-", signif(max(it$FDR), 3), "\n")
write.csv(it[order(-abs(Zscore))][1:min(5, nrow(it))], row.names = FALSE)

cat("\n=== FinnGen variants in the same window (whole file is already local) ===\n")
fg <- fread("codex_output/mr/finngen_R12_DM_NEPHROPATHY.gz",
            select = c("#chrom","pos","ref","alt","rsids","pval","beta","sebeta","af_alt"),
            showProgress = FALSE)
setnames(fg, c("#chrom","ref","alt","rsids","pval","beta","sebeta","af_alt"),
         c("chr","ref","alt","SNP","p","beta","se","af_alt"))
chr_i <- unique(it$SNPChr)[1]
lo <- min(it$SNPPos) - 5e5; hi <- max(it$SNPPos) + 5e5
win <- fg[chr == chr_i & pos >= lo & pos <= hi]
cat("FinnGen variants in the +/-500kb window:", nrow(win), "\n")
overlap <- intersect(it$SNP, win$SNP)
cat("SNPs present in BOTH datasets:", length(overlap), "\n")
cat("\nconclusion: a basic coloc.abf run needs both traits to cover the same SNPs.\n")
cat("eQTLGen full cis-eQTL data would give thousands of SNPs per gene; the significant-only\n")
cat("file gives", nrow(it), "for ITPR3.\n")
