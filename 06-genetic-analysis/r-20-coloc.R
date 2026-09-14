# Step 18: colocalisation of the cis-eQTL and the DKD GWAS signal
# coloc.abf needs only summary statistics for the same SNPs in both datasets, which we already
# have locally (eQTLGen significant cis-eQTLs + FinnGen genome-wide summary statistics).
suppressPackageStartupMessages({ library(data.table) })
if (!requireNamespace("coloc", quietly = TRUE)) {
  options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  install.packages("coloc", update = FALSE)
}
suppressPackageStartupMessages(library(coloc))
log_con <- file("codex_output/logs/r-20-coloc.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 18: colocalisation ===\n")
cat("date:", format(Sys.time()), "\n\n")

# ---------- load exposure --------------------------------------------------
cat("reading eQTLGen cis-eQTLs...\n")
eq <- fread("codex_output/mr/eqtlgen-cis-eQTL-FDR0.05.txt.gz",
            select = c("SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele",
                       "Zscore","GeneSymbol","NrSamples"),
            showProgress = FALSE)
setnames(eq, c("SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele","Zscore","GeneSymbol","NrSamples"),
         c("SNP","chr","pos","ea_exp","oa_exp","z_exp","gene","n_exp"))

genes <- read.csv("codex_output/results/ml-consensus-genes.csv", stringsAsFactors = FALSE)$gene
genes <- c(genes, "ITPR3")
eq <- eq[gene %in% genes]
cat("exposure SNPs for the candidate genes:", nrow(eq), "\n")

# ---------- load outcome ---------------------------------------------------
cat("reading FinnGen summary statistics...\n")
fg <- fread("codex_output/mr/finngen_R12_DM_NEPHROPATHY.gz",
            select = c("#chrom","pos","ref","alt","rsids","pval","beta","sebeta","af_alt"),
            showProgress = FALSE)
setnames(fg, c("#chrom","ref","alt","rsids","pval","beta","sebeta","af_alt"),
         c("chr","ref","alt","SNP","p_out","b_out","se_out","af_alt"))
cat("FinnGen variants:", nrow(fg), "\n\n")

N_CASE <- 5579; N_CONTROL <- 90951
S_PROP <- N_CASE / (N_CASE + N_CONTROL)

# ---------- run coloc per gene ---------------------------------------------
results <- list()
detail <- list()
for (g in sort(unique(eq$gene))) {
  e <- eq[gene == g]
  e <- merge(e, fg, by = "SNP", suffixes = c("", "_fg"))
  # keep SNPs whose alleles can be aligned between the two datasets
  same <- (e$ea_exp == e$alt & e$oa_exp == e$ref) | (e$ea_exp == e$ref & e$oa_exp == e$alt)
  e <- e[same]
  if (nrow(e) < 50) {
    cat(sprintf("%-10s skipped (only %d alignable SNPs)\n", g, nrow(e)))
    next
  }
  # eQTLGen reports Z-scores; convert to beta / varbeta using the allele frequency
  e[, maf := pmin(af_alt, 1 - af_alt)]
  e <- e[maf > 0.01 & maf < 0.99]
  e[, b_exp := z_exp / sqrt(2 * maf * (1 - maf) * (n_exp + z_exp^2))]
  e[, vb_exp := 1 / (2 * maf * (1 - maf) * (n_exp + z_exp^2))]

  d1 <- list(beta = e$b_exp, varbeta = e$vb_exp, snp = e$SNP,
             MAF = e$maf, N = as.integer(median(e$n_exp)), type = "quant")
  d2 <- list(beta = e$b_out, varbeta = e$se_out^2, snp = e$SNP,
             MAF = e$maf, N = N_CASE + N_CONTROL, s = S_PROP, type = "cc")
  res <- tryCatch(coloc.abf(d1, d2), error = function(err) NULL)
  if (is.null(res)) { cat(sprintf("%-10s coloc failed\n", g)); next }
  pp <- res$summary
  results[[g]] <- data.frame(gene = g, n_snp = nrow(e),
                             PP_H0 = round(pp["PP.H0.abf"], 3),
                             PP_H1 = round(pp["PP.H1.abf"], 3),
                             PP_H2 = round(pp["PP.H2.abf"], 3),
                             PP_H3 = round(pp["PP.H3.abf"], 3),
                             PP_H4 = round(pp["PP.H4.abf"], 3))
  top <- res$results[which.max(res$results$SNP.PP.H4), ]
  detail[[g]] <- data.frame(gene = g, top_snp = top$snp,
                            PP_H4_top_snp = round(top$SNP.PP.H4, 3))
  cat(sprintf("%-10s n=%5d  PP.H4 = %.3f  (top SNP %s, %.3f)\n",
              g, nrow(e), pp["PP.H4.abf"], top$snp, top$SNP.PP.H4))
}

tab <- do.call(rbind, results)
det <- do.call(rbind, detail)
tab <- merge(tab, det, by = "gene")
tab <- tab[order(-tab$PP_H4), ]
cat("\n=== colocalisation results ===\n")
print(tab, row.names = FALSE)
write.csv(tab, "codex_output/results/coloc-results.csv", row.names = FALSE)

cat("\n=== interpretation ===\n")
cat("PP.H4 > 0.8  : strong evidence for a single shared causal variant (supports the MR)\n")
cat("PP.H4 > 0.5  : suggestive\n")
cat("PP.H3 high   : two distinct causal variants (MR result should be questioned)\n")
cat("\ndone\n")
sink()
close(log_con)
