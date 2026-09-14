# Step 6b: differential expression analysis on the training cohort (GSE96804)

suppressPackageStartupMessages({
  library(limma)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

dir.create("codex_output/results", showWarnings = FALSE, recursive = TRUE)
dir.create("codex_output/figures", showWarnings = FALSE, recursive = TRUE)

log_con <- file("codex_output/logs/r-03-deg.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 6b: differential expression (GSE96804) ===\n")
cat("date:", format(Sys.time()), "\n\n")

expr <- readRDS("codex_output/data/GSE96804-expr.rds")
grp <- read.csv("codex_output/data/GSE96804-group.csv", stringsAsFactors = FALSE)
stopifnot(identical(colnames(expr), grp$sample))
group <- factor(grp$group, levels = c("Control", "DN"))
cat("dimensions:", dim(expr), "\n")
print(table(group))

# filtering: drop genes in the lowest expression quartile
q75 <- apply(expr, 1, function(x) quantile(x, 0.75))
keep <- q75 > quantile(q75, 0.25)
# NOTE: the HTA 2.0 array reports many non-coding transcripts. Restrict the analysis to
# protein-coding genes so the results are biologically interpretable.
ann_gc <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL"),
                                columns = "GENETYPE", keytype = "SYMBOL")
pc_genes <- unique(ann_gc$SYMBOL[!is.na(ann_gc$GENETYPE) & ann_gc$GENETYPE == "protein-coding"])
keep_pc <- rownames(expr) %in% pc_genes
cat("\nprotein-coding genes measured:", sum(keep_pc), "of", nrow(expr), "\n")
cat("genes retained after expression filter:", sum(keep & keep_pc), "\n")
expr_f <- expr[keep & keep_pc, , drop = FALSE]
expr_f <- normalizeBetweenArrays(expr_f, method = "quantile")

# limma
design <- model.matrix(~ group)
fit <- lmFit(expr_f, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
res$gene <- rownames(res)
res <- res[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
write.csv(res, "codex_output/results/GSE96804-limma-all.csv", row.names = FALSE)

cat("\n=== top 20 genes ===\n")
print(head(res, 20), row.names = FALSE)

n_up   <- sum(res$adj.P.Val < 0.05 & res$logFC >  0.5)
n_down <- sum(res$adj.P.Val < 0.05 & res$logFC < -0.5)
cat(sprintf("\nDEGs (adj.P<0.05, |log2FC|>0.5): up=%d down=%d total=%d\n", n_up, n_down, n_up + n_down))
cat(sprintf("DEGs (adj.P<0.05 only): %d\n", sum(res$adj.P.Val < 0.05)))

# senescence gene sets
sets_wanted <- c("SAUL_SEN_MAYO", "REACTOME_CELLULAR_SENESCENCE",
                 "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
                 "GOBP_CELLULAR_SENESCENCE", "FRIDMAN_SENESCENCE_UP",
                 "RODWELL_AGING_KIDNEY_UP")
msig <- msigdbr(species = "Homo sapiens")
sen <- unique(msig[msig$gs_name %in% sets_wanted, c("gs_name", "gene_symbol")])
cat("\n=== senescence gene sets ===\n")
print(table(sen$gs_name))

universe <- rownames(expr_f)
sen_in_data <- sen[sen$gene_symbol %in% universe, ]
cat("senescence genes measured:", length(unique(sen_in_data$gene_symbol)),
    "of", length(unique(sen$gene_symbol)), "\n")
write.csv(unique(sen_in_data), "codex_output/results/senescence-genesets-in-data.csv", row.names = FALSE)

# senescence-related DEGs
deg <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 0.5, ]
sen_deg <- deg[deg$gene %in% unique(sen_in_data$gene_symbol), ]
sen_deg <- sen_deg[order(sen_deg$adj.P.Val), ]
write.csv(sen_deg, "codex_output/results/GSE96804-senescence-DEGs.csv", row.names = FALSE)
cat(sprintf("\nsenescence-related DEGs: %d (up=%d, down=%d)\n",
            nrow(sen_deg), sum(sen_deg$logFC > 0), sum(sen_deg$logFC < 0)))
cat("\ntop 30 senescence-related DEGs:\n")
print(head(sen_deg[, c("gene", "logFC", "adj.P.Val")], 30), row.names = FALSE)

# volcano plot
res$sig <- "NS"
res$sig[res$adj.P.Val < 0.05 & res$logFC >  0.5] <- "Up in DN"
res$sig[res$adj.P.Val < 0.05 & res$logFC < -0.5] <- "Down in DN"
res$sig <- factor(res$sig, levels = c("Down in DN", "NS", "Up in DN"))
res$is_sen <- res$gene %in% unique(sen_in_data$gene_symbol)

top_lab <- rbind(head(res[res$sig == "Up in DN", ], 15),
                 head(res[res$sig == "Down in DN", ], 10))

p <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = sig), size = 1.2, alpha = 0.7) +
  geom_point(data = subset(res, is_sen & sig != "NS"), shape = 21,
             color = "black", size = 2, stroke = 0.3, fill = NA) +
  scale_color_manual(values = c("Down in DN" = "#2C7BB6", "NS" = "grey80", "Up in DN" = "#D7191C")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_text_repel(data = top_lab, aes(label = gene), size = 2.6, max.overlaps = 30,
                  segment.size = 0.2, min.segment.length = 0) +
  labs(x = expression(log[2]~"fold change (DN vs Control)"),
       y = expression(-log[10]~"adjusted P"),
       color = NULL,
       title = "GSE96804: DN vs Control",
       subtitle = "black circles = senescence-related genes") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave("codex_output/figures/fig-volcano-GSE96804.pdf", p, width = 7, height = 6)
ggsave("codex_output/figures/fig-volcano-GSE96804.png", p, width = 7, height = 6, dpi = 300)
cat("\nvolcano saved\n")

# heatmap of top senescence DEGs
hm_genes <- head(sen_deg$gene, 40)
hm_genes <- hm_genes[hm_genes %in% rownames(expr_f)]
mat <- expr_f[hm_genes, , drop = FALSE]
mat <- t(scale(t(mat)))
mat[mat >  2] <-  2
mat[mat < -2] <- -2
ann <- data.frame(Group = group, row.names = colnames(mat))
ph <- function() pheatmap(mat, annotation_col = ann, show_colnames = FALSE,
                          cluster_cols = TRUE, cluster_rows = TRUE,
                          fontsize_row = 8, fontsize = 9,
                          color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
                          main = "Top senescence-related DEGs (GSE96804)")
pdf("codex_output/figures/fig-heatmap-senescence-DEGs.pdf", width = 8, height = 9); ph(); dev.off()
png("codex_output/figures/fig-heatmap-senescence-DEGs.png", width = 8, height = 9,
    units = "in", res = 300); ph(); dev.off()
cat("heatmap saved\n")

saveRDS(list(expr_f = expr_f, group = group, res = res, sen_deg = sen_deg,
             sen_sets = unique(sen_in_data)),
        "codex_output/data/GSE96804-deg.rds")
cat("\ndone\n")
sink()
close(log_con)
