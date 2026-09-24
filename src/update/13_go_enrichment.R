#!/usr/bin/env Rscript
# =============================================================================
# STEP 13 (R) - Gene-ontology enrichment: local equivalent of the paper's
#               Metascape panels.
#
#   Fig S8A  GO:BP of ALL DUBR-associated genes (GREAT, step 09)
#   Fig 5F   GO:BP of the direct DUBR-associated genes (GREAT genes that are
#            also DEGs; 130 here / 129 in the paper). The paper draws a Sankey
#            term -> TF plot; we draw a dot plot and list the TFs per term.
#   Fig 6F   GO:BP of the DEGs bound by HES1 within 5 kb of the TSS (step 11)
#
# HOW  clusterProfiler::enrichGO does an over-representation test (hypergeometric,
#      Benjamini-Hochberg FDR, same statistics as Metascape's GO module).
#      Metascape additionally clusters redundant terms; enrichGO's simplify()
#      (semantic similarity) plays that role here.
#      Thresholds: Metascape REPORTS terms at raw P < 0.01 and the paper's
#      Fig 5F/6F colour scales are -log10(P) (~3-4), i.e. not FDR-filtered. We
#      use the same P < 0.01 reporting cutoff and keep the BH FDR in the table
#      so you can see how many terms would survive a stricter filter (for the
#      130-gene list: few - small lists rarely pass FDR with GO:BP).
#      Metascape also mixes in Reactome/KEGG ("Signaling by NOTCH", "CMYB
#      pathway" in Fig 5F are Reactome/PID terms); here GO:BP only.
#      Metascape itself is a website: paste the gene lists written by steps 09
#      and 11 (fig5E_direct_DUBR_associated_genes.tsv, fig6E_HES1_bound_DEGs.txt)
#      into https://metascape.org if you want the exact figures.
#
# RUN  sbatch 13_go_enrichment.sh   (wrapper activates chirp-fig and calls this file)
# =============================================================================
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
  library(dplyr); library(ggplot2)
})

out     <- file.path("/tscc/lustre/ddn/scratch", Sys.getenv("USER"), "chirp_seq_analysis/dubr")
D       <- out
fig_dir <- file.path(D, "figures")
go_dir  <- file.path(D, "go_enrichment"); dir.create(go_dir, showWarnings = FALSE)

# ---- gene lists produced by steps 09 and 11 -----------------------------------
assoc  <- read.delim(file.path(fig_dir, "fig5D_GREAT_region_gene_associations.tsv"))
direct <- read.delim(file.path(fig_dir, "fig5E_direct_DUBR_associated_genes.tsv"))
hes1   <- readLines(file.path(fig_dir, "fig6E_HES1_bound_DEGs.txt"))
gene_sets <- list(
  figS8A_DUBR_associated_genes = unique(assoc$gene),
  fig5F_direct_DUBR_targets    = unique(direct$GeneID),
  fig6F_HES1_bound_DEGs        = unique(hes1)
)
# Background ("universe"): all genes GREAT could have assigned = all symbols in
# org.Hs.eg.db that GREAT knows. Simplest defensible choice: all annotated
# protein-coding-ish symbols in org.Hs.eg.db (Metascape's default is similar).
universe <- unique(na.omit(AnnotationDbi::mapIds(org.Hs.eg.db,
              keys(org.Hs.eg.db, "ENTREZID"), "SYMBOL", "ENTREZID")))

# ---- TF annotation for the Fig 5F "term -> TF" view -----------------------------
# GO:0003700 "DNA-binding transcription factor activity" as a practical TF list.
tf_symbols <- unique(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = "GO:0003700",
                 keytype = "GOALL", columns = "SYMBOL")$SYMBOL))

run_go <- function(name, genes) {
  message(sprintf("=== %s : %d genes", name, length(genes)))
  ego <- enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  universe = universe, pAdjustMethod = "BH",
                  # NB: clusterProfiler applies pvalueCutoff to the ADJUSTED p. Keep
                  # everything here and filter on raw P < 0.01 below (Metascape default).
                  pvalueCutoff = 1, qvalueCutoff = 1,
                  minGSSize = 10, maxGSSize = 500)
  if (is.null(ego) || nrow(ego) == 0) { message("  no enriched terms"); return(invisible(NULL)) }
  # keep only P < 0.01 BEFORE simplify(): simplify() computes pairwise semantic
  # similarity between terms and is very slow on thousands of them.
  ego@result <- ego@result[ego@result$pvalue < 0.01, ]
  if (nrow(ego@result) == 0) { message("  no terms at P < 0.01"); return(invisible(NULL)) }
  ego <- simplify(ego, cutoff = 0.7, by = "pvalue", select_fun = min)     # collapse redundant terms
  res <- as.data.frame(ego) %>% arrange(pvalue) %>%
    mutate(TFs_in_term = sapply(strsplit(geneID, "/"), function(g) paste(intersect(g, tf_symbols), collapse = ",")))
  if (nrow(res) == 0) { message("  no terms at P < 0.01"); return(invisible(NULL)) }
  write.table(res, file.path(go_dir, paste0(name, "_GO_BP.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("  %d terms at P<0.01 after simplify() (%d with FDR<0.05); top 5:", nrow(res), sum(res$p.adjust < 0.05)))
  print(head(res[, c("Description", "Count", "pvalue", "p.adjust")], 5), row.names = FALSE)

  # Bubble plot in the paper's style: terms on y, -log10(P) colour, gene count size
  top <- res %>% head(15) %>% mutate(Description = factor(Description, levels = rev(Description)))
  p <- ggplot(top, aes(x = -log10(pvalue), y = Description, size = Count, colour = -log10(pvalue))) +
    geom_point() + scale_colour_gradient(low = "#f4a582", high = "#b2182b", name = "-log10 P") +
    scale_size(range = c(2, 8), name = "Gene number") +
    labs(x = "-log10(P)", y = NULL, title = sprintf("%s  (n = %d genes)", name, length(genes)),
         subtitle = sprintf("GO:BP over-representation; %d of %d terms also FDR < 0.05", sum(res$p.adjust < 0.05), nrow(res))) +
    theme_classic(base_size = 11)
  ggsave(file.path(go_dir, paste0(name, "_GO_BP_dotplot.png")), p, width = 8, height = 5.5, dpi = 200)
  invisible(res)
}

results <- Map(run_go, names(gene_sets), gene_sets)

# ---- Fig 5F companion table: which direct targets are TFs, and in which terms
tf_direct <- intersect(gene_sets$fig5F_direct_DUBR_targets, tf_symbols)
message(sprintf("\nTFs among the direct DUBR targets (%d): %s", length(tf_direct), paste(sort(tf_direct), collapse = ", ")))
message("paper Fig 5F highlights: BTAF1, CTCFL, TAF12, RUNX1, NCOR1, HES1, ATF5, TAF1")
if (!is.null(results$fig5F_direct_DUBR_targets)) {
  tf_terms <- results$fig5F_direct_DUBR_targets %>% filter(TFs_in_term != "") %>%
    dplyr::select(Description, pvalue, p.adjust, Count, TFs_in_term) %>% arrange(pvalue)
  write.table(tf_terms, file.path(go_dir, "fig5F_terms_with_TFs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  print(head(tf_terms, 10), row.names = FALSE)
}
message("done. results in ", go_dir)
