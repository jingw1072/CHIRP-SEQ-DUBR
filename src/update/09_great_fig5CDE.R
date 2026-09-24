#!/usr/bin/env Rscript
# =============================================================================
# STEP 9 (R) - Fig 5C, 5D, 5E, EVEN/ODD scatter, and Supp. Fig S8C.
#
#   Fig 5C  Venn of EVEN / ODD peaks (counts from step 05)
#   ----    EVEN vs ODD signal scatter + Pearson r (the paper's specificity QC;
#           the authors did this with ggpubr::ggscatter - see their R script)
#   Fig 5D  GREAT "two nearest genes" peak->gene assignment; bar chart of
#           region-gene associations binned by distance to TSS
#   Fig 5E  Venn: DUBR-associated genes (GREAT)  vs  DEGs shared by dEX1 & dEX2
#           (Tables S2/S3)  ->  "direct DUBR-associated genes" (paper: 129)
#   Fig S8C Heatmap of RNA abundance (z-score) of those direct targets in
#           WT / dEX1 / dEX2 (approximated from the GEO RNA-seq bigWigs)
#   + BED files of peaks linked to up- / down-regulated targets (for S8D, step 11)
#
# WHY R HERE
#   GREAT is the one piece of this figure that has no Python equivalent:
#   rGREAT (Bioconductor) re-implements the GREAT algorithm locally and can
#   also submit to the GREAT web server (what the authors used, hg19). Venn /
#   scatter / heatmap are equally easy in either language; keeping them in the
#   same script as GREAT avoids passing gene lists back and forth.
#
# RUN   sbatch 09_great_fig5CDE.sh   (the .sh wrapper activates chirp-fig and calls this file)
# =============================================================================
suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer); library(rGREAT)
  library(readxl); library(dplyr); library(tidyr)
  library(ggplot2); library(ggpubr); library(eulerr); library(ComplexHeatmap); library(circlize)
})

# ---- paths -------------------------------------------------------------------
# TSCC layout: everything under the scratch "dubr" folder (see README section 6)
out      <- file.path("/tscc/lustre/ddn/scratch", Sys.getenv("USER"), "chirp_seq_analysis/dubr")
D        <- out
peak_dir <- file.path(D, "peaks")
fig_dir  <- file.path(D, "figures")
trk      <- file.path(D, "tracks")
xlsx     <- file.path(D, "paper/Table_S1-S6.xlsx")     # copy from doc/dubr_paper/ (scp) - see README
refgene  <- file.path(D, "ref/annot/refGene.txt.gz")   # downloaded by step 07
dir.create(fig_dir, showWarnings = FALSE)

# GREAT parameters (GREAT web defaults for the "two nearest genes" rule)
GREAT_RULE     <- "twoClosest"
GREAT_MAX_EXT  <- 1000000      # 1000 kb max extension
# TRUE = submit to great.stanford.edu (the web tool the authors used). With
# rGREAT 2.8 this currently fails because the server redirects http->https;
# the local mode with tss_source = "GREAT:hg19" uses GREAT's own TSS table
# and gives near-identical assignments, so it is the default here.
USE_ONLINE     <- FALSE

# =============================================================================
# 1. Fig 5C - EVEN / ODD Venn
# =============================================================================
vc <- read.delim(file.path(peak_dir, "fig5C_venn_counts.tsv"))
v  <- setNames(vc$count, vc$set)
# eulerr wants disjoint parts. "common" = EVEN peaks that have an ODD partner.
fit5c <- euler(c(EVEN = v[["EVEN_only"]], ODD = v[["ODD_only"]], "EVEN&ODD" = v[["common_EVEN_side"]]))
png(file.path(fig_dir, "fig5C_venn_EVEN_ODD.png"), width = 1000, height = 800, res = 200)
print(plot(fit5c, quantities = TRUE, fills = c("#6fa8c9", "#e0a070"),
           main = list(label = sprintf("Fig 5C: %d | %d | %d   (paper 783 | 2031 | 1253)",
                                       v[["EVEN_only"]], v[["common_EVEN_side"]], v[["ODD_only"]]), fontsize = 9)))
dev.off()

# =============================================================================
# 2. EVEN vs ODD signal concordance scatter (Pearson) - from step 06 table
# =============================================================================
sig <- read.delim(file.path(peak_dir, "DUBR_peaks_EVEN_ODD_signal.tsv"))
sig <- sig %>% mutate(log2EVEN = log2(EVEN_cpm + 0.01), log2ODD = log2(ODD_cpm + 0.01),
                      set = ifelse(pass_0.5_2 == 1, "final (0.5<=EVEN/ODD<=2)", "dropped"))
p <- ggscatter(sig, x = "log2ODD", y = "log2EVEN", color = "set", palette = c("grey60", "firebrick"),
               size = 0.8, alpha = 0.6, cor.coef = TRUE, cor.method = "pearson",
               xlab = "ODD  log2(CPM)", ylab = "EVEN  log2(CPM)",
               title = sprintf("EVEN vs ODD signal at %d confident peaks (%d kept)", nrow(sig), sum(sig$pass_0.5_2)))
ggsave(file.path(fig_dir, "fig5_EVEN_ODD_scatter_pearson.png"), p, width = 5.5, height = 5, dpi = 200)

# =============================================================================
# 3. Fig 5D - GREAT peak -> gene assignment ("two nearest genes", 1000 kb)
# =============================================================================
peaks <- import(file.path(peak_dir, "DUBR_peaks_final.bed"), format = "BED")
names(peaks) <- peaks$name
message(sprintf("DUBR final peaks: %d", length(peaks)))

great_mode <- NA
assoc <- NULL
if (USE_ONLINE) {
  # Online GREAT (v4.0.4, hg19) = the exact service the authors used.
  assoc <- tryCatch({
    job <- submitGreatJob(peaks, species = "hg19", rule = GREAT_RULE,
                          adv_twoDistance = GREAT_MAX_EXT / 1000, request_interval = 10)
    a <- getRegionGeneAssociations(job)
    great_mode <- "online GREAT 4.0.4 (hg19)"
    pdf(file.path(fig_dir, "fig5D_GREAT_online_association_plots.pdf"), width = 12, height = 4)
    plotRegionGeneAssociations(job); dev.off()
    a
  }, error = function(e) { message("online GREAT failed: ", conditionMessage(e)); NULL })
}
if (is.null(assoc)) {
  # Local re-implementation. "GREAT:hg19" = the TSS list used by the GREAT web
  # server (downloaded once by rGREAT); falls back to UCSC knownGene TxDb.
  tss_used <- "GREAT:hg19"
  res <- tryCatch(
    great(peaks, "GO:BP", tss_source = tss_used, mode = GREAT_RULE, extension = GREAT_MAX_EXT, verbose = FALSE),
    error = function(e) {
      message("GREAT:hg19 TSS unavailable (", conditionMessage(e), "); using TxDb.Hsapiens.UCSC.hg19.knownGene")
      tss_used <<- "TxDb.Hsapiens.UCSC.hg19.knownGene"
      great(peaks, "GO:BP", tss_source = tss_used, mode = GREAT_RULE, extension = GREAT_MAX_EXT, verbose = FALSE)
    })
  great_mode <- sprintf("local rGREAT (%s)", tss_used)
  assoc <- getRegionGeneAssociations(res)
  pdf(file.path(fig_dir, "fig5D_GREAT_local_association_plots.pdf"), width = 12, height = 4)
  plotRegionGeneAssociations(res); dev.off()
}
message("GREAT mode: ", great_mode)

# Long table: one row per region-gene association
ra <- as.data.frame(assoc) %>%
  mutate(peak = names(assoc)) %>%
  dplyr::select(peak, seqnames, start, end, annotated_genes, dist_to_TSS) %>%
  unnest(cols = c(annotated_genes, dist_to_TSS)) %>%
  rename(gene = annotated_genes)
# Entrez -> symbol if the TSS source returned numeric IDs (TxDb fallback)
if (all(grepl("^[0-9]+$", ra$gene))) {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  ra$gene <- AnnotationDbi::mapIds(org.Hs.eg.db, ra$gene, "SYMBOL", "ENTREZID")
  ra <- ra %>% filter(!is.na(gene))
}
write.table(ra, file.path(fig_dir, "fig5D_GREAT_region_gene_associations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
dubr_genes <- unique(ra$gene)
message(sprintf("region-gene associations: %d ; DUBR-associated genes: %d (paper Fig 5D: 2936 assoc., Fig 5E: 2807 genes)",
                nrow(ra), length(dubr_genes)))

# Fig 5D bar chart with the paper's distance bins (kb, signed: - = upstream of TSS)
brks <- c(-Inf, -500, -50, -5, 0, 5, 50, 500, Inf) * 1000
labs <- c("< -500", "-500 to -50", "-50 to -5", "-5 to 0", "0 to 5", "5 to 50", "50 to 500", "> 500")
ra$bin <- cut(ra$dist_to_TSS, breaks = brks, labels = labs, right = FALSE)
d5 <- ra %>% count(bin) %>% mutate(pct = 100 * n / sum(n))
prom_pct <- 100 * mean(abs(ra$dist_to_TSS) <= 5000)
p5d <- ggplot(d5, aes(bin, pct)) + geom_col(fill = "#3b6fb6") +
  labs(x = "Distance to TSS (kb)", y = "Region-gene associations (%)",
       title = sprintf("Fig 5D  %d associations, %d genes  [%s]", nrow(ra), length(dubr_genes), great_mode),
       subtitle = sprintf("%.1f%% within +/-5 kb of a TSS (paper: ~10%% at promoters)", prom_pct)) +
  theme_classic(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "fig5D_GREAT_distance_to_TSS.png"), p5d, width = 6, height = 4.2, dpi = 200)

# =============================================================================
# 4. Fig 5E - DUBR-associated genes  vs  DEGs shared by dEX1 and dEX2
# =============================================================================
read_deg <- function(sheet) read_xlsx(xlsx, sheet = sheet) %>%
  dplyr::select(GeneID, logFC, FDR, state) %>% mutate(GeneID = as.character(GeneID))
s2 <- read_deg("Table S2")    # WT vs dEX1
s3 <- read_deg("Table S3")    # WT vs dEX2
deg <- inner_join(s2, s3, by = "GeneID", suffix = c("_EX1", "_EX2")) %>%
  filter(state_EX1 == state_EX2) %>%                       # concordant direction
  mutate(direction = ifelse(state_EX1 == "upregulated", "up", "down"))
message(sprintf("shared DEGs: %d (up %d, down %d)  [paper Fig 3B: 1216 = 673 up + 543 down]",
                nrow(deg), sum(deg$direction == "up"), sum(deg$direction == "down")))
write.table(deg, file.path(fig_dir, "shared_DEGs_dEX1_dEX2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Table S4 (differential H3K27ac, dEX1 vs WT) -> gained / reduced BEDs for Fig 6D (step 11)
s4 <- read_xlsx(xlsx, sheet = "Table S4") %>% mutate(start = as.integer(start), end = as.integer(end))
for (st in c("gained", "reduced")) {
  b <- s4 %>% filter(state == st) %>% transmute(chr, start, end, name = sprintf("%s_%d", st, row_number()), fold)
  write.table(b, file.path(peak_dir, sprintf("TableS4_H3K27ac_%s.bed", st)), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  message(sprintf("Table S4 %s H3K27ac sites: %d", st, nrow(b)))
}

# GREAT only knows genes in its own annotation; the paper's 1087 (vs 1216) is
# the DEG set restricted to that universe. Recover the universe when we can
# (local mode: res@extended_tss holds Entrez IDs -> map to symbols).
universe <- tryCatch({
  if (exists("res")) {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    ids <- unique(as.character(res@extended_tss$gene_id))
    sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, ids, "SYMBOL", "ENTREZID"))
    unique(c(na.omit(sym), dubr_genes))          # dubr_genes are symbols already
  } else NULL
}, error = function(e) { message("universe mapping failed: ", conditionMessage(e)); NULL })
message(sprintf("GREAT gene universe: %s genes", if (is.null(universe)) "unknown" else length(universe)))
deg_u <- if (!is.null(universe)) deg %>% filter(GeneID %in% universe) else deg
message(sprintf("DEGs in GREAT gene universe: %d  (paper: 1087)%s", nrow(deg_u),
                if (is.null(universe)) "  [universe unknown for online mode; using all DEGs]" else ""))

direct <- deg_u %>% filter(GeneID %in% dubr_genes)
message(sprintf("*** direct DUBR-associated genes: %d (up %d, down %d)   [paper Fig 5E: 129; S8C: 58 up, 71 down]",
                nrow(direct), sum(direct$direction == "up"), sum(direct$direction == "down")))

# NB: eulerr splits set names on "&", so no "&" inside a set label
fit5e <- euler(c("DUBR-associated genes" = length(dubr_genes) - nrow(direct),
                 "DEGs (dEX1, dEX2)"     = nrow(deg_u) - nrow(direct),
                 "DUBR-associated genes&DEGs (dEX1, dEX2)" = nrow(direct)))
png(file.path(fig_dir, "fig5E_venn_DUBRgenes_vs_DEGs.png"), width = 1200, height = 900, res = 200)
print(plot(fit5e, quantities = TRUE, fills = c("#6fa8c9", "#e0a070"),
           main = list(label = sprintf("Fig 5E: %d vs %d -> %d   (paper 2807 vs 1087 -> 129)",
                                       length(dubr_genes), nrow(deg_u), nrow(direct)), fontsize = 9)))
dev.off()

# gene table + the peaks behind each direct target
direct_full <- direct %>% left_join(ra %>% group_by(gene) %>%
                  summarise(peaks = paste(peak, collapse = ","), dist_to_TSS = paste(dist_to_TSS, collapse = ","), .groups = "drop"),
                by = c("GeneID" = "gene")) %>% arrange(direction, FDR_EX1)
write.table(direct_full, file.path(fig_dir, "fig5E_direct_DUBR_associated_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
labelled <- c("HES1", "CSF1", "ACVRL1", "ATF5", "TAF12", "TAF1", "RUNX1", "NCOR1", "CTCFL", "CBFA2T3", "FOSB", "CD33")
message("paper-highlighted genes recovered: ", paste(intersect(labelled, direct$GeneID), collapse = ", "))

# BEDs of DUBR peaks linked to up / down direct targets (used for Fig S8D in step 11)
for (dir_ in c("up", "down")) {
  pk <- ra %>% filter(gene %in% direct$GeneID[direct$direction == dir_]) %>%
    distinct(peak, .keep_all = TRUE) %>% dplyr::select(seqnames, start, end, peak)
  pk$start <- pk$start - 1                                  # GRanges (1-based) -> BED (0-based)
  write.table(pk, file.path(peak_dir, sprintf("DUBR_peaks_direct_%sregulated.bed", dir_)),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# =============================================================================
# 5. Fig S8C - RNA abundance heatmap of the direct targets (WT / dEX1 / dEX2)
#    GEO has no count matrix, only bigWigs -> approximate gene expression as
#    mean exon coverage per replicate (a stand-in for CPM; fine for a z-score).
# =============================================================================
rg <- read.delim(refgene, header = FALSE, stringsAsFactors = FALSE)
rg <- rg[rg$V13 %in% direct$GeneID, c("V3", "V4", "V10", "V11", "V13")]
exons <- do.call(rbind, lapply(seq_len(nrow(rg)), function(i) {
  s <- as.integer(strsplit(rg$V10[i], ",")[[1]]); e <- as.integer(strsplit(rg$V11[i], ",")[[1]])
  data.frame(chr = rg$V3[i], start = s + 1, end = e, gene = rg$V13[i])
}))
exons_gr <- reduce(split(GRanges(exons$chr, IRanges(exons$start, exons$end)), exons$gene))
exon_bw_mean <- function(bwfile, grl) {
  cov <- import(bwfile, format = "BigWig", selection = BigWigSelection(unlist(grl)))
  sapply(grl, function(g) {
    h <- findOverlaps(g, cov)
    if (length(h) == 0) return(0)
    ov <- pintersect(g[queryHits(h)], cov[subjectHits(h)])
    sum(width(ov) * cov$score[subjectHits(h)]) / sum(width(g))
  })
}
rna_files <- list.files(file.path(trk, "geo_rna"), pattern = "bigwig$", full.names = TRUE)
rna_files <- rna_files[order(match(sub(".*_RNA_(WT|EX1|EX2)_.*", "\\1", rna_files), c("WT", "EX1", "EX2")))]
expr <- sapply(rna_files, exon_bw_mean, grl = exons_gr)
colnames(expr) <- sub(".*_RNA_", "", sub(".bigwig$", "", colnames(expr)))
lz  <- t(scale(t(log2(expr + 0.01))))                       # z-score per gene
lz  <- lz[rowSums(is.na(lz)) == 0, , drop = FALSE]
dir_vec <- setNames(direct$direction, direct$GeneID)[rownames(lz)]
cond    <- sub("_[0-9]$", "", colnames(lz))
show    <- intersect(labelled, rownames(lz))
ha_top  <- HeatmapAnnotation(condition = cond, col = list(condition = c(WT = "grey50", EX1 = "firebrick", EX2 = "steelblue")))
ha_row  <- rowAnnotation(mark = anno_mark(at = match(show, rownames(lz)), labels = show))
png(file.path(fig_dir, "figS8C_RNA_heatmap_direct_targets.png"), width = 1400, height = 1800, res = 200)
draw(Heatmap(lz, name = "RNA z-score\n(log2 exon cov.)", col = colorRamp2(c(-2, 0, 2), c("#2166ac", "#f7f7f7", "#d6604d")),
             row_split = factor(dir_vec, levels = c("up", "down")), cluster_columns = FALSE, column_split = factor(cond, levels = c("WT", "EX1", "EX2")),
             show_row_names = FALSE, top_annotation = ha_top, right_annotation = ha_row,
             column_title = sprintf("Fig S8C  direct DUBR targets: up %d / down %d (bigWig-derived expression)",
                                    sum(dir_vec == "up"), sum(dir_vec == "down"))))
dev.off()
write.table(data.frame(gene = rownames(expr), direction = dir_vec[rownames(expr)], round(expr, 3)),
            file.path(fig_dir, "figS8C_RNA_exon_coverage_direct_targets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("done. figures in ", fig_dir)
