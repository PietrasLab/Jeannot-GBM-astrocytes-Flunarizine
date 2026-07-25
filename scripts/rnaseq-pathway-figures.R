# ================================================================
# Title: TGF-beta pathway and GSEA figures (manuscript revision)
# Project: Jeannot et al., GBM astrocytes / Flunarizine
# Description:
#   Pathway-level analyses of this project's own RNA-seq, added
#   during revision to complement the PROGENy analysis in
#   scripts/rnaseq-figures.R (Figures 3F/3G).
#
#   Figures produced (results/RevFig_*):
#     1. PROGENy TGF-beta activity per sample, by condition
#     2. GSEA Hallmark dotplot, Flunarizine effect on IR background
#     3. GSEA GO:BP top-20 barplot, Flunarizine effect on IR background
#     4. Heatmap of PROGENy TGF-beta pathway genes
#     5. TGF-beta ligands TGFB1/2/3 expression
#
# Design:
#   16 samples, 4 conditions (nonIR/IR x DMSO/Flunarizine) x 4
#   replicates. Each replicate contains all four conditions, i.e. a
#   randomized block design. Every test below blocks on replicate,
#   matching the DESeq2 design ~ Rep + Condition.
#
# Statistics notes:
#   - Gene-level significance is taken from DESeq2 (Rep-blocked,
#     genome-wide FDR), NOT recomputed with ad hoc t-tests on
#     transformed values. An unpaired t-test on VST values ignores
#     the blocking and disagrees with the model: for TGFB1 it gives
#     p = 0.15 ("ns") where DESeq2 gives padj = 0.0016.
#   - Sample-level scores are compared with PAIRED t-tests blocked on
#     replicate, BH-adjusted across the three contrasts shown.
#     Wilcoxon is not usable here: with 4 pairs its smallest
#     attainable two-sided p is 0.125, so it can never be significant.
#
# Usage:  Rscript scripts/rnaseq-pathway-figures.R
# ================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(decoupleR)
  library(msigdbr)
  library(fgsea)
  library(ComplexHeatmap)
  library(circlize)
  library(ggpubr)
})

set.seed(42)

OUT_PREFIX <- "RevFig"
PROGENY_TOP <- 500   # matches scripts/rnaseq-figures.R
TOP_N_GENES <- 100   # PROGENy TGF-beta genes carried into the heatmap

# ---- Project root ---------------------------------------------
SELF <- file.path("scripts", "rnaseq-pathway-figures.R")

find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, SELF))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the project root: no ", SELF, " at or above\n  ",
           normalizePath(start, mustWork = FALSE))
    }
    path <- parent
  }
}

project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
project_root <- if (nzchar(project_root)) normalizePath(path.expand(project_root), mustWork = TRUE) else find_project_root()

out_dir <- file.path(project_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- function(...) file.path(out_dir, paste0(OUT_PREFIX, "_", ...))

save_fig <- function(plot, stem, width, height) {
  ggsave(out_path(stem, ".pdf"), plot, width = width, height = height)
  ggsave(out_path(stem, ".png"), plot, width = width, height = height, dpi = 300)
  invisible(NULL)
}

# ---- Load data ------------------------------------------------
dds <- readRDS(file.path(project_root, "data", "rnaseq", "deseq2.filtered.rds"))

fdata <- as.data.frame(rowData(dds))
cdata <- as.data.frame(colData(dds))
stopifnot(identical(as.character(colnames(dds)), as.character(cdata$sample_id)))

# Group the axis logically: vehicle then drug, within each irradiation
# state. The stored factor order interleaves them.
COND_LEVELS <- c("nonIR_DMSO", "nonIR_Flunarizine", "IR_DMSO", "IR_Flunarizine")
COND_COLS <- c(
  "nonIR_DMSO"        = "#E8E8E8",
  "nonIR_Flunarizine" = "#A8D5BA",
  "IR_DMSO"           = "#FFB347",
  "IR_Flunarizine"    = "#FF6961"
)
cdata$Condition <- factor(as.character(cdata$Condition), levels = COND_LEVELS)

message("Samples per condition:")
print(table(cdata$Condition))

# ---- DESeq2 contrasts -----------------------------------------
# Each contrast is re-fitted on its own 8-sample subset, blocked on
# replicate, exactly as in the original sandbox analysis.
run_contrast <- function(conds, num, den) {
  d <- dds[, colData(dds)$Condition %in% conds]
  colData(d) <- droplevels(colData(d))
  design(d) <- ~ Rep + Condition
  d <- DESeq(d, quiet = TRUE)
  results(d, contrast = c("Condition", num, den)) %>%
    as.data.frame() %>%
    mutate(ENSG = rownames(.)) %>%
    left_join(fdata, by = c("ENSG" = "ensembl_gene_id")) %>%
    filter(!is.na(stat))
}

message("\nFitting DESeq2 contrasts ...")
res_ir  <- run_contrast(c("nonIR_DMSO", "IR_DMSO"), "IR_DMSO", "nonIR_DMSO")
res_fln <- run_contrast(c("IR_DMSO", "IR_Flunarizine"), "IR_Flunarizine", "IR_DMSO")

# ---- PROGENy TGF-beta gene set --------------------------------
net <- get_progeny(organism = "human", top = PROGENY_TOP)

tgfb_pathway <- net %>% filter(source == "TGFb") %>% arrange(desc(weight))
tgfb_top_pos <- tgfb_pathway %>% filter(weight > 0) %>% slice_head(n = TOP_N_GENES)

write.csv(tgfb_pathway, out_path("PROGENy_TGFB_signature.csv"), row.names = FALSE)

# get_progeny() queries OmniPath at run time, so the network is not
# pinned by the package version the way MSigDB is by msigdbr. Save the
# whole thing, not just the TGF-beta rows: run_mlm is multivariate and
# fits all pathways jointly, so the TGF-beta score depends on every
# row here. Without this the figure cannot be reproduced if OmniPath
# changes.
write.csv(net, out_path("PROGENy_network_full.csv"), row.names = FALSE)
message("PROGENy network: ", nrow(net), " rows across ",
        dplyr::n_distinct(net$source), " pathways")

tgfb_ligands <- c("TGFB1", "TGFB2", "TGFB3")
genes_of_interest <- fdata$ensembl_gene_id[
  fdata$external_gene_name_ensembl %in% unique(c(tgfb_top_pos$target, tgfb_ligands))
]

message("\nPROGENy TGF-beta signature: ", nrow(tgfb_pathway), " genes (",
        sum(tgfb_pathway$weight > 0), " positive, ", sum(tgfb_pathway$weight < 0), " negative)")
message("Genes carried into the heatmap: ", length(genes_of_interest))

# ---- Per-sample PROGENy activity ------------------------------
# run_mlm on z-scored VST, matching the method used for Figures
# 3F/3G in scripts/rnaseq-figures.R (which applies it to the
# contrast statistic rather than per sample).
vst_mat <- assay(vst(dds))
sym_mat <- vst_mat
rownames(sym_mat) <- fdata$external_gene_name_ensembl[match(rownames(sym_mat), fdata$ensembl_gene_id)]
sym_mat <- sym_mat[!is.na(rownames(sym_mat)), ]
sym_mat_z <- t(scale(t(sym_mat)))

progeny_sample <- run_mlm(
  mat = sym_mat_z, net = net,
  .source = "source", .target = "target", .mor = "weight", minsize = 5
) %>%
  filter(source == "TGFb") %>%
  select(sample_id = condition, PROGENy_TGFB_score = score) %>%
  left_join(cdata %>% select(sample_id, Condition, Rep), by = "sample_id") %>%
  mutate(Condition = factor(as.character(Condition), levels = COND_LEVELS)) %>%
  arrange(Condition, Rep)

write.csv(progeny_sample, out_path("PROGENy_TGFB_per_sample.csv"), row.names = FALSE)

# ---- Paired tests on per-sample scores ------------------------
# Blocked on replicate, matching the DESeq2 design.
paired_t <- function(df, a, b) {
  x <- df %>% filter(Condition == a) %>% arrange(Rep) %>% pull(PROGENy_TGFB_score)
  y <- df %>% filter(Condition == b) %>% arrange(Rep) %>% pull(PROGENy_TGFB_score)
  stopifnot(length(x) == length(y), length(x) > 1)
  tt <- t.test(y, x, paired = TRUE)
  tibble(group1 = a, group2 = b, n_pairs = length(x),
         mean_diff = unname(tt$estimate), p_raw = tt$p.value)
}

CONTRASTS <- list(
  c("nonIR_DMSO", "IR_DMSO"),               # effect of irradiation
  c("IR_DMSO", "IR_Flunarizine"),           # effect of drug on IR background
  c("nonIR_DMSO", "nonIR_Flunarizine")      # effect of drug without IR
)

progeny_stats <- map_dfr(CONTRASTS, ~ paired_t(progeny_sample, .x[1], .x[2])) %>%
  mutate(
    p_adj = p.adjust(p_raw, method = "BH"),
    test  = "Paired t-test on per-sample PROGENy score, blocked on replicate, BH-adjusted"
  )

write.csv(progeny_stats, out_path("PROGENy_TGFB_stats.csv"), row.names = FALSE)
print(progeny_stats)

signif_label <- function(p) {
  as.character(cut(p, c(-Inf, 1e-4, 1e-3, 1e-2, 5e-2, Inf),
                   labels = c("****", "***", "**", "*", "ns")))
}

# ---- Figure 1: per-sample PROGENy TGF-beta activity -----------
y_rng  <- range(progeny_sample$PROGENy_TGFB_score, na.rm = TRUE)
y_step <- diff(y_rng) * 0.12

progeny_annot <- progeny_stats %>%
  mutate(label = signif_label(p_adj),
         y.position = max(y_rng) + y_step * seq_len(n()))

p_progeny <- ggplot(progeny_sample, aes(Condition, PROGENy_TGFB_score, fill = Condition)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 2.4, alpha = 0.8, colour = "grey20") +
  stat_pvalue_manual(progeny_annot, label = "label", tip.length = 0.01,
                     label.size = 4, bracket.size = 0.4) +
  scale_fill_manual(values = COND_COLS) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(
    title = "TGF-beta pathway activity (PROGENy)",
    subtitle = "Per-sample activity inferred with decoupleR run_mlm",
    x = NULL, y = "PROGENy TGF-beta activity score",
    caption = "Paired t-test blocked on replicate, Benjamini-Hochberg adjusted; n = 4 per condition"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1))

save_fig(p_progeny, "TGFB_PROGENy_per_sample", width = 6.5, height = 5.5)

# ---- Figure 5: TGF-beta ligands -------------------------------
# Significance comes from DESeq2, not from a test recomputed on the
# plotted VST values.
ligand_stats <- bind_rows(
  res_ir  %>% filter(Symbol %in% tgfb_ligands) %>%
    transmute(Symbol, log2FoldChange, padj, group1 = "nonIR_DMSO", group2 = "IR_DMSO",
              contrast = "IR vs nonIR (DMSO background)"),
  res_fln %>% filter(Symbol %in% tgfb_ligands) %>%
    transmute(Symbol, log2FoldChange, padj, group1 = "IR_DMSO", group2 = "IR_Flunarizine",
              contrast = "Flunarizine vs DMSO (IR background)")
) %>%
  mutate(test = "DESeq2 Wald test, design ~ Rep + Condition, BH-adjusted genome-wide")

write.csv(ligand_stats, out_path("TGFB_ligands_stats.csv"), row.names = FALSE)
print(ligand_stats)

ligand_expr <- as.data.frame(t(vst_mat[fdata$ensembl_gene_id[match(tgfb_ligands, fdata$Symbol)], ])) %>%
  setNames(tgfb_ligands) %>%
  rownames_to_column("sample_id") %>%
  left_join(cdata %>% select(sample_id, Condition, Rep), by = "sample_id") %>%
  pivot_longer(all_of(tgfb_ligands), names_to = "Symbol", values_to = "vst") %>%
  mutate(Condition = factor(as.character(Condition), levels = COND_LEVELS))

write.csv(ligand_expr, out_path("TGFB_ligands_expression.csv"), row.names = FALSE)

ligand_annot <- ligand_stats %>%
  group_by(Symbol) %>%
  group_modify(function(d, key) {
    rng <- range(ligand_expr$vst[ligand_expr$Symbol == key$Symbol])
    d$label <- signif_label(d$padj)
    d$y.position <- max(rng) + diff(rng) * c(0.10, 0.24)[seq_len(nrow(d))]
    d
  }) %>%
  ungroup()

p_ligands <- ggplot(ligand_expr, aes(Condition, vst, fill = Condition)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.8, colour = "grey20") +
  stat_pvalue_manual(ligand_annot, label = "label", tip.length = 0.01,
                     label.size = 3.6, bracket.size = 0.35) +
  facet_wrap(~ Symbol, scales = "free_y") +
  scale_fill_manual(values = COND_COLS) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(
    title = "TGF-beta ligand expression",
    subtitle = "VST-normalised counts; significance from the DESeq2 model",
    x = NULL, y = "VST-normalised expression",
    caption = "DESeq2 Wald test, design ~ Rep + Condition, BH-adjusted genome-wide"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1))

save_fig(p_ligands, "TGFB_ligands", width = 11, height = 4.8)

# ---- Figure 4: TGF-beta pathway heatmap -----------------------
mat <- vst_mat[genes_of_interest, ]
rownames(mat) <- fdata$Symbol[match(genes_of_interest, fdata$ensembl_gene_id)]
mat_z <- t(scale(t(mat)))

# Rows ordered by the magnitude of the irradiation effect.
row_ord <- res_ir %>%
  filter(Symbol %in% rownames(mat_z)) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  pull(Symbol)
mat_z <- mat_z[row_ord[row_ord %in% rownames(mat_z)], ]

write.csv(as.data.frame(mat_z) %>% rownames_to_column("Symbol"),
          out_path("TGFB_heatmap_zscores.csv"), row.names = FALSE)

col_annot <- HeatmapAnnotation(
  Treatment   = cdata$Treatment,
  Irradiation = cdata$Irradiation,
  col = list(
    Treatment   = c("DMSO" = "grey80", "Flunarizine" = "darkgreen"),
    Irradiation = c("nonIR" = "grey80", "IR" = "orange")
  )
)

ht <- Heatmap(
  mat_z, name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  top_annotation = col_annot,
  cluster_rows = FALSE, cluster_columns = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8)
)

pdf(out_path("TGFB_pathway_heatmap.pdf"), width = 9, height = 13); draw(ht); dev.off()
png(out_path("TGFB_pathway_heatmap.png"), width = 9, height = 13, units = "in", res = 300)
draw(ht); dev.off()

# ---- GSEA -----------------------------------------------------
# fgsea's multilevel algorithm (no nperm) - the nperm interface is
# the deprecated fgseaSimple and floors p-values at the permutation
# resolution. scripts/rnaseq-figures.R does not use fgsea, so this
# does not conflict with the published analysis.
fetch_msigdb <- function(collection, sub = NULL) {
  args <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% args) {
    if (is.null(sub)) msigdbr(species = "Homo sapiens", collection = collection)
    else msigdbr(species = "Homo sapiens", collection = collection, subcollection = sub)
  } else {
    if (is.null(sub)) msigdbr(species = "Homo sapiens", category = collection)
    else msigdbr(species = "Homo sapiens", category = collection, subcategory = sub)
  }
}

ranked <- res_fln$stat
names(ranked) <- res_fln$external_gene_name_ensembl
ranked <- sort(ranked[!is.na(names(ranked))], decreasing = TRUE)
stopifnot(!any(duplicated(names(ranked))))

run_gsea <- function(sets) {
  set.seed(42)
  fgsea(pathways = sets, stats = ranked, minSize = 15, maxSize = 500) %>%
    as_tibble() %>%
    mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";"),
           direction = ifelse(NES > 0, "Up with Flunarizine", "Down with Flunarizine"),
           significant = padj < 0.05) %>%
    arrange(desc(abs(NES)))
}

message("\nRunning GSEA (Hallmark) ...")
gsea_h <- run_gsea(split(x = fetch_msigdb("H")$gene_symbol, f = fetch_msigdb("H")$gs_name))
write.csv(gsea_h, out_path("GSEA_Hallmark_FlunarizineEffect.csv"), row.names = FALSE)

message("Running GSEA (GO:BP) ...")
gobp <- fetch_msigdb("C5", "GO:BP")
gsea_gobp <- run_gsea(split(x = gobp$gene_symbol, f = gobp$gs_name))
write.csv(gsea_gobp, out_path("GSEA_GOBP_FlunarizineEffect.csv"), row.names = FALSE)

message("Significant (padj<0.05): Hallmark ", sum(gsea_h$padj < 0.05),
        ", GO:BP ", sum(gsea_gobp$padj < 0.05))

# Sentence case for pathway labels, restoring acronyms that a naive
# title-case would mangle (DNA -> Dna).
ACRONYMS <- c("DNA","RNA","MRNA","TRNA","RRNA","ATP","GTP","NADH","TGF","NFKB","TNFA",
              "MHC","UV","IL2","IL6","STAT3","STAT5","MYC","E2F","G2M","MTORC1","MTOR",
              "P53","KRAS","WNT","PI3K","AKT","JAK","IFN","ROS","EMT","UPR","GTPASE")
pretty_pathway <- function(x, prefix) {
  s <- gsub(prefix, "", x)
  s <- gsub("_", " ", s)
  s <- paste0(toupper(substr(s, 1, 1)), tolower(substr(s, 2, nchar(s))))
  # ignore.case matters: sentence-casing capitalises the first word, so
  # "DNA REPAIR" becomes "Dna repair" and a lowercase pattern misses it.
  for (a in ACRONYMS) s <- gsub(paste0("\\b", a, "\\b"), a, s, ignore.case = TRUE)
  s
}

# ---- Figure 2: Hallmark dotplot -------------------------------
# All Hallmark pathways significant at FDR 5%. Single, stated rule -
# the earlier version mixed collections with a keyword filter.
dot_h <- gsea_h %>%
  filter(padj < 0.05) %>%
  mutate(pathway_label = pretty_pathway(pathway, "^HALLMARK_"))

p_dot <- ggplot(dot_h, aes(NES, reorder(pathway_label, NES))) +
  geom_vline(xintercept = 0, colour = "grey70", linetype = "dashed") +
  geom_point(aes(size = -log10(padj), colour = NES)) +
  scale_colour_gradient2(low = "#155F83", mid = "whitesmoke", high = "#FFA319", midpoint = 0) +
  scale_size_continuous(range = c(2, 7)) +
  labs(
    title = "GSEA: Hallmark pathways",
    subtitle = "Flunarizine vs DMSO, IR background",
    x = "Normalized enrichment score (NES)", y = NULL,
    size = "-log10(padj)", colour = "NES",
    caption = sprintf("All %d Hallmark pathways significant at FDR < 0.05", nrow(dot_h))
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 8))

save_fig(p_dot, "GSEA_Hallmark_dotplot", width = 8.5, height = 9)

# ---- Figure 3: GO:BP top-20 barplot ---------------------------
# Significant at FDR 5%, then the 20 largest |NES| in either
# direction, drawn in signed NES order.
bar_gobp <- gsea_gobp %>%
  filter(padj < 0.05) %>%
  slice_max(abs(NES), n = 20) %>%
  mutate(pathway_label = pretty_pathway(pathway, "^GOBP_"))

write.csv(bar_gobp, out_path("GSEA_GOBP_top20.csv"), row.names = FALSE)

p_bar <- ggplot(bar_gobp, aes(NES, reorder(pathway_label, NES), fill = NES)) +
  geom_col() +
  scale_fill_gradient2(low = "#155F83", mid = "whitesmoke", high = "#FFA319", midpoint = 0) +
  labs(
    title = "GSEA: GO biological process (top 20)",
    subtitle = "Flunarizine vs DMSO, IR background",
    x = "Normalized enrichment score (NES)", y = NULL,
    caption = "Pathways with FDR < 0.05, ranked by absolute NES"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 8),
        legend.position = "none")

save_fig(p_bar, "GSEA_GOBP_top20_barplot", width = 9, height = 8)

# ---- Session info ---------------------------------------------
writeLines(
  c("scripts/rnaseq-pathway-figures.R",
    paste0("PROGENy top-n: ", PROGENY_TOP, "; heatmap genes: ", length(genes_of_interest)),
    "GSEA: fgsea multilevel (no nperm), seed 42",
    "", capture.output(sessionInfo())),
  file.path(out_dir, "sessionInfo-rnaseq-pathway-figures.txt")
)

message("\nDone. Outputs written to: ", out_dir)
