# ================================================================
# Title: External data figures - in vivo TGF-beta activity in
#        reactive astrocytes (murine GBM, 10x Flex scRNA-seq)
# Project: Jeannot et al., GBM astrocytes / Flunarizine
# Description:
#   Self-contained analysis of an EXTERNAL, independently published
#   dataset, used as in vivo support for the in vitro TGF-beta
#   findings in this manuscript. No data or results from this
#   project's own RNA-seq are used as input, so the two lines of
#   evidence remain independent.
#
#   Comparisons (both on per-sample medians, not per-cell):
#     1) Astrocyte R vs Astrocyte TE/NT, paired within sample
#     2) Astrocyte R across Healthy / Primary / Recurrent
#
# External data source
# --------------------
#   Lindgren D, Rosberg R, Smolag KI, Johansson E, et al.
#   "Processed single-cell and spatial transcriptomics Seurat
#    objects for murine glioblastoma" (v1.0, 2025-12-13).
#   Zenodo. https://doi.org/10.5281/zenodo.17922854
#   Record: https://zenodo.org/records/17922854
#   Licence: CC-BY-4.0
#
#   Required file: seurat_flex_filtered_v1.0.rds  (520 MB)
#
#   The .rds is far too large to version-control, so it is NOT in
#   this repository. Download it from the record above and either
#     (a) place it at  data/external/seurat_flex_filtered_v1.0.rds
#     (b) or point the FLEX_SEURAT_RDS environment variable at it.
#
# Usage:  Rscript scripts/external-data-figures.R
# Outputs: results/ExtData_TGFB_*
# ================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(msigdbr)
  library(ggpubr)
  library(ggbeeswarm)
  library(patchwork)
})

ZENODO_DOI    <- "10.5281/zenodo.17922854"
ZENODO_URL    <- "https://zenodo.org/records/17922854"
SEURAT_FILE   <- "seurat_flex_filtered_v1.0.rds"
OUT_PREFIX    <- "ExtData_TGFB"

# AddModuleScore draws its control gene sets with sample(). Seurat
# seeds that internally (seed = 1), which is why the scores are
# stable, but the default is Seurat's to change - pin it here so an
# upstream change cannot silently shift the published numbers.
MODULE_SCORE_SEED <- 1
set.seed(MODULE_SCORE_SEED)

# ---- Resolve project root -------------------------------------
# Walk up from the working directory looking for this script, so the
# run works from the repo root or from scripts/ alike. The marker is
# this file rather than the .Rproj, which is git-ignored and so is
# absent in a fresh clone.
SELF <- file.path("scripts", "external-data-figures.R")

find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, SELF))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the project root: no ", SELF, " found at or above\n  ",
           normalizePath(start, mustWork = FALSE),
           "\nRun this from the repository, or set PROJECT_ROOT.")
    }
    path <- parent
  }
}

project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
project_root <- if (nzchar(project_root)) normalizePath(path.expand(project_root), mustWork = TRUE) else find_project_root()
out_dir      <- file.path(project_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_path <- function(...) file.path(out_dir, paste0(OUT_PREFIX, "_", ...))

# Save every figure as vector PDF (for the manuscript) and 300 dpi
# PNG (for slides and quick review) from the same plot object.
save_fig <- function(plot, stem, width, height) {
  ggsave(out_path(stem, ".pdf"), plot, width = width, height = height)
  ggsave(out_path(stem, ".png"), plot, width = width, height = height, dpi = 300)
  invisible(NULL)
}

# ---- Locate the external dataset ------------------------------
resolve_seurat_path <- function() {
  from_env <- Sys.getenv("FLEX_SEURAT_RDS", unset = "")
  if (nzchar(from_env)) {
    p <- path.expand(from_env)
    if (!file.exists(p)) {
      stop("FLEX_SEURAT_RDS is set to a file that does not exist:\n  ", p)
    }
    return(p)
  }

  p <- file.path(project_root, "data", "external", SEURAT_FILE)
  if (file.exists(p)) return(p)

  stop(
    "Could not find the external dataset '", SEURAT_FILE, "'.\n\n",
    "Download it from ", ZENODO_URL, " (DOI ", ZENODO_DOI, ") and either:\n",
    "  - place it at ", p, "\n",
    "  - or set FLEX_SEURAT_RDS to its location.\n"
  )
}

seurat_path <- resolve_seurat_path()
message("Loading external Seurat object:\n  ", seurat_path)
obj <- readRDS(seurat_path)

required_cols <- c("Level_3AC", "sample_type", "sample_id")
missing_cols  <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

# AddModuleScore needs log-normalised values, not raw counts.
norm_mat <- tryCatch(GetAssayData(obj, layer = "data"), error = function(e) NULL)
if (is.null(norm_mat) || ncol(norm_mat) == 0 || nrow(norm_mat) == 0) {
  message("No normalised data in the active assay; running NormalizeData()")
  obj <- NormalizeData(obj)
}

# ---- Subset astrocytes and define groups ----------------------
astro_keep <- c("Astrocyte TE", "Astrocyte NT", "Astrocyte R")
obj_astro  <- subset(obj, subset = Level_3AC %in% astro_keep)

obj_astro$astro_group <- factor(
  ifelse(obj_astro$Level_3AC %in% c("Astrocyte TE", "Astrocyte NT"),
         "Astrocyte TE/NT", "Astrocyte R"),
  levels = c("Astrocyte TE/NT", "Astrocyte R")
)
obj_astro$sample_type <- factor(
  as.character(obj_astro$sample_type),
  levels = c("Healthy", "Primary", "Recurrent")
)

message("\nCells in astrocyte subset:")
print(table(obj_astro$Level_3AC))

# ---- Hallmark TGF-beta gene set (mouse) -----------------------
# msigdbr renamed the `category` argument to `collection` in v10.
fetch_hallmark <- function() {
  args <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% args) {
    msigdbr(species = "Mus musculus", collection = "H")
  } else {
    msigdbr(species = "Mus musculus", category = "H")
  }
}

tgfb_set <- fetch_hallmark() %>%
  filter(gs_name == "HALLMARK_TGF_BETA_SIGNALING") %>%
  pull(gene_symbol) %>%
  unique()

tgfb_present <- intersect(tgfb_set, rownames(obj_astro))
message("\nHALLMARK_TGF_BETA_SIGNALING genes: ", length(tgfb_set),
        " total, ", length(tgfb_present), " present in the object")

if (length(tgfb_present) < 10) {
  stop("Too few TGF-beta hallmark genes found in the object (", length(tgfb_present), ").")
}

write.csv(data.frame(gene = tgfb_present),
          out_path("hallmark_genes_present.csv"), row.names = FALSE)

# ---- Per-cell module score ------------------------------------
obj_astro <- AddModuleScore(obj_astro, features = list(tgfb_present),
                            name = "Hallmark_TGFB", seed = MODULE_SCORE_SEED)

cell_scores <- obj_astro@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  transmute(
    cell_id,
    sample_id   = as.character(sample_id),
    sample_type = as.character(sample_type),
    Level_3AC   = as.character(Level_3AC),
    astro_group = as.character(astro_group),
    tgfb_score  = .data[["Hallmark_TGFB1"]]
  )

write.csv(cell_scores, out_path("cell_scores.csv"), row.names = FALSE)

# ---- Collapse to per-sample medians ---------------------------
# Sample, not cell, is the unit of replication here. Testing on
# cells would inflate n by ~4 orders of magnitude.
sample_scores_group <- cell_scores %>%
  group_by(sample_id, sample_type, astro_group) %>%
  summarise(
    n_cells           = n(),
    tgfb_score_median = median(tgfb_score, na.rm = TRUE),
    tgfb_score_mean   = mean(tgfb_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(astro_group = factor(astro_group, levels = c("Astrocyte TE/NT", "Astrocyte R")))

write.csv(sample_scores_group, out_path("sample_scores_by_astro_group.csv"), row.names = FALSE)

sample_scores_acr <- sample_scores_group %>%
  filter(astro_group == "Astrocyte R") %>%
  select(sample_id, sample_type, n_cells, tgfb_score_median, tgfb_score_mean) %>%
  mutate(sample_type = factor(sample_type, levels = c("Healthy", "Primary", "Recurrent")))

write.csv(sample_scores_acr, out_path("sample_scores_AstrocyteR.csv"), row.names = FALSE)

# ---- Stats 1: Astrocyte R vs TE/NT, paired within sample ------
paired_group <- sample_scores_group %>%
  select(sample_id, astro_group, tgfb_score_median) %>%
  tidyr::pivot_wider(names_from = astro_group, values_from = tgfb_score_median) %>%
  filter(!is.na(`Astrocyte R`) & !is.na(`Astrocyte TE/NT`)) %>%
  mutate(delta_R_minus_TENT = `Astrocyte R` - `Astrocyte TE/NT`)

write.csv(paired_group, out_path("paired_AstroR_vs_TENT.csv"), row.names = FALSE)

if (nrow(paired_group) < 3) stop("Too few paired samples for a meaningful test.")

# Exact tests throughout: n is small (5-6 per stage, 14 pairs) and
# there are no ties, so the normal approximation buys nothing and is
# not trustworthy at this size.
wt1 <- wilcox.test(paired_group$`Astrocyte R`,
                   paired_group$`Astrocyte TE/NT`,
                   paired = TRUE, exact = TRUE)

stats1 <- data.frame(
  comparison                = "Astrocyte R vs Astrocyte TE/NT (paired per-sample medians)",
  test                      = "Wilcoxon signed-rank (paired, exact)",
  n_pairs                   = nrow(paired_group),
  n_higher_in_AstroR        = sum(paired_group$delta_R_minus_TENT > 0),
  p_value                   = wt1$p.value,
  median_delta_R_minus_TENT = median(paired_group$delta_R_minus_TENT, na.rm = TRUE)
)
write.csv(stats1, out_path("stats_AstroR_vs_TENT.csv"), row.names = FALSE)

# ---- Stats 2: Astrocyte R across disease stage ----------------
pairwise_stage <- pairwise.wilcox.test(
  x = sample_scores_acr$tgfb_score_median,
  g = sample_scores_acr$sample_type,
  p.adjust.method = "BH",
  exact = TRUE
)

# Keep the raw p alongside the adjusted one so the correction is
# auditable rather than implied.
raw_p <- function(a, b) {
  wilcox.test(
    sample_scores_acr$tgfb_score_median[sample_scores_acr$sample_type == a],
    sample_scores_acr$tgfb_score_median[sample_scores_acr$sample_type == b],
    exact = TRUE
  )$p.value
}

pw_df <- as.data.frame(as.table(pairwise_stage$p.value)) %>%
  filter(!is.na(Freq)) %>%
  rename(group1 = Var1, group2 = Var2, p_adj = Freq) %>%
  mutate(
    group1 = as.character(group1),
    group2 = as.character(group2),
    p_raw  = map2_dbl(group1, group2, raw_p),
    test   = "Exact Wilcoxon rank-sum on per-sample medians, BH-adjusted"
  ) %>%
  select(group1, group2, p_raw, p_adj, test)

write.csv(pw_df, out_path("stats_AstroR_stage_pairwise.csv"), row.names = FALSE)

# ---- Plot 1: Astrocyte R vs TE/NT -----------------------------
paired_long <- paired_group %>%
  pivot_longer(cols = c("Astrocyte TE/NT", "Astrocyte R"),
               names_to = "astro_group", values_to = "score") %>%
  mutate(astro_group = factor(astro_group, levels = c("Astrocyte TE/NT", "Astrocyte R")))

p1 <- ggplot(sample_scores_group,
             aes(x = astro_group, y = tgfb_score_median, fill = astro_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.55) +
  geom_line(data = paired_long,
            aes(x = astro_group, y = score, group = sample_id),
            inherit.aes = FALSE, color = "grey55", alpha = 0.5) +
  geom_beeswarm(size = 2.2, cex = 2.2, color = "grey20") +
  scale_fill_manual(values = c("Astrocyte TE/NT" = "#9ECAE1", "Astrocyte R" = "#EF3B2C")) +
  labs(
    title    = "In vivo TGF-beta hallmark score",
    subtitle = "Astrocyte TE/NT vs Astrocyte R (per-sample medians)",
    x = NULL,
    y = "Hallmark TGF-beta module score",
    caption  = sprintf("Exact paired Wilcoxon signed-rank, n = %d samples, p = %s",
                       nrow(paired_group), format.pval(wt1$p.value, digits = 2))
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

# p1 and p2 are saved only as the combined panel below. To export a
# half-panel on its own, call save_fig(p1, "AstroR_vs_TENT_boxplot", 5.8, 5).

# ---- Plot 2: Astrocyte R across stages ------------------------
# Annotate with the BH-adjusted p values computed above. Letting
# ggpubr's stat_compare_means() draw the brackets instead would
# silently recompute UNCORRECTED p values, so the stars on the
# figure would disagree with the exported stats table.
signif_label <- function(p) {
  cut(p,
      breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
      labels = c("****", "***", "**", "*", "ns"))
}

y_max  <- max(sample_scores_acr$tgfb_score_median, na.rm = TRUE)
y_min  <- min(sample_scores_acr$tgfb_score_median, na.rm = TRUE)
y_step <- (y_max - y_min) * 0.12

bracket_order <- tibble::tribble(
  ~group1,    ~group2,     ~rank,
  "Healthy",  "Primary",   1,
  "Primary",  "Recurrent", 2,
  "Healthy",  "Recurrent", 3
)

stage_annot <- bracket_order %>%
  left_join(
    bind_rows(pw_df, pw_df %>% rename(group1 = group2, group2 = group1)),
    by = c("group1", "group2")
  ) %>%
  mutate(
    label    = as.character(signif_label(p_adj)),
    y.position = y_max + y_step * rank
  )

stopifnot(!any(is.na(stage_annot$p_adj)))

p2 <- ggplot(sample_scores_acr,
             aes(x = sample_type, y = tgfb_score_median, fill = sample_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.55) +
  geom_beeswarm(size = 2.2, cex = 2.2, color = "grey20") +
  stat_pvalue_manual(stage_annot, label = "label", tip.length = 0.01,
                     label.size = 4, bracket.size = 0.4) +
  scale_fill_manual(values = c("Healthy" = "#74C476", "Primary" = "#FB6A4A", "Recurrent" = "#CB181D")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(
    title    = "In vivo TGF-beta hallmark score in Astrocyte R",
    subtitle = "Per-sample medians across disease stages",
    x = NULL,
    y = "Hallmark TGF-beta module score",
    caption  = "Exact Wilcoxon rank-sum on per-sample medians, Benjamini-Hochberg adjusted"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))


# ---- Combined panel -------------------------------------------
p_combined <- p1 + p2 + plot_layout(ncol = 2)
save_fig(p_combined, "invivo_combined_panel", width = 12, height = 5)

# ---- Provenance and session info ------------------------------
writeLines(
  c(
    "External dataset used by scripts/external-data-figures.R",
    "",
    paste0("Zenodo record : ", ZENODO_URL),
    paste0("DOI           : ", ZENODO_DOI),
    paste0("File          : ", SEURAT_FILE),
    "Licence       : CC-BY-4.0",
    "",
    paste0("Resolved local path at runtime: ", seurat_path),
    paste0("AddModuleScore seed: ", MODULE_SCORE_SEED),
    "",
    "--- sessionInfo() ---",
    capture.output(sessionInfo())
  ),
  file.path(out_dir, "sessionInfo-external-data-figures.txt")
)

message("\nDone. Outputs written to: ", out_dir)
