# ================================================================
# Title: Flunarizine Transcriptome Analysis – DESeq2 and decoupleR
# Author: [Your Name]
# Date: [Insert Date]
# Description: Full pipeline for normalization, DE analysis, PCA,
#              transcriptional activity inference (decoupleR), and
#              heatmap/boxplot visualization of RNA-seq data. This
#              script is part of the Flunarizine project manuscript.
#              It includes PCA plots, GSEA signature heatmaps, and
#              differential expression analysis. Data is loaded
#              from preprocessed DESeq2 objects and metadata files
#              available in the `data/` directory.
#              Results are saved in the `results/` directory.
# ================================================================

# ---- Load libraries ----
library(DESeq2)
library(tidyverse)
suppressMessages(library(ggplot2))
library(cli)
library(glue)
library(ComplexHeatmap)
library(msigdbr)
library(decoupleR)

# ---- Set project paths and parameters ----
project.path <- "~/Projects/Flunarizine-Submission-Git-v1/"
setwd(project.path)

analysis.name <- "FlnProject-Manuscript"
analysis.version <- "v1"

deseq.raw.file <- "./data/deseq2.dds.RData"
fdata.file     <- "./data/fdata_full.rds"
pdata.file     <- "./data/pdata.rds"

# ---- Load and preprocess DESeq2 object ----
load(deseq.raw.file)
fdata.full <- readRDS(fdata.file)
pdata <- readRDS(pdata.file)

# Match and append sample metadata
u <- match(pdata$sample_fastq, colnames(dds))
dds <- dds[,u]
stopifnot(identical(colnames(dds), pdata$sample_fastq))

# # ---- save raw data matrix for Annotare Data submission
# # Extract raw counts matrix
raw_counts <- counts(dds)
head(raw_counts)
write.table(raw_counts,
    file = "data/raw_counts_matrix.txt",
    row.names = TRUE, sep = "\t", quote = FALSE)


# add coldata
colData(dds) <- cbind(colData(dds), DataFrame(pdata))
colnames(dds) <- dds$sample_id

# Append gene annotations from biomaRt
rowData(dds) <- fdata.full
rowData(dds)$Symbol <- rowData(dds)$external_gene_name_ensembl

# ---- Filter and clean feature data ----
cli::cli_alert_info("Filtering on gene transcript types")
dds <- dds[rowData(dds)$gene_biotype %in% c("protein_coding"), ]

# Remove genes with no symbol or duplicated symbols
dds <- dds[!rowData(dds)$Symbol == "", ]
dds <- dds[!duplicated(rowData(dds)$Symbol), ]

# Remove genes with zero counts
dds <- dds[rowSums(counts(dds)) > 1, ]

# add designd and run deseq
pdata <- colData(dds)
# design(dds) <- ~ Rep + Irradiation + Treatment + Irradiation:Treatment
design(dds) <- ~ Rep + Irradiation + Treatment
dds <- DESeq(dds)
# Save cleaned object
saveRDS(dds, file = file.path("./data/deseq2.filtered.rds"))





# ---- PCA Plot: (FIGURE S7C) -----
dds <- readRDS(file = file.path("./data/deseq2.filtered.rds"))
pdata <- colData(dds)
dds <- vst(dds) # vst transform for PCA
pcaData <- plotPCA(dds, intgroup=c("Rep", "Treatment", "Irradiation","Condition"), returnData=TRUE, ntop=5000)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# plot tp pdf
pdf(file = file.path("./results", "PCAplot.pdf"))
ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Irradiation)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()
dev.off()
# save csv
write.csv(pcaData, file = file.path("./results", "PCAdata.csv"), row.names = FALSE)



# ---- GSEA Signature Heatmaps (FIGURE S7A and S7B) ----
msig_list_h <- msigdbr(species = "Homo sapiens")
gene_sig_names <- c(
    "GOBP_POSITIVE_REGULATION_OF_CALCIUM_MEDIATED_SIGNALING",
    "GOBP_CALMODULIN_DEPENDENT_KINASE_SIGNALING_PATHWAY"
)

gene_sig_list <- list()
for (sig.i in gene_sig_names) {
    df.i <- msig_list_h[msig_list_h$gs_name == sig.i, ]
    if (nrow(df.i) > 0) {
        cli::cli_alert("found signature: {.var {sig.i}} ")
        cli::cli_alert("n.genes: {nrow(df.i)}")
        gene_sig_list[[sig.i]] <- df.i
    }
}

for (sig.x in names(gene_sig_list)) {
    cli::cli_alert(paste0("Calculating scores for signature: ", sig.x))

    df.i <- gene_sig_list[[sig.x]]
    df.ensgs <- df.i$ensembl_gene

    fdata.i <- rowData(dds[rownames(dds) %in% df.ensgs, ])
    rownames(fdata.i) <- fdata.i$Symbol

    mat.i <- assay(dds[rownames(dds) %in% df.ensgs, ])
    rownames(mat.i) <- fdata.i$Symbol

    u <- order(apply(mat.i, 1, var), decreasing = TRUE)
    mat.i <- mat.i[u, ]
    mat.i.rownorm <- mat.i - apply(mat.i, 1, median)

    heatmap1 <- ComplexHeatmap::Heatmap(
        mat.i.rownorm,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        name = "mat.i"
    )

    filename <- file.path("./results", glue("Heatmap_{sig.x}.pdf"))
    pdf(filename, width = 8, height = 16)
    plot(heatmap1)
    dev.off()

    # save the matrix data values
    write.csv(mat.i.rownorm, file = file.path("./results", glue("Heatmap_{sig.x}.csv")), row.names = TRUE)
}




# ---- Subset: IR_DMSO vs IR_Flunarizine (n = 8) ----
# Investigate the differential expression between Flunarizine treatement on IR bacckground
dds <- readRDS(file = file.path("./data/deseq2.filtered.rds"))
dds.ir <- dds[,colData(dds)$Condition %in% c("IR_DMSO", "IR_Flunarizine")]
colData(dds.ir) <- droplevels(colData(dds.ir))
design(dds.ir) <- ~ Rep + Condition
dds.ir <- estimateSizeFactors(dds.ir)
dds.ir <- DESeq(dds.ir)
rowData(dds.ir)


# ---- Differential expression: IR_Flunarizine vs IR_DMSO ----
# ---- decoupleR activity inference on IR subset ----
dds.ir_res <- results(dds.ir, contrast = c("Condition", "IR_Flunarizine", "IR_DMSO")) %>%
    as.data.frame(.) %>%
    mutate(ENSG = rownames(.)) %>%
    left_join(as.data.frame(rowData(dds.ir)), by = c("ENSG" = "ensembl_gene_id")) %>%
    filter(!is.na(stat)) %>%
    arrange(desc(stat))
dds.ir_res_df <- as.matrix(select(dds.ir_res, stat))
rownames(dds.ir_res_df) <- dds.ir_res$Symbol

net <- get_progeny(organism = "human", top = 500)

contrast_acts <- run_mlm(
    mat = dds.ir_res_df,
    net = net,
    .source = 'source',
    .target = 'target',
    .mor = 'weight',
    minsize = 5
)

write.csv(contrast_acts, "./results/contrast_acts_FlunarizineEffect-on-IR.csv", row.names = FALSE)

# ---- Plot decoupleR results: Figure 3G ----
p <- ggplot(contrast_acts, aes(x = reorder(source, score), y = score, fill = score)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient2(low = "#155F83E5", high = "#FFA319E5", mid = "whitesmoke", midpoint = 0) +
    theme_minimal() +
    xlab("Pathways") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

pdf("./results/Fig_3G_contrast_acts_FlunarizineEffect-on-IR.pdf", width = 11, height = 12)
print(p)
dev.off()


# ---- Heatmap of rlog values, top/bottom DEGs: Figure 3E (IR subset) ----
ngenes <- 32
topgenes <- dds.ir_res$ENSG[1:ngenes]
bottomgenes <- rev(dds.ir_res$ENSG)[1:ngenes]

mat <- assay(rlog(dds.ir))[c(topgenes, bottomgenes), ]
mat.rownorm <- t(apply(mat, 1, function(x) (x - min(x)) / (max(x) - min(x))))
rownames(mat.rownorm) <- c(dds.ir_res$Symbol[1:ngenes], rev(dds.ir_res$Symbol)[1:ngenes])

heatmap1 <- ComplexHeatmap::Heatmap(
    mat.rownorm,
    col = c("#155F83E5", "white", "#FFA319E5"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    name = "Expression"
)

pdf("./results/heatmap_top_bottom_32genes_FlunarizineEffect-on-IR.pdf", width = 11, height = 18)
heatmap1
dev.off()

write.csv(as.data.frame(mat.rownorm),
    "./results/heatmap_top_bottom_32genes_FlunarizineEffect-on-IR.csv",
    row.names = TRUE)



# ---- Subset: IR_DMSO vs DMSO (n = 8) ----
# Investigate the differential expression between IR vs nonIR (DMSO)
dds <- readRDS(file = file.path("./data/deseq2.filtered.rds"))
pdata
dds.fln <- dds[,colData(dds)$Condition %in% c("nonIR_DMSO", "IR_DMSO")]
colData(dds.fln) <- droplevels(colData(dds.fln))
design(dds.fln) <- ~ Rep + Condition
dds.fln <- estimateSizeFactors(dds.fln)
dds.fln <- DESeq(dds.fln)
rowData(dds.fln)


# ---- Differential expression: DMSO IR vs DMSO CTRL ----
# ---- decoupleR activity inference on non-IR subset ----
dds.fln_res <- results(dds.fln, contrast = c("Condition", "IR_DMSO", "nonIR_DMSO")) %>%
    as.data.frame(.) %>%
    mutate(ENSG = rownames(.)) %>%
    left_join(as.data.frame(rowData(dds.fln)), by = c("ENSG" = "ensembl_gene_id")) %>%
    filter(!is.na(stat)) %>%
    arrange(desc(stat))
dds.fln_res_df <- as.matrix(select(dds.fln_res, stat))
rownames(dds.fln_res_df) <- dds.fln_res$Symbol

net <- get_progeny(organism = "human", top = 500)

contrast_acts <- run_mlm(
    mat = dds.fln_res_df,
    net = net,
    .source = 'source',
    .target = 'target',
    .mor = 'weight',
    minsize = 5
)

write.csv(contrast_acts, "./results/contrast_acts_IR-effect-DMSO.csv", row.names = FALSE)

# ---- Plot decoupleR results: Figure 3F ----
p <- ggplot(contrast_acts, aes(x = reorder(source, score), y = score, fill = score)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient2(low = "#155F83E5", high = "#FFA319E5", mid = "whitesmoke", midpoint = 0) +
    theme_minimal() +
    xlab("Pathways") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

pdf("./results/Fig_3F_contrast_acts_IR-effect-DMSO.pdf", width = 11, height = 12)
print(p)
dev.off()


# ---- Save session info ----
writeLines(capture.output(sessionInfo()), "./results/sessionInfo.txt")