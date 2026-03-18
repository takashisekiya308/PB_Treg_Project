#!/usr/bin/env Rscript
#The R script to determine euclidean distances between M0 macrophage and IFN-g-polarized cells.

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
})

# ---- settings ----
samplesheet_path <- "samplesheet.txt"  # <- a samplesheet compatible with DESeq2
outdir <- "euclid_out"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- read samplesheet ----
meta <- read.table(samplesheet_path, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("sample","group","path") %in% colnames(meta)))
stopifnot(all(file.exists(meta$path)))

# ---- confirm original groups ----
original_groups <- sort(unique(meta$group))
cat("Original groups:\n")
print(original_groups)

# ---- merge groups as requested ----
meta$group_merged <- meta$group
meta$group_merged[meta$group %in% c("IFNg_Treg_421","IFNg_Treg_647")]   <- "IFNg_Treg"
meta$group_merged[meta$group %in% c("IFNg_Tconv_421","IFNg_Tconv_647")] <- "IFNg_Tconv"

meta$group_merged <- factor(meta$group_merged,
                            levels = c("M0_plus_AC","IFNg","IFNg_Treg","IFNg_Tconv"))

cat("\nMerged groups:\n")
print(table(meta$group_merged, useNA = "ifany"))
stopifnot(all(!is.na(meta$group_merged)))

rownames(meta) <- meta$sample
files <- setNames(meta$path, meta$sample)

# ---- tximport (RSEM genes.results) ----
txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

# DESeq2 requires length > 0
len <- txi$length
bad <- !is.finite(len) | len <= 0
bad_genes <- rowSums(bad) > 0
txi$counts    <- txi$counts[!bad_genes, , drop = FALSE]
txi$abundance <- txi$abundance[!bad_genes, , drop = FALSE]
txi$length    <- txi$length[!bad_genes, , drop = FALSE]

# ---- DESeq2 + VST ----
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ group_merged)

# low-count filter
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)   # VST_normalization

mat <- assay(vsd)  # genes x samples
stopifnot(all(colnames(mat) == rownames(meta)))

# ---- group mean vectors (centroids) ----
groups <- levels(meta$group_merged)
group_means <- sapply(groups, function(g) {
  cols <- rownames(meta)[meta$group_merged == g]
  rowMeans(mat[, cols, drop = FALSE])
})
# group_means: genes x groups

# ---- Euclidean distances from M0_plus_AC to others ----
m0 <- group_means[, "M0_plus_AC"]
dist_tbl <- data.frame(
  group = c("IFNg","IFNg_Treg","IFNg_Tconv"),
  euclidean_distance = c(
    sqrt(sum((m0 - group_means[, "IFNg"])^2)),
    sqrt(sum((m0 - group_means[, "IFNg_Treg"])^2)),
    sqrt(sum((m0 - group_means[, "IFNg_Tconv"])^2))
  ),
  stringsAsFactors = FALSE
)

# Also output full group-to-group distance matrix
dmat <- as.matrix(dist(t(group_means), method = "euclidean"))

write.table(dist_tbl,
            file = file.path(outdir, "euclidean_from_M0_plus_AC.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(dmat,
            file = file.path(outdir, "euclidean_group_distance_matrix.tsv"),
            sep = "\t", quote = FALSE)

cat("\nSaved:\n",
    file.path(outdir, "euclidean_from_M0_plus_AC.tsv"), "\n",
    file.path(outdir, "euclidean_group_distance_matrix.tsv"), "\n")