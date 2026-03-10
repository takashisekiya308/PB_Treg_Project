#!/usr/bin/env Rscript

# =========================
# DESeq2 pipeline (4 groups, triplicate)
# Inputs: samplesheet.txt (tab-delimited; columns: sample, group, path)
# Outputs:
#  1) PCA png
#  2) group-specific high/low DEG lists (strict: significant & consistent vs ALL other groups)
#  3) top500 DEG heatmap png (clustering shown)
# =========================

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)  
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(matrixStats)
})

# ---- user settings ----
samplesheet_path <- "samplesheet.txt"  
outdir <- "deseq_out"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# threshold for "specific" DEG
padj_cut <- 0.1
lfc_cut  <- 0.58  # |log2FC| >= 1 
min_support <- 2

# ---- read samplesheet ----
meta <- read.table(samplesheet_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("sample","group","path") %in% colnames(meta)))
stopifnot(all(file.exists(meta$path)))

meta$group <- factor(meta$group, levels = unique(meta$group)) 
rownames(meta) <- meta$sample

files <- setNames(meta$path, meta$sample)

# ---- import counts (RSEM genes.results) ----
txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

# ---- fix: drop genes with non-positive / invalid lengths (required by DESeq2) ----
len <- txi$length
bad <- !is.finite(len) | len <= 0
bad_genes <- rowSums(bad) > 0

cat("Genes with length<=0 or invalid in >=1 sample: ", sum(bad_genes), "\n")

txi$counts    <- txi$counts[!bad_genes, , drop = FALSE]
txi$abundance <- txi$abundance[!bad_genes, , drop = FALSE]
txi$length    <- txi$length[!bad_genes, , drop = FALSE]

# =========================
# Output "TPM-like" matrix (RSEM abundance from tximport)
# =========================
tpm_like <- txi$abundance 
tpm_like_df <- data.frame(gene_id = rownames(tpm_like), as.data.frame(tpm_like), check.names = FALSE)

tpm_like_df <- add_gene_symbol(tpm_like_df, "gene_id")

write.table(tpm_like_df,
            file = file.path(outdir, "RSEM_TPM_like_matrix.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ group)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds <- DESeq(dds)

# =========================
# Output normalized expression matrices
# =========================

# 1) DESeq2 normalized counts (size-factor normalized)
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- data.frame(gene_id = rownames(norm_counts), as.data.frame(norm_counts), check.names = FALSE)

norm_counts_df <- add_gene_symbol(norm_counts_df, "gene_id")

write.table(norm_counts_df,
            file = file.path(outdir, "DESeq2_normalized_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 2) VST matrix
# 既に vsd を作っているならそれを使ってOK
vst_mat <- assay(vsd)
vst_df <- data.frame(gene_id = rownames(vst_mat), as.data.frame(vst_mat), check.names = FALSE)
vst_df <- add_gene_symbol(vst_df, "gene_id")

write.table(vst_df,
            file = file.path(outdir, "DESeq2_VST_matrix.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- VST for PCA / heatmap ----
vsd <- vst(dds, blind = FALSE)

# =========================
# 1) PCA plot (png)
# =========================
pca <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

p_pca <- ggplot(pca, aes(x = PC1, y = PC2, color = group, label = name)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw(base_size = 13) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  ggtitle("PCA (VST)")

ggsave(filename = file.path(outdir, "PCA.png"), plot = p_pca, width = 6.5, height = 5.5, dpi = 300)


# ---- ENSEMBL -> SYMBOL mapper ----
add_gene_symbol <- function(df, ensembl_col = "gene_id") {
  if (!ensembl_col %in% colnames(df)) return(df)

  ens <- sub("\\..*$", "", as.character(df[[ensembl_col]]))

  sym <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ens,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  df$gene_symbol <- unname(sym[ens])
  gi <- match(ensembl_col, colnames(df))
  df <- df[, c(seq_len(gi), ncol(df), setdiff(seq_len(ncol(df)-1), seq_len(gi))), drop = FALSE]
  df
}


# =========================
# Helper: results for pairwise contrasts
# =========================
groups <- levels(meta$group)

get_res <- function(g1, g2) {
  # g1 vs g2
  res <- results(dds, contrast = c("group", g1, g2))
  res <- as.data.frame(res)
  res$gene_id <- rownames(res)
  res
}

# =========================
# 2) Group-specific high/low DEG lists
# Definition (strict):
#   "G specific HIGH" = (G vs each other group) all padj<cut AND log2FC>lfc_cut
#   "G specific LOW"  = (G vs each other group) all padj<cut AND log2FC<-lfc_cut
# =========================
spec_dir <- file.path(outdir, "group_specific_DEG")
dir.create(spec_dir, showWarnings = FALSE, recursive = TRUE)

for (g in groups) {
  others <- setdiff(groups, g)
  res_list <- lapply(others, function(o) get_res(g, o))
  names(res_list) <- others

  # Build matrices for padj and lfc aligned by gene_id
  all_genes <- Reduce(intersect, lapply(res_list, function(x) x$gene_id))

  padj_mat <- sapply(res_list, function(x) {
    x <- x[match(all_genes, x$gene_id), ]
    x$padj
  })
  lfc_mat <- sapply(res_list, function(x) {
    x <- x[match(all_genes, x$gene_id), ]
    x$log2FoldChange
  })

  # strict criteria across all pairwise
# ---- relaxed criteria: supported by >= min_support comparisons ----
ok_high <- (!is.na(padj_mat) & padj_mat < padj_cut & !is.na(lfc_mat) & lfc_mat >= lfc_cut)
ok_low  <- (!is.na(padj_mat) & padj_mat < padj_cut & !is.na(lfc_mat) & lfc_mat <= -lfc_cut)

high_idx <- rowSums(ok_high) >= min_support
low_idx  <- rowSums(ok_low)  >= min_support

  high_genes <- all_genes[high_idx]
  low_genes  <- all_genes[low_idx]

  # Summarize with min padj and min lfc across comparisons (for reporting)
  make_summary <- function(sel_genes) {
    if (length(sel_genes) == 0) return(data.frame())
    padj_sub <- padj_mat[match(sel_genes, all_genes), , drop = FALSE]
    lfc_sub  <- lfc_mat[match(sel_genes, all_genes), , drop = FALSE]
    data.frame(
      gene_id = sel_genes,
      min_padj = rowMins(padj_sub, na.rm = TRUE),
      min_abs_lfc = apply(abs(lfc_sub), 1, min, na.rm = TRUE),
      mean_lfc = rowMeans(lfc_sub, na.rm = TRUE),
      stringsAsFactors = FALSE
    )[order(rowMins(padj_sub, na.rm = TRUE)), ]
  }

  high_df <- make_summary(high_genes)
  low_df  <- make_summary(low_genes)
  high_df <- add_gene_symbol(high_df, "gene_id")
  low_df  <- add_gene_symbol(low_df,  "gene_id")
  write.table(high_df,
              file = file.path(spec_dir, paste0("DEG_specific_HIGH_", g, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  write.table(low_df,
              file = file.path(spec_dir, paste0("DEG_specific_LOW_", g, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# =========================
# 3) Top500 DEG heatmap (with clustering)
# Approach:
#  - Overall group effect by LRT (reduced ~ 1)
#  - Take top500 by padj
#  - Plot heatmap on VST values (row-scaled)
# =========================
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
res_lrt <- results(dds_lrt)
res_lrt_df <- as.data.frame(res_lrt)
res_lrt_df$gene_id <- rownames(res_lrt_df)
res_lrt_df <- res_lrt_df[order(res_lrt_df$padj), ]
res_lrt_df <- add_gene_symbol(res_lrt_df, "gene_id")

topN <- 500
top_genes <- na.omit(res_lrt_df$gene_id)[1:min(topN, sum(!is.na(res_lrt_df$gene_id)))]

mat <- assay(vsd)[top_genes, ]
# row z-score
mat_z <- t(scale(t(mat)))

ann_col <- data.frame(group = meta$group)
rownames(ann_col) <- rownames(meta)

png(filename = file.path(outdir, "Heatmap_top500_DEG.png"), width = 1600, height = 1800, res = 200)
pheatmap(
  mat_z,
  annotation_col = ann_col,
  show_colnames = TRUE,
  show_rownames = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 500 DEG (LRT) - VST row Z-score"
)
dev.off()

# Save full LRT table too (handy)
write.table(res_lrt_df,
            file = file.path(outdir, "DESeq2_LRT_all_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Save session info
writeLines(capture.output(sessionInfo()), con = file.path(outdir, "sessionInfo.txt"))

message("Done. Outputs in: ", normalizePath(outdir))