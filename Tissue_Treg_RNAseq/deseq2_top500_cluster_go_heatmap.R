#!/usr/bin/env Rscript
# Independent script:
# 1) Take top 500 DEGs across groups (DESeq2 LRT) from samplesheet_RNAseq_Treg.txt
# 2) Cluster those genes and extract clusters
# 3) Run GO BP enrichment for each cluster
# 4) Output a heatmap PNG for top 500 with cluster regions demarcated
#
# Input samplesheet (tab-delimited) columns:
#   sample   group   path
# where 'path' points to quant files (e.g., RSEM *.genes.results)
#
# Usage:
#   Rscript deseq2_top500_cluster_go_heatmap.R samplesheet_RNAseq_Treg.txt [outdir]
#
# Outputs (in outdir):
#   Heatmap_DEG_top500_clusters.png
#   Top500_DEG_LRT.tsv
#   Cluster_assignments_top500.tsv
#   Cluster_<ID>_genes.txt
#   GO_cluster<ID>_top20.tsv
#   GO_cluster<ID>_top20.png
#
# Notes:
# - DEGs: DESeq2 LRT (design ~ group vs ~1), top 500 by padj
# - Clustering: hierarchical on VST values (row-scaled), dynamic tree cut if available;
#   otherwise fallback to cutree(k=6).
# - GO: clusterProfiler::enrichGO (BP) on each cluster gene set.

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---- Package helpers ---------------------------------------------------------
ensure_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
ensure_bioc <- function(pkg) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

ensure_bioc("DESeq2")
ensure_cran("pheatmap")
ensure_bioc("clusterProfiler")
ensure_bioc("AnnotationDbi")
ensure_cran("dynamicTreeCut")
ensure_cran("stringr")
ensure_cran("scales") # optional but we attempt to use it

suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(dynamicTreeCut)
  library(grid)
  library(stringr)
  library(scales)
})

# ---- Args -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
samplesheet_path <- if (length(args) >= 1) args[1] else "samplesheet_RNAseq_Treg.txt"
outdir <- if (length(args) >= 2) args[2] else "."

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(samplesheet_path)) stop("Sample sheet not found: ", samplesheet_path)

ss <- read.delim(samplesheet_path, stringsAsFactors = FALSE, check.names = FALSE)
req <- c("sample","group","path")
miss <- setdiff(req, colnames(ss))
if (length(miss) > 0) stop("Missing column(s) in sample sheet: ", paste(miss, collapse=", "))

# ---- Preflight checks --------------------------------------------------------
missing_files <- which(is.na(ss$path) | !file.exists(ss$path))
if (length(missing_files) > 0) {
  bad <- data.frame(sample=ss$sample[missing_files], path=ss$path[missing_files], stringsAsFactors=FALSE)
  stop("Some quant files do not exist (or are NA). Examples:\n",
       paste(capture.output(print(head(bad, 10), row.names=FALSE)), collapse="\n"),
       if (nrow(bad) > 10) sprintf("\n... plus %d more", nrow(bad)-10) else "")
}
finfo <- file.info(ss$path)
bad_size <- which(is.na(finfo$size) | finfo$size <= 0)
if (length(bad_size) > 0) {
  bad <- data.frame(sample=ss$sample[bad_size], path=ss$path[bad_size], size=finfo$size[bad_size], stringsAsFactors=FALSE)
  stop("Some quant files are empty/unreadable (size NA or <=0). Examples:\n",
       paste(capture.output(print(head(bad, 10), row.names=FALSE)), collapse="\n"),
       if (nrow(bad) > 10) sprintf("\n... plus %d more", nrow(bad)-10) else "")
}

# ---- Read counts (RSEM genes.results compatible) -----------------------------
read_counts_one <- function(file) {
  df <- read.delim(file, stringsAsFactors=FALSE, check.names=FALSE)

  gene_col <- if ("gene_id" %in% colnames(df)) "gene_id" else if ("gene" %in% colnames(df)) "gene" else colnames(df)[1]
  count_col <- NULL
  for (cand in c("expected_count","counts","count","num_reads")) {
    if (cand %in% colnames(df)) { count_col <- cand; break }
  }
  if (is.null(count_col)) {
    stop("Could not find a count column in: ", file, "\nAvailable columns: ", paste(colnames(df), collapse=", "))
  }

  symbol_col <- NULL
  for (cand in c("gene_name","gene_symbol","symbol","gene")) {
    if (cand %in% colnames(df) && cand != gene_col) { symbol_col <- cand; break }
  }

  out <- data.frame(gene_id=df[[gene_col]], count=df[[count_col]], stringsAsFactors=FALSE)
  if (!is.null(symbol_col)) out$gene_symbol <- df[[symbol_col]]
  out
}

message("Reading quant files...")
first <- read_counts_one(ss$path[1])
genes <- first$gene_id

counts <- matrix(0, nrow=length(genes), ncol=nrow(ss))
rownames(counts) <- genes
colnames(counts) <- ss$sample

gene_symbol <- if ("gene_symbol" %in% colnames(first)) first$gene_symbol else rep(NA_character_, length(genes))
names(gene_symbol) <- genes

for (i in seq_len(nrow(ss))) {
  df <- read_counts_one(ss$path[i])
  idx <- match(genes, df$gene_id)
  v <- df$count[idx]
  v[is.na(v)] <- 0
  counts[, i] <- v
  if (any(is.na(gene_symbol)) && "gene_symbol" %in% colnames(df)) {
    tmp <- df$gene_symbol[idx]
    gene_symbol[is.na(gene_symbol)] <- tmp[is.na(gene_symbol)]
  }
}
counts <- round(counts)

# ---- Symbol mapping fallback (mouse/human) -----------------------------------
detect_org <- function(gene_ids) {
  if (any(grepl("^ENSMUSG", gene_ids))) return("mouse")
  if (any(grepl("^ENSG", gene_ids))) return("human")
  return("unknown")
}
org <- detect_org(genes)

if (all(is.na(gene_symbol))) {
  if (org == "mouse") {
    ensure_bioc("org.Mm.eg.db")
    suppressPackageStartupMessages({ library(org.Mm.eg.db) })
    gene_symbol <- mapIds(org.Mm.eg.db, keys=genes, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
    gene_symbol <- as.character(gene_symbol); names(gene_symbol) <- genes
  } else if (org == "human") {
    ensure_bioc("org.Hs.eg.db")
    suppressPackageStartupMessages({ library(org.Hs.eg.db) })
    gene_symbol <- mapIds(org.Hs.eg.db, keys=genes, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
    gene_symbol <- as.character(gene_symbol); names(gene_symbol) <- genes
  } else {
    gene_symbol <- genes
    names(gene_symbol) <- genes
    warning("Could not infer organism from gene IDs, and gene_symbol column not found. Using gene_id as label.")
  }
}
gene_symbol[is.na(gene_symbol) | !nzchar(gene_symbol)] <- names(gene_symbol)[is.na(gene_symbol) | !nzchar(gene_symbol)]

# ---- DESeq2 LRT and top500 ---------------------------------------------------
coldata <- data.frame(group=factor(ss$group), row.names=ss$sample, stringsAsFactors=FALSE)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ group)

dds_lrt <- DESeq(dds, test="LRT", reduced= ~ 1)
res_lrt <- results(dds_lrt)
res_df <- as.data.frame(res_lrt)
res_df$gene_id <- rownames(res_df)
res_df$gene <- gene_symbol[res_df$gene_id]
res_df <- res_df[order(res_df$padj, na.last=TRUE), , drop=FALSE]
res_df <- res_df[!is.na(res_df$padj), , drop=FALSE]

topN <- 500
top_deg <- head(res_df, topN)

write.table(top_deg, file.path(outdir, "Top500_DEG_LRT.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
message("Saved: ", file.path(outdir, "Top500_DEG_LRT.tsv"))

# ---- VST and matrix for clustering/heatmap -----------------------------------
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
mat <- assay(vsd)[top_deg$gene_id, , drop=FALSE]
rownames(mat) <- top_deg$gene

# row-scale for clustering/heat display
mat_z <- t(scale(t(mat)))
mat_z[is.na(mat_z)] <- 0

# ---- 1) Gene clustering and cluster extraction -------------------------------
# Hierarchical clustering
d_rows <- as.dist(1 - cor(t(mat_z), method="pearson"))
hc <- hclust(d_rows, method="average")

# Use dynamic tree cut if possible, else fallback to fixed k
clusters <- NULL
if (requireNamespace("dynamicTreeCut", quietly = TRUE)) {
  clusters <- dynamicTreeCut::cutreeDynamic(
    dendro = hc,
    distM = as.matrix(d_rows),
    deepSplit = 2,
    pamRespectsDendro = TRUE,
    minClusterSize = 15
  )
  clusters[clusters == 0] <- NA
  if (all(is.na(clusters))) clusters <- NULL
}

if (is.null(clusters)) {
  k <- 6
  clusters <- cutree(hc, k = k)
  message("dynamicTreeCut not used; fallback to cutree(k=", k, ")")
} else {
  message("dynamicTreeCut used. Found clusters: ", paste(sort(unique(na.omit(clusters))), collapse=", "))
}

# Ensure clusters ordered by appearance in dendrogram
ord <- hc$order
labs <- rownames(mat_z)[ord]
cl_ord <- clusters[match(labs, rownames(mat_z))]
uniq <- unique(cl_ord)
uniq <- uniq[!is.na(uniq)]
cluster_levels <- uniq
clusters_f <- factor(clusters, levels = cluster_levels)
clusters_num <- as.integer(as.character(clusters_f))
names(clusters_num) <- rownames(mat_z)

# Save cluster assignments
cl_out <- data.frame(
  gene = rownames(mat_z),
  cluster = clusters_num,
  stringsAsFactors = FALSE
)
cl_out <- cl_out[order(cl_out$cluster, cl_out$gene), ]
write.table(cl_out, file.path(outdir, "Cluster_assignments_top500.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
message("Saved: ", file.path(outdir, "Cluster_assignments_top500.tsv"))

# Write per-cluster gene lists
for (cid in sort(unique(na.omit(clusters_num)))) {
  genes_c <- cl_out$gene[cl_out$cluster == cid]
  writeLines(genes_c, con = file.path(outdir, paste0("Cluster_", cid, "_genes.txt")))
}

# ---- 2) GO enrichment per cluster --------------------------------------------
# Pick OrgDb
detect_org_from_ids <- function(gene_ids) {
  if (any(grepl("^ENSMUSG", gene_ids))) return("mouse")
  if (any(grepl("^ENSG", gene_ids))) return("human")
  return("unknown")
}
org2 <- detect_org_from_ids(top_deg$gene_id)

if (org2 == "mouse") {
  ensure_bioc("org.Mm.eg.db")
  suppressPackageStartupMessages({ library(org.Mm.eg.db) })
  orgdb <- org.Mm.eg.db
} else if (org2 == "human") {
  ensure_bioc("org.Hs.eg.db")
  suppressPackageStartupMessages({ library(org.Hs.eg.db) })
  orgdb <- org.Hs.eg.db
} else {
  ensure_bioc("org.Mm.eg.db")
  suppressPackageStartupMessages({ library(org.Mm.eg.db) })
  orgdb <- org.Mm.eg.db
  warning("Organism unknown; defaulting to org.Mm.eg.db for GO mapping.")
}

sym_to_entrez <- function(sym_vec) {
  sym_vec <- unique(na.omit(as.character(sym_vec)))
  sym_vec <- sym_vec[nzchar(sym_vec)]
  entrez <- mapIds(orgdb, keys=sym_vec, column="ENTREZID", keytype="SYMBOL", multiVals="first")
  unique(na.omit(as.character(entrez)))
}

run_go <- function(symbols, cid) {
  symbols <- unique(na.omit(as.character(symbols)))
  symbols <- symbols[nzchar(symbols)]

  entrez <- sym_to_entrez(symbols)
  mapped_n <- length(entrez)
  total_n <- length(symbols)

  # Prepare output paths up front
  tsv <- file.path(outdir, paste0("GO_cluster", cid, "_top20.tsv"))
  pngf <- file.path(outdir, paste0("GO_cluster", cid, "_top20.png"))

  if (mapped_n < 10) {
    # Always emit files, even if GO can't run
    msg <- sprintf("Cluster %s: too few mapped genes for GO (mapped ENTREZ n=%d / symbols n=%d).", cid, mapped_n, total_n)
    warning(msg)
    write.table(
      data.frame(message=msg, stringsAsFactors=FALSE),
      tsv, sep="	", quote=FALSE, row.names=FALSE
    )
    p <- ggplot(data.frame(x=1, y=1, label=msg), aes(x, y, label=label)) +
      geom_text(size=5) +
      theme_void() +
      labs(title=paste0("GO BP enrichment — cluster ", cid))
    ggsave(pngf, p, width=13, height=7.5, dpi=200, limitsize = FALSE)
    message("Saved: ", tsv)
    message("Saved: ", pngf)
    return(invisible(NULL))
  }

  ego <- enrichGO(
    gene          = entrez,
    OrgDb         = orgdb,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  ego_df <- as.data.frame(ego)

  if (is.null(ego_df) || nrow(ego_df) == 0) {
    msg <- sprintf("Cluster %s: no significant GO BP terms (mapped ENTREZ n=%d / symbols n=%d).", cid, mapped_n, total_n)
    warning(msg)
    write.table(
      data.frame(message=msg, stringsAsFactors=FALSE),
      tsv, sep="	", quote=FALSE, row.names=FALSE
    )
    p <- ggplot(data.frame(x=1, y=1, label=msg), aes(x, y, label=label)) +
      geom_text(size=5) +
      theme_void() +
      labs(title=paste0("GO BP enrichment — cluster ", cid))
    ggsave(pngf, p, width=13, height=7.5, dpi=200, limitsize = FALSE)
    message("Saved: ", tsv)
    message("Saved: ", pngf)
    return(invisible(NULL))
  }

  ego_df <- ego_df[order(ego_df$p.adjust, ego_df$pvalue), , drop=FALSE]
  # Optional readability filter: drop very long GO terms (>=9 words)
  ego_df$word_n <- sapply(strsplit(as.character(ego_df$Description), "\s+"), length)
  ego_df_f <- ego_df[ego_df$word_n <= 8, , drop=FALSE]
  if (nrow(ego_df_f) >= 5) {
    ego_df <- ego_df_f
  }
  top <- head(ego_df, 15)
  top$mlog10 <- -log10(pmax(top$p.adjust, 1e-300))
  top$mlog10_plot <- (2/3) * top$mlog10
  top$Description_wrapped <- stringr::str_wrap(as.character(top$Description), width = 40)
top$Description_wrapped <- factor(top$Description_wrapped, levels = rev(unique(top$Description_wrapped)))

# Save TSV
  write.table(top, tsv, sep="	", quote=FALSE, row.names=FALSE)

  # Barplot
  local_max <- max(top$mlog10_plot, na.rm=TRUE)
  int_breaks <- seq(0, ceiling(local_max), by=1)

  p <- ggplot(top, aes(x=Description_wrapped, y=mlog10_plot)) +
    geom_col() +
    coord_flip() +
    scale_y_continuous(breaks = int_breaks, labels = int_breaks) +
    theme_classic(base_size=22) +
    theme(
      text = element_text(family = "Arial"),
      axis.text.y = element_text(size = 28, lineheight = 0.85),
      plot.margin = margin(15, 15, 15, 70)
    ) +
    labs(
      title=paste0("GO BP enrichment — cluster ", cid),
      x=NULL,
      y="scaled -log10(adj p)"
    )

  plot_h <-plot_h <- max(10, 0.75 * nrow(top) + 5)
  plot_w <- 16
  ggsave(pngf, p, width=plot_w, height=plot_h, dpi=200, limitsize = FALSE)

  message("Saved: ", tsv)
  message("Saved: ", pngf)
  invisible(list(go=ego, top=top))
}

for (cid in sort(unique(na.omit(clusters_num)))) {
  symbols <- cl_out$gene[cl_out$cluster == cid]
  run_go(symbols, cid)
}

# ---- 3) Heatmap (top500) with cluster blocks --------------------------------
# Determine row order and gaps between clusters in dendrogram order
row_order <- hc$order
ordered_genes <- rownames(mat_z)[row_order]
ordered_clusters <- clusters_num[ordered_genes]

# compute gaps row indices where cluster changes
gaps <- which(ordered_clusters[-1] != ordered_clusters[-length(ordered_clusters)])

# Row annotation for cluster id
ann_row <- data.frame(cluster = factor(ordered_clusters))
rownames(ann_row) <- ordered_genes
# ---- Cluster color palette (high-contrast) -----------------------------------
cluster_ids <- sort(unique(ordered_clusters))

cluster_colors <- c(
  "1"  = "black",
  "2"  = "red",
  "3"  = "gold",
  "4"  = "blue",
  "5"  = "darkgreen",
  "6"  = "orange",
  "7"  = "purple",
  "8"  = "cyan",
  "9"  = "brown",
  "10" = "deeppink"
)

# Safety: if cluster count exceeds prepared colors, append rainbow colors
if (length(cluster_ids) > length(cluster_colors)) {
  extra_cols <- grDevices::rainbow(length(cluster_ids) - length(cluster_colors))
  names(extra_cols) <- setdiff(as.character(cluster_ids), names(cluster_colors))
  cluster_colors <- c(cluster_colors, extra_cols)
}

cluster_palette <- cluster_colors[as.character(cluster_ids)]


# Column annotation (group)
group_levels <- levels(coldata$group)
group_pal <- setNames(grDevices::hcl.colors(length(group_levels), "Dark 3"), group_levels)

# Build pheatmap grob first, then size device from grob to avoid cropping
ph <- pheatmap(
  mat_z,
  cluster_rows = hc,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = data.frame(group=coldata$group, row.names=rownames(coldata)),
  annotation_row = ann_row,
  annotation_colors = list(
    group = group_pal,
    cluster = cluster_palette
  ),
  gaps_row = gaps,
  fontsize_row = 6,
  silent = TRUE
)

g <- ph$gtable
w_in <- convertWidth(sum(g$widths), unitTo="in", valueOnly=TRUE) * 1.03
h_in <- convertHeight(sum(g$heights), unitTo="in", valueOnly=TRUE) * 1.03
res <- 160

heat_png <- file.path(outdir, "Heatmap_DEG_top500_clusters.png")
png(heat_png, width = w_in * res, height = h_in * res, res = res)
grid.newpage(); grid.draw(g); dev.off()
message("Saved: ", heat_png, " (computed size ", round(w_in,2), "x", round(h_in,2), " inches @", res, "dpi)")

message("Done.")
