#!/usr/bin/env Rscript
#An R script used for generating PCA plots of Treg-phagocytosed BMDMs. 
#usage: move to a folder which contains samplesheet.txt for DESeq2, and source this R file.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
samplesheet_path <- if (length(args) >= 1) args[1] else "samplesheet.txt"
outdir <- if (length(args) >= 2) args[2] else "."

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(samplesheet_path)) stop("Sample sheet not found: ", samplesheet_path)

ss <- read.delim(samplesheet_path, stringsAsFactors = FALSE, check.names = FALSE)
req <- c("sample","group","path")
miss <- setdiff(req, colnames(ss))
if (length(miss) > 0) stop("Missing column(s) in sample sheet: ", paste(miss, collapse=", "))

# Preflight file existence
missing_files <- which(is.na(ss$path) | !file.exists(ss$path))
if (length(missing_files) > 0) {
  bad <- data.frame(sample=ss$sample[missing_files], path=ss$path[missing_files], stringsAsFactors=FALSE)
  stop("Some quant files do not exist (or are NA). Examples:\n",
       paste(capture.output(print(head(bad, 10), row.names=FALSE)), collapse="\n"),
       if (nrow(bad) > 10) sprintf("\n... plus %d more", nrow(bad)-10) else "")
}

# ---- Read counts (RSEM genes.results-compatible) -----------------------------
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

  out <- data.frame(
    gene_id = df[[gene_col]],
    count   = df[[count_col]],
    stringsAsFactors=FALSE
  )
  out
}

message("Reading quant files...")
first <- read_counts_one(ss$path[1])
genes <- first$gene_id

counts <- matrix(0, nrow=length(genes), ncol=nrow(ss))
rownames(counts) <- genes
colnames(counts) <- ss$sample

for (i in seq_len(nrow(ss))) {
  df <- read_counts_one(ss$path[i])
  idx <- match(genes, df$gene_id)
  v <- df$count[idx]
  v[is.na(v)] <- 0
  counts[, i] <- v
}
counts <- round(counts)

# ---- DESeq2 VST + PCA -------------------------------------------------------
coldata <- data.frame(group=ss$group, row.names=ss$sample, stringsAsFactors=FALSE)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ 1)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
pca <- plotPCA(vsd, intgroup="group", returnData=TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

# ---- Color mapping (requested scheme) ----------------------------------------
# - M0_plus_AC: black
# - IFNg_Treg_421 & IFNg_Treg_647: red family
# - IFNg_Tconv_421 & IFNg_Tconv_647: blue family
# - IFNg: clearly distinct from red/blue/black (dark green)
requested_pal <- c(
  "M0_plus_AC"     = "black",
  "IFNg_Treg_421"  = "firebrick2",
  "IFNg_Treg_647"  = "indianred3",
  "IFNg_Tconv_421" = "dodgerblue3",
  "IFNg_Tconv_647" = "deepskyblue3",
  "IFNg"           = "darkgreen"
)

groups <- sort(unique(pca$group))
pal <- requested_pal

# If there are any unexpected groups, assign additional distinct colors
extra <- setdiff(groups, names(pal))
if (length(extra) > 0) {
  extra_cols <- grDevices::hcl.colors(length(extra), palette="Dark 3")
  names(extra_cols) <- extra
  pal <- c(pal, extra_cols)
}
pal <- pal[names(pal) %in% groups]

p <- ggplot(pca, aes(x=PC1, y=PC2, color=group)) +
  geom_point(size=4, alpha=0.95) +
  scale_color_manual(values=pal, drop=FALSE) +
  theme_classic(base_size=14) +
  labs(
    title="PCA (VST, DESeq2)",
    x=paste0("PC1: ", percentVar[1], "% variance"),
    y=paste0("PC2: ", percentVar[2], "% variance"),
    color="group"
  ) +
  theme(plot.title = element_text(face="bold"))

out_png <- file.path(outdir, "PCA_by_group.png")
ggsave(out_png, p, width=8.5, height=6.5, dpi=200)
message("Saved: ", out_png)

# ---- Second PCA plot: aspect reflects variance explained (robust) -------------
# Goal: make the *physical axis lengths* roughly follow PC2:PC1 = percentVar[2]:percentVar[1],
# without over-squashing the panel.
#
# With coord_fixed(ratio=R): (panel_height / panel_width) = R * (y_range / x_range)
# So to target panel_height/panel_width = percentVar[2]/percentVar[1], set:
#   R = (percentVar[2]/percentVar[1]) * (x_range / y_range)

x_rng <- diff(range(pca$PC1, na.rm=TRUE))
y_rng <- diff(range(pca$PC2, na.rm=TRUE))
if (!is.finite(x_rng) || x_rng == 0) x_rng <- 1
if (!is.finite(y_rng) || y_rng == 0) y_rng <- 1

target_hw <- percentVar[2] / percentVar[1]     # desired (height/width)
R <- target_hw * (x_rng / y_rng)

# Add a bit of padding so points never touch borders
pad_mult <- 0.06

p_aspect <- ggplot(pca, aes(x=PC1, y=PC2, color=group)) +
  geom_point(size=4, alpha=0.95) +
  scale_color_manual(values=pal, drop=FALSE) +
  scale_x_continuous(expand = expansion(mult = pad_mult)) +
  scale_y_continuous(expand = expansion(mult = pad_mult)) +
  coord_fixed(ratio = R) +
  theme_classic(base_size=14) +
  labs(
    title=paste0("PCA (VST, DESeq2) — axis length ~ PC2:PC1 = ", percentVar[2], ":", percentVar[1]),
    x=paste0("PC1: ", percentVar[1], "% variance"),
    y=paste0("PC2: ", percentVar[2], "% variance"),
    color="group"
  ) +
  theme(plot.title = element_text(face="bold"))

out_png2 <- file.path(outdir, "PCA_by_group_aspect_by_variance.png")
ggsave(out_png2, p_aspect, width=10.5, height=6.5, dpi=200)
message("Saved: ", out_png2)
