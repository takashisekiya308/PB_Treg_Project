#!/usr/bin/env Rscript

# DiffBind sample sheet -> PCA plot
# - Color: tissue (PB forced to red)
# - Shape: Treg = filled circle, TN = open circle
#
# Usage:
#   Rscript diffbind_pca_tissue_celltype.R samplesheet.csv PCA.png
#
# Notes:
#   - The sample sheet must contain at least: Condition, SampleID, Replicate, bamReads, Peaks, PeakCaller
#   - Condition must look like "LN_TN" or "SP_Treg" (tissue_celltype)

suppressPackageStartupMessages({
  library(DiffBind)
  library(ggplot2)
  library(DESeq2)
})

args <- commandArgs(trailingOnly=TRUE)
samplesheet_path <- if (length(args) >= 1) args[1] else "samplesheet.csv"
out_png          <- if (length(args) >= 2) args[2] else "PCA_tissue_celltype.png"

if (!file.exists(samplesheet_path)) {
  stop("Sample sheet not found: ", samplesheet_path)
}

# --- Read sample sheet --------------------------------------------------------
ss <- read.csv(samplesheet_path, stringsAsFactors=FALSE, check.names=FALSE)

required_cols <- c("Condition","SampleID","Replicate","bamReads","Peaks","PeakCaller")
missing <- setdiff(required_cols, colnames(ss))
if (length(missing) > 0) {
  stop("Missing required column(s) in sample sheet: ", paste(missing, collapse=", "))
}

# Parse Condition -> Tissue / CellType
# Accept both "Sp_TN" and "SP_TN" (standardize Tissue to upper case)
parts <- strsplit(ss$Condition, "_", fixed=TRUE)
if (any(vapply(parts, length, integer(1)) < 2)) {
  stop("Condition must be 'TISSUE_CELLTYPE' like 'LN_TN' or 'SP_Treg'. Bad rows: ",
       paste(which(vapply(parts, length, integer(1)) < 2), collapse=", "))
}
ss$Tissue   <- toupper(vapply(parts, `[[`, character(1), 1))
ss$CellType <- vapply(parts, `[[`, character(1), 2)

# Standardize cell type labels a bit
ss$CellType <- ifelse(tolower(ss$CellType) %in% c("treg","tregs"), "Treg", ss$CellType)
ss$CellType <- ifelse(tolower(ss$CellType) %in% c("tn","tnaive","naive"), "TN", ss$CellType)

if (!all(ss$CellType %in% c("TN","Treg"))) {
  warning("CellType contains values other than TN/Treg: ",
          paste(sort(unique(ss$CellType)), collapse=", "),
          "\nShapes will be assigned only for TN and Treg.")
}

# --- Build DiffBind object & count -------------------------------------------
# This uses peaks provided in the sample sheet. For ATAC-seq BED peaks, PeakCaller='bed' is fine.
dba_obj <- dba(sampleSheet = ss)

# Counting can take time depending on BAM size / peakset size.
dba_obj <- dba.count(dba_obj, bUseSummarizeOverlaps=TRUE)

# --- Extract count matrix -----------------------------------------------------
# DiffBind's dba.getMatrix() sometimes returns a data.frame that includes genomic coordinates
# as leading columns (e.g., CHR/START/END). We strip those and keep only sample columns.

counts_raw <- NULL
counts_raw <- tryCatch({
  dba.getMatrix(dba_obj, bNormalized=FALSE)
}, error=function(e) NULL)

if (is.null(counts_raw)) {
  # DiffBind stores a counts/binding matrix internally as $binding in many versions
  if (!is.null(dba_obj$binding)) {
    counts_raw <- dba_obj$binding
  } else {
    stop("Could not extract count matrix. dba.getMatrix failed and dba_obj$binding is NULL.")
  }
}

# Prepare metadata keyed by SampleID
meta <- ss
rownames(meta) <- meta$SampleID

# If counts_raw is a data.frame with CHR/START/END (or similar), keep only sample columns
if (is.data.frame(counts_raw)) {
  sample_cols <- intersect(colnames(counts_raw), meta$SampleID)
  if (length(sample_cols) == 0) {
    stop("Counts matrix columns and SampleID mismatch.
",
         "Counts columns: ", paste(head(colnames(counts_raw)), collapse=", "), "...
",
         "SampleIDs: ", paste(head(meta$SampleID), collapse=", "), "...")
  }

  # Set rownames from genomic coordinates if present
  chr_candidates   <- c("CHR","CHROM","chr","seqnames","SEQNAME")
  start_candidates <- c("START","start")
  end_candidates   <- c("END","end")

  chr_col   <- intersect(chr_candidates, colnames(counts_raw))
  start_col <- intersect(start_candidates, colnames(counts_raw))
  end_col   <- intersect(end_candidates, colnames(counts_raw))

  if (length(chr_col) > 0 && length(start_col) > 0 && length(end_col) > 0) {
    rownames(counts_raw) <- paste(counts_raw[[chr_col[1]]], counts_raw[[start_col[1]]], counts_raw[[end_col[1]]], sep=":")
  }

  counts <- as.matrix(counts_raw[, sample_cols, drop=FALSE])
} else {
  counts <- as.matrix(counts_raw)
  if (is.null(colnames(counts))) {
    colnames(counts) <- dba_obj$samples$SampleID
  }
}

# Ensure columns are exactly sample IDs (drop any non-sample columns if present)
sample_cols <- intersect(colnames(counts), meta$SampleID)
if (length(sample_cols) == 0) {
  stop("Counts matrix columns and SampleID mismatch.
",
       "Counts columns: ", paste(head(colnames(counts)), collapse=", "), "...
",
       "SampleIDs: ", paste(head(meta$SampleID), collapse=", "), "...")
}
counts <- counts[, sample_cols, drop=FALSE]

# Reorder sample metadata to match matrix columns
meta <- meta[colnames(counts), , drop=FALSE]

# Basic sanity: numeric matrix
mode(counts) <- "numeric"

# Ensure matrix with rownames and colnames
counts <- as.matrix(counts)
if (is.null(colnames(counts))) {
  # Try to recover colnames from DiffBind samples
  colnames(counts) <- dba_obj$samples$SampleID
}

# Reorder sample metadata to match matrix columns
meta <- ss
rownames(meta) <- meta$SampleID
if (!all(colnames(counts) %in% rownames(meta))) {
  stop("Counts matrix columns and SampleID mismatch.\n",
       "Counts columns: ", paste(head(colnames(counts)), collapse=", "), "...\n",
       "SampleIDs: ", paste(head(meta$SampleID), collapse=", "), "...")
}
meta <- meta[colnames(counts), , drop=FALSE]

# Basic filtering to reduce noise (optional but usually helps PCA)
keep <- rowSums(counts) > 10
counts_f <- counts[keep, , drop=FALSE]
if (nrow(counts_f) < 100) {
  warning("After filtering, very few peaks remain (n=", nrow(counts_f), "). PCA may be unstable.")
}

# --- VST transform (DESeq2) & PCA --------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(counts_f),
                              colData   = meta,
                              design    = ~ 1)

vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
mat <- assay(vsd)

pca <- prcomp(t(mat), center=TRUE, scale.=FALSE)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

pca_df <- data.frame(
  SampleID = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Tissue = meta$Tissue,
  CellType = meta$CellType,
  stringsAsFactors=FALSE
)

# --- Colors (PB forced to red) -----------------------------------------------
tissues <- sort(unique(pca_df$Tissue))
# Nice qualitative palette from base R
pal <- setNames(grDevices::hcl.colors(length(tissues), palette="Dark 3"), tissues)
if ("PB" %in% tissues) pal["PB"] <- "red"

# Shapes: TN open circle, Treg filled circle
shape_map <- c("TN"=1, "Treg"=16)

# --- Plot --------------------------------------------------------------------
p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Tissue, shape=CellType)) +
  geom_point(size=4, stroke=1.2) +
  scale_color_manual(values=pal, drop=FALSE) +
  scale_shape_manual(values=shape_map, drop=FALSE) +
  labs(
    title = "PCA (DiffBind binding matrix; VST-transformed)",
    x = sprintf("PC1 (%.1f%%)", 100*var_expl[1]),
    y = sprintf("PC2 (%.1f%%)", 100*var_expl[2]),
    color = "Tissue",
    shape = "Cell type"
  ) +
  theme_classic(base_size=14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face="bold")
  )

ggsave(out_png, plot=p, width=8.5, height=6.5, dpi=200)
message("Saved PCA plot: ", out_png)
