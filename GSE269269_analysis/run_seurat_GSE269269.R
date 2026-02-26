#An R script for integration of GSE269269 data

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

library(future)
options(future.globals.maxSize = 64 * 1024^3)
plan("sequential") 

root <- "."                   
outdir <- "seurat_out"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- samples autodetect ----------
sample_dirs <- list.dirs(root, full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[basename(sample_dirs) %in% c(paste0("PR", 1:5), paste0("NPR", 1:5))]

stopifnot(length(sample_dirs) > 0)

# ---------- helper ----------
read_one <- function(d) {
  sample <- basename(d)
  condition <- ifelse(grepl("^PR", sample), "PR", "NPR")

  # ---- find files with prefixes (e.g. GSMxxxx_NPR1_*) ----
  bc <- list.files(d, pattern = "barcodes\\.tsv\\.gz$", full.names = TRUE)
  ft <- list.files(d, pattern = "features\\.tsv\\.gz$", full.names = TRUE)
  mt <- list.files(d, pattern = "matrix\\.mtx\\.gz$", full.names = TRUE)

  if (length(bc) != 1 || length(ft) != 1 || length(mt) != 1) {
    stop("Could not uniquely find 10x files in: ", d,
         "\nbarcodes candidates: ", paste(bc, collapse = ", "),
         "\nfeatures candidates: ", paste(ft, collapse = ", "),
         "\nmatrix candidates: ", paste(mt, collapse = ", "))
  }

  # ---- ReadMtx handles arbitrary filenames ----
  mat <- ReadMtx(
    mtx = mt[1],
    features = ft[1],
    cells = bc[1],
    feature.column = 2   # features.tsv.gz: col2 = gene symbol
  )

  obj <- CreateSeuratObject(mat, project = "GSE269269", min.cells = 3, min.features = 200)
  obj$sample <- sample
  obj$condition <- condition

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
}

objs <- lapply(sample_dirs, read_one)
names(objs) <- basename(sample_dirs)

# ---------- QC plots (before filtering) ----------
qc_pdf <- file.path(outdir, "QC_violin_before_filtering.pdf")
pdf(qc_pdf, width = 11, height = 6)
print(VlnPlot(merge(objs[[1]], objs[-1]), features = c("nFeature_RNA","nCount_RNA","percent.mt"), group.by = "sample", ncol = 3))
dev.off()

# ---------- filtering (adjust if needed) ----------
filter_one <- function(obj) {
  subset(obj,
         subset = nFeature_RNA >= 200 &
                  nCount_RNA >= 500 &
                  percent.mt <= 25)
}
objs <- lapply(objs, filter_one)

# ---------- SCTransform integration (robust default) ----------
objs <- lapply(objs, function(x) SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes=TRUE))

features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)
objs <- PrepSCTIntegration(object.list = objs, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = objs, normalization.method = "SCT",
                                  anchor.features = features, verbose = FALSE)
integ <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)

# ---------- dimensional reduction / clustering ----------
integ <- RunPCA(integ, verbose = FALSE)
integ <- RunUMAP(integ, dims = 1:30, verbose = FALSE)
integ <- FindNeighbors(integ, dims = 1:30, verbose = FALSE)
integ <- FindClusters(integ, resolution = 0.5, verbose = FALSE)

# ---------- plots ----------
p1 <- DimPlot(integ, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("Clusters")
p2 <- DimPlot(integ, reduction = "umap", group.by = "condition") + ggtitle("PR vs NPR")
p3 <- DimPlot(integ, reduction = "umap", group.by = "sample") + ggtitle("Samples")
ggsave(file.path(outdir, "UMAP_clusters.png"), p1, width = 7, height = 6, dpi = 300)
ggsave(file.path(outdir, "UMAP_condition.png"), p2, width = 7, height = 6, dpi = 300)
ggsave(file.path(outdir, "UMAP_samples.png"), p3, width = 9, height = 6, dpi = 300)

# ---------- marker genes per cluster ----------
DefaultAssay(integ) <- "RNA"
integ <- JoinLayers(integ, assay = "RNA")

integ <- NormalizeData(integ, verbose = FALSE)
integ <- FindVariableFeatures(integ, verbose = FALSE)
integ <- ScaleData(integ, verbose = FALSE)

markers <- FindAllMarkers(integ, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

markers <- FindAllMarkers(integ, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers, file = file.path(outdir, "Markers_by_cluster.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- PR vs NPR differential (global) ----------
Idents(integ) <- "condition"
deg_global <- FindMarkers(integ, ident.1 = "PR", ident.2 = "NPR", min.pct = 0.1, logfc.threshold = 0.25)
deg_global$gene <- rownames(deg_global)
write.table(deg_global, file = file.path(outdir, "DEG_PR_vs_NPR_global.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- PR vs NPR differential within each cluster ----------
Idents(integ) <- "seurat_clusters"
clusters <- levels(Idents(integ))

deg_list <- list()
for (cl in clusters) {
  sub <- subset(integ, idents = cl)
  if (length(unique(sub$condition)) < 2) next

  Idents(sub) <- "condition"
  deg <- FindMarkers(sub, ident.1 = "PR", ident.2 = "NPR", min.pct = 0.1, logfc.threshold = 0.25)
  deg$gene <- rownames(deg)
  deg$cluster <- cl
  deg_list[[cl]] <- deg
}
deg_by_cluster <- do.call(rbind, deg_list)
write.table(deg_by_cluster, file = file.path(outdir, "DEG_PR_vs_NPR_by_cluster.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- save Seurat object ----------
saveRDS(integ, file = file.path(outdir, "GSE269269_integrated_seurat.rds"))

message("Done. Outputs in: ", normalizePath(outdir))