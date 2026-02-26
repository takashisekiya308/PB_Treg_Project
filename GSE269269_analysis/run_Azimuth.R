#An R script to label integrated GSE269269 data with Azimuth

library(Seurat)
library(Azimuth)

obj <- readRDS("seurat_out/GSE269269_integrated_seurat.rds")
DefaultAssay(obj) <- "RNA"

obj <- JoinLayers(obj, assay = "RNA")

# Azimuth
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- RunAzimuth(obj, reference = "pbmcref")

pred_cols <- grep("^predicted", colnames(obj@meta.data), value = TRUE)
print(pred_cols)

if ("predicted.celltype.l2" %in% pred_cols) print(table(obj@meta.data[["predicted.celltype.l2"]]))
if ("predicted.celltype.l3" %in% pred_cols) print(table(obj@meta.data[["predicted.celltype.l3"]]))

saveRDS(obj, file = "seurat_out/GSE269269_azimuth_pbmc.rds")