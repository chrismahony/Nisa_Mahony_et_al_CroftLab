analysis_sc_mm2 <- readRDS("/rds/projects/c/croftap-stia-atac/CM_multiome/data_from_Ilya/analysis_sc_mm2.rds")
counts<-analysis_sc_mm2[["exprs_norm"]]
metadata<-analysis_sc_mm2[["meta_data"]]
row.names(metadata) <- metadata$cell_id
metadata <- subset(metadata, select = -cell_id)

analysis_sc_mm2_seurat<-CreateSeuratObject(counts = counts, meta.data = metadata)
analysis_sc_mm2_seurat <- FindVariableFeatures(analysis_sc_mm2_seurat, selection.method = "vst", nfeatures = 2000)
analysis_sc_mm2_seurat <- ScaleData(analysis_sc_mm2_seurat, verbose = FALSE)
analysis_sc_mm2_seurat <- RunPCA(analysis_sc_mm2_seurat, npcs = 30, verbose = FALSE)
analysis_sc_mm2_seurat <- RunUMAP(analysis_sc_mm2_seurat, reduction = "pca", dims = 1:30)
DimPlot(analysis_sc_mm2_seurat, reduction = "umap", group.by = "label")

Idents(analysis_sc_mm2_seurat)<-'label_cm'
All.markers_fig3_notch<-FindAllMarkers(analysis_sc_mm2_seurat, only.pos = T)
