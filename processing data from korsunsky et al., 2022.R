counts <- readRDS("/rds/projects/c/croftap-stia-atac/CM_multiome/data_from_Ilya/med_fibroblasts/exprs_raw.rds")
meta_data <- readRDS("/rds/projects/c/croftap-stia-atac/CM_multiome/data_from_Ilya/med_fibroblasts/meta_data.rds")

row.names(meta_data) <- meta_data$CellID
meta_data <- subset(meta_data, select = -CellID)

med_fibros<-CreateSeuratObject(counts = counts, meta.data = meta_data)
med_fibros<-NormalizeData(med_fibros)
med_fibros <- FindVariableFeatures(med_fibros, selection.method = "vst", nfeatures = 2000)
med_fibros <- ScaleData(med_fibros, verbose = FALSE)
med_fibros <- RunPCA(med_fibros, npcs = 30, verbose = FALSE)
med_fibros <- RunUMAP(med_fibros, reduction = "pca", dims = 1:30)
DimPlot(med_fibros, reduction = "umap", group.by = "Tissue")

Idents(med_fibros) <- 'Cluster_name'
all_markeers_med_firbos <- FindAllMarkers(med_fibros, only.pos = T)
write.csv(all_markeers_med_firbos, "/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/paper_figures/med_fibs_markers.csv")
