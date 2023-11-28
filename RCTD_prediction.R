#read in notch raw data n(from authros)
analysis_sc_mm2 <- readRDS("/rds/projects/c/croftap-lab-data-2021/CM/analysis_sc_mm2.rds")
metadata_notch<-analysis_sc_mm2[["meta_data"]]
metadata_notch<-as.data.frame(metadata_notch)
library(tidyverse)
metadata_notch<-metadata_notch %>% remove_rownames %>% column_to_rownames(var="cell_id")
counts_notch<-analysis_sc_mm2[["exprs_raw"]]
head(counts_notch)

Notch <- CreateSeuratObject(counts = counts_notch, min.cells = 3, min.features = 200, meta.data = metadata_notch)
Notch$dataset<-'Notch'
Notch$type_cm<-Notch$cell_type_nice
current.sample.ids <- c( "Endothelial cells"   ,             "NOTCH- fibroblasts"    ,          "NOTCH+ fibroblasts"      ,         "proliferating fibroblasts")
new.sample.ids <- c("Endothelial cells"   ,             "Fibroblast"    ,          "Fibroblast"      ,         "Fibroblast")
Notch@meta.data[["type_cm"]] <- plyr::mapvalues(x = Notch@meta.data[["type_cm"]], from = current.sample.ids, to = new.sample.ids)
table(Notch@meta.data[["type_cm"]])


amp <- readRDS("~/croftap-stia-atac-path/Visium/NotchData/amp.rds")
amp$dataset<-'amp'

Amp_notch_merge<-merge(x=amp, y=Notch)

#make RCTD ref
counts<-Amp_notch_merge@assays[["RNA"]]@counts
head(counts)


metaData<-Amp_notch_merge@meta.data
cell_types <- metaData$type_cm; names(cell_types) <- rownames(metaData)
head(cell_types)
cell_types <- as.factor(cell_types)
reference <- Reference(counts, cell_types)


load("~/croftap-stia-atac-path/Visium/NotchData/merging.RData")
rm(list= ls()[!(ls() %in% c('merging','reference'))])
#making RCTD object
spatialcounts<-merging@assays[["Spatial"]]@counts
slice1<-merging@images[["slice1"]]
slice1.coordinates<-slice1@coordinates
slice1.coordinates = subset(slice1.coordinates, select = -c(1:3) )
library(tibble)
slice1.coordinates <- tibble::rownames_to_column(slice1.coordinates, "barcodes")

slice1.1<-merging@images[["slice1.1"]]
slice1.1.coordinates<-slice1.1@coordinates
slice1.1.coordinates = subset(slice1.1.coordinates, select = -c(1:3) )
library(tibble)
slice1.1.coordinates <- tibble::rownames_to_column(slice1.1.coordinates, "barcodes")

slice1.2<-merging@images[["slice1.2"]]
slice1.2.coordinates<-slice1.2@coordinates
slice1.2.coordinates = subset(slice1.2.coordinates, select = -c(1:3) )
library(tibble)
slice1.2.coordinates <- tibble::rownames_to_column(slice1.2.coordinates, "barcodes")

slice1.3<-merging@images[["slice1.3"]]
slice1.3.coordinates<-slice1.3@coordinates
slice1.3.coordinates = subset(slice1.3.coordinates, select = -c(1:3) )
library(tibble)
slice1.3.coordinates <- tibble::rownames_to_column(slice1.3.coordinates, "barcodes")

slice1.4<-merging@images[["slice1.4"]]
slice1.4.coordinates<-slice1.4@coordinates
slice1.4.coordinates = subset(slice1.4.coordinates, select = -c(1:3) )
library(tibble)
slice1.4.coordinates <- tibble::rownames_to_column(slice1.4.coordinates, "barcodes")

slice1.5<-merging@images[["slice1.5"]]
slice1.5.coordinates<-slice1.5@coordinates
slice1.5.coordinates = subset(slice1.5.coordinates, select = -c(1:3) )
library(tibble)
slice1.5.coordinates <- tibble::rownames_to_column(slice1.5.coordinates, "barcodes")

slice1.6<-merging@images[["slice1.6"]]
slice1.6.coordinates<-slice1.6@coordinates
slice1.6.coordinates = subset(slice1.6.coordinates, select = -c(1:3) )
library(tibble)
slice1.6.coordinates <- tibble::rownames_to_column(slice1.6.coordinates, "barcodes")

slice1.7<-merging@images[["slice1.7"]]
slice1.7.coordinates<-slice1.7@coordinates
slice1.7.coordinates = subset(slice1.7.coordinates, select = -c(1:3) )
library(tibble)
slice1.7.coordinates <- tibble::rownames_to_column(slice1.7.coordinates, "barcodes")

slice1.8<-merging@images[["slice1.8"]]
slice1.8.coordinates<-slice1.8@coordinates
slice1.8.coordinates = subset(slice1.8.coordinates, select = -c(1:3) )
library(tibble)
slice1.8.coordinates <- tibble::rownames_to_column(slice1.8.coordinates, "barcodes")

slice1.9<-merging@images[["slice1.9"]]
slice1.9.coordinates<-slice1.9@coordinates
slice1.9.coordinates = subset(slice1.9.coordinates, select = -c(1:3) )
library(tibble)
slice1.9.coordinates <- tibble::rownames_to_column(slice1.9.coordinates, "barcodes")

slice1.10<-merging@images[["slice1.10"]]
slice1.10.coordinates<-slice1.10@coordinates
slice1.10.coordinates = subset(slice1.10.coordinates, select = -c(1:3) )
library(tibble)
slice1.10.coordinates <- tibble::rownames_to_column(slice1.10.coordinates, "barcodes")

slice1.11<-merging@images[["slice1.11"]]
slice1.11.coordinates<-slice1.11@coordinates
slice1.11.coordinates = subset(slice1.11.coordinates, select = -c(1:3) )
library(tibble)
slice1.11.coordinates <- tibble::rownames_to_column(slice1.11.coordinates, "barcodes")

slice1.12<-merging@images[["slice1.12"]]
slice1.12.coordinates<-slice1.12@coordinates
slice1.12.coordinates = subset(slice1.12.coordinates, select = -c(1:3) )
library(tibble)
slice1.12.coordinates <- tibble::rownames_to_column(slice1.12.coordinates, "barcodes")

slice1.13<-merging@images[["slice1.13"]]
slice1.13.coordinates<-slice1.13@coordinates
slice1.13.coordinates = subset(slice1.13.coordinates, select = -c(1:3) )
library(tibble)
slice1.13.coordinates <- tibble::rownames_to_column(slice1.13.coordinates, "barcodes")

slice1.14<-merging@images[["slice1.14"]]
slice1.14.coordinates<-slice1.14@coordinates
slice1.14.coordinates = subset(slice1.14.coordinates, select = -c(1:3) )
library(tibble)
slice1.14.coordinates <- tibble::rownames_to_column(slice1.14.coordinates, "barcodes")

slice1.15<-merging@images[["slice1.15"]]
slice1.15.coordinates<-slice1.15@coordinates
slice1.15.coordinates = subset(slice1.15.coordinates, select = -c(1:3) )
library(tibble)
slice1.15.coordinates <- tibble::rownames_to_column(slice1.15.coordinates, "barcodes")

slice1.16<-merging@images[["slice1.16"]]
slice1.16.coordinates<-slice1.16@coordinates
slice1.16.coordinates = subset(slice1.16.coordinates, select = -c(1:3) )
library(tibble)
slice1.16.coordinates <- tibble::rownames_to_column(slice1.16.coordinates, "barcodes")

slice1.17<-merging@images[["slice1.17"]]
slice1.17.coordinates<-slice1.17@coordinates
slice1.17.coordinates = subset(slice1.17.coordinates, select = -c(1:3) )
library(tibble)
slice1.17.coordinates <- tibble::rownames_to_column(slice1.17.coordinates, "barcodes")

slice1.18<-merging@images[["slice1.18"]]
slice1.18.coordinates<-slice1.18@coordinates
slice1.18.coordinates = subset(slice1.18.coordinates, select = -c(1:3) )
library(tibble)
slice1.18.coordinates <- tibble::rownames_to_column(slice1.18.coordinates, "barcodes")

slice1.19<-merging@images[["slice1.19"]]
slice1.19.coordinates<-slice1.19@coordinates
slice1.19.coordinates = subset(slice1.19.coordinates, select = -c(1:3) )
library(tibble)
slice1.19.coordinates <- tibble::rownames_to_column(slice1.19.coordinates, "barcodes")

slice1.20<-merging@images[["slice1.20"]]
slice1.20.coordinates<-slice1.20@coordinates
slice1.20.coordinates = subset(slice1.20.coordinates, select = -c(1:3) )
library(tibble)
slice1.20.coordinates <- tibble::rownames_to_column(slice1.20.coordinates, "barcodes")

slice1.21<-merging@images[["slice1.21"]]
slice1.21.coordinates<-slice1.21@coordinates
slice1.21.coordinates = subset(slice1.21.coordinates, select = -c(1:3) )
library(tibble)
slice1.21.coordinates <- tibble::rownames_to_column(slice1.21.coordinates, "barcodes")

slice1.22<-merging@images[["slice1.22"]]
slice1.22.coordinates<-slice1.22@coordinates
slice1.22.coordinates = subset(slice1.22.coordinates, select = -c(1:3) )
library(tibble)
slice1.22.coordinates <- tibble::rownames_to_column(slice1.22.coordinates, "barcodes")

slice1.23<-merging@images[["slice1.23"]]
slice1.23.coordinates<-slice1.23@coordinates
slice1.23.coordinates = subset(slice1.23.coordinates, select = -c(1:3) )
library(tibble)
slice1.23.coordinates <- tibble::rownames_to_column(slice1.23.coordinates, "barcodes")

slice1.24<-merging@images[["slice1.24"]]
slice1.24.coordinates<-slice1.24@coordinates
slice1.24.coordinates = subset(slice1.24.coordinates, select = -c(1:3) )
library(tibble)
slice1.24.coordinates <- tibble::rownames_to_column(slice1.24.coordinates, "barcodes")

slice1.25<-merging@images[["slice1.25"]]
slice1.25.coordinates<-slice1.25@coordinates
slice1.25.coordinates = subset(slice1.25.coordinates, select = -c(1:3) )
library(tibble)
slice1.25.coordinates <- tibble::rownames_to_column(slice1.25.coordinates, "barcodes")

slice1.26<-merging@images[["slice1.26"]]
slice1.26.coordinates<-slice1.26@coordinates
slice1.26.coordinates = subset(slice1.26.coordinates, select = -c(1:3) )
library(tibble)
slice1.26.coordinates <- tibble::rownames_to_column(slice1.26.coordinates, "barcodes")

slice1.27<-merging@images[["slice1.27"]]
slice1.27.coordinates<-slice1.27@coordinates
slice1.27.coordinates = subset(slice1.27.coordinates, select = -c(1:3) )
library(tibble)
slice1.27.coordinates <- tibble::rownames_to_column(slice1.27.coordinates, "barcodes")

coords<-rbind(slice1.coordinates,slice1.1.coordinates,slice1.2.coordinates,slice1.3.coordinates,slice1.4.coordinates,slice1.5.coordinates,slice1.6.coordinates,slice1.7.coordinates,slice1.8.coordinates,slice1.9.coordinates,slice1.10.coordinates,slice1.11.coordinates,slice1.12.coordinates,slice1.13.coordinates,slice1.14.coordinates,slice1.15.coordinates,slice1.16.coordinates,slice1.17.coordinates,slice1.18.coordinates,slice1.19.coordinates,slice1.20.coordinates,slice1.21.coordinates,slice1.22.coordinates,slice1.23.coordinates,slice1.24.coordinates,slice1.25.coordinates,slice1.26.coordinates,slice1.27.coordinates)
colnames(coords) <- c("barcodes","xcoord", "ycoord")
rownames(coords) <- coords$barcodes; coords$barcodes <- NULL
puck <- SpatialRNA(coords, spatialcounts)
print(dim(puck@counts))
hist(log(puck@nUMI,2))


print(head(puck@coords)) 
barcodes <- colnames(puck@counts) 
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 

                     
myRCTD <- create.RCTD(puck, reference, max_cores = 40)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
results <- myRCTD@results
norm_weights = normalize_weights(results$weights) 
cell_type_names <- myRCTD@cell_type_info$info[[2]]
spatialRNA <- myRCTD@spatialRNA
resultsdir <- '/rds/projects/c/croftap-stia-atac/Visium/RCTD_analysis/'
plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights)
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
                     results$results_df) 
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir) 

doublets <- results$results_df[results$results_df$spot_class == "doublet_certain",]
plot_doublets(spatialRNA, doublets, resultsdir, cell_type_names) 


plot_doublets_type(spatialRNA, doublets, resultsdir, cell_type_names)
doub_occur <- table(doublets$second_type, doublets$first_type) 
plot_doub_occur_stack(doub_occur, resultsdir, cell_type_names) 

                     
RCTD.annotations<-myRCTD@results[["results_df"]]
RCTD.annotations_meta <- subset(RCTD.annotations, select = c(first_type))
table(row.names(RCTD.annotations_meta) %in% colnames(merging@assays$Spatial))
table(colnames(merging@assays$Spatial) %in% row.names(RCTD.annotations_meta))
cellstoremove <- colnames(merging@assays$Spatial)[!(colnames(merging@assays$Spatial) %in% row.names(RCTD.annotations_meta))]
merging_subset <- merging[,!(colnames(merging@assays$Spatial) %in% cellstoremove)]
length(colnames(merging_subset))

                     
merging_subset<-AddMetaData(merging_subset, RCTD.annotations_meta)

merging_subset$cm_clusters<-merging_subset@active.ident

RCTDpredictions <- table(merging_subset$cm_clusters, merging_subset$first_type)
RCTDpredictions <- RCTDpredictions/rowSums(RCTDpredictions)  # normalize for number of cells in each cell type
RCTDpredictions <- as.data.frame(RCTDpredictions)
p3 <- ggplot(RCTDpredictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("Seurat annoatated") + ylab("Predicted (RCTD)") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p3


merging_subset$first_type_nolining<-merging_subset$first_type
Idents(merging_subset)<-'first_type_nolining'
levels(merging_subset)

current.sample.ids <-c("B cell"     ,       "Endothelial cells" ,"Fibroblast"  ,      "linning"     ,      "Monocyte"   ,       "T cell" )
new.sample.ids <-c("B cell"     ,       "Endothelial cells" ,"Fibroblast"  ,      "Fibroblast"     ,      "Monocyte"   ,       "T cell" )
merging_subset@meta.data[["first_type_nolining"]] <- plyr::mapvalues(x = merging_subset@meta.data[["first_type_nolining"]], from = current.sample.ids, to = new.sample.ids)

RCTDpredictions <- table(merging_subset$cm_clusters, merging_subset$first_type_nolining)
RCTDpredictions <- RCTDpredictions/rowSums(RCTDpredictions)  # normalize for number of cells in each cell type
RCTDpredictions <- as.data.frame(RCTDpredictions)
p3 <- ggplot(RCTDpredictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("Seurat annoatated") + ylab("Predicted (RCTD)") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p3

Idents(merging_subset)<-'first_type'
levels(merging_subset)
merging_subset$first_type_firbomerge<-merging_subset$first_type
current.sample.ids <-c("B cell"      ,      "Endothelial cells", "Fibroblast"  ,      "linning"   ,        "Monocyte"     ,     "T cell")
new.sample.ids <-c("B cell"      ,      "Endothelial cells", "Fibroblast"  ,      "Fibroblast"   ,        "Monocyte"     ,     "T cell")
merging_subset@meta.data[["first_type_firbomerge"]] <- plyr::mapvalues(x = merging_subset@meta.data[["first_type_firbomerge"]], from = current.sample.ids, to = new.sample.ids)
table(merging_subset$first_type_firbomerge)
saveRDS(merging_subset, "/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/RCTD_analysis/merging_subset.rds")

DimPlot(merging_subset, group.by = "cm_clusters")
DimPlot(merging_subset, group.by = "first_type", label=T, label.box = T)
