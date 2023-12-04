#Patritcia can you add the script you ran ofr building the initial visiuim object?
#Reading in data (DONE), creating seuart objects (DONE), merging data preprocessing (DONE), clustering etc.
#Thanks!

#Load data####

S1datA <- Load10X_Spatial("./SampleA/outs")
S1datA@meta.data$orig.ident <- gsub("SeuratProject", "Slide1SampleA", S1datA@meta.data$orig.ident)
S1datB <- Load10X_Spatial("./SampleB/outs")
S1datB@meta.data$orig.ident <- gsub("SeuratProject", "Slide1SampleB", S1datB@meta.data$orig.ident)
S1datC <- Load10X_Spatial("./SampleC/outs")
S1datC@meta.data$orig.ident <- gsub("SeuratProject", "Slide1SampleC", S1datC@meta.data$orig.ident)
S1datD <- Load10X_Spatial("./SampleD/outs")
S1datD@meta.data$orig.ident <- gsub("SeuratProject", "Slide1SampleD", S1datD@meta.data$orig.ident)
S2datA <- Load10X_Spatial("./Slide2SampleA/outs")
S2datA@meta.data$orig.ident <- gsub("SeuratProject", "Slide2SampleA", S2datA@meta.data$orig.ident)
S2datB <- Load10X_Spatial("./Slide2SampleB/outs")
S2datB@meta.data$orig.ident <- gsub("SeuratProject", "Slide2SampleB", S2datB@meta.data$orig.ident)
S2datC <- Load10X_Spatial("./Slide2SampleC/outs")
S2datC@meta.data$orig.ident <- gsub("SeuratProject", "Slide2SampleC", S2datC@meta.data$orig.ident)
S2datD <- Load10X_Spatial("./Slide2SampleD/outs")
S2datD@meta.data$orig.ident <- gsub("SeuratProject", "Slide2SampleD", S2datD@meta.data$orig.ident)
S3datA <- Load10X_Spatial("./Slide3SampleA/outs")
S3datA@meta.data$orig.ident <- gsub("SeuratProject", "Slide3SampleA", S3datA@meta.data$orig.ident)
S3datB <- Load10X_Spatial("./Slide3SampleB/outs")
S3datB@meta.data$orig.ident <- gsub("SeuratProject", "Slide3SampleB", S3datB@meta.data$orig.ident)
S3datC <- Load10X_Spatial("./Slide3SampleC/outs")
S3datC@meta.data$orig.ident <- gsub("SeuratProject", "Slide3SampleC", S3datC@meta.data$orig.ident)
S3datD <- Load10X_Spatial("./Slide3SampleD/outs")
S3datD@meta.data$orig.ident <- gsub("SeuratProject", "Slide3SampleD", S3datD@meta.data$orig.ident)
S4datA <- Load10X_Spatial("./Slide4SampleA/outs")
S4datA@meta.data$orig.ident <- gsub("SeuratProject", "Slide4SampleA", S4datA@meta.data$orig.ident)
S4datB <- Load10X_Spatial("./Slide4SampleB/outs")
S4datB@meta.data$orig.ident <- gsub("SeuratProject", "Slide4SampleB", S4datB@meta.data$orig.ident)
S4datC <- Load10X_Spatial("./Slide4SampleC/outs")
S4datC@meta.data$orig.ident <- gsub("SeuratProject", "Slide4SampleC", S4datC@meta.data$orig.ident)
S4datD <- Load10X_Spatial("./Slide4SampleD/outs")
S4datD@meta.data$orig.ident <- gsub("SeuratProject", "Slide4SampleD", S4datD@meta.data$orig.ident)
S5datA <- Load10X_Spatial("./Slide5SampleA/outs")
S5datA@meta.data$orig.ident <- gsub("SeuratProject", "Slide5SampleA", S5datA@meta.data$orig.ident)
S5datB <- Load10X_Spatial("./Slide5SampleB/outs")
S5datB@meta.data$orig.ident <- gsub("SeuratProject", "Slide5SampleB", S5datB@meta.data$orig.ident)
S5datC <- Load10X_Spatial("./Slide5SampleC/outs")
S5datC@meta.data$orig.ident <- gsub("SeuratProject", "Slide5SampleC", S5datC@meta.data$orig.ident)
S5datD <- Load10X_Spatial("./Slide5SampleD/outs")
S5datD@meta.data$orig.ident <- gsub("SeuratProject", "Slide5SampleD", S5datD@meta.data$orig.ident)
S6datA <- Load10X_Spatial("./Slide6SampleA/outs")
S6datA@meta.data$orig.ident <- gsub("SeuratProject", "Slide6SampleA", S6datA@meta.data$orig.ident)
S6datB <- Load10X_Spatial("./Slide6SampleB/outs")
S6datB@meta.data$orig.ident <- gsub("SeuratProject", "Slide6SampleB", S6datB@meta.data$orig.ident)
S6datC <- Load10X_Spatial("./Slide6SampleC/outs")
S6datC@meta.data$orig.ident <- gsub("SeuratProject", "Slide6SampleC", S6datC@meta.data$orig.ident)
S6datD <- Load10X_Spatial("./Slide6SampleD/outs")
S6datD@meta.data$orig.ident <- gsub("SeuratProject", "Slide6SampleD", S6datD@meta.data$orig.ident)
S7datA <- Load10X_Spatial("./Slide7SampleA/outs")
S7datA@meta.data$orig.ident <- gsub("SeuratProject", "Slide7SampleA", S7datA@meta.data$orig.ident)
S7datB <- Load10X_Spatial("./Slide7SampleB/outs")
S7datB@meta.data$orig.ident <- gsub("SeuratProject", "Slide7SampleB", S7datB@meta.data$orig.ident)
S7datC <- Load10X_Spatial("./Slide7SampleC/outs")
S7datC@meta.data$orig.ident <- gsub("SeuratProject", "Slide7SampleC", S7datC@meta.data$orig.ident)
S7datD <- Load10X_Spatial("./Slide7SampleD/outs")
S7datD@meta.data$orig.ident <- gsub("SeuratProject", "Slide7SampleD", S7datD@meta.data$orig.ident)
S8datA <- Load10X_Spatial("./Slide8SampleA/outs")
S8datA@meta.data$orig.ident <- gsub("SeuratProject", "Slide8SampleA", S8datA@meta.data$orig.ident)
S8datB <- Load10X_Spatial("./Slide8SampleB/outs")
S8datB@meta.data$orig.ident <- gsub("SeuratProject", "Slide8SampleB", S8datB@meta.data$orig.ident)
S8datC <- Load10X_Spatial("./Slide8SampleC/outs")
S8datC@meta.data$orig.ident <- gsub("SeuratProject", "Slide8SampleC", S8datC@meta.data$orig.ident)

#QC and normalisation####

#_S1A####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S1datA, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S1datA, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
S1datA[["percent.mt"]] <- PercentageFeatureSet(object = S1datA, pattern = "^MT-")
S1datA <- SCTransform(S1datA, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(S1datA, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"))
S1datA.2 <- S1datA
S1datA.2$group <- 1

# Visualize QC metrics as a violin plot
VlnPlot(S1datA.2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3, group.by = "group")

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(S1datA, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(S1datA, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2

#_S1B####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S1datB, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S1datB, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
S1datB[["percent.mt"]] <- PercentageFeatureSet(S1datB, pattern = "^MT-")
S1datB <- SCTransform(S1datB, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S1datB, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(S1datB, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(S1datB, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2

#_S1C####
plot1 <- VlnPlot(S1datC, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S1datC, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S1datC[["percent.mt"]] <- PercentageFeatureSet(S1datC, pattern = "^MT-")
S1datC <- SCTransform(S1datC, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S1datC, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(S1datC, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(S1datC, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2



#_S1D####
plot1 <- VlnPlot(S1datD, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S1datD, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S1datD[["percent.mt"]] <- PercentageFeatureSet(S1datD, pattern = "^MT-")
S1datD <- SCTransform(S1datD, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S1datD, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(S1datD, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(S1datD, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2





#_S2A####
plot1 <- VlnPlot(S2datA, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S2datA, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S2datA[["percent.mt"]] <- PercentageFeatureSet(S2datA, pattern = "^MT-")
S2datA <- SCTransform(S2datA, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S2datA, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(S2datA, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(S2datA, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2



#_S2B####
plot1 <- VlnPlot(S2datB, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S2datB, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S2datB[["percent.mt"]] <- PercentageFeatureSet(S2datB, pattern = "^MT-")
S2datB <- SCTransform(S2datB, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S2datB, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(S2datB, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(S2datB, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2



#_S2C####
plot1 <- VlnPlot(S2datC, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S2datC, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S2datC[["percent.mt"]] <- PercentageFeatureSet(S2datC, pattern = "^MT-")
S2datC <- SCTransform(S2datC, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S2datC, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(S2datC, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(S2datC, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2



#_S2D####
plot1 <- VlnPlot(S2datD, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S2datD, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S2datD[["percent.mt"]] <- PercentageFeatureSet(S2datD, pattern = "^MT-")
S2datD <- SCTransform(S2datD, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S2datD, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(S2datD, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(S2datD, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2



#_S3A####
plot1 <- VlnPlot(S3datA, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S3datA, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S3datA[["percent.mt"]] <- PercentageFeatureSet(S3datA, pattern = "^MT-")
S3datA <- SCTransform(S3datA, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S3datA, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(S3datA, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(S3datA, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2



#_S3B####
plot1 <- VlnPlot(S3datB, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S3datB, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S3datB[["percent.mt"]] <- PercentageFeatureSet(S3datB, pattern = "^MT-")
S3datB <- SCTransform(S3datB, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S3datB, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(S3datB, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(S3datB, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2



#_S3C####
plot1 <- VlnPlot(S3datC, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S3datC, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S3datC[["percent.mt"]] <- PercentageFeatureSet(S3datC, pattern = "^MT-")
S3datC <- SCTransform(S3datC, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S3datC, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(S3datC, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(S3datC, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2



#_S3D####
plot1 <- VlnPlot(S3datD, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S3datD, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S3datD[["percent.mt"]] <- PercentageFeatureSet(S3datD, pattern = "^MT-")
S3datD <- SCTransform(S3datD, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S3datD, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(S3datD, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(S3datD, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2


#_S4A####

#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S4datA, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S4datA, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S4datA[["percent.mt"]] <- PercentageFeatureSet(object = S4datA, pattern = "^MT-")
#SCTransform will select variable genes and normalize in one step.
S4datA <- SCTransform(S4datA, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S4datA, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(S4datA, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(S4datA, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2


#_S4B####

#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S4datB, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S4datB, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S4datB[["percent.mt"]] <- PercentageFeatureSet(object = S4datB, pattern = "^MT-")
#SCTransform will select variable genes and normalize in one step.

S4datB <- SCTransform(S4datB, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S4datB, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(S4datB, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(S4datB, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2


#_S4C####

#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S4datC, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S4datC, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S4datC[["percent.mt"]] <- PercentageFeatureSet(object = S4datC, pattern = "^MT-")
#SCTransform will select variable genes and normalize in one step.
S4datC <- SCTransform(S4datC, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S4datC, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(S4datC, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(S4datC, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2

#_S4D####

#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S4datD, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S4datD, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S4datD[["percent.mt"]] <- PercentageFeatureSet(object = S4datD, pattern = "^MT-")
#SCTransform will select variable genes and normalize in one step.
S4datD <- SCTransform(S4datD, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S4datD, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(S4datD, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(S4datD, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2






#_S5A####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S5datA, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S5datA, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S5datA[["percent.mt"]] <- PercentageFeatureSet(S5datA, pattern = "^MT-")
S5datA <- SCTransform(S5datA, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S5datA, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)


#_S5B####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S5datB, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S5datB, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S5datB[["percent.mt"]] <- PercentageFeatureSet(S5datB, pattern = "^MT-")
S5datB <- SCTransform(S5datB, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S5datB, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)


#_S5C####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S5datC, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S5datC, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S5datC[["percent.mt"]] <- PercentageFeatureSet(S5datC, pattern = "^MT-")
S5datC <- SCTransform(S5datC, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S5datC, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)


#_S5D####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S5datD, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S5datD, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S5datD[["percent.mt"]] <- PercentageFeatureSet(S5datD, pattern = "^MT-")
S5datD <- SCTransform(S5datD, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S5datD, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)


#_S6A ####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S6datA, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S6datA, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S6datA[["percent.mt"]] <- PercentageFeatureSet(S6datA, pattern = "^MT-")
S6datA <- SCTransform(S6datA, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S6datA, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)


#_S6B ####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S6datB, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S6datB, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S6datB[["percent.mt"]] <- PercentageFeatureSet(S6datB, pattern = "^MT-")
S6datB <- SCTransform(S6datB, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S6datB, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)


#_S6C ####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S6datC, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S6datC, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S6datC[["percent.mt"]] <- PercentageFeatureSet(S6datC, pattern = "^MT-")
S6datC <- SCTransform(S6datC, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S6datC, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)


#_S6D ####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S6datD, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S6datD, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S6datD[["percent.mt"]] <- PercentageFeatureSet(S6datD, pattern = "^MT-")
S6datD <- SCTransform(S6datD, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S6datD, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)


#_S7A ####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S7datA, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S7datA, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S7datA[["percent.mt"]] <- PercentageFeatureSet(S7datA, pattern = "^MT-")
S7datA <- SCTransform(S7datA, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S7datA, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)


#_S7B ####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S7datB, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S7datB, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S7datB[["percent.mt"]] <- PercentageFeatureSet(S7datB, pattern = "^MT-")
S7datB <- SCTransform(S7datB, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S7datB, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)

#_S7C ####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S7datC, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S7datC, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S7datC[["percent.mt"]] <- PercentageFeatureSet(S7datC, pattern = "^MT-")
S7datC <- SCTransform(S7datC, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S7datC, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)

#_S7D ####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S7datD, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S7datD, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S7datD[["percent.mt"]] <- PercentageFeatureSet(S7datD, pattern = "^MT-")
S7datD <- SCTransform(S7datD, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S7datD, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)

#_S8A ####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S8datA, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S8datA, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S8datA[["percent.mt"]] <- PercentageFeatureSet(S8datA, pattern = "^MT-")
S8datA <- SCTransform(S8datA, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S8datA, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)

#_S8B ####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S8datB, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S8datB, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S8datB[["percent.mt"]] <- PercentageFeatureSet(S8datB, pattern = "^MT-")
S8datB <- SCTransform(S8datB, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S8datB, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)

#_S8C ####
#.........................................QUALITY CONTROL...............................................
plot1 <- VlnPlot(S8datC, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S8datC, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

S8datC[["percent.mt"]] <- PercentageFeatureSet(S8datC, pattern = "^MT-")
S8datC <- SCTransform(S8datC, assay = "Spatial", verbose = FALSE)

# Visualize QC metrics as a violin plot
VlnPlot(S8datC, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)



#DOWNSTREAM ANALYSIS####
#Dimensionality reduction and clustering# to run on the SCT assay!!!

#_S1A####

#NO HARMONY
S1datA <- RunPCA(S1datA, assay = "SCT", verbose = FALSE)
S1datA <- FindNeighbors(S1datA, reduction = "pca", dims = 1:30)
S1datA <- FindClusters(S1datA, verbose = FALSE, resolution = c(0.2, 
                                                                   0.3,
                                                                   0.4,
                                                                   0.5,
                                                                   0.6,
                                                                   0.7,
                                                                   0.8,
                                                                   0.9,
                                                                   1.0))

clustree(S1datA)
S1datA <- FindClusters(S1datA, verbose = FALSE, resolution = 0.4) #SET
S1datA <- RunUMAP(S1datA, reduction = "pca", dims = 1:30)
DimPlot(S1datA, pt.size = 1.5)

#HARMONY

S1datAH <- RunPCA(S1datA, assay = "SCT", verbose = FALSE)
S1datAH <- FindNeighbors(S1datAH, reduction = "pca", dims = 1:30)
S1datAH <- FindClusters(S1datAH, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))

clustree(S1datAH)
S1datAH <- FindClusters(S1datAH, verbose = FALSE, resolution = 0.4) #SET
S1datAH <- RunHarmony(S1datAH, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S1datAH <- RunUMAP(S1datAH, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S1datAH)
DimPlot(S1datAH, reduction="HarmonyUMAP", pt.size = 1.5)


#find all markers and plot them in heatmap
S1datAdegs <- FindAllMarkers(S1datAH, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S1datAdegs$cluster)) {
  genes <- c(genes,
             head(S1datAdegs$gene[grep(ind, 
                                       S1datAdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S1datAH, genes)
DotPlot(S1datAH, features = genes) + RotatedAxis()

#_S1B####

#NO HARMONY
S1datB <- RunPCA(S1datB, assay = "SCT", verbose = FALSE)
S1datB <- FindNeighbors(S1datB, reduction = "pca", dims = 1:30)
S1datB <- FindClusters(S1datB, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S1datB)
S1datB <- FindClusters(S1datB, verbose = FALSE, resolution = 0.5)
S1datB <- RunUMAP(S1datB, reduction = "pca", dims = 1:30)
DimPlot(S1datB)

#HARMONY
S1datBH <- RunPCA(S1datB, assay = "SCT", verbose = FALSE)
S1datBH <- FindNeighbors(S1datBH, reduction = "pca", dims = 1:30)
S1datBH <- FindClusters(S1datBH, verbose = FALSE, resolution = c(0.2, 
                                                                 0.3,
                                                                 0.4,
                                                                 0.5,
                                                                 0.6,
                                                                 0.7,
                                                                 0.8,
                                                                 0.9,
                                                                 1.0))

clustree(S1datBH)
S1datBH <- FindClusters(S1datBH, verbose = FALSE, resolution = 0.5) #SET
S1datBH <- RunHarmony(S1datBH, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S1datBH <- RunUMAP(S1datBH, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S1datBH)
DimPlot(S1datBH, reduction="HarmonyUMAP", pt.size = 1.5)



S1datBdegs <- FindAllMarkers(S1datBH, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S1datBdegs$cluster)) {
  genes <- c(genes,
             head(S1datBdegs$gene[grep(ind, 
                                       S1datBdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S1datBH, genes)

#_S1C####

#NO HARMONY
S1datC <- RunPCA(S1datC, assay = "SCT", verbose = FALSE)
S1datC <- FindNeighbors(S1datC, reduction = "pca", dims = 1:30)
S1datC <- FindClusters(S1datC, verbose = FALSE, resolution = c(0.2, 
                                                           0.3,
                                                           0.4,
                                                           0.5,
                                                           0.6,
                                                           0.7,
                                                           0.8,
                                                           0.9,
                                                           1.0))
clustree(S1datC)
S1datC <- FindClusters(S1datC, verbose = FALSE, resolution = 0.6) #Set
S1datC <- RunUMAP(S1datC, reduction = "pca", dims = 1:30)
DimPlot(S1datC, pt.size = 1.5)

S1datCdegs <- FindAllMarkers(S1datC, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S1datCdegs$cluster)) {
  genes <- c(genes,
             head(S1datCdegs$gene[grep(ind, 
                                       S1datCdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S1datC, genes)
DotPlot(S1datC, features = c())

#HARMONY
S1datCH <- RunPCA(S1datC, assay = "SCT", verbose = FALSE)
S1datCH <- FindNeighbors(S1datCH, reduction = "pca", dims = 1:30)
S1datCH <- FindClusters(S1datCH, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S1datCH)
S1datCH <- FindClusters(S1datCH, verbose = FALSE, resolution = 0.6) #Set
S1datCH <- RunUMAP(S1datCH, reduction = "pca", dims = 1:30)
DimPlot(S1datCH, pt.size = 1.5)

S1datCdegs <- FindAllMarkers(S1datCH, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S1datCdegs$cluster)) {
  genes <- c(genes,
             head(S1datCdegs$gene[grep(ind, 
                                       S1datCdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S1datCH, genes)


#_S1D####
S1datD <- RunPCA(S1datD, assay = "SCT", verbose = FALSE)
S1datD <- FindNeighbors(S1datD, reduction = "pca", dims = 1:30)
S1datD <- FindClusters(S1datD, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S1datD)
S1datD <- FindClusters(S1datD, verbose = FALSE, resolution = 0.8)
S1datD <- RunUMAP(S1datD, reduction = "pca", dims = 1:30)
DimPlot(S1datD)

S1datCdegs <- FindAllMarkers(S1datD, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S1datDdegs$cluster)) {
  genes <- c(genes,
             head(S1datDdegs$gene[grep(ind, 
                                       S1datDdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S1datD, genes)
DotPlot(S1datD, features = genes) + RotatedAxis()

    #NOT RELEVANT - TOO SMALL


#_S2A####

#NO HARMONY
S2datA <- RunPCA(S2datA, assay = "SCT", verbose = FALSE)
S2datA <- FindNeighbors(S2datA, reduction = "pca", dims = 1:30)
S2datA <- FindClusters(S2datA, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S2datA)
S2datA <- FindClusters(S2datA, verbose = FALSE, resolution = 0.6)
S2datA <- RunUMAP(S2datA, reduction = "pca", dims = 1:30)
DimPlot(S2datA, pt.size = 1.5)

S2datAdegs <- FindAllMarkers(S2datA, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S2datAdegs$cluster)) {
  genes <- c(genes,
             head(S2datAdegs$gene[grep(ind, 
                                       S2datAdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S2datA, genes)
DotPlot(S2datA, features = genes) + RotatedAxis()


#HARMONY
S2datAH <- RunPCA(S2datA, assay = "SCT", verbose = FALSE)
S2datAH <- FindNeighbors(S2datAH, reduction = "pca", dims = 1:30)
S2datAH <- FindClusters(S2datAH, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S2datAH)
S2datAH <- FindClusters(S2datAH, verbose = FALSE, resolution = 0.6) #set
S2datAH <- RunUMAP(S2datAH, reduction = "pca", dims = 1:30)
DimPlot(S2datAH, pt.size = 1.5)

S2datAdegs <- FindAllMarkers(S2datAH, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S2datAdegs$cluster)) {
  genes <- c(genes,
             head(S2datAdegs$gene[grep(ind, 
                                       S2datAdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S2datAH, genes)
DotPlot(S2datAH, features = genes) + RotatedAxis()

#_S2B####

#NO HARMONY
S2datB <- RunPCA(S2datB, assay = "SCT", verbose = FALSE)
S2datB <- FindNeighbors(S2datB, reduction = "pca", dims = 1:30)
S2datB <- FindClusters(S2datB, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S2datB)
S2datB <- FindClusters(S2datB, verbose = FALSE, resolution = 0.7)
S2datB <- RunUMAP(S2datB, reduction = "pca", dims = 1:30)
DimPlot(S2datB)
   
#HARMONY
S2datBH <- RunPCA(S2datB, assay = "SCT", verbose = FALSE)
S2datBH <- FindNeighbors(S2datBH, reduction = "pca", dims = 1:30)
S2datBH <- FindClusters(S2datBH, verbose = FALSE, resolution = c(0.2, 
                                                                 0.3,
                                                                 0.4,
                                                                 0.5,
                                                                 0.6,
                                                                 0.7,
                                                                 0.8,
                                                                 0.9,
                                                                 1.0))

clustree(S2datBH)
S2datBH <- FindClusters(S2datBH, verbose = FALSE, resolution = 0.5) #SET
S2datBH <- RunHarmony(S2datBH, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S2datBH <- RunUMAP(S2datBH, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S2datBH)
DimPlot(S2datBH, reduction="HarmonyUMAP", pt.size = 1.5)

#Find all markers

S2datBdegs <- FindAllMarkers(S2datBH, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S2datBdegs$cluster)) {
  genes <- c(genes,
             head(S2datBdegs$gene[grep(ind, 
                                       S2datBdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S2datBH, genes)

#Find all markers

S2datBdegs <- FindAllMarkers(S2datBH, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S2datBdegs$cluster)) {
  genes <- c(genes,
             head(S2datBdegs$gene[grep(ind, 
                                       S2datBdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S2datBH, genes)


#_S2C####

#NO HARMONY
S2datC <- RunPCA(S2datC, assay = "SCT", verbose = FALSE)
S2datC <- FindNeighbors(S2datC, reduction = "pca", dims = 1:30)
S2datC <- FindClusters(S2datC, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S2datC)
S2datC <- FindClusters(S2datC, verbose = FALSE, resolution = 0.5)
S2datC <- RunUMAP(S2datC, reduction = "pca", dims = 1:30)
DimPlot(S2datC)

#HARMONY
S2datCH <- RunPCA(S2datC, assay = "SCT", verbose = FALSE)
S2datCH  <- FindNeighbors(S2datCH, reduction = "pca", dims = 1:30)
S2datCH <- FindClusters(S2datCH, verbose = FALSE, resolution = c(0.2, 
                                                                 0.3,
                                                                 0.4,
                                                                 0.5,
                                                                 0.6,
                                                                 0.7,
                                                                 0.8,
                                                                 0.9,
                                                                 1.0))

clustree(S2datCH)
S2datCH <- FindClusters(S2datCH, verbose = FALSE, resolution = 0.5) #SET
S2datCH <- RunHarmony(S2datCH, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S2datCH <- RunUMAP(S2datCH, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S2datCH)
DimPlot(S2datCH, reduction="HarmonyUMAP", pt.size = 1.5)


#.....

S2datCdegs <- FindAllMarkers(S2datCH, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S2datCdegs$cluster)) {
  genes <- c(genes,
             head(S2datCdegs$gene[grep(ind, 
                                       S2datCdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S2datCH, genes)


#_S2D####

#NO HARMONY
S2datD <- RunPCA(S2datD, assay = "SCT", verbose = FALSE)
S2datD <- FindNeighbors(S2datD, reduction = "pca", dims = 1:30)
S2datD <- FindClusters(S2datD, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S2datD)
S2datD <- FindClusters(S2datD, verbose = FALSE, resolution = 0.6)
S2datD <- RunUMAP(S2datD, reduction = "pca", dims = 1:30)
DimPlot(S2datD, pt.size = 1.5)

S2datDdegs <- FindAllMarkers(S2datD, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S2datDdegs$cluster)) {
  genes <- c(genes,
             head(S2datDdegs$gene[grep(ind, 
                                       S2datDdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S2datD, genes)
DotPlot(S2datD, features = genes) + RotatedAxis()

#HARMONY

S2datDH <- RunPCA(S2datD, assay = "SCT", verbose = FALSE)
S2datDH <- FindNeighbors(S2datDH, reduction = "pca", dims = 1:30)
S2datDH <- FindClusters(S2datDH, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S2datDH)
S2datDH <- FindClusters(S2datDH, verbose = FALSE, resolution = 0.6) #SET
S2datDH <- RunUMAP(S2datDH, reduction = "pca", dims = 1:30)
DimPlot(S2datDH, pt.size = 1.5)

S2datDdegs <- FindAllMarkers(S2datDH, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S2datDdegs$cluster)) {
  genes <- c(genes,
             head(S2datDdegs$gene[grep(ind, 
                                       S2datDdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S2datDH, genes)
DotPlot(S2datDH, features = genes) + RotatedAxis()

#_S3A####

#NO HARMONY
S3datA <- RunPCA(S3datA, assay = "SCT", verbose = FALSE)
S3datA <- FindNeighbors(S3datA, reduction = "pca", dims = 1:30)
S3datA <- FindClusters(S3datA, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S3datA)
S3datA <- FindClusters(S3datA, verbose = FALSE, resolution = 0.5)
S3datA <- RunUMAP(S3datA, reduction = "pca", dims = 1:30)
DimPlot(S3datA, pt.size = 1.5)


S3datAdegs <- FindAllMarkers(S3datA, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S3datAdegs$cluster)) {
  genes <- c(genes,
             head(S3datAdegs$gene[grep(ind, 
                                       S3datAdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S3datA, genes)
DotPlot(S3datA, features = genes) + RotatedAxis()


#HARMONY

S3datAH <- RunPCA(S3datA, assay = "SCT", verbose = FALSE)
S3datAH <- FindNeighbors(S3datAH, reduction = "pca", dims = 1:30)
S3datAH <- FindClusters(S3datAH, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S3datAH)
S3datAH <- FindClusters(S3datAH, verbose = FALSE, resolution = 0.5) #SET
S3datAH <- RunUMAP(S3datAH, reduction = "pca", dims = 1:30)
DimPlot(S3datAH, pt.size = 1.5)


S3datAdegs <- FindAllMarkers(S3datAH, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S3datAdegs$cluster)) {
  genes <- c(genes,
             head(S3datAdegs$gene[grep(ind, 
                                       S3datAdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S3datAH, genes)
DotPlot(S3datAH, features = genes) + RotatedAxis()

#_S3B####

#NO HARMONY
S3datB <- RunPCA(S3datB, assay = "SCT", verbose = FALSE)
S3datB <- FindNeighbors(S3datB, reduction = "pca", dims = 1:30)
S3datB <- FindClusters(S3datB, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S3datB)
S3datB <- FindClusters(S3datB, verbose = FALSE, resolution = 0.6)
S3datB <- RunUMAP(S3datB, reduction = "pca", dims = 1:30)
DimPlot(S3datB)

#HARMONY
S3datBH <- RunPCA(S3datB, assay = "SCT", verbose = FALSE)
S3datBH <- FindNeighbors(S3datBH, reduction = "pca", dims = 1:30)
S3datBH <- FindClusters(S3datBH, verbose = FALSE, resolution = c(0.2, 
                                                                 0.3,
                                                                 0.4,
                                                                 0.5,
                                                                 0.6,
                                                                 0.7,
                                                                 0.8,
                                                                 0.9,
                                                                 1.0))

clustree(S3datBH)
S3datBH <- FindClusters(S3datBH, verbose = FALSE, resolution = 0.5) #SET
S3datBH <- RunHarmony(S3datBH, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S3datBH <- RunUMAP(S3datBH, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S3datBH)
DimPlot(S3datBH, reduction="HarmonyUMAP", pt.size = 1.5)

S3datBdegs <- FindAllMarkers(S3datB, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S3datBdegs$cluster)) {
  genes <- c(genes,
             head(S3datBdegs$gene[grep(ind, 
                                       S3datBdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S3datB, genes)


#_S3C####

#NO HARMONY
S3datC <- RunPCA(S3datC, assay = "SCT", verbose = FALSE)
S3datC <- FindNeighbors(S3datC, reduction = "pca", dims = 1:30)
S3datC <- FindClusters(S3datC, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S3datC)
S3datC <- FindClusters(S3datC, verbose = FALSE, resolution = 0.9)
S3datC <- RunUMAP(S3datC, reduction = "pca", dims = 1:30)
DimPlot(S3datC)


#HARMONY
S3datCH <- RunPCA(S3datC, assay = "SCT", verbose = FALSE)
S3datCH <- FindNeighbors(S3datCH, reduction = "pca", dims = 1:30)
S3datCH <- FindClusters(S3datCH, verbose = FALSE, resolution = c(0.2, 
                                                                 0.3,
                                                                 0.4,
                                                                 0.5,
                                                                 0.6,
                                                                 0.7,
                                                                 0.8,
                                                                 0.9,
                                                                 1.0))

clustree(S3datCH)
S3datCH <- FindClusters(S3datCH, verbose = FALSE, resolution = 0.9) #SET
S3datCH <- RunHarmony(S3datCH, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S3datCH <- RunUMAP(S3datCH, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S3datCH)
DimPlot(S3datCH, reduction="HarmonyUMAP", pt.size = 1.5)





S3datCdegs <- FindAllMarkers(S3datCH, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S3datCdegs$cluster)) {
  genes <- c(genes,
             head(S3datCdegs$gene[grep(ind, 
                                       S3datCdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S3datCH, genes)


#_S3D####

#NO HARMONY
S3datD <- RunPCA(S3datD, assay = "SCT", verbose = FALSE)
S3datD <- FindNeighbors(S3datD, reduction = "pca", dims = 1:30)
S3datD <- FindClusters(S3datD, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S3datD)
S3datD <- FindClusters(S3datD, verbose = FALSE, resolution = 0.6)
S3datD <- RunUMAP(S3datD, reduction = "pca", dims = 1:30)
DimPlot(S3datD, pt.size = 2)
FeaturePlot(S3datD, c(SC_MallMar))


#HARMONY
S3datDH <- RunPCA(S3datD, assay = "SCT", verbose = FALSE)
S3datDH <- FindNeighbors(S3datDH, reduction = "pca", dims = 1:30)
S3datDH <- FindClusters(S3datDH, verbose = FALSE, resolution = c(0.2, 
                                                                 0.3,
                                                                 0.4,
                                                                 0.5,
                                                                 0.6,
                                                                 0.7,
                                                                 0.8,
                                                                 0.9,
                                                                 1.0))

clustree(S3datDH)
S3datDH <- FindClusters(S3datDH, verbose = FALSE, resolution = 0.6) #SET
S3datDH <- RunHarmony(S3datDH, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S3datDH <- RunUMAP(S3datDH, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S3datDH)
DimPlot(S3datDH, reduction="HarmonyUMAP", pt.size = 1.5)

S3datDdegs <- FindAllMarkers(S3datDH, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S3datDdegs$cluster)) {
  genes <- c(genes,
             head(S3datDdegs$gene[grep(ind, 
                                       S3datDdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S3datDH, genes)


#_S4A####

S4datA <- RunPCA(S4datA, assay = "SCT", verbose = FALSE)
S4datA <- FindNeighbors(S4datA, reduction = "pca", dims = 1:30)
S4datA <- FindClusters(S4datA, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S4datA)
S4datA <- FindClusters(S4datA, verbose = FALSE, resolution = 0.7)
S4datA <- RunUMAP(S4datA, reduction = "pca", dims = 1:30)
DimPlot(S4datA, pt.size = 1.5)

S4datAdegs <- FindAllMarkers(S4datA, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S4datAdegs$cluster)) {
  genes <- c(genes,
             head(S4datAdegs$gene[grep(ind, 
                                       S4datAdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S4datA, genes)
DotPlot(S4datA, features = genes) + RotatedAxis()

#_S4B####

S4datB <- RunPCA(S4datB, assay = "SCT", verbose = FALSE)
S4datB <- FindNeighbors(S4datB, reduction = "pca", dims = 1:30)
S4datB <- FindClusters(S4datB, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S4datB)
S4datB <- FindClusters(S4datB, verbose = FALSE, resolution = 0.6)
S4datB <- RunUMAP(S4datB, reduction = "pca", dims = 1:30)
DimPlot(S4datB, pt.size = 1.5)


#gene expression visualisation 
FeaturePlot(S4datB, c(SC_AFibs), pt.size = 0.5, label.size = 0.5)


p1 <- DimPlot(S4datB, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(S4datB, label = TRUE, label.size = 3, )
p1 + p2

p3 <- SpatialFeaturePlot(S4datB, "PRG4")
p2 + p3

S4datBdegs <- FindAllMarkers(S4datB, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S4datBdegs$cluster)) {
  genes <- c(genes,
             head(S4datBdegs$gene[grep(ind, 
                                       S4datBdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S4datB, genes)
DotPlot(S4datB, features = genes) + RotatedAxis()


#_S4C####


S4datC <- RunPCA(S4datC, assay = "SCT", verbose = FALSE)
S4datC <- FindNeighbors(S4datC, reduction = "pca", dims = 1:30)
S4datC <- FindClusters(S4datC, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S4datC)
S4datC <- FindClusters(S4datC, verbose = FALSE, resolution = 0.7)
S4datC <- RunUMAP(S4datC, reduction = "pca", dims = 1:30)
DimPlot(S4datC, pt.size = 1.5)

S4datCdegs <- FindAllMarkers(S4datC, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S4datCdegs$cluster)) {
  genes <- c(genes,
             head(S4datCdegs$gene[grep(ind, 
                                       S4datCdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S4datC, genes)
DotPlot(S4datC, features = genes) + RotatedAxis()


#_S4D####


S4datD <- RunPCA(S4datD, assay = "SCT", verbose = FALSE)
S4datD <- FindNeighbors(S4datD, reduction = "pca", dims = 1:30)
S4datD <- FindClusters(S4datD, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S4datD)
S4datD <- FindClusters(S4datD, verbose = FALSE, resolution = 1)
S4datD <- RunUMAP(S4datD, reduction = "pca", dims = 1:30)
DimPlot(S4datD, pt.size = 1.5)

S4datDdegs <- FindAllMarkers(S4datD, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S4datDdegs$cluster)) {
  genes <- c(genes,
             head(S4datDdegs$gene[grep(ind, 
                                       S4datDdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S4datD, genes)
DotPlot(S4datD, features = genes) + RotatedAxis()


#_S5A####

#NO HARMONY
S5datA <- RunPCA(S5datA, assay = "SCT", verbose = FALSE)
S5datA <- FindNeighbors(S5datA, reduction = "pca", dims = 1:30)
S5datA <- FindClusters(S5datA, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S5datA)
S5datA <- FindClusters(S5datA, verbose = FALSE, resolution = 0.6)
S5datA <- RunUMAP(S5datA, reduction = "pca", dims = 1:30)
DimPlot(S5datA, pt.size = 1.5)

#HARMONY
S5datAH <- RunPCA(S5datA, assay = "SCT", verbose = FALSE)
S5datAH <- FindNeighbors(S5datAH, reduction = "pca", dims = 1:30)
S5datAH <- FindClusters(S5datAH, verbose = FALSE, resolution = c(0.2, 
                                                                 0.3,
                                                                 0.4,
                                                                 0.5,
                                                                 0.6,
                                                                 0.7,
                                                                 0.8,
                                                                 0.9,
                                                                 1.0))
clustree(S5datAH)
S5datAH <- FindClusters(S5datAH, verbose = FALSE, resolution = 0.6) #set
S5datAH <- RunHarmony(S5datAH, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S5datAH <- RunUMAP(S5datAH, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S5datAH, reduction="HarmonyUMAP", pt.size = 1.5)


S5datAdegs <- FindAllMarkers(S5datA, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S5datAdegs$cluster)) {
  genes <- c(genes,
             head(S5datAdegs$gene[grep(ind, 
                                       S5datAdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S5datA, genes)
DotPlot(S5datA, features = genes) + RotatedAxis()


#_S5B####


#NO HARMONY
S5datB <- RunPCA(S5datB, assay = "SCT", verbose = FALSE)
S5datB <- FindNeighbors(S5datB, reduction = "pca", dims = 1:30)
S5datB <- FindClusters(S5datB, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S5datB)
S5datB <- FindClusters(S5datB, verbose = FALSE, resolution = 0.6)
S5datB <- RunUMAP(S5datB, reduction = "pca", dims = 1:30)
DimPlot(S5datB, pt.size = 1.5)
FeaturePlot(S5datB, c(SC_MallMar))


#HARMONY
S5datBH <- RunPCA(S5datB, assay = "SCT", verbose = FALSE)
S5datBH <- FindNeighbors(S5datBH, reduction = "pca", dims = 1:30)
S5datBH <- FindClusters(S5datBH, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S5datBH)
S5datBH <- FindClusters(S5datBH, verbose = FALSE, resolution = 0.6) #set
S5datBH <- RunHarmony(S5datBH, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S5datBH <- RunUMAP(S5datBH, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S5datBH, reduction="HarmonyUMAP", pt.size = 1.5)

S5datBdegs <- FindAllMarkers(S5datBH, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S5datBdegs$cluster)) {
  genes <- c(genes,
             head(S5datBdegs$gene[grep(ind, 
                                       S5datBdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S5datBH, genes)
DotPlot(S5datBH, features = genes) + RotatedAxis()


#_S5C OA####


#NO HARMONY
S5datC <- RunPCA(S5datC, assay = "SCT", verbose = FALSE)
S5datC <- FindNeighbors(S5datC, reduction = "pca", dims = 1:30)
S5datC <- FindClusters(S5datC, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S5datC)
S5datC <- FindClusters(S5datC, verbose = FALSE, resolution = 0.5)
S5datC <- RunUMAP(S5datC, reduction = "pca", dims = 1:30)
DimPlot(S5datC, pt.size = 1.5)
FeaturePlot(S5datC, c(SC_MallMar))

#HARMONY
S5datCH <- RunPCA(S5datC, assay = "SCT", verbose = FALSE)
S5datCH <- FindNeighbors(S5datCH, reduction = "pca", dims = 1:30)
S5datCH <- FindClusters(S5datCH, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S5datCH)
S5datCH <- FindClusters(S5datCH, verbose = FALSE, resolution = 0.5) #set
S5datCH <- RunHarmony(S5datCH, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S5datCH <- RunUMAP(S5datCH, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S5datCH, reduction="HarmonyUMAP", pt.size = 1.5)

S5datCdegs <- FindAllMarkers(S5datCH, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S5datCdegs$cluster)) {
  genes <- c(genes,
             head(S5datCdegs$gene[grep(ind, 
                                       S5datCdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S5datCH, genes)
DotPlot(S5datCH, features = genes) + RotatedAxis()


#_S5D####

#NO HARMONY
S5datD <- RunPCA(S5datD, assay = "SCT", verbose = FALSE)
S5datD <- FindNeighbors(S5datD, reduction = "pca", dims = 1:30)
S5datD <- FindClusters(S5datD, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S5datD)
S5datD <- FindClusters(S5datD, verbose = FALSE, resolution = 0.4)
S5datD <- RunUMAP(S5datD, reduction = "pca", dims = 1:30)
DimPlot(S5datD, pt.size = 1.5)
FeaturePlot(S5datD, c(SC_MallMar))


#HARMONY
S5datDH <- RunPCA(S5datD, assay = "SCT", verbose = FALSE)
S5datDH <- FindNeighbors(S5datDH, reduction = "pca", dims = 1:30)
S5datDH <- FindClusters(S5datDH, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S5datDH)
S5datDH <- FindClusters(S5datDH, verbose = FALSE, resolution = 0.4) #set
S5datDH <- RunHarmony(S5datDH, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S5datDH <- RunUMAP(S5datDH, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S5datDH, reduction="HarmonyUMAP", pt.size = 1.5)

S5datDdegs <- FindAllMarkers(S5datDH, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S5datDdegs$cluster)) {
  genes <- c(genes,
             head(S5datDdegs$gene[grep(ind, 
                                       S5datDdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S5datDH, genes)
DotPlot(S5datDH, features = genes) + RotatedAxis()


#_S6A####

#HARMONY
S6datA <- RunPCA(S6datA, assay = "SCT", verbose = FALSE)
S6datA <- FindNeighbors(S6datA, reduction = "pca", dims = 1:30)
S6datA <- FindClusters(S6datA, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S6datA)
S6datA <- FindClusters(S6datA, verbose = FALSE, resolution = 0.6) #set
S6datA <- RunHarmony(S6datA, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S6datA <- RunUMAP(S6datA, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S6datA, reduction="HarmonyUMAP", pt.size = 1.5)

S6datAdegs <- FindAllMarkers(S6datA, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S6datAdegs$cluster)) {
  genes <- c(genes,
             head(S6datAdegs$gene[grep(ind, 
                                       S6datAdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S6datA, genes)
DotPlot(S6datA, features = genes) + RotatedAxis()


#_S6B####

#HARMONY
S6datB <- RunPCA(S6datB, assay = "SCT", verbose = FALSE)
S6datB <- FindNeighbors(S6datB, reduction = "pca", dims = 1:30)
S6datB <- FindClusters(S6datB, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S6datB)
S6datB <- FindClusters(S6datB, verbose = FALSE, resolution = 0.9) #set
S6datB <- RunHarmony(S6datB, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S6datB <- RunUMAP(S6datB, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S6datB, reduction="HarmonyUMAP", pt.size = 1.5)

S6datBdegs <- FindAllMarkers(S6datB, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S6datBdegs$cluster)) {
  genes <- c(genes,
             head(S6datBdegs$gene[grep(ind, 
                                       S6datBdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S6datB, genes)
DotPlot(S6datB, features = genes) + RotatedAxis()


#_S6C####

#HARMONY
S6datC <- RunPCA(S6datC, assay = "SCT", verbose = FALSE)
S6datC <- FindNeighbors(S6datC, reduction = "pca", dims = 1:30)
S6datC <- FindClusters(S6datC, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S6datC)
S6datC <- FindClusters(S6datC, verbose = FALSE, resolution = 0.7) #set
S6datC <- RunHarmony(S6datC, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S6datC <- RunUMAP(S6datC, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S6datC, reduction="HarmonyUMAP", pt.size = 1.5)


S6datCdegs <- FindAllMarkers(S6datC, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S6datCdegs$cluster)) {
  genes <- c(genes,
             head(S6datCdegs$gene[grep(ind, 
                                       S6datCdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S6datC, genes)
DotPlot(S6datC, features = genes) + RotatedAxis()

#_S6D####

#HARMONY
S6datD <- RunPCA(S6datD, assay = "SCT", verbose = FALSE)
S6datD <- FindNeighbors(S6datD, reduction = "pca", dims = 1:30)
S6datD <- FindClusters(S6datD, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S6datD)
S6datD <- FindClusters(S6datD, verbose = FALSE, resolution = 0.8) #set
S6datD <- RunHarmony(S6datD, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S6datD <- RunUMAP(S6datD, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S6datC, reduction="HarmonyUMAP", pt.size = 1.5)

S6datDdegs <- FindAllMarkers(S6datD, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S6datDdegs$cluster)) {
  genes <- c(genes,
             head(S6datDdegs$gene[grep(ind, 
                                       S6datDdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S6datD, genes)
DotPlot(S6datD, features = genes) + RotatedAxis()

#_S7A####

#HARMONY
S7datA <- RunPCA(S7datA, assay = "SCT", verbose = FALSE)
S7datA <- FindNeighbors(S7datA, reduction = "pca", dims = 1:30)
S7datA <- FindClusters(S7datA, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S7datA)
S7datA <- FindClusters(S7datA, verbose = FALSE, resolution = 0.7) #set
S7datA <- RunHarmony(S7datA, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S7datA <- RunUMAP(S7datA, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S7datA, reduction="HarmonyUMAP", pt.size = 1.5)

S7datAdegs <- FindAllMarkers(S7datA, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S7datAdegs$cluster)) {
  genes <- c(genes,
             head(S7datAdegs$gene[grep(ind, 
                                       S7datAdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S7datA, genes)
DotPlot(S7datA, features = genes) + RotatedAxis()

#_S7B####

#HARMONY
S7datB <- RunPCA(S7datB, assay = "SCT", verbose = FALSE)
S7datB <- FindNeighbors(S7datB, reduction = "pca", dims = 1:30)
S7datB <- FindClusters(S7datB, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S7datB)
S7datB <- FindClusters(S7datB, verbose = FALSE, resolution = 0.7) #set
S7datB <- RunHarmony(S7datB, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S7datB <- RunUMAP(S7datB, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S7datB, reduction="HarmonyUMAP", pt.size = 1.5)

S7datBdegs <- FindAllMarkers(S7datB, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S7datBdegs$cluster)) {
  genes <- c(genes,
             head(S7datBdegs$gene[grep(ind, 
                                       S7datBdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S7datB, genes)
DotPlot(S7datB, features = genes) + RotatedAxis()


#_S7C####

#HARMONY
S7datC <- RunPCA(S7datC, assay = "SCT", verbose = FALSE)
S7datC <- FindNeighbors(S7datC, reduction = "pca", dims = 1:30)
S7datC <- FindClusters(S7datC, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S7datC)
S7datC <- FindClusters(S7datC, verbose = FALSE, resolution = 0.5) #set
S7datC <- RunHarmony(S7datC, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S7datC <- RunUMAP(S7datC, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S7datC, reduction="HarmonyUMAP", pt.size = 1.5)


S7datCdegs <- FindAllMarkers(S7datC, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S7datCdegs$cluster)) {
  genes <- c(genes,
             head(S7datCdegs$gene[grep(ind, 
                                       S7datCdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S7datC, genes)
DotPlot(S7datC, features = genes) + RotatedAxis()


#_S7D####

#HARMONY
S7datD <- RunPCA(S7datD, assay = "SCT", verbose = FALSE)
S7datD <- FindNeighbors(S7datD, reduction = "pca", dims = 1:30)
S7datD <- FindClusters(S7datD, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S7datD)
S7datD <- FindClusters(S7datD, verbose = FALSE, resolution = 0.6) #set
S7datD <- RunHarmony(S7datD, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S7datD <- RunUMAP(S7datD, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S7datD, reduction="HarmonyUMAP", pt.size = 1.5)


S7datDdegs <- FindAllMarkers(S7datD, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S7datDdegs$cluster)) {
  genes <- c(genes,
             head(S7datDdegs$gene[grep(ind, 
                                       S7datDdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S7datD, genes)
DotPlot(S7datD, features = genes) + RotatedAxis()

SpatialDimPlot(S7datD)


#_S8A####

#HARMONY
S8datA <- RunPCA(S8datA, assay = "SCT", verbose = FALSE)
S8datA <- FindNeighbors(S8datA, reduction = "pca", dims = 1:30)
S8datA <- FindClusters(S8datA, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S8datA)
S8datA <- FindClusters(S8datA, verbose = FALSE, resolution = 0.6) #set
S8datA <- RunHarmony(S8datA, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S8datA <- RunUMAP(S8datA, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S8datA, reduction="HarmonyUMAP", pt.size = 1.5)


S8datAdegs <- FindAllMarkers(S8datA, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S8datAdegs$cluster)) {
  genes <- c(genes,
             head(S8datAdegs$gene[grep(ind, 
                                       S8datAdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S8datA, genes)
DotPlot(S8datA, features = genes) + RotatedAxis()


#_S8B####

#HARMONY
S8datB <- RunPCA(S8datB, assay = "SCT", verbose = FALSE)
S8datB <- FindNeighbors(S8datB, reduction = "pca", dims = 1:30)
S8datB <- FindClusters(S8datB, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S8datB)
S8datB <- FindClusters(S8datB, verbose = FALSE, resolution = 0.5) #set
S8datB <- RunHarmony(S8datB, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S8datB <- RunUMAP(S8datB, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S8datB, reduction="HarmonyUMAP", pt.size = 1.5)


S8datBdegs <- FindAllMarkers(S8datB, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S8datBdegs$cluster)) {
  genes <- c(genes,
             head(S8datBdegs$gene[grep(ind, 
                                       S8datBdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S8datB, genes)
DotPlot(S8datB, features = genes) + RotatedAxis()

#_S8C####

#HARMONY
S8datC <- RunPCA(S8datC, assay = "SCT", verbose = FALSE)
S8datC <- FindNeighbors(S8datC, reduction = "pca", dims = 1:30)
S8datC <- FindClusters(S8datC, verbose = FALSE, resolution = c(0.2, 
                                                               0.3,
                                                               0.4,
                                                               0.5,
                                                               0.6,
                                                               0.7,
                                                               0.8,
                                                               0.9,
                                                               1.0))
clustree(S8datC)
S8datC <- FindClusters(S8datC, verbose = FALSE, resolution = 0.6) #set
S8datC <- RunHarmony(S8datC, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
S8datC <- RunUMAP(S8datC, reduction = "harmony", dims = 1:30, reduction.key="HarmUMAP_", reduction.name= "HarmonyUMAP" )
DimPlot(S8datC, reduction="HarmonyUMAP", pt.size = 1.5)


S8datCdegs <- FindAllMarkers(S8datC, test.use = "MAST")

genes <- vector(mode= "character")
for(ind in levels(S8datCdegs$cluster)) {
  genes <- c(genes,
             head(S8datCdegs$gene[grep(ind, 
                                       S8datCdegs$cluster)],
                  n=10))
}
genes <- unique(genes)
DoHeatmap(S8datC, genes)
DotPlot(S8datC, features = genes) + RotatedAxis()


#MERGING####

merging <- merge(x = S1datAH, c(S1datBH, 
                                S2datAH, S2datCH, S2datDH, 
                                S3datAH, S3datBH, S3datCH, S3datDH, 
                                S4datA, S4datB, S4datC, S4datD, 
                                S5datAH, S5datBH, S5datCH, S5datDH,
                                S6datA, S6datB, S6datC, S6datD, 
                                S7datA, S7datB, S7datC, S7datD,
                                S8datA, S8datB, S8datC))


current.sample.ids <- c("Slide1SampleA", "Slide1SampleB", "Slide2SampleA", "Slide2SampleC", "Slide2SampleD", "Slide3SampleA", "Slide3SampleB", "Slide3SampleC", "Slide3SampleD", "Slide4SampleA", "Slide4SampleB", 
                        "Slide4SampleC", "Slide4SampleD", "Slide5SampleA", "Slide5SampleB", "Slide5SampleC", "Slide5SampleD", 
                        "Slide6SampleA", "Slide6SampleB", "Slide6SampleC", "Slide6SampleD", "Slide7SampleA", "Slide7SampleB", "Slide7SampleC", "Slide7SampleD", "Slide8SampleA", "Slide8SampleB", "Slide8SampleC")

new.sample.ids <- c("RA","RA","RA","RA","RA","RA","RA","RA","RA","RA","RA",
                    "OA","OA","OA","OA","OA","OA",
                    "RA","RA","RA","RA","RA","RA","RA","RA","RA","RA","RA")

merging$condition <- merging$orig.ident
merging@meta.data[["condition"]] <- plyr::mapvalues(x = merging@meta.data[["condition"]], from = current.sample.ids, to = new.sample.ids)



current.sample.ids <- c("Slide1SampleA", "Slide1SampleB", "Slide2SampleA", "Slide2SampleC", "Slide2SampleD", "Slide3SampleA", "Slide3SampleB", "Slide3SampleC", "Slide3SampleD", "Slide4SampleA", "Slide4SampleB", 
                        "Slide4SampleC", "Slide4SampleD", "Slide5SampleA", "Slide5SampleB", "Slide5SampleC", "Slide5SampleD", 
                        "Slide6SampleA", "Slide6SampleB", "Slide6SampleC", "Slide6SampleD", "Slide7SampleA", "Slide7SampleB", "Slide7SampleC", "Slide7SampleD", "Slide8SampleA", "Slide8SampleB", "Slide8SampleC")


new.sample.ids <- c("Non-Resolving","Non-Resolving","Non-Resolving","Non-Resolving","Non-Resolving","Non-Resolving","Non-Resolving","Non-Resolving","Non-Resolving","Non-Resolving","Non-Resolving",
                    "OA","OA","OA","OA","OA","OA",
                    "Non-Resolving","Resolving","Non-Resolving","Resolving","Non-Resolving","Resolving","Non-Resolving","Resolving","Non-Resolving","Non-Resolving","Resolving")

merging$resvsnon <- merging$orig.ident
merging@meta.data[["resvsnon"]] <- plyr::mapvalues(x = merging@meta.data[["resvsnon"]], from = current.sample.ids, to = new.sample.ids)


save.image("merging.RData")


Idents(merging)<- "pathotype"
levels(merging)<-c("Diffuse", "Lymphoid", "PI", "OA")
merging$pathotype<-merging@active.ident


#QC 

merging[["percent.mt"]] <- PercentageFeatureSet(merging, pattern = "^MT-")
VlnPlot(merging, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
VlnPlot(merging, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), group.by = "orig.ident")
VlnPlot(merging, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), group.by = "condition")


#HARMONY
DefaultAssay (merging) <- "SCT"

VariableFeatures(merging) <- c(VariableFeatures(S1datAH),
                               VariableFeatures(S1datBH),
                               VariableFeatures(S2datAH),
                               VariableFeatures(S2datCH),
                               VariableFeatures(S2datDH),
                               VariableFeatures(S3datAH),
                               VariableFeatures(S3datBH),
                               VariableFeatures(S3datCH),
                               VariableFeatures(S3datAH),
                               VariableFeatures(S3datDH),
                               VariableFeatures(S4datA),
                               VariableFeatures(S4datB),
                               VariableFeatures(S4datC),
                               VariableFeatures(S4datD),
                               VariableFeatures(S5datAH),
                               VariableFeatures(S5datBH),
                               VariableFeatures(S5datCH),
                               VariableFeatures(S5datDH),
                               VariableFeatures(S6datA),
                               VariableFeatures(S6datB),
                               VariableFeatures(S6datC),
                               VariableFeatures(S6datD),
                               VariableFeatures(S7datA),
                               VariableFeatures(S7datB),
                               VariableFeatures(S7datC),
                               VariableFeatures(S7datD),
                               VariableFeatures(S8datA),
                               VariableFeatures(S8datB),
                               VariableFeatures(S8datC))

merging <- RunPCA(merging, verbose = FALSE)  
merging <- RunHarmony(merging, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 20)
#RunUMAP
merging <- RunUMAP(merging, reduction = "harmony", dims = 1:30, reduction.key="HarmonisedUMAP_", reduction.name= "HarmonyUMAP")
merging <- FindNeighbors(merging, dims = 1:30, reduction = "harmony")
merging <- FindClusters(merging, 
                            verbose = FALSE, 
                            resolution = c(0.1,0.2,0.3,0.4, 0.5, 0.6, 0.7, 0.8,0.9))
clustree(merging)
merging <- FindClusters(merging, verbose = FALSE, resolution = 0.2)
DimPlot(merging, pt.size = 1.5, reduction = "HarmonyUMAP", label = TRUE, label.size = 8, label.box = T, cols = palette)
DimPlot(merging, pt.size = 1.5, reduction = "HarmonyUMAP", group.by = "pathotype", split.by = "pathotype")
DimPlot(merging, pt.size = 1.5, reduction = "HarmonyUMAP", group.by = "orig.ident")

test <- DimPlot(merging, pt.size = 1.5, reduction = "HarmonyUMAP", group.by = "cluster_ids")
saveRDS(test, file="test.rds")

DoHeatmap(merging, genes)

global_markers <- mergingAllHMarkers %>%
  group_by(cluster)
write.csv(global_markers, file = "mergingAllHMarkers_0.2.csv")

top40_markers <- mergingAllHMarkers %>%
  group_by(cluster) %>%
  top_n(n = 40, wt = avg_log2FC)
write.csv(top40_markers, file = "mergingAllHMarkerstop40markers_0.2.csv")

top8_markers <- mergingAllHMarkers %>%
  group_by(cluster) %>%
  top_n(n = 8, wt = avg_log2FC)
write.csv(top8_markers, file = "mergingAllHMarkerstop8markers_0.2.csv")


DoHeatmap(merging, genes)

lymphoidUMAP <- merging[,grepl("Slide3SampleA", merging$orig.ident, ignore.case=TRUE)]
DimPlot(lymphoidUMAP, pt.size = 1.5, reduction = "HarmonyUMAP", label = TRUE, label.size = 6, label.box = T, repel = T, cols = palette)

PIUMAP <- merging[,grepl("Slide1SampleB", merging$orig.ident, ignore.case=TRUE)]
DimPlot(PIUMAP, pt.size = 1.5, reduction = "HarmonyUMAP",  label = TRUE, label.size = 6, label.box = T, repel = T, cols = palette)

diffuseUMAP <- merging[,grepl("Slide8SampleC", merging$orig.ident, ignore.case=TRUE)]
DimPlot(diffuseUMAP, pt.size = 1.5, reduction = "HarmonyUMAP",  label = TRUE, label.size = 6, label.box = T, repel = T, cols = palette)

OAUMAP <- merging[,grepl("Slide5SampleB", merging$orig.ident, ignore.case=TRUE)]
DimPlot(OAUMAP, pt.size = 1.5, reduction = "HarmonyUMAP", group.by = "cluster_ids", label = TRUE, label.size = 6, label.box = T, repel = T)


#RuntSNE
merging <- RunTSNE(merging, reduction = "harmony", dims = 1:30, reduction.key="HarmonisedtSNE_", reduction.name= "HarmonytSNE")
DimPlot(merging, pt.size = 1.5, reduction = "HarmonytSNE")
DimPlot(merging, pt.size = 1.5, reduction = "HarmonytSNE", group.by = "condition", split.by = "condition")

mergingAllHMarkers <- FindAllMarkers(merging)

genes <- vector(mode= "character")
for(ind in levels(mergingAllHMarkers$cluster)) {
  genes <- c(genes,
             head(mergingAllHMarkers$gene[grep(ind, 
                                               mergingAllHMarkers$cluster)],
                  n=10))
}
genes <- unique(genes)

DoHeatmap(merging, genes)

global_markers <- mergingAllHMarkers %>%
  group_by(cluster)
write.csv(global_markers, file = "mergingAllHMarkers_0.2.csv")

top40_markers <- mergingAllHMarkers %>%
  group_by(cluster) %>%
  top_n(n = 40, wt = avg_log2FC)
write.csv(top40_markers, file = "mergingAllHMarkerstop40markers_0.2.csv")

top8_markers <- mergingAllHMarkers %>%
  group_by(cluster) %>%
  top_n(n = 8, wt = avg_log2FC)
write.csv(top8_markers, file = "mergingAllHMarkerstop8markers_0.2.csv")

DoHeatmap(merging, genes)


#palette
palette <- c("T Cell Rich Niche"="#0C5BB0FF", "Erythrocytes"="#EE0011FF", "B Cell Rich Niche"="#15983DFF", "COMP+ Fibroblast Niche"="#EC579AFF", "APOD+ GAS5+ FABP4+"="#FA6B09FF", "Vascular Niche"="#149BEDFF", "Lining Layer Cells"= "#A1C720FF")


