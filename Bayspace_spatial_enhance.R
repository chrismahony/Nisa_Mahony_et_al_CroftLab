#can you add the code you ran to do the spatial enhance? Thanks!

#Convert to SCE

merging.sce = Seurat::DietSeurat(merging, graphs = "pca") #slim down Seurat obj prior to conversion
sce = as.SingleCellExperiment(merging.sce) #convert seurat to SCE
coord1 = merging@images$slice1@coordinates #replace 1.2, 1.3, 1.4 with actual names of your slices
coord2 = merging@images$slice1.1@coordinates
coord3 = merging@images$slice1.2@coordinates 
coord4 = merging@images$slice1.3@coordinates
coord5 = merging@images$slice1.4@coordinates
coord6 = merging@images$slice1.5@coordinates
coord7 = merging@images$slice1.6@coordinates
coord8 = merging@images$slice1.7@coordinates
coord9 = merging@images$slice1.8@coordinates
coord10 = merging@images$slice1.9@coordinates
coord11 = merging@images$slice1.10@coordinates
coord12 = merging@images$slice1.11@coordinates
coord13 = merging@images$slice1.12@coordinates
coord14 = merging@images$slice1.13@coordinates
coord15 = merging@images$slice1.14@coordinates
coord16 = merging@images$slice1.15@coordinates
coord17 = merging@images$slice1.16@coordinates
coord18 = merging@images$slice1.17@coordinates
coord19 = merging@images$slice1.18@coordinates
coord20 = merging@images$slice1.19@coordinates
coord21 = merging@images$slice1.20@coordinates
coord22 = merging@images$slice1.21@coordinates
coord23 = merging@images$slice1.22@coordinates
coord24 = merging@images$slice1.23@coordinates
coord25 = merging@images$slice1.24@coordinates
coord26 = merging@images$slice1.25@coordinates
coord27 = merging@images$slice1.26@coordinates
coord28 = merging@images$slice1.27@coordinates

coords = rbind(coord1, coord2, coord3, coord4, coord5, coord6, coord7, coord8, coord9, coord10, 
coord11, coord12, coord13, coord14, coord15, coord16, coord17, coord18, coord19, coord20, coord21, 
coord22, coord23, coord24, coord25, coord26, coord27, coord28)

colData(sce) = cbind(colData(sce), coords)


#BayesSpace Workflow

sce1 = sce[ , sce$orig.ident == "Slide1SampleA"] #replace sample with the name of the column containing sample IDs
sce1.1 = sce[ , sce$orig.ident == "Slide1SampleB"]
sce1.2 = sce[ , sce$orig.ident == "Slide2SampleA"]
sce1.3 = sce[ , sce$orig.ident == "Slide2SampleC"]
sce1.4 = sce[ , sce$orig.ident == "Slide2SampleD"]
sce1.5 = sce[ , sce$orig.ident == "Slide3SampleA"]
sce1.6 = sce[ , sce$orig.ident == "Slide3SampleB"]
sce1.7 = sce[ , sce$orig.ident == "Slide3SampleC"]
sce1.8 = sce[ , sce$orig.ident == "Slide3SampleD"]
sce1.9 = sce[ , sce$orig.ident == "Slide4SampleA"]
sce1.10 = sce[ , sce$orig.ident == "Slide4SampleB"]
sce1.11 = sce[ , sce$orig.ident == "Slide4SampleC"]
sce1.12 = sce[ , sce$orig.ident == "Slide4SampleD"]
sce1.13 = sce[ , sce$orig.ident == "Slide5SampleA"]
sce1.14 = sce[ , sce$orig.ident == "Slide5SampleB"]
sce1.15 = sce[ , sce$orig.ident == "Slide5SampleC"]
sce1.16 = sce[ , sce$orig.ident == "Slide5SampleD"]
sce1.17 = sce[ , sce$orig.ident == "Slide6SampleA"]
sce1.18 = sce[ , sce$orig.ident == "Slide6SampleB"]
sce1.19 = sce[ , sce$orig.ident == "Slide6SampleC"]
sce1.20 = sce[ , sce$orig.ident == "Slide6SampleD"]
sce1.21 = sce[ , sce$orig.ident == "Slide7SampleA"]
sce1.22 = sce[ , sce$orig.ident == "Slide7SampleB"]
sce1.23 = sce[ , sce$orig.ident == "Slide7SampleC"]
sce1.24 = sce[ , sce$orig.ident == "Slide7SampleD"]
sce1.25 = sce[ , sce$orig.ident == "Slide8SampleA"]
sce1.26 = sce[ , sce$orig.ident == "Slide8SampleB"]
sce1.27 = sce[ , sce$orig.ident == "Slide8SampleC"]



sce1 = spatialPreprocess(sce1, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1 = spatialCluster(sce1, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1e <- spatialEnhance(sce1, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)

sce1.1 = spatialPreprocess(sce1.1, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.1 = spatialCluster(sce1.1, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.1e <- spatialEnhance(sce1.1, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.2 = spatialPreprocess(sce1.2, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.2 = spatialCluster(sce1.2, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.2e <- spatialEnhance(sce1.2, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.3 = spatialPreprocess(sce1.3, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.3 = spatialCluster(sce1.3, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.3e <- spatialEnhance(sce1.3, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.4 = spatialPreprocess(sce1.4, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.4 = spatialCluster(sce1.4, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.4e <- spatialEnhance(sce1.4, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.5 = spatialPreprocess(sce1.5, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.5 = spatialCluster(sce1.5, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.5e <- spatialEnhance(sce1.5, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)

sce1.6 = spatialPreprocess(sce1.6, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.6 = spatialCluster(sce1.6, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.6e <- spatialEnhance(sce1.6, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.7 = spatialPreprocess(sce1.7, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.7 = spatialCluster(sce1.7, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.7e <- spatialEnhance(sce1.7, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.8 = spatialPreprocess(sce1.8, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.8 = spatialCluster(sce1.8, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.8e <- spatialEnhance(sce1.8, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.9 = spatialPreprocess(sce1.9, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.9 = spatialCluster(sce1.9, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.9e <- spatialEnhance(sce1.9, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.10 = spatialPreprocess(sce1.10, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.10 = spatialCluster(sce1.10, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.10e <- spatialEnhance(sce1.10, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.11 = spatialPreprocess(sce1.11, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.11 = spatialCluster(sce1.11, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.11e <- spatialEnhance(sce1.11, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.12 = spatialPreprocess(sce1.12, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.12 = spatialCluster(sce1.12, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.12e <- spatialEnhance(sce1.12, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.13 = spatialPreprocess(sce1.13, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.13 = spatialCluster(sce1.13, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.13e <- spatialEnhance(sce1.13, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.14 = spatialPreprocess(sce1.14, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.14 = spatialCluster(sce1.14, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.14e <- spatialEnhance(sce1.14, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.15 = spatialPreprocess(sce1.15, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.15 = spatialCluster(sce1.15, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.15e <- spatialEnhance(sce1.15, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.16 = spatialPreprocess(sce1.16, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.16 = spatialCluster(sce1.16, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.16e <- spatialEnhance(sce1.16, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.17 = spatialPreprocess(sce1.17, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.17 = spatialCluster(sce1.17, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.17e <- spatialEnhance(sce1.17, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.18 = spatialPreprocess(sce1.18, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.18 = spatialCluster(sce1.18, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.18e <- spatialEnhance(sce1.18, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.19 = spatialPreprocess(sce1.19, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.19 = spatialCluster(sce1.19, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.19e <- spatialEnhance(sce1.19, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.20 = spatialPreprocess(sce1.20, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.20 = spatialCluster(sce1.20, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.20e <- spatialEnhance(sce1.20, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


sce1.21 = spatialPreprocess(sce1.21, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.21 = spatialCluster(sce1.21, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.21e <- spatialEnhance(sce1.21, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)

sce1.22 = spatialPreprocess(sce1.22, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.22 = spatialCluster(sce1.22, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.22e <- spatialEnhance(sce1.22, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)

sce1.23 = spatialPreprocess(sce1.23, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.23 = spatialCluster(sce1.23, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.23e <- spatialEnhance(sce1.23, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)

sce1.24 = spatialPreprocess(sce1.24, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.24 = spatialCluster(sce1.24, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.24e <- spatialEnhance(sce1.24, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)

sce1.25 = spatialPreprocess(sce1.25, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.25 = spatialCluster(sce1.25, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.25e <- spatialEnhance(sce1.25, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)

sce1.26 = spatialPreprocess(sce1.26, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.26 = spatialCluster(sce1.26, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.26e <- spatialEnhance(sce1.26, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)

sce1.27 = spatialPreprocess(sce1.27, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
sce1.27 = spatialCluster(sce1.27, nrep = 50000, burn.in = 100, q = 7) #quickly cluster via BayesSpace q = number of clusters

sce1.27e <- spatialEnhance(sce1.27, q=7, platform="Visium",
                          nrep=50000, gamma=2, 
                          verbose=TRUE, save.chain=TRUE,
                          jitter_scale=3.5, jitter_prior=0.3)


## Using the same 2,000 HVGs previously computed for PCA


set.seed(101)
dec <- scran::modelGeneVar(sce1)
top <- scran::getTopHVGs(dec, n = 20000)

sce1e_de <- enhanceFeatures(sce1e, sce1, 
                            model="xgboost",
                            feature_names=top,
                            nrounds=0)

## DE Using the same 2,000 HVGs previously computed for PCA, excluding ribosomal
hvgs <- top[grep("^RP[LS]", top, invert=TRUE)]

sce1e_de <- enhanceFeatures(sce1e, sce1, 
                              model="xgboost",
                              feature_names=top,
                              nrounds=0)
## Seurat DE workflow

## Convert SCE to Seurat object and use BayesSpace cluster as identifier
sobj <- Seurat::CreateSeuratObject(counts=logcounts(sce1e_de),
                                   assay='Spatial',
                                   meta.data=as.data.frame(colData(sce1e_de)))

sobj <- Seurat::SetIdent(sobj, value = "spatial.cluster")

#DO THIS FOR ALL SAMPLES


counts_S1a<-seurat1A@assays[["Spatial"]]@counts
counts_S1a[is.na(counts_S1a)] <- 0
seurat1A <- CreateSeuratObject(counts = counts_S1a, min.cells =0, min.features =0, assay= "Spatial")
seurat1A$sample <-'seurat1A'

counts_S1b<-seurat1B@assays[["Spatial"]]@counts
counts_S1b[is.na(counts_S1b)] <- 0
seurat1B <- CreateSeuratObject(counts = counts_S1b, min.cells =0, min.features =0, assay= "Spatial")
seurat1B$sample <-'seurat1B'

counts_S2a<-seurat2A@assays[["Spatial"]]@counts
counts_S2a[is.na(counts_S2a)] <- 0
seurat2A <- CreateSeuratObject(counts = counts_S2a, min.cells =0, min.features =0, assay= "Spatial")
seurat2A$sample <-'seurat2A'

counts_S2c<-seurat2C@assays[["Spatial"]]@counts
counts_S2c[is.na(counts_S2c)] <- 0
seurat2C <- CreateSeuratObject(counts = counts_S2c, min.cells =0, min.features =0, assay= "Spatial")
seurat2C$sample <-'seurat2C'


counts_S2d<-seurat2D@assays[["Spatial"]]@counts
counts_S2d[is.na(counts_S2d)] <- 0
seurat2D <- CreateSeuratObject(counts = counts_S2d, min.cells =0, min.features =0, assay= "Spatial")
seurat2D$sample <-'seurat2D'


counts_S3a<-seurat3A@assays[["Spatial"]]@counts
counts_S3a[is.na(counts_S3a)] <- 0
seurat3A <- CreateSeuratObject(counts = counts_S3a, min.cells =0, min.features =0, assay= "Spatial")
seurat3A$sample <-'seurat3A'


counts_S3b<-seurat3B@assays[["Spatial"]]@counts
counts_S3b[is.na(counts_S3b)] <- 0
seurat3B <- CreateSeuratObject(counts = counts_S3b, min.cells =0, min.features =0, assay= "Spatial")
seurat3B$sample <-'seurat3B'


counts_S3c<-seurat3C@assays[["Spatial"]]@counts
counts_S3c[is.na(counts_S3c)] <- 0
seurat3C <- CreateSeuratObject(counts = counts_S3c, min.cells =0, min.features =0, assay= "Spatial")
seurat3C$sample <-'seurat3C'

counts_S3d<-seurat3D@assays[["Spatial"]]@counts
counts_S3d[is.na(counts_S3d)] <- 0
seurat3D <- CreateSeuratObject(counts = counts_S3d, min.cells =0, min.features =0, assay= "Spatial")
seurat3D$sample <-'seurat3D'

counts_S4a<-seurat4A@assays[["Spatial"]]@counts
counts_S4a[is.na(counts_S4a)] <- 0
seurat4A <- CreateSeuratObject(counts = counts_S4a, min.cells =0, min.features =0, assay= "Spatial")
seurat4A$sample <-'seurat4A'

counts_S4b<-seurat4B@assays[["Spatial"]]@counts
counts_S4b[is.na(counts_S4b)] <- 0
seurat4B <- CreateSeuratObject(counts = counts_S4b, min.cells =0, min.features =0, assay= "Spatial")
seurat4B$sample <-'seurat4B'

counts_S4c<-seurat4C@assays[["Spatial"]]@counts
counts_S4c[is.na(counts_S4c)] <- 0
seurat4C <- CreateSeuratObject(counts = counts_S4c, min.cells =0, min.features =0, assay= "Spatial")
seurat4C$sample <-'seurat4C'

counts_S4d<-seurat4D@assays[["Spatial"]]@counts
counts_S4d[is.na(counts_S4d)] <- 0
seurat4D <- CreateSeuratObject(counts = counts_S4d, min.cells =0, min.features =0, assay= "Spatial")
seurat4D$sample <-'seurat4D'

counts_S5a<-seurat5A@assays[["Spatial"]]@counts
counts_S5a[is.na(counts_S5a)] <- 0
seurat5A <- CreateSeuratObject(counts = counts_S5a, min.cells =0, min.features =0, assay= "Spatial")
seurat5A$sample <-'seurat5A'

counts_S5b<-seurat5B@assays[["Spatial"]]@counts
counts_S5b[is.na(counts_S5b)] <- 0
seurat5B <- CreateSeuratObject(counts = counts_S5b, min.cells =0, min.features =0, assay= "Spatial")
seurat5B$sample <-'seurat5B'

counts_S5c<-seurat5C@assays[["Spatial"]]@counts
counts_S5c[is.na(counts_S5c)] <- 0
seurat5C <- CreateSeuratObject(counts = counts_S5c, min.cells =0, min.features =0, assay= "Spatial")
seurat5C$sample <-'seurat5C'

counts_S5d<-seurat5D@assays[["Spatial"]]@counts
counts_S5d[is.na(counts_S5d)] <- 0
seurat5D <- CreateSeuratObject(counts = counts_S5d, min.cells =0, min.features =0, assay= "Spatial")
seurat5D$sample <-'seurat5D'

counts_S6a<-seurat6A@assays[["Spatial"]]@counts
counts_S6a[is.na(counts_S6a)] <- 0
seurat6A <- CreateSeuratObject(counts = counts_S6a, min.cells =0, min.features =0, assay= "Spatial")
seurat6A$sample <-'seurat6A'

counts_S6b<-seurat6B@assays[["Spatial"]]@counts
counts_S6b[is.na(counts_S6b)] <- 0
seurat6B <- CreateSeuratObject(counts = counts_S6b, min.cells =0, min.features =0, assay= "Spatial")
seurat6B$sample <-'seurat6B'

counts_S6c<-seurat6C@assays[["Spatial"]]@counts
counts_S6c[is.na(counts_S6c)] <- 0
seurat6C <- CreateSeuratObject(counts = counts_S6c, min.cells =0, min.features =0, assay= "Spatial")
seurat6C$sample <-'seurat6C'

counts_S6d<-seurat6D@assays[["Spatial"]]@counts
counts_S6d[is.na(counts_S6d)] <- 0
seurat6D <- CreateSeuratObject(counts = counts_S6d, min.cells =0, min.features =0, assay= "Spatial")
seurat6D$sample <-'seurat6D'

counts_S7a<-seurat7A@assays[["Spatial"]]@counts
counts_S7a[is.na(counts_S7a)] <- 0
seurat7A <- CreateSeuratObject(counts = counts_S7a, min.cells =0, min.features =0, assay= "Spatial")
seurat7A$sample <-'seurat7A'

counts_S7b<-seurat7B@assays[["Spatial"]]@counts
counts_S7b[is.na(counts_S7b)] <- 0
seurat7B <- CreateSeuratObject(counts = counts_S7b, min.cells =0, min.features =0, assay= "Spatial")
seurat7B$sample <-'seurat7B'

counts_S7c<-seurat7C@assays[["Spatial"]]@counts
counts_S7c[is.na(counts_S7c)] <- 0
seurat7C <- CreateSeuratObject(counts = counts_S7c, min.cells =0, min.features =0, assay= "Spatial")
seurat7C$sample <-'seurat7C'

counts_S7d<-seurat7D@assays[["Spatial"]]@counts
counts_S7d[is.na(counts_S7d)] <- 0
seurat7D <- CreateSeuratObject(counts = counts_S7d, min.cells =0, min.features =0, assay= "Spatial")
seurat7D$sample <-'seurat7D'

counts_S8a<-seurat8A@assays[["Spatial"]]@counts
counts_S8a[is.na(counts_S8a)] <- 0
seurat8A <- CreateSeuratObject(counts = counts_S8a, min.cells =0, min.features =0, assay= "Spatial")
seurat8A$sample <-'seurat8A'

counts_S8b<-seurat8B@assays[["Spatial"]]@counts
counts_S8b[is.na(counts_S8b)] <- 0
seurat8B <- CreateSeuratObject(counts = counts_S8b, min.cells =0, min.features =0, assay= "Spatial")
seurat8B$sample <-'seurat8B'

counts_S8c<-seurat8C@assays[["Spatial"]]@counts
counts_S8c[is.na(counts_S8c)] <- 0
seurat8C <- CreateSeuratObject(counts = counts_S8c, min.cells =0, min.features =0, assay= "Spatial")
seurat8C$sample <-'seurat8C'

save.image("BS2S_counts.RData")

#_S1A

seurat1A <- SCTransform(seurat1A, assay = "Spatial", verbose = FALSE)
seurat1A <- RunPCA(seurat1A, assay = "SCT", verbose = FALSE)
seurat1A <- FindNeighbors(seurat1A, reduction = "pca", dims = 1:30)
seurat1A <- RunUMAP(seurat1A, reduction = "pca", dims = 1:30)
seurat1A_DE <- FindAllMarkers(seurat1A, test.use = "MAST")
seurat1A <- FindVariableFeatures(seurat1A, assay = "SCT", features = VariableFeatures(seurat1A), 
                                          selection.method = "vst")

#_S1B
seurat1B <- SCTransform(seurat1B, assay = "Spatial", verbose = FALSE)
seurat1B <- RunPCA(seurat1B, assay = "SCT", verbose = FALSE)
seurat1B <- FindNeighbors(seurat1B, reduction = "pca", dims = 1:30)
seurat1B <- RunUMAP(seurat1B, reduction = "pca", dims = 1:30)
seurat1B_DE <- FindAllMarkers(seurat1B, test.use = "MAST")
seurat1B <- FindVariableFeatures(seurat1B, assay = "SCT", features = VariableFeatures(seurat1B), 
                                 selection.method = "vst")

#_S2A

seurat2A <- SCTransform(seurat2A, assay = "Spatial", verbose = FALSE)
seurat2A <- RunPCA(seurat2A, assay = "SCT", verbose = FALSE)
seurat2A <- FindNeighbors(seurat2A, reduction = "pca", dims = 1:30)
seurat2A <- RunUMAP(seurat2A, reduction = "pca", dims = 1:30)
seurat2A_DE <- FindAllMarkers(seurat2A, test.use = "MAST")
seurat2A <- FindVariableFeatures(seurat2A, assay = "SCT", features = VariableFeatures(seurat2A), 
                                 selection.method = "vst")


#_S2C

seurat2C <- SCTransform(seurat2C, assay = "Spatial", verbose = FALSE)
seurat2C <- RunPCA(seurat2C, assay = "SCT", verbose = FALSE)
seurat2C <- FindNeighbors(seurat2C, reduction = "pca", dims = 1:30)
seurat2C <- RunUMAP(seurat2C, reduction = "pca", dims = 1:30)
seurat2C_DE <- FindAllMarkers(seurat2C, test.use = "MAST")
seurat2C <- FindVariableFeatures(seurat2C, assay = "SCT", features = VariableFeatures(seurat2C), 
                                 selection.method = "vst")


#_S2D

seurat2D <- SCTransform(seurat2D, assay = "Spatial", verbose = FALSE)
seurat2D <- RunPCA(seurat2D, assay = "SCT", verbose = FALSE)
seurat2D <- FindNeighbors(seurat2D, reduction = "pca", dims = 1:30)
seurat2D <- RunUMAP(seurat2D, reduction = "pca", dims = 1:30)
seurat2D_DE <- FindAllMarkers(seurat2D, test.use = "MAST")
seurat2D <- FindVariableFeatures(seurat2D, assay = "SCT", features = VariableFeatures(seurat2D), 
                                 selection.method = "vst")


#_S3A

seurat3A <- SCTransform(seurat3A, assay = "Spatial", verbose = FALSE)
seurat3A <- RunPCA(seurat3A, assay = "SCT", verbose = FALSE)
seurat3A <- FindNeighbors(seurat3A, reduction = "pca", dims = 1:30)
seurat3A <- RunUMAP(seurat3A, reduction = "pca", dims = 1:30)
seurat3A_DE <- FindAllMarkers(seurat3A, test.use = "MAST")
seurat3A <- FindVariableFeatures(seurat3A, assay = "SCT", features = VariableFeatures(seurat3A), 
                                 selection.method = "vst")

#_S3B

seurat3B <- SCTransform(seurat3B, assay = "Spatial", verbose = FALSE)
seurat3B <- RunPCA(seurat3B, assay = "SCT", verbose = FALSE)
seurat3B <- FindNeighbors(seurat3B, reduction = "pca", dims = 1:30)
seurat3B <- RunUMAP(seurat3B, reduction = "pca", dims = 1:30)
seurat3B_DE <- FindAllMarkers(seurat3B, test.use = "MAST")
seurat3B <- FindVariableFeatures(seurat3B, assay = "SCT", features = VariableFeatures(seurat3B), 
                                 selection.method = "vst")

#_S3C

seurat3C <- SCTransform(seurat3C, assay = "Spatial", verbose = FALSE)
seurat3C <- RunPCA(seurat3C, assay = "SCT", verbose = FALSE)
seurat3C <- FindNeighbors(seurat3C, reduction = "pca", dims = 1:30)
seurat3C <- RunUMAP(seurat3C, reduction = "pca", dims = 1:30)
seurat3C_DE <- FindAllMarkers(seurat3C, test.use = "MAST")
seurat3C <- FindVariableFeatures(seurat3C, assay = "SCT", features = VariableFeatures(seurat3C), 
                                 selection.method = "vst")

#_S3D

seurat3D <- SCTransform(seurat3D, assay = "Spatial", verbose = FALSE)
seurat3D <- RunPCA(seurat3D, assay = "SCT", verbose = FALSE)
seurat3D <- FindNeighbors(seurat3D, reduction = "pca", dims = 1:30)
seurat3D <- RunUMAP(seurat3D, reduction = "pca", dims = 1:30)
seurat3D_DE <- FindAllMarkers(seurat3D, test.use = "MAST")
seurat3D <- FindVariableFeatures(seurat3D, assay = "SCT", features = VariableFeatures(seurat3D), 
                                 selection.method = "vst")


#_S4A

seurat4A <- SCTransform(seurat4A, assay = "Spatial", verbose = FALSE)
seurat4A <- RunPCA(seurat4A, assay = "SCT", verbose = FALSE)
seurat4A <- FindNeighbors(seurat4A, reduction = "pca", dims = 1:30)
seurat4A <- RunUMAP(seurat4A, reduction = "pca", dims = 1:30)
seurat4A_DE <- FindAllMarkers(seurat4A, test.use = "MAST")
seurat4A <- FindVariableFeatures(seurat4A, assay = "SCT", features = VariableFeatures(seurat4A), 
                                 selection.method = "vst")


#_S4B

seurat4B <- SCTransform(seurat4B, assay = "Spatial", verbose = FALSE)
seurat4B <- RunPCA(seurat4B, assay = "SCT", verbose = FALSE)
seurat4B <- FindNeighbors(seurat4B, reduction = "pca", dims = 1:30)
seurat4B <- RunUMAP(seurat4B, reduction = "pca", dims = 1:30)
seurat4B_DE <- FindAllMarkers(seurat4B, test.use = "MAST")
seurat4B <- FindVariableFeatures(seurat4B, assay = "SCT", features = VariableFeatures(seurat4B), 
                                 selection.method = "vst")

#_S4C

seurat4C <- SCTransform(seurat4C, assay = "Spatial", verbose = FALSE)
seurat4C <- RunPCA(seurat4C, assay = "SCT", verbose = FALSE)
seurat4C <- FindNeighbors(seurat4C, reduction = "pca", dims = 1:30)
seurat4C <- RunUMAP(seurat4C, reduction = "pca", dims = 1:30)
seurat4C_DE <- FindAllMarkers(seurat4C, test.use = "MAST")
seurat4C <- FindVariableFeatures(seurat4C, assay = "SCT", features = VariableFeatures(seurat4C), 
                                 selection.method = "vst")

#_S4D

seurat4D <- SCTransform(seurat4D, assay = "Spatial", verbose = FALSE)
seurat4D <- RunPCA(seurat4D, assay = "SCT", verbose = FALSE)
seurat4D <- FindNeighbors(seurat4D, reduction = "pca", dims = 1:30)
seurat4D <- RunUMAP(seurat4D, reduction = "pca", dims = 1:30)
seurat4D_DE <- FindAllMarkers(seurat4D, test.use = "MAST")
seurat4D <- FindVariableFeatures(seurat4D, assay = "SCT", features = VariableFeatures(seurat4D), 
                                 selection.method = "vst")


#_S5A

seurat5A <- SCTransform(seurat5A, assay = "Spatial", verbose = FALSE)
seurat5A <- RunPCA(seurat5A, assay = "SCT", verbose = FALSE)
seurat5A <- FindNeighbors(seurat5A, reduction = "pca", dims = 1:30)
seurat5A <- RunUMAP(seurat5A, reduction = "pca", dims = 1:30)
seurat5A_DE <- FindAllMarkers(seurat5A, test.use = "MAST")
seurat5A <- FindVariableFeatures(seurat5A, assay = "SCT", features = VariableFeatures(seurat5A), 
                                 selection.method = "vst")

#_S5B

seurat5B <- SCTransform(seurat5B, assay = "Spatial", verbose = FALSE)
seurat5B <- RunPCA(seurat5B, assay = "SCT", verbose = FALSE)
seurat5B <- FindNeighbors(seurat5B, reduction = "pca", dims = 1:30)
seurat5B <- RunUMAP(seurat5B, reduction = "pca", dims = 1:30)
seurat5B_DE <- FindAllMarkers(seurat5B, test.use = "MAST")
seurat5B <- FindVariableFeatures(seurat5B, assay = "SCT", features = VariableFeatures(seurat5B), 
                                 selection.method = "vst")


#_S5C

seurat5C <- SCTransform(seurat5C, assay = "Spatial", verbose = FALSE)
seurat5C <- RunPCA(seurat5C, assay = "SCT", verbose = FALSE)
seurat5C <- FindNeighbors(seurat5C, reduction = "pca", dims = 1:30)
seurat5C <- RunUMAP(seurat5C, reduction = "pca", dims = 1:30)
seurat5C_DE <- FindAllMarkers(seurat5C, test.use = "MAST")
seurat5C <- FindVariableFeatures(seurat5C, assay = "SCT", features = VariableFeatures(seurat5C), 
                                 selection.method = "vst")


#_S5D

seurat5D <- SCTransform(seurat5D, assay = "Spatial", verbose = FALSE)
seurat5D <- RunPCA(seurat5D, assay = "SCT", verbose = FALSE)
seurat5D <- FindNeighbors(seurat5D, reduction = "pca", dims = 1:30)
seurat5D <- RunUMAP(seurat5D, reduction = "pca", dims = 1:30)
seurat5D_DE <- FindAllMarkers(seurat5D, test.use = "MAST")
seurat5D <- FindVariableFeatures(seurat5D, assay = "SCT", features = VariableFeatures(seurat5D), 
                                 selection.method = "vst")

#_S6A

seurat6A <- SCTransform(seurat6A, assay = "Spatial", verbose = FALSE)
seurat6A <- RunPCA(seurat6A, assay = "SCT", verbose = FALSE)
seurat6A <- FindNeighbors(seurat6A, reduction = "pca", dims = 1:30)
seurat6A <- RunUMAP(seurat6A, reduction = "pca", dims = 1:30)
seurat6A_DE <- FindAllMarkers(seurat6A, test.use = "MAST")
seurat6A <- FindVariableFeatures(seurat6A, assay = "SCT", features = VariableFeatures(seurat6A), 
                                 selection.method = "vst")

#_S6B

seurat6B <- SCTransform(seurat6B, assay = "Spatial", verbose = FALSE)
seurat6B <- RunPCA(seurat6B, assay = "SCT", verbose = FALSE)
seurat6B <- FindNeighbors(seurat6B, reduction = "pca", dims = 1:30)
seurat6B <- RunUMAP(seurat6B, reduction = "pca", dims = 1:30)
seurat6B_DE <- FindAllMarkers(seurat6B, test.use = "MAST")
seurat6B <- FindVariableFeatures(seurat6B, assay = "SCT", features = VariableFeatures(seurat6B), 
                                 selection.method = "vst")


#_S6C

seurat6C <- SCTransform(seurat6C, assay = "Spatial", verbose = FALSE)
seurat6C <- RunPCA(seurat6C, assay = "SCT", verbose = FALSE)
seurat6C <- FindNeighbors(seurat6C, reduction = "pca", dims = 1:30)
seurat6C <- RunUMAP(seurat6C, reduction = "pca", dims = 1:30)
seurat6C_DE <- FindAllMarkers(seurat6C, test.use = "MAST")
seurat6C <- FindVariableFeatures(seurat6C, assay = "SCT", features = VariableFeatures(seurat6C), 
                                 selection.method = "vst")


#_S6D

seurat6D <- SCTransform(seurat6D, assay = "Spatial", verbose = FALSE)
seurat6D <- RunPCA(seurat6D, assay = "SCT", verbose = FALSE)
seurat6D <- FindNeighbors(seurat6D, reduction = "pca", dims = 1:30)
seurat6D <- RunUMAP(seurat6D, reduction = "pca", dims = 1:30)
seurat6D_DE <- FindAllMarkers(seurat6D, test.use = "MAST")
seurat6D <- FindVariableFeatures(seurat6D, assay = "SCT", features = VariableFeatures(seurat6D), 
                                 selection.method = "vst")


#_S7A

seurat7A <- SCTransform(seurat7A, assay = "Spatial", verbose = FALSE)
seurat7A <- RunPCA(seurat7A, assay = "SCT", verbose = FALSE)
seurat7A <- FindNeighbors(seurat7A, reduction = "pca", dims = 1:30)
seurat7A <- RunUMAP(seurat7A, reduction = "pca", dims = 1:30)
seurat7A_DE <- FindAllMarkers(seurat7A, test.use = "MAST")
seurat7A <- FindVariableFeatures(seurat7A, assay = "SCT", features = VariableFeatures(seurat7A), 
                                 selection.method = "vst")

#_S7B

seurat7B <- SCTransform(seurat7B, assay = "Spatial", verbose = FALSE)
seurat7B <- RunPCA(seurat7B, assay = "SCT", verbose = FALSE)
seurat7B <- FindNeighbors(seurat7B, reduction = "pca", dims = 1:30)
seurat7B <- RunUMAP(seurat7B, reduction = "pca", dims = 1:30)
seurat7B_DE <- FindAllMarkers(seurat7B, test.use = "MAST")
seurat7B <- FindVariableFeatures(seurat7B, assay = "SCT", features = VariableFeatures(seurat7B), 
                                 selection.method = "vst")


#_S7C

seurat7C <- SCTransform(seurat7C, assay = "Spatial", verbose = FALSE)
seurat7C <- RunPCA(seurat7C, assay = "SCT", verbose = FALSE)
seurat7C <- FindNeighbors(seurat7C, reduction = "pca", dims = 1:30)
seurat7C <- RunUMAP(seurat7C, reduction = "pca", dims = 1:30)
seurat7C_DE <- FindAllMarkers(seurat7C, test.use = "MAST")
seurat7C <- FindVariableFeatures(seurat7C, assay = "SCT", features = VariableFeatures(seurat7C), 
                                 selection.method = "vst")


#_S7D

seurat7D <- SCTransform(seurat7D, assay = "Spatial", verbose = FALSE)
seurat7D <- RunPCA(seurat7D, assay = "SCT", verbose = FALSE)
seurat7D <- FindNeighbors(seurat7D, reduction = "pca", dims = 1:30)
seurat7D <- RunUMAP(seurat7D, reduction = "pca", dims = 1:30)
seurat7D_DE <- FindAllMarkers(seurat7D, test.use = "MAST")
seurat7D <- FindVariableFeatures(seurat7D, assay = "SCT", features = VariableFeatures(seurat7D), 
                                 selection.method = "vst")

#_S8A

seurat8A <- SCTransform(seurat8A, assay = "Spatial", verbose = FALSE)
seurat8A <- RunPCA(seurat8A, assay = "SCT", verbose = FALSE)
seurat8A <- FindNeighbors(seurat8A, reduction = "pca", dims = 1:30)
seurat8A <- RunUMAP(seurat8A, reduction = "pca", dims = 1:30)
seurat8A_DE <- FindAllMarkers(seurat8A, test.use = "MAST")
seurat8A <- FindVariableFeatures(seurat8A, assay = "SCT", features = VariableFeatures(seurat8A), 
                                 selection.method = "vst")

#_S8B

seurat8B <- SCTransform(seurat8B, assay = "Spatial", verbose = FALSE)
seurat8B <- RunPCA(seurat8B, assay = "SCT", verbose = FALSE)
seurat8B <- FindNeighbors(seurat8B, reduction = "pca", dims = 1:30)
seurat8B <- RunUMAP(seurat8B, reduction = "pca", dims = 1:30)
seurat8B_DE <- FindAllMarkers(seurat8B, test.use = "MAST")
seurat8B <- FindVariableFeatures(seurat8B, assay = "SCT", features = VariableFeatures(seurat8B), 
                                 selection.method = "vst")

#_S8C

seurat8C <- SCTransform(seurat8C, assay = "Spatial", verbose = FALSE)
seurat8C <- RunPCA(seurat8C, assay = "SCT", verbose = FALSE)
seurat8C <- FindNeighbors(seurat8C, reduction = "pca", dims = 1:30)
seurat8C <- RunUMAP(seurat8C, reduction = "pca", dims = 1:30)
seurat8C_DE <- FindAllMarkers(seurat8C, test.use = "MAST")
seurat8C <- FindVariableFeatures(seurat8C, assay = "SCT", features = VariableFeatures(seurat8C), 
                                 selection.method = "vst")




bsmerge <- merge(x = seurat1A, c(seurat1B))
bsmerge <- merge(x = bsmerge, c(seurat2A))
bsmerge <- merge(x = bsmerge, c(seurat2C))
bsmerge <- merge(x = bsmerge, c(seurat2D))
bsmerge <- merge(x = bsmerge, c(seurat3A))
bsmerge <- merge(x = bsmerge, c(seurat3B))
bsmerge <- merge(x = bsmerge, c(seurat3C))
bsmerge <- merge(x = bsmerge, c(seurat3D))
bsmerge <- merge(x = bsmerge, c(seurat4A))
bsmerge <- merge(x = bsmerge, c(seurat4B))
bsmerge <- merge(x = bsmerge, c(seurat4C))
bsmerge <- merge(x = bsmerge, c(seurat4D))
bsmerge <- merge(x = bsmerge, c(seurat5A))
bsmerge <- merge(x = bsmerge, c(seurat5B))
bsmerge <- merge(x = bsmerge, c(seurat5C))
bsmerge <- merge(x = bsmerge, c(seurat5D))
bsmerge <- merge(x = bsmerge, c(seurat6A))
bsmerge <- merge(x = bsmerge, c(seurat6B))
bsmerge <- merge(x = bsmerge, c(seurat6C))
bsmerge <- merge(x = bsmerge, c(seurat6D))
bsmerge <- merge(x = bsmerge, c(seurat7A))
bsmerge <- merge(x = bsmerge, c(seurat7B))
bsmerge <- merge(x = bsmerge, c(seurat7C))
bsmerge <- merge(x = bsmerge, c(seurat7D))
bsmerge <- merge(x = bsmerge, c(seurat8A))
bsmerge <- merge(x = bsmerge, c(seurat8B))
bsmerge <- merge(x = bsmerge, c(seurat8C)
