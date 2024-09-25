

```{r}
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Mmusculus.v79)
#library(EnsDb.Mmusculus.v75)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(cowplot)
library(patchwork) #latest version is required!
library(TFBSTools)
library(JASPAR2020)
library(gsfisher)
library(EnhancedVolcano)
library(monocle)
setwd("/rds/projects/c/croftap-stia-atac/CM_multiome/STIA_andATAC/")
options(bitmapType='cairo')
library(CellChat)
#load()
```

```{r}
#setwd("/rds/projects/c/croftap-celldive01/Panel 2/segmentation")
results <- dir("./", pattern = "*txt", 
    full.names = TRUE)

results <- sub('./', '', results)
results<-paste("/rds/projects/c/croftap-celldive01/Panel_2/segmentation (adult)", results, sep="")

data = list()
for (i in 1:length(results)) {
    data[[i]] <- read.delim(results[[i]])
      
}
library("stringr")   
names <- dir("./", pattern = "*txt", 
    full.names = TRUE)
names <- sub('.//', '', names)


library(dplyr)

for (i in 1:length(results)) {
  data[[i]] <- data[[i]] %>% select(contains('.mean'))
  data[[i]] <- data[[i]] %>% select(contains('cell'))
  rownames(data[[i]])<-seq_along(data[[i]][,1])
}
  
names(data) <- names

  
data_x_y = list()
for (i in 1:length(results)) {
    data_x_y[[i]] <- read.delim(results[[i]])
    rownames(data_x_y[[i]])<-seq_along(data_x_y[[i]][,1])
    data_x_y[[i]] <- data_x_y[[i]] %>% select(c('Centroid.X.µm', 'Centroid.Y.µm'))
}

names(data_x_y) <- names

library(sctransform)

for (i in 1:length(results)) {
data[[i]]<-as.data.frame(t(data[[i]]))
data[[i]] <- data[[i]] %>% replace(is.na(.), 0)
data[[i]]<-CreateSeuratObject(data[[i]])
data[[i]]<-SCTransform(data[[i]])
data[[i]]<-RunPCA(data[[i]])
data[[i]]<-RunUMAP(data[[i]], dims = 1:20)
}



```

```{r}
conditons <- c("eRA", "eRA","eRA", "EstRA", "EstRA","EstRA","OA", "OA", "OA", "Res", "Res", "Res" )
for (i in 1:length(results)) {
data[[i]]$orig.ident<-names[[i]];
data[[i]]$condition<-conditons[[i]]
}

all <- merge(x = data[[1]], y = data[-1])
int_featuers <- SelectIntegrationFeatures(data)
VariableFeatures(all) <-int_featuers
#all<-FindVariableFeatures(all)
all<-ScaleData(all)
all<-RunPCA(all)
all<-RunUMAP(all, dims = 1:20)
DimPlot(all, group.by="orig.ident")
all <- RunHarmony.Seurat_CM(all, group.by.vars = "orig.ident")
all <- RunUMAP(all, reduction="harmony", dims=1:20)
DimPlot(all, group.by = "orig.ident")


all<-FindNeighbors(all, dims = 1:20, reduction = "harmony")
all<-FindClusters(all, resolution = c(0.3, 0.4, 0.5, 0.6))
all<-FindClusters(all, resolution = c(0.05, 0.1, 0.2))
all<-FindClusters(all, resolution = c(0.075))

DimPlot(all, group.by = "SCT_snn_res.0.4")


Idents(all)<-'SCT_snn_res.0.05'
DotPlot(all, features =rownames(all)) +RotatedAxis()

DefaultAssay(all) <- 'SCT'
Idents(all)<-'SCT_snn_res.0.1'
DotPlot(all, features =rownames(all)) +RotatedAxis()
Idents(all)<-'SCT_snn_res.0.2'
DotPlot(all, features =rownames(all)) +RotatedAxis()
Idents(all)<-'SCT_snn_res.0.3'
DotPlot(all, features =rownames(all)) +RotatedAxis()
Idents(all)<-'SCT_snn_res.0.05'
DotPlot(all, features =rownames(all)) +RotatedAxis()

#all<-FindClusters(all, resolution = c(0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3))


#res 0.1 used

Idents(merging_vis) <- 'named_clusters'
VlnPlot(merging_vis, features = "POSTN")
VlnPlot(merging_vis, features = "SPARC")
DimPlot(merging_vis)
ncol(merging_vis)

all$named_niches <- all@meta.data[["SCT_snn_res.0.1"]]
Idents(all) <- 'named_niches'
levels(all)



#TODO-split out cluster 2 into aggr (CD20, CD3 high), thenrename all and check
current.sample.ids <- c("0"  ,"1",  "10" ,"11", "12", "2", "3" , "4",  "5" , "6",  "7" , "8" , "9")
new.sample.ids <- c("Fibrotic_lymph"  ,"LL",  "?" ,"Lymphocytic1", "Interstitial", "2", "vascular" , "Interstitial",  "Interstitial" , "vascular",  "Lymphocytic2" , "??" , "Fibrotic_Col")

all@meta.data[["named_niches"]] <- plyr::mapvalues(x = all@meta.data[["named_niches"]], from = current.sample.ids, to = new.sample.ids)

Idents(all) <- 'named_niches'

DotPlot(all, features =rownames(all)) +RotatedAxis()


```

```{r}
#plot one section and zoom in

#best section#
names <- names(data_x_y)
xy_list <- list()
for (i in 1:length(data_x_y)) {
xy_list[[i]]<-data_x_y[[names[[i]]]];
}


all_meta <- all@meta.data
s_obj_meta <- list()
ggplot_ls <- list()
for (i in 1:length(names)) {
s_obj_meta[[i]] <- all_meta[all_meta$orig.ident == names[[i]],];
xy_list[[i]]$named <- s_obj_meta[[i]]$named_niches;
cols <- ArchR::paletteDiscrete(s_obj_meta[[i]][, "named_niches"]);
ggplot_ls[[i]] <- ggplot(xy_list[[i]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.1)+theme_classic()+scale_color_manual(values = c(cols));
print(ggplot_ls[[i]])
}


for (i in 1:length(names)) {
print(ggplot_ls[[i]])
}

for (i in 1:length(names)) {
print(table(xy_list[[i]]$named))
}
  
  
cols2=c("grey", "grey","red","grey","grey","grey","grey","green","grey")

ggplot(xy_list[[1]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.5)+theme_classic()+scale_color_manual(values = c(cols))+xlim(3000,4500)+ ylim(2000,3500)

```





```{r}

Idents(all) <- 'orig.ident'
eRA_BX020 <- subset(all, idents="eRA_BX020.txt")

Idents(eRA_BX020) <- 'named_niches'
eRA_BX020_LL <- subset(eRA_BX020, idents="LL")

eRA_BX020_LL<-ScaleData(eRA_BX020_LL)
eRA_BX020_LL<-RunPCA(eRA_BX020_LL)
eRA_BX020_LL<-RunUMAP(eRA_BX020_LL, dims = 1:20)
DimPlot(eRA_BX020_LL, group.by="named_niches")

eRA_BX020_LL<-FindNeighbors(eRA_BX020_LL, dims = 1:20)
eRA_BX020_LL<-FindClusters(eRA_BX020_LL, resolution = c(0.3, 0.4))

eRA_BX020_LL<-FindClusters(eRA_BX020_LL, resolution = c(0.1, 0.05))
eRA_BX020_LL<-FindClusters(eRA_BX020_LL, resolution = c(0.2, 0.5))

DimPlot(eRA_BX020_LL, group.by="SCT_snn_res.0.3")

Idents(eRA_BX020_LL)<-'SCT_snn_res.0.3'
DotPlot(eRA_BX020_LL, features =rownames(eRA_BX020_LL)) +RotatedAxis()


eRA_BX020_LL$named_niches <- eRA_BX020_LL@meta.data[["SCT_snn_res.0.3"]]
Idents(eRA_BX020_LL) <- 'named_niches'
levels(eRA_BX020_LL)
current.sample.ids <- c( "0",  "1" , "2" , "3" , "4" , "5" , "6" , "7",  "8" , "9" , "10")
new.sample.ids <- c( "Fibrotic_lymph",  "LL" , "Lymphocytic1" , "vascular" , "Fibrotic_lymph" , "LL" , "LL" , "LL",  "Fibrotic" , "LL" , "vascular")
eRA_BX020_LL@meta.data[["named_niches"]] <- plyr::mapvalues(x = eRA_BX020_LL@meta.data[["named_niches"]], from = current.sample.ids, to = new.sample.ids)

eRA_BX020_LL_meta <- eRA_BX020_LL@meta.data
all_meta <- all@meta.data
all_meta_f <- all_meta[!rownames(all_meta) %in% rownames(eRA_BX020_LL_meta),]
all_meta_f <- rbind(all_meta_f, eRA_BX020_LL_meta)
all <- AddMetaData(all, all_meta_f)


Idents(all) <- 'named_niches'
levels(all)
current.sample.ids <- c( "Fibrotic_lymph" ,"LL"     ,        "?"      ,        "Lymphocytic1"  , "Interstitial" ,  "2"     ,         "vascular"  ,     "Lymphocytic2",   "??",  "Fibrotic_Col" ,  "Fibrotic" )

new.sample.ids <- c(  "Fibrotic", "LL"    ,         "?"         ,     "Lymphocytic1" ,  "Adipose" ,  "2"  ,     "vascular",       "Lymphocytic2",   "??" ,    "Fibrotic", "RBC_contam")
all@meta.data[["named_niches"]] <- plyr::mapvalues(x = all@meta.data[["named_niches"]], from = current.sample.ids, to = new.sample.ids)


all_meta <- all@meta.data
s_obj_meta <- list()
ggplot_ls <- list()
for (i in 1:length(names)) {
s_obj_meta[[i]] <- all_meta[all_meta$orig.ident == names[[i]],];
xy_list[[i]]$named <- s_obj_meta[[i]]$named_niches;
cols <- ArchR::paletteDiscrete(s_obj_meta[[i]][, "named_niches"]);
ggplot_ls[[i]] <- ggplot(xy_list[[i]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.1)+theme_classic()+scale_color_manual(values = c(cols));
print(ggplot_ls[[i]])
}

cols <- ArchR::paletteDiscrete(all@meta.data[, "named_niches"]);


cols2=c("grey", "grey","red","grey","grey","green","grey","grey","grey","yellow","grey")
ggplot(xy_list[[6]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.3)+theme_classic()+scale_color_manual(values = c(cols))+xlim(5000,6500)+ylim(4000,6000)


```

```{r}

#check these two are OK for lymphocytic1 and that if there are any other aggr?

Idents(all) <- 'orig.ident'
levels(all)

rename <- subset(all, idents=c("EstRA_BX240.txt", "OA_JRP127.txt"))

Idents(rename) <- 'named_niches'
levels(rename)
current.sample.ids <- c( "Fibrotic"   ,  "LL"   ,        "?"      ,      "Lymphocytic1", "Adipose"  ,    "2"       ,     "vascular"  ,   "Lymphocytic2", "??" ,          "RBC_contam"   )

new.sample.ids <- c("Fibrotic"   ,  "LL"   ,        "?"      ,      "Lymphocytic1", "Adipose"  ,    "Lymphocytic1"       ,     "vascular"  ,   "Lymphocytic2", "??" ,          "RBC_contam"   )
rename@meta.data[["named_niches"]] <- plyr::mapvalues(x = rename@meta.data[["named_niches"]], from = current.sample.ids, to = new.sample.ids)

rename_meta <- rename@meta.data
all_meta <- all@meta.data
all_meta_f <- all_meta[!rownames(all_meta) %in% rownames(rename_meta),]
all_meta_f <- rbind(all_meta_f, rename_meta)
all <- AddMetaData(all, all_meta_f)


all_meta <- all@meta.data
s_obj_meta <- list()
ggplot_ls <- list()
for (i in 1:length(names)) {
s_obj_meta[[i]] <- all_meta[all_meta$orig.ident == names[[i]],];
xy_list[[i]]$named <- s_obj_meta[[i]]$named_niches;
cols <- ArchR::paletteDiscrete(s_obj_meta[[i]][, "named_niches"]);
ggplot_ls[[i]] <- ggplot(xy_list[[i]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.1)+theme_classic()+scale_color_manual(values = c(cols));
print(ggplot_ls[[i]])
}


ggplot(xy_list[[1]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=1)+theme_classic()+scale_color_manual(values = c(cols))+xlim(100, 1500)+ylim(2500,3500)

#figure 3k
ggplot(xy_list[[6]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.3)+theme_classic()+scale_color_manual(values = c(cols))+xlim(250, 1500)+ylim(2300,4250)


ggplot(xy_list[[6]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.3)+theme_classic()+scale_color_manual(values = c("Lymphocytic1"="#0C5BB0FF", "2"="#EE0011FF", "Lymphocytic2"="#15983DFF", "Fibrotic"="#EC579AFF", "Adipose"="#FA6B09FF", "vascular"="#149BEDFF", "LL"= "#A1C720FF", "??"="#C06CAB"))+xlim(250, 1500)+ylim(2300,4250)


c("Lymphocytic1"="#0C5BB0FF", "2"="#EE0011FF", "Lymphocytic2"="#15983DFF", "Fibrotic"="#EC579AFF", "Adipose"="#FA6B09FF", "vascular"="#149BEDFF", "LL"= "#A1C720FF", "??"="#C06CAB")

ggplot(xy_list[[6]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.3)+theme_classic()+scale_color_manual(values = c(cols))+xlim(4750, 6500)+ylim(3000,5000)

ggplot(xy_list[[6]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.3)+theme_classic()+scale_color_manual(values = c("Lymphocytic1"="#0C5BB0FF", "2"="#EE0011FF", "Lymphocytic2"="#15983DFF", "Fibrotic"="#EC579AFF", "Adipose"="#FA6B09FF", "vascular"="#149BEDFF", "LL"= "#A1C720FF", "??"="#C06CAB"))+xlim(4750, 6500)+ylim(3000,5000)

```


```{r}
cols <- ArchR::paletteDiscrete(s_obj_meta[[i]][, "named_niches"])
cols2=c("grey","grey","grey","grey","grey","grey", "red")
 ggplot(xy_list[[6]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.1)+theme_classic()+scale_color_manual(values = c(cols))+xlim(3000, 7300)+ylim(1200,5500)
 
 
 ggplot(xy_list[[6]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.3)+theme_classic()+scale_color_manual(values = c(cols))
 
 
 ggplot(xy_list[[6]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.2)+theme_classic()+scale_color_manual(values = c(cols))
 
 #figure 
 palette <- c("2"="#", "??"="#", "Lymphocytic2"="#", "Fibrotic"="#D51F26", "Adipose"="#89288F", "vascular"="#", "LL"= "#")
 
 

   
   

 ggplot(xy_list[[6]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=1)+theme_classic()+scale_color_manual(values =  c("Lymphocytic1"="#0C5BB0FF", "2"="#EE0011FF", "Lymphocytic2"="#15983DFF", "Fibrotic"="#EC579AFF", "Adipose"="#FA6B09FF", "vascular"="#149BEDFF", "LL"= "#A1C720FF", "??"="#C06CAB"))+xlim(4900, 5750)+ylim(4500,5800)
 
 
 
 #lymphocyic 
 ggplot(xy_list[[6]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=1)+theme_classic()+scale_color_manual(values = c(cols))+xlim(5400, 5750)+ylim(4800,5600)
 
```


```{r}

DimPlot(merging_vis, group.by = "orig.ident")

DimPlot(merging_vis, split.by = "diseasestate")

```

```{r}

all$named_niches_final <- all$named_niches
Idents(all) <- 'named_niches_final'
levels(all)

current.sample.ids <- c( "Fibrotic"  ,   "LL"      ,     "?"    ,        "Lymphocytic1", "Adipose"   ,   "2"      ,      "vascular"  ,   "Lymphocytic2",
  "??"  ,         "RBC_contam"    )

new.sample.ids <- c("Fibrotic"  ,   "LL"      ,     "low_quality"    ,        "Lymphocytic1", "Adipose"   ,   "low_quality"      ,      "vascular"  ,   "Lymphocytic2",
  "low_quality"  ,         "RBC_contam"    )
all@meta.data[["named_niches_final"]] <- plyr::mapvalues(x = all@meta.data[["named_niches_final"]], from = current.sample.ids, to = new.sample.ids)


all_meta <- all@meta.data
s_obj_meta <- list()
ggplot_ls <- list()
for (i in 1:length(names)) {
s_obj_meta[[i]] <- all_meta[all_meta$orig.ident == names[[i]],];
xy_list[[i]]$named <- s_obj_meta[[i]]$named_niches_final;
cols <- ArchR::paletteDiscrete(s_obj_meta[[i]][, "named_niches_final"]);
ggplot_ls[[i]] <- ggplot(xy_list[[i]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.1)+theme_classic()+scale_color_manual(values = c(cols));
print(ggplot_ls[[i]])
}


```
```{r}
#proximity enrichemnt analysis

#one cell type
coloc_one_type = function(index_type, adj, y, nperm = 100, max_dist=30, compartments=NULL, verbose=TRUE) {
    if (verbose) message(index_type)
    types = unique(y)
    i_index = which(y == index_type)
    i_shuffle = setdiff(seq_len(length(y)), i_index)
    X = adj[i_index, ] %*% Matrix::sparse.model.matrix(~0+y) %>% as.matrix()
    colnames(X) = gsub('^y', '', colnames(X))
    freq = (colSums(X) / nrow(X))[types]
    freq_perm = map(seq_len(nperm), function(i) {
        set.seed(i)
        yperm = y
        if (is.null(compartments)) {
            yperm[i_shuffle] = sample(y[i_shuffle])
        } else {
            ## shuffle inside compartments, to preserve total composition within compartment
            .x = split(i_shuffle, compartments[i_shuffle]) %>%
                map(function(.i) {
                    ## CAUTION: if .i is a single number, sample will interpret it as 1:.i
                    if (length(.i) == 1) {
                        message('No shuffling is taking place, check code')
                        res = .i
                    } else {
                        res = sample(.i) ## shuffle non-index cells inside hub
                    }
                    names(res) = .i
                    return(res)
                }) %>%
                reduce(c)
            yperm[as.integer(names(.x))] <- y[.x]
        }
        X = adj[i_index, ] %*% Matrix::sparse.model.matrix(~0+yperm) %>% as.matrix() #%>% prop.table(1)
        colnames(X) = gsub('^yperm', '', colnames(X))
        (colSums(X) / nrow(X))[types]
    }) %>%
        purrr::reduce(rbind2)
    stats = tibble(
        type = types,
        freq,
        zscore = (freq - apply(freq_perm, 2, mean)) / apply(freq_perm, 2, sd),
        pval = exp(log(2) + (pnorm(-abs(zscore), log.p = TRUE, lower.tail = TRUE))), ## one-tailed
        fdr = p.adjust(pval)
    ) %>%
        cbind(dplyr::rename(data.frame(t(apply(freq_perm, 2, quantile, c(.025, .975)))), q025 = `X2.5.`, q975 = `X97.5.`)) %>% ## 95% CI
        subset(type != index_type) %>%
        dplyr::mutate(index_type = index_type) %>%
        dplyr::select(index_type, type, everything()) %>%
        arrange(fdr)
    return(stats)
}H



#loop through all cell types
coloc_all_types = function(index_types, coords, y, nperm = 100, nsteps=1, max_dist=30, compartments=NULL, parallel=TRUE, verbose=TRUE) {
    if (parallel & length(index_types) > 1) {
        plan(multicore)
    } else {
        plan(sequential)
    }
    ## Define neighbors
    ## NOTE: max_dist only refers to directly adjacent neighbors
    adj = spatula::getSpatialNeighbors(coords, return_weights = TRUE)
    adj@x[adj@x > max_dist] = 0
    adj = Matrix::drop0(adj)
    adj@x = rep(1, length(adj@x))
    ## If nsteps>1, consider not only your adjacent neighbors
    ##   but also your neighbor’s neighbors etc.
    if (nsteps > 1) {
        adj = adj + Matrix::Diagonal(n = nrow(adj)) ## add self
        for (iter in seq_len(nsteps - 1)) {
            adj = adj %*% adj
        }
        ## Ignore weights. Only care if cell is a neighbor or not
        adj@x = rep(1, length(adj@x))
        ## Remove self as neighbor
        adj = adj - Matrix::Diagonal(n = nrow(adj))
        adj = Matrix::drop0(adj)
    }
    index_types %>%
        future_map(coloc_one_type, adj, y, nperm, max_dist, compartments, verbose, .options = furrr::furrr_options(seed = 1)) %>%
        rbindlist()  %>%
        identity
}


library(devtools)
#install_github("korsunskylab/spatula")
library(spatula)
library(furrr)

xy_1 <- xy_list[[1]] %>%  as.data.frame() %>% select("Centroid.X.µm", "Centroid.Y.µm", "named")


#add max of 10000
off_sets <- (rep(1:((length(xy_list)-1)))*10000)
off_sets <- c(0, off_sets)
xy_list_coloc <- xy_list

for (i in 1:length(xy_list_coloc)){
  xy_list_coloc[[i]]$Centroid.X.µm <- xy_list_coloc[[i]]$Centroid.X.µm + off_sets[[i]]
}

library(data.table)
xy_list_coloc_list <- rbindlist(xy_list_coloc, use.names=TRUE)

# repeat but check named column is coreect
coloc_res_coarse<-tryCatch({
    coloc_res_coarse = coloc_all_types(
        index_type = unique(xy_list_coloc_list$named),
        coords = xy_list_coloc_list[, c("Centroid.X.µm", "Centroid.Y.µm")],
        y = xy_list_coloc_list$named,
        compartments = NULL,
        max_dist = 40,
        nperm = 1000,
        parallel = TRUE
    )}, error = function(e) {
    message('Failed run on: ', xy_list[[i]])
    return(NULL)
    }
)



names(xy_list_coloc) <- names(data_x_y)




library(dplyr)
library(data.table)
library(tidyverse)


plt_df<-coloc_res_coarse %>%
    subset(pval < 0.05) %>%
    dplyr::select(index_type, type, zscore) %>%
    spread(type, zscore, fill = 0) %>%
    column_to_rownames('index_type') %>%
    as.matrix

plt_df[-c(3,8),-c(3,8)] %>%
    Heatmap(col = circlize::colorRamp2(c(-20, 0, 20), c('blue', 'white', 'red')), border = T, cluster_columns = F, cluster_rows = F)

plt_df_use <- plt_df[-c(3,4,8),-c(3,4,8)]

plt_df_use2 <- plt_df_use+t(plt_df_use)


plt_df_use2 %>%
    Heatmap(col = circlize::colorRamp2(c(-20, 0, 20), c('blue', 'white', 'red')), border = T, cluster_columns = F, cluster_rows = F)

```



```{r}
Idents(all) <- 'named_niches'
rownames(all)
levels(all)

VlnPlot(all, idents = c("Fibrotic_lymph", "LL"  ,  "Lymphocytic1" ,  "Interstitial",   "vascular"   ,    "Lymphocytic2", "Fibrotic_Col"), features = c("Cell..CD20.mean"   ,    "Cell..CD3.mean"   ,     "Cell..CLU.mean"    ,    "Cell..COL1.mean"  ,    
"Cell..COL6A1.mean" ,    "Cell..COMP.mean"  ,    "Cell..MMP3.mean"   ,    "Cell..PDPN.mean"    ,   "Cell..POSTN.mean"   ), stack=T, pt.size = 0)+NoLegend()

DotPlot(all, features = c(  "Cell..COL1.mean"  ,  "Cell..COMP.mean"  ,   "Cell..POSTN.mean"   ))+NoLegend()


DefaultAssay(all) <- "SCT"
dotplot<-DotPlot(all, features = c("Cell..COL1.mean", "Cell..COL6A1.mean", "Cell..POSTN.mean",  "Cell..MMP3.mean", "Cell..PDPN.mean","Cell..CD68.mean", "Cell..CD3.mean", "Cell..CD20.mean", "Cell..FABP4.mean",  "Cell..SMA.mean" ,  "Cell..CD31.mean" ,  "Cell..CD146.mean",  "Cell..CD138.mean"), idents = c("Fibrotic" ,    "LL"    ,          "Lymphocytic1", "Adipose"  ,           "vascular"     ))+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

dotplot<-dotplot$data

dotplot<-dotplot %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

library(ComplexHeatmap)

Heatmap(dotplot, cluster_columns = F, cluster_rows = F)





Idents(merging) <- 'named_clusters'
levels(merging)
DefaultAssay(merging) <- "SCT"
levels(merging) <- c("COMP+ Fibroblast Niche",  "Lining Layer Cells" ,"T Cell Rich Niche", "APOD+ GAS5+ FABP4+","Vascular Niche","B Cell Rich Niche" , "Erythrocytes")

dotplot<-DotPlot(merging,features=c("COL1A1", "COL6A1", "POSTN", "MMP3", "PDPN", "CD68", "MS4A1", "FABP4", "ACTA2", "PECAM1", "MCAM", "SDC1"), idents= levels(merging)[-7])+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

dotplot<-dotplot$data

dotplot<-dotplot %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

library(ComplexHeatmap)

Heatmap(dotplot, cluster_columns = F, cluster_rows = F)



DefaultAssay(all) <- "SCT"
DotPlot(all, features = rownames(all), idents = c("Fibrotic" ,    "LL"    ,          "Lymphocytic1", "Adipose"  ,           "vascular"   ,  "Lymphocytic2"  ))+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

```

```{r}
Idents(all) <- 'orig.ident'
levels(all)
all$group <- all$orig.ident
current.sample.ids <- c( "eRA_BX020.txt" ,  "eRA_BX031.txt" ,  "eRA_BX115.txt"  , "EstRA_BX127.txt", "EstRA_BX195.txt", "EstRA_BX240.txt", "OA_JRP127.txt" ,  "OA_JRP146.txt" , "OA_JRP149.txt" ,  "Res_BX026.txt"  , "Res_BX028.txt",   "Res_BX202.txt"    )

new.sample.ids <- c( "early" ,  "early" ,  "early"  , "established", "established", "established", "OA" ,  "OA" , "OA" ,  "resolving"  , "resolving",   "resolving"   )
all@meta.data[["group"]] <- plyr::mapvalues(x = all@meta.data[["group"]], from = current.sample.ids, to = new.sample.ids)

pt <- table(all$named_niches, all$group)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
pt <- pt/rowSums(pt)

cols <- as.data.frame(ArchR::paletteDiscrete(all@meta.data[, "named_niches"]))
colnames(cols) <- "colors"

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1))  + 
        geom_bar(stat = 'identity', position = position_fill())+
        theme( axis.ticks.y = element_blank()) +
        coord_flip()+  scale_fill_manual(values = cols$colors)+
        theme(legend.position = 'bottom', panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank()) + 
        labs(fill = 'Tissue-Defined Cluster', y = 'Cluster Frequency in Sample') + 
        guides(fill = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 2)), alpha = FALSE)

```
```{r}

pt <- table(merging_vis$named_clusters, merging_vis$diseasestate)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

cols <- as.data.frame(ArchR::paletteDiscrete(all@meta.data[, "named_niches"]))
colnames(cols) <- "colors"

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1))  + 
        geom_bar(stat = 'identity', position = position_fill())+
        theme( axis.ticks.y = element_blank()) +
        coord_flip()+  scale_fill_manual(values = cols$colors)+
        theme(legend.position = 'bottom', panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank()) + 
        labs(fill = 'Tissue-Defined Cluster', y = 'Cluster Frequency in Sample') + 
        guides(fill = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 2)), alpha = FALSE)
```
```{r}

all_meta <- all@meta.data
all_metaC_vascular <- subset(all_meta, named_niches == "vascular")
numbers <- 1:12
coord <- xy_list
coord_vasc <- list()
for ( i in 1:length(xy_list)){
  coord[[i]]$cells <- rownames(coord[[i]]);
  coord[[i]]$cells <- paste(coord[[i]]$cells, numbers[[i]], sep="_");
coord_vasc[[i]] <- coord[[i]][(coord[[i]]$cells %in% rownames(all_metaC_vascular)),];
coord[[i]]$dist1 <- NA;
coord[[i]]$dist.x1 <- NA;
coord[[i]]$dist.y1 <- NA;
}


library(sp)


for (i in 1:length(coord)){
for ( spot in 1:nrow(coord[[i]])){
  dists <- spDistsN1(as.matrix(coord_vasc[[i]][, 2:1]), pt=as.numeric(coord[[i]][spot, c( "Centroid.X.µm", "Centroid.Y.µm")]))
  coord[[i]]$dist1[spot] <- min(dists)
  coord[[i]]$dist.x1[spot] <- coord_vasc[[i]][which(dists == min(dists))[1], 1]
  coord[[i]]$dist.y1[spot] <- coord_vasc[[i]][which(dists == min(dists))[1], 2]
}
}



```


```{r}

library(data.table)
coord_norm <- coord
for (i in 1:length(coord)){
coord_norm[[i]]$dist1 <- as.numeric(coord[[i]]$dist1) / as.numeric(max(coord[[i]]$dist1))
}


rownames(all)
DotPlot(all, features = c("Cell..COL1A1.mean", "Cell..COL3A1.mean", "Cell..COL6A1.mean"))

DefaultAssay(all) <- 'RNA'
all <- AddModuleScore(all, features = list("Cell..COL1.mean", "Cell..COL3A1.mean", "Cell..COL6A1.mean"), name="perivascualr_mod", ctrl=1)

expression_genes <- all@assays[["SCT"]]@data[c("Cell..COL1.mean", "Cell..COL4A1.mean","Cell..COL6A1.mean"),] %>% t() %>% as.data.frame()


library(data.table)

coord_norm <- rbindlist(coord_norm)


df <- data.frame(Dist=coord_norm$dist1, 
				Score=as.numeric(expression_genes$Cell..COL1.mean) / as.numeric(max(expression_genes$Cell..COL1.mean)), #scaling cell type signature to max here 
				CellType="ALL")
df$cluster <- all$named_niches

df %>% filter(cluster == "Fibrotic") %>%   ggplot(aes(Dist, Score, color=cluster, lty=cluster)) + geom_smooth(alpha=.1) + labs(x="Distance (towards lumen)", y="Cell Type Signal") + theme_light(base_size = 16)+theme_classic()+ggtitle("COL1")


expression_genes$mean <- rowMeans(expression_genes)


df2 <- data.frame(Dist=coord_norm$dist1, 
				Score=as.numeric(expression_genes$mean) / as.numeric(max(expression_genes$mean)), #scaling cell type signature to max here 
				CellType="ALL")
df2$cluster <- all$named_niches

df2 %>% filter(cluster == "Fibrotic") %>%   ggplot(aes(Dist, Score, color=cluster, lty=cluster)) + geom_smooth(alpha=.1) + labs(x="Distance (towards lumen)", y="Cell Type Signal") + theme_light(base_size = 16)+theme_classic()+ggtitle("COL1")
```
```{r}
coord_norm_to_vasc <- coord_norm

coord <- rbindlist(coord)
qunat <- quantile(coord$dist1,0.75)
coord_norm_to_vasc_perivasc <-  coord %>% filter(dist1 > qunat)




coord_norm_to_vasc$Annotations_new<-'far'
coord_norm_to_vasc$Annotations_new[coord_norm_to_vasc$dist1 > -0.1 & coord_norm_to_vasc$dist1 < 0.00999] <- "Vasc"
coord_norm_to_vasc$Annotations_new[coord_norm_to_vasc$dist1 > 0.002 & coord_norm_to_vasc$dist1 < 0.05] <- "periVasc"
coord_norm_to_vasc$Annotations_new[coord_norm_to_vasc$dist1 > 0.05 & coord_norm_to_vasc$dist1 < 0.3] <- "int1"
coord_norm_to_vasc$Annotations_new[coord_norm_to_vasc$dist1 > 0.3 ]<-'far'
coord_norm_to_vasc_f <- coord_norm_to_vasc %>% select(dist1, Annotations_new) %>%  as.data.frame()
rownames(coord_norm_to_vasc_f) <- coord_norm_to_vasc$cells

coord_norm_to_vasc_f$dist1 <- NULL
all <- AddMetaData(all, metadata=coord_norm_to_vasc_f)





all_meta <- all@meta.data
s_obj_meta <- list()
ggplot_ls <- list()
for (i in 1:length(names)) {
s_obj_meta[[i]] <- all_meta[all_meta$orig.ident == names[[i]],];
xy_list[[i]]$Annotations_new <- s_obj_meta[[i]]$Annotations_new;
cols <- ArchR::paletteDiscrete(s_obj_meta[[i]][, "Annotations_new"]);
ggplot_ls[[i]] <- ggplot(xy_list[[i]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=Annotations_new)) +
    geom_point(size=0.1)+theme_classic()+scale_color_manual(values = c(cols));
print(ggplot_ls[[i]])
}


all_meta <- all@meta.data
s_obj_meta <- list()
ggplot_ls <- list()
for (i in 1:length(names)) {
s_obj_meta[[i]] <- all_meta[all_meta$orig.ident == names[[i]],];
xy_list[[i]]$named_niches <- s_obj_meta[[i]]$named_niches;
cols <- ArchR::paletteDiscrete(s_obj_meta[[i]][, "named_niches"]);
ggplot_ls[[i]] <- ggplot(xy_list[[i]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named_niches)) +
    geom_point(size=0.1)+theme_classic()+scale_color_manual(values = c(cols));
print(ggplot_ls[[i]])
}




```


```{r}
all_meta <- all@meta.data
all_metaC_lymph <- subset(all_meta, named_niches == "Lymphocytic1")
numbers <- 1:12
coord_lymph_main <- xy_list
coord_lymph <- list()
for ( i in 1:length(xy_list)){
  coord_lymph_main[[i]]$cells <- rownames(coord_lymph_main[[i]]);
  coord_lymph_main[[i]]$cells <- paste(coord_lymph_main[[i]]$cells, numbers[[i]], sep="_");
coord_lymph[[i]] <- coord_lymph_main[[i]][(coord_lymph_main[[i]]$cells %in% rownames(all_metaC_lymph)),];
coord_lymph_main[[i]]$dist1 <- NA;
coord_lymph_main[[i]]$dist.x1 <- NA;
coord_lymph_main[[i]]$dist.y1 <- NA;
}


library(sp)


for (i in 1:length(coord_lymph_main)){
for ( spot in 1:nrow(coord_lymph_main[[i]])){
  dists <- spDistsN1(as.matrix(coord_lymph_main[[i]][, 2:1]), pt=as.numeric(coord_lymph_main[[i]][spot, c( "Centroid.X.µm", "Centroid.Y.µm")]))
  coord_lymph_main[[i]]$dist1[spot] <- min(dists)
  coord_lymph_main[[i]]$dist.x1[spot] <- coord_lymph[[i]][which(dists == min(dists))[1], 1]
  coord_lymph_main[[i]]$dist.y1[spot] <- coord_lymph[[i]][which(dists == min(dists))[1], 2]
}
}


library(data.table)
coord_lymph_main_norm <- coord_lymph_main
for (i in 1:length(coord)){
coord_lymph_main_norm[[i]]$dist1 <- as.numeric(coord_lymph_main_norm[[i]]$dist1) / as.numeric(max(coord_lymph_main_norm[[i]]$dist1))
}

coord_lymph_main_norm <- rbindlist(coord_lymph_main_norm,fill=TRUE)


#run from here after sbatch complete



coord_lymph_main_norm$aggr_prox<-'far'
coord_lymph_main_norm$aggr_prox[coord_lymph_main_norm$dist1 > -0.1 & coord_lymph_main_norm$dist1 < 0.009999] <- "aggr"
coord_lymph_main_norm$aggr_prox[coord_lymph_main_norm$dist1 > 0.002 & coord_lymph_main_norm$dist1 < 0.05] <- "periAggr"
coord_lymph_main_norm$aggr_prox[coord_lymph_main_norm$dist1 > 0.05 & coord_lymph_main_norm$dist1 < 0.3] <- "int1"
coord_lymph_main_norm$aggr_prox[coord_lymph_main_norm$dist1 > 0.3 ]<-'far'
coord_lymph_main_norm_f <- coord_lymph_main_norm %>% select(dist1, aggr_prox) %>%  as.data.frame()
rownames(coord_lymph_main_norm_f) <- coord_lymph_main_norm$cells

coord_lymph_main_norm_f$dist1 <- NULL

all <- AddMetaData(all, metadata = coord_lymph_main_norm_f)

all_meta <- all@meta.data
s_obj_meta <- list()
ggplot_ls <- list()
for (i in 1:length(names)) {
s_obj_meta[[i]] <- all_meta[all_meta$orig.ident == names[[i]],];
xy_list[[i]]$aggr_prox <- s_obj_meta[[i]]$aggr_prox;
cols <- ArchR::paletteDiscrete(s_obj_meta[[i]][, "aggr_prox"]);
ggplot_ls[[i]] <- ggplot(xy_list[[i]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=aggr_prox)) +
    geom_point(size=0.1)+theme_classic()+scale_color_manual(values = c(cols))+ggtitle("aggr_prox");
print(ggplot_ls[[i]])
}

all_meta <- all@meta.data
s_obj_meta <- list()
ggplot_ls <- list()
for (i in 1:length(names)) {
s_obj_meta[[i]] <- all_meta[all_meta$orig.ident == names[[i]],];
xy_list[[i]]$Annotations_new <- s_obj_meta[[i]]$Annotations_new;
cols <- ArchR::paletteDiscrete(s_obj_meta[[i]][, "Annotations_new"]);
ggplot_ls[[i]] <- ggplot(xy_list[[i]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=Annotations_new)) +
    geom_point(size=0.1)+theme_classic()+scale_color_manual(values = c(cols))+ggtitle("Annotations_new");
print(ggplot_ls[[i]])
}


all_meta <- all@meta.data
s_obj_meta <- list()
ggplot_ls <- list()
for (i in 1:length(names)) {
s_obj_meta[[i]] <- all_meta[all_meta$orig.ident == names[[i]],];
xy_list[[i]]$named_niches <- s_obj_meta[[i]]$named_niches;
cols <- ArchR::paletteDiscrete(s_obj_meta[[i]][, "named_niches"]);
ggplot_ls[[i]] <- ggplot(xy_list[[i]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named_niches)) +
    geom_point(size=0.1)+theme_classic()+scale_color_manual(values = c(cols))+ggtitle("named_niches");
print(ggplot_ls[[i]])
}




cols <- ArchR::paletteDiscrete(all@meta.data[, "named_niches_proximity"])

ggplot(xy_list[[1]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named_niches_proximity)) +
    geom_point(size=1)+theme_classic()+scale_color_manual(values = c(cols))



```

```{r}

all$aggr_vasc_prox <- paste(all$aggr_prox, all$Annotations_new, sep="_")
table(all$aggr_vasc_prox)

all$aggr_vasc_prox_named_niches <- paste(all$aggr_vasc_prox, all$named_niches, sep="_")


pt2 <- table(all$aggr_vasc_prox, all$named_niches)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)

cols <- as.data.frame(ArchR::paletteDiscrete(all@meta.data[, "aggr_vasc_prox"]))
colnames(cols)<-"colors"

pt2 %>% filter(Var1 == "aggr_periVasc") %>%  ggplot(aes(x = Var1, y = Freq, fill = Var2))  + 
        geom_bar(stat = 'identity', position = position_fill())+
        theme(axis.text.y= element_blank(), axis.ticks.y = element_blank()) +
        coord_flip() +  scale_fill_manual(values = cols$colors)+
        theme(legend.position = 'bottom', panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank()) + 
        labs(fill = 'Tissue-Defined Cluster', y = 'Cluster Frequency in Sample') + 
        guides(fill = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 2)), alpha = FALSE)


pt2 %>% filter(Var1 == "far_periVasc") %>%  ggplot(aes(x = Var1, y = Freq, fill = Var2))  + 
        geom_bar(stat = 'identity', position = position_fill())+
        theme(axis.text.y= element_blank(), axis.ticks.y = element_blank()) +
        coord_flip() +  scale_fill_manual(values = cols$colors)+
        theme(legend.position = 'bottom', panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank()) + 
        labs(fill = 'Tissue-Defined Cluster', y = 'Cluster Frequency in Sample') + 
        guides(fill = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 2)), alpha = FALSE)


```

```{r}


library(tidyr)
#all$sample_condition<-paste(all$aggr_vasc_prox_named_niches, obj$condition, sep=".") #condition would be resting/peak/resolved etc

pt <- table(all$orig.ident, all$aggr_vasc_prox_named_niches)
pt <- as.data.frame(pt)

pt_f <- rbind(pt[pt$Var2 == "far_periVasc_Fibrotic",], pt[pt$Var2 == "periAggr_periVasc_Fibrotic",])


cluster1 <- ggplot(pt_f, aes(x=Var2, y=Freq)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("Cluster1")+theme_classic()+RotatedAxis()
cluster1
res.aov_cluster1 <- aov(Freq ~ Var2, data = pt_f)
summary(res.aov_cluster1)
stats_cluster1 <- TukeyHSD(res.aov_cluster1)
stats_cluster1

```






```{r}
#distance calculations (spatstat)

library(spatstat)
library(splitstackshape)
all_split <- SplitObject(all, split.by = "orig.ident")

data_x_y_dist <- data_x_y
ln_list <- list()
ANN_list <- list()
nnda_list <- list()
dist_mean <- list()

for (i in 1:length(data_x_y_dist)){
  data_x_y_dist[[i]]$class <- all_split[[i]]$named_niches
colnames(data_x_y_dist[[i]])<-c("x", "y", "class")
xlim = range(data_x_y_dist[[i]]$x)
ylim = range(data_x_y_dist[[i]]$y)
ln_list[[i]] = with(data_x_y_dist[[i]], ppp(x = x, y = y, marks = class, xrange = xlim, yrange = ylim))
 ANN_list[[i]] <- apply(nndist(ln_list[[i]], k=1:100),2,FUN=mean, by=marks(vascular))
plot(ANN_list[[i]] ~ eval(1:100), type="b", main=NULL, las=1)
#d<-nndist(ln_list[[i]])
nnda_list[[i]]<-nndist(ln_list[[i]], by=marks(ln_list[[i]]))
dist_mean[[i]]<-aggregate(nnda_list[[i]], by=list(from=marks(ln_list[[i]])), mean)
}

library(data.table)
for (i in 1:length(nnda_list)){
nnda_list[[i]] <- as.data.frame(nnda_list[[i]])
}

nnda_list_all <- rbindlist(nnda_list)



hist(nnda_list_all$Fibrotic)
median(nnda_list_all$Fibrotic)
max(nnda_list_all$Fibrotic)
min(nnda_list_all$Fibrotic)
range(nnda_list_all$Fibrotic)


ggplot(nnda_list_all, aes(x = Fibrotic)) +
  geom_histogram()+xlim(0,100)

ggplot(nnda_list_all, aes(x = Lymphocytic1)) +
  geom_histogram()+xlim(0,100)

ggplot(nnda_list_all, aes(x = vascular)) +
  geom_histogram()+xlim(0,100)


nnda_list_all$cell_id <- paste("c", colnames(all), sep=".")
quant <- quantile((nnda_list_all$Fibrotic), 0.10)
nnda_list_all_perivasc <-  nnda_list_all %>% filter_all(any_vars(. < 8))
nnda_list_all_perivasc <- cSplit(nnda_list_all_perivasc, splitCols = "cell_id", sep="." )

data_x_y_dist_lymph <- data_x_y
ln_list_lymph <- list()
ANN_list_lymph <- list()
nnda_list_lymph <- list()
dist_mean_lymph <- list()
for (i in 1:length(data_x_y_dist_lymph)){
  data_x_y_dist_lymph[[i]]$class <- all_split[[i]]$named_niches
colnames(data_x_y_dist_lymph[[i]])<-c("x", "y", "class")
xlim = range(data_x_y_dist_lymph[[i]]$x)
ylim = range(data_x_y_dist_lymph[[i]]$y)
ln_list_lymph[[i]] = with(data_x_y_dist_lymph[[i]], ppp(x = x, y = y, marks = class, xrange = xlim, yrange = ylim))
 ANN_list_lymph[[i]] <- apply(nndist(ln_list_lymph[[i]], k=1:100),2,FUN=mean, by=marks(Lymphocytic1))
plot(ANN_list_lymph[[i]] ~ eval(1:100), type="b", main=NULL, las=1)
#d<-nndist(ln_list[[i]])
nnda_list_lymph[[i]]<-nndist(ln_list_lymph[[i]], by=marks(ln_list_lymph[[i]]))
dist_mean_lymph[[i]]<-aggregate(nnda_list_lymph[[i]], by=list(from=marks(ln_list_lymph[[i]])), mean)
}

library(data.table)
for (i in 1:length(nnda_list_lymph)){
nnda_list_lymph[[i]] <- as.data.frame(nnda_list_lymph[[i]])
}

nnda_list_lymph_all <- rbindlist(nnda_list_lymph)
nnda_list_lymph_all$cell_id <- paste("c", colnames(all), sep=".")
quant <- quantile((nnda_list_lymph_all$Fibrotic), 0.010)
nnda_list_all_peri_aggr <-  nnda_list_lymph_all %>% filter_all(any_vars(. < 8))
nnda_list_all_peri_aggr <- cSplit(nnda_list_all_peri_aggr, splitCols = "cell_id", sep="." )

nnda_list_all
nnda_list_lymph_all



#mead/disy


dist_mean_vasc <- list()
for (i in 1:length(dist_mean)){
 dist_mean_vasc[[i]]<-dist_mean[[i]] %>% 
  select(-`?`, -Adipose, -`2`, -vascular, -Lymphocytic2, -`??`, -RBC_contam) 
dist_mean_vasc[[i]] <- dist_mean_vasc[[i]][dist_mean_vasc[[i]]$from == 'vascular', ]
}



dist_mean_vasc_all <- dist_mean_vasc %>%  rbindlist()

data_long1 <- melt(dist_mean_vasc_all,                                 # Apply melt function
                  id.vars = c("from"))


wisteria <- c("grey65", "burlywood3", "khaki2", "plum1", "lightcyan2", "cornflowerblue", "slateblue3")
 
  

  
  
dist_mean_lymph <- list()
for (i in 1:length(dist_mean)){
 dist_mean_lymph[[i]]<-dist_mean[[i]] %>% 
  select(-`?`, -Adipose, -`2`, -Lymphocytic1, -Lymphocytic2, -`??`, -RBC_contam) 
dist_mean_lymph[[i]] <- dist_mean_lymph[[i]][dist_mean_lymph[[i]]$from == 'Lymphocytic1', ]
}



dist_mean_lymph_all <- dist_mean_lymph %>%  rbindlist()

data_long1_lymph <- melt(dist_mean_lymph_all,                                 # Apply melt function
                  id.vars = c("from"))


wisteria <- c("grey65", "burlywood3", "khaki2", "plum1", "lightcyan2", "cornflowerblue", "slateblue3")




ggplot(data_long1, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill = variable),            #You can make a box plot too!
               alpha = 0.8, width = 0.7) +      
  geom_point(aes(fill = variable), shape = 21, color = "black", alpha = 0.8,
             position = position_jitter(width = 0.1, seed = 666))+
  scale_fill_manual(values = wisteria[c(3, 1, 6)]) +
  labs(x = "niche",
       y = "mean nnd/sample from vascular") +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_line(size = 1.2),
        text = element_text(size = 12, color = "black", face = "bold"),
        axis.text = element_text(size = 12, color = "black", face = "bold")
        )

 
  ggplot(data_long1_lymph, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill = variable),            #You can make a box plot too!
               alpha = 0.8, width = 0.7) +      
  geom_point(aes(fill = variable), shape = 21, color = "black", alpha = 0.8,
             position = position_jitter(width = 0.1, seed = 666))+
  scale_fill_manual(values = wisteria[c(3, 1, 6)]) +
  labs(x = "niche",
       y = "mean nnd/samplefrom lymphocytic1") +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_line(size = 1.2),
        text = element_text(size = 12, color = "black", face = "bold"),
        axis.text = element_text(size = 12, color = "black", face = "bold")  )
  







```


```{r}
#int <- (1:length(data_x_y))*1000
#int2 <- (1:length(data_x_y))
#data_x_y_niches <- data_x_y
#for (i in 1:length(data_x_y)){
#  data_x_y_niches[[i]]$Centroid.X.µm <- data_x_y_niches[[i]]$Centroid.X.µm +int[[i]];
#  data_x_y_niches[[i]]$cell_id <- paste(rownames(data_x_y_niches[[i]]), int2[[i]], sep="_")
#}
#library(data.table) 
#data_x_y_niches <- rbindlist(data_x_y_niches)
#data_x_y_niches <- data_x_y_niches %>%  as.data.frame()
#rownames(data_x_y_niches) <- data_x_y_niches$cell_id
#data_x_y_niches$cell_id <- NULL
```


```{r}
all_obj_split <- SplitObject(all, split.by = "orig.ident")

data_x_y_niches <- data_x_y
library(dbscan)
eps <- 50
nn_list <- list()
nn_df_list <- list()
cluster_ids <- list()
nn_count <- list()
nn_mat <- list()
k_means_res <- list()
k_means_id <- list()
k_means_df <- list()
nn_objs <- list()

for (i in 1:length(data_x_y_niches)){
  nn_list[[i]] <- frNN(x= data_x_y_niches[[i]] %>% as.matrix(), eps = eps)
  nn_df_list[[i]]<- nn_list[[i]]$id %>%
 stack()
cluster_ids[[i]] <- all_obj_split[[i]]$named_niches_final %>% unname()
nn_df_list[[i]]$cluster_id<- cluster_ids[[i]][nn_df_list[[i]]$values]
nn_df_list[[i]]$cluster_id<- factor(nn_df_list[[i]]$cluster_id)
nn_count[[i]]<- nn_df_list[[i]] %>%
  group_by(ind) %>%
  dplyr::count(cluster_id, .drop=F)
nn_count[[i]]<- nn_count[[i]] %>%
  tidyr::pivot_wider(names_from = cluster_id, values_from = n)
nn_mat[[i]]<- nn_count[[i]][,-1] %>% as.matrix()
rownames(nn_mat[[i]])<- nn_count[[i]]$ind
k_means_res[[i]]<- kmeans(nn_mat[[i]], centers = 6)
k_means_id[[i]]<- k_means_res[[i]]$cluster %>%
tibble::enframe(name = "cell_id", value = "kmeans_cluster")
k_means_df[[i]]<- as.data.frame(k_means_id[[i]])
rownames(k_means_df[[i]])<- k_means_id[[i]]$cell_id
  nn_objs[[i]] <- CreateSeuratObject(counts = t(nn_mat[[i]]),  min.features = 1)
  DefaultAssay(nn_objs[[i]]) <- 'RNA'
  nn_objs[[i]] <- ScaleData(nn_objs[[i]])
  nn_objs[[i]] <- RunPCA(nn_objs[[i]], npcs = 10, features = rownames(nn_objs[[i]]))
  nn_objs[[i]] <- FindNeighbors(nn_objs[[i]], reduction = "pca", dims = 1:5)
  nn_objs[[i]] <- RunUMAP(nn_objs[[i]], dims = 1:5)
  nn_objs[[i]] <- FindClusters(nn_objs[[i]], resolution = c(0.02))
}


names(nn_objs) <- names(data_x_y_niches)
names_orig <- names(data_x_y_niches)

for (i in 1:length(nn_objs)){
nn_objs[[i]]$orig.ident <- names[[i]] 
}

library(harmony)
all_niches_merged <- merge(x=nn_objs[[1]], y=nn_objs[-1])
all_niches_merged <-ScaleData(all_niches_merged)
all_niches_merged <-FindVariableFeatures(all_niches_merged)
all_niches_merged<-RunPCA(all_niches_merged)
all_niches_merged<-RunUMAP(all_niches_merged, dims = 1:5)
#re-run from here
DimPlot(all_niches_merged, group.by="orig.ident", raster=FALSE)
all_niches_merged <- RunHarmony(all_niches_merged, group.by.vars = "orig.ident")
all_niches_merged <- RunUMAP(all_niches_merged, reduction="harmony", dims=1:5)
DimPlot(all_niches_merged, group.by = "orig.ident", raster = F)

all_niches_merged <- FindNeighbors(all_niches_merged, reduction="harmony", dims = 1:5)
all_niches_merged <- FindClusters(all_niches_merged, resolution = c(0.02, 0.01, 0.05))
all_niches_merged <- FindClusters(all_niches_merged, resolution = c(0.03))

DimPlot(all_niches_merged, raster=F, group.by = "RNA_snn_res.0.01")
all_niches_merged_slit <- SplitObject(all_niches_merged, split.by = "orig.ident")


library(data.table)
ggplot_niches <- list()
for (i in 1:length(all_niches_merged_slit)){
data_x_y_niches[[i]] <- data_x_y_niches[[i]][rownames(data_x_y_niches[[i]]) %in% colnames(nn_objs[[i]]),]
data_x_y_niches[[i]]$niches <- all_niches_merged_slit[[i]]$RNA_snn_res.0.05
data_x_y_niches[[i]]$cell_id <- colnames(all_niches_merged_slit[[i]])
ggplot_niches[[i]] <- ggplot(data_x_y_niches[[i]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=niches)) +
    geom_point(size=0.1)+theme_classic()+ggtitle("niches_spatial") 
print(ggplot_niches[[i]])  
}
table(all_niches_merged$RNA_snn_res.0.03)


#remove clusters with <1000 cells
Idents(all_niches_merged) <- 'RNA_snn_res.0.03'
all_niches_merged_fil <- subset(all_niches_merged, idents=c("0", "2", "1", "3", "4", "5"))
all_niches_merged_slit <- SplitObject(all_niches_merged_fil, split.by = "orig.ident")

table(all_niches_merged_fil$RNA_snn_res.0.01, all_niches_merged_fil$orig.ident)


#plot res0.03 looks good
library(data.table)
data_x_y_niches_f <- data_x_y_niches
ggplot_niches <- list()
smaple_id <- c(1:12)
for (i in 1:length(all_niches_merged_slit)){
rownames(data_x_y_niches_f[[i]]) <- paste(rownames(data_x_y_niches_f[[i]]), smaple_id[[i]], sep="_")
data_x_y_niches_f[[i]] <- data_x_y_niches_f[[i]][rownames(data_x_y_niches_f[[i]]) %in% colnames(all_niches_merged_slit[[i]]),]
data_x_y_niches_f[[i]]$niches <- all_niches_merged_slit[[i]]$RNA_snn_res.0.03
data_x_y_niches_f[[i]]$cell_id <- colnames(all_niches_merged_slit[[i]])
ggplot_niches[[i]] <- ggplot(data_x_y_niches_f[[i]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=niches)) +
    geom_point(size=0.1)+theme_classic()+ggtitle("niches_spatial") 
print(ggplot_niches[[i]])  
}



cell_fun = function(j, i, x, y, width, height, fill) {
                grid::grid.rect(x = x, y = y, width = width *0.99, 
                                height = height *0.99,
                                gp = grid::gpar(col = "grey", 
                                                fill = fill, lty = 1, lwd = 0.5))
}

col_fun=circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))


avg_abun<- AverageExpression(
  all_niches_merged_fil,
  assays = NULL,
  features = rownames(all_niches_merged_fil),
  return.seurat = FALSE,
  group.by = "RNA_snn_res.0.03")
```
```{r}
 ggplot(data_x_y_niches_f[[6]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=niches)) +
    geom_point(size=0.8)+theme_classic()+ggtitle("niches_spatial")+xlim(4750,6400)+ylim(4000,6000)


for (i in 1:length(data_x_y_niches_f)){
  print(ggplot_niches[[i]])
 }

cols3=c("grey", "grey", "grey", "red", "grey", "blue", "grey", "grey")
all_meta <- all@meta.data
s_obj_meta <- list()
ggplot_ls <- list()
for (i in 1:length(names)) {
s_obj_meta[[i]] <- all_meta[all_meta$orig.ident == names[[i]],];
xy_list[[i]]$named_niches <- s_obj_meta[[i]]$named_niches_final;
cols <- ArchR::paletteDiscrete(s_obj_meta[[i]][, "named_niches_final"]);
ggplot_ls[[i]] <- ggplot(xy_list[[i]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named_niches)) +
    geom_point(size=0.1)+theme_classic()+scale_color_manual(values = c(cols3))+ggtitle("named_niches");
print(ggplot_ls[[i]])
}



cols_n <- ArchR::paletteDiscrete(data_x_y_niches_f[[i]][, "niches"]);

#early RA
ggplot(xy_list[[1]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named_niches)) +
    geom_point(size=0.6)+theme_classic()+scale_color_manual(values = c(cols3))+xlim(3000,4500)+ylim(1800,3200)

ggplot(data_x_y_niches_f[[1]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=niches)) +
    geom_point(size=0.6)+theme_classic()+xlim(3000,4500)+ylim(1800,3200)


#est RA
ggplot(xy_list[[6]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named_niches)) +
    geom_point(size=0.6)+theme_classic()+scale_color_manual(values = c(cols3))+xlim(4000,6000)+ylim(2200,3500)

ggplot(data_x_y_niches_f[[6]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=niches)) +
    geom_point(size=0.6)+theme_classic()+xlim(4000,6000)+ylim(2200,3500)


#OA
ggplot(xy_list[[9]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named_niches)) +
    geom_point(size=0.6)+theme_classic()+scale_color_manual(values = c(cols3))+xlim(0,3000)+ylim(1000,5000)

ggplot(data_x_y_niches_f[[9]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=niches)) +
    geom_point(size=0.6)+theme_classic()+xlim(0,3000)+ylim(1000,5000)


#res
ggplot(xy_list[[12]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named_niches)) +
    geom_point(size=0.6)+theme_classic()+scale_color_manual(values = c(cols3))+xlim(0,2500)+ylim(4000,6500)

ggplot(data_x_y_niches_f[[12]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=niches)) +
    geom_point(size=0.6)+theme_classic()+xlim(0,2500)+ylim(4000,6500)

```


```{r}
Heatmap(t(scale(t(avg_abun$RNA))),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        rect_gp = grid::gpar(type = "none"),
        cell_fun = cell_fun,
        col = col_fun,
        column_names_rot = 45)
```


```{r}
library(data.table)
data_x_y_niches_all <- rbindlist(data_x_y_niches_f)
all_niches_f <- all[,colnames(all) %in% data_x_y_niches_all$cell_id]
all_niches_f$spatial_niches <- data_x_y_niches_all$niches
Idents(all_niches_f) <- "spatial_niches"
levels(all_niches_f)
DefaultAssay(all_niches_f) <- 'SCT'
Idents(all_niches_f) <- 'spatial_niches'
DotPlot(all_niches_f, features = rownames(all_niches_f))+RotatedAxis()
DotPlot(all_niches_f, features = c("Cell..COL6A1.mean","Cell..COMP.mean", "Cell..CD31.mean", "Cell..CD3.mean"))+RotatedAxis()


Idents(all_niches_f) <- 'spatial_niches'
DefaultAssay(all_niches_f) <- "SCT"
levels(all_niches_f)<-c( "1", "0", "2", "3", "4", "5")

DotPlot(all_niches_f, features = c("Cell..COL1.mean", "Cell..COL6A1.mean", "Cell..POSTN.mean",  "Cell..MMP3.mean", "Cell..PDPN.mean","Cell..CD68.mean", "Cell..CD3.mean", "Cell..CD20.mean", "Cell..FABP4.mean",  "Cell..SMA.mean" ,  "Cell..CD31.mean" ,  "Cell..CD146.mean",  "Cell..CD138.mean"), idents=levels(all_niches_f)[-c(4,6)])+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))


DefaultAssay(all_niches_f) <- "SCT"
DotPlot(all_niches_f, features = c("Cell..COL1.mean", "Cell..COL6A1.mean", "Cell..COL4A1.mean"), idents = c("0", "1", "4"))+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))


#instead of GEX look at avg 'expr' in celll numbers or neighbours to define niches, then rename/sample and combine

rownames(all_niches_f)


```

```{r}

pt2 <- table(all_niches_f$named_niches_final, all_niches_f$spatial_niches)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
pt2$Var2 <- as.integer(pt2$Var2)
pt2 <- pt2 %>% filter(Var2 < 6)

cols <- as.data.frame(ArchR::paletteDiscrete(all_niches_f@meta.data[, "named_niches"]))
colnames(cols)<-"colors"

ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1))  + 
        geom_bar(stat = 'identity', position = position_fill())+
        theme( axis.ticks.y = element_blank()) +
        coord_flip() +  scale_fill_manual(values = cols$colors)+
        theme(legend.position = 'bottom', panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank()) + 
        labs(fill = 'Tissue-Defined Cluster', y = 'Cluster Frequency in Sample') + 
        guides(fill = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 2)), alpha = FALSE)


```


```{r}


#how these niches chnage over disease?

#scProprtiontest across pathotypes, how the neches cngae between two?

Idents(all_niches_f) <- 'spatial_niches'
all_niches_f2 <- subset(all_niches_f, idents=levels(all_niches_f)[-c(4,6)])
levels(all_niches_f2)
test <- sc_utils(all_niches_f2)
Idents(all_niches_f2) <- 'condition'
levels <- (levels(all_niches_f2)[-3])

library(scProportionTest)
prop.test_1 <- permutation_test(test, cluster_identity = "spatial_niches", sample_1="Res", sample_2="eRA", sample_identity="condition", n_permutations=10000)
permutation_plot(prop.test_1, FDR_threshold = 0.01, log2FD_threshold = 0.5, order_clusters = F)+
  geom_point(size=5,color = c( "grey", "red", "grey", "red"))
##
library(scProportionTest)
prop.test_2 <- permutation_test(test, cluster_identity = "spatial_niches", sample_1="eRA", sample_2="EstRA", sample_identity="condition", n_permutations=10000)
permutation_plot(prop.test_2, FDR_threshold = 0.01, log2FD_threshold = 0.5, order_clusters = F)+
  geom_point(size=5,color = c( "red", "grey", "red", "red"))

library(scProportionTest)
prop.test_3 <- permutation_test(test, cluster_identity = "spatial_niches", sample_1="OA", sample_2="EstRA", sample_identity="condition", n_permutations=10000)
permutation_plot(prop.test_3, FDR_threshold = 0.01, log2FD_threshold = 0.5, order_clusters = F)+
 geom_point(size=5,color = c("grey", "grey", "grey", "red"))


library(scProportionTest)
prop.test_4 <- permutation_test(test, cluster_identity = "spatial_niches", sample_1="Res", sample_2="OA", sample_identity="condition", n_permutations=10000)
permutation_plot(prop.test_4, FDR_threshold = 0.01, log2FD_threshold = 0.5, order_clusters = F)#+
  #geom_point(size=4,color = c("grey", "grey", "red", "red", "red"))


library(scProportionTest)
prop.test_5 <- permutation_test(test, cluster_identity = "spatial_niches", sample_1="eRA", sample_2="OA", sample_identity="condition", n_permutations=10000)
permutation_plot(prop.test_5, FDR_threshold = 0.01, log2FD_threshold = 0.5, order_clusters = F)


pt2 <- table(all_niches_f2$spatial_niches, all_niches_f2$condition) %>% as.data.frame()
pt2$Var1 <- as.character(pt2$Var1)

cols <- as.data.frame(ArchR::paletteDiscrete(all_niches_f@meta.data[, "spatial_niches"]))
colnames(cols)<-"colors"

pt2 <- pt2 %>% filter(Freq > 0)

#DimPlot()+ggtitle

ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1))  + 
        geom_bar(stat = 'identity', position = position_fill())+
        theme() +
        coord_flip() +  scale_fill_manual(values = cols$colors)+
        theme(legend.position = 'bottom', panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank()) + 
        labs(fill = 'Tissue-Defined Cluster', y = 'Cluster Frequency in Sample') + 
        guides(fill = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 2)), alpha = FALSE)

library(scProportionTest)
for (i in 1:length(prop.test_list)){
  print(permutation_plot(prop.test_list[[i]], FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = F))}



```

```{r}
table(all_niches_f$spatial_niches, all_niches_f$condition)
```

```{r}
pt2 <- table(all_niches_f$spatial_niches, all_niches_f$orig.ident) %>% as.data.frame()
pt2$Var1 <- as.character(pt2$Var1)

cols <- as.data.frame(ArchR::paletteDiscrete(all_niches_f@meta.data[, "spatial_niches"]))
colnames(cols)<-"colors"

pt2 <- pt2 %>% filter(Freq > 0)


ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1))  + 
        geom_bar(stat = 'identity', position = position_fill())+
        theme() +
        coord_flip() +  scale_fill_manual(values = cols$colors)+
        theme(legend.position = 'bottom', panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank()) + 
        labs(fill = 'Tissue-Defined Cluster', y = 'Cluster Frequency in Sample') + 
        guides(fill = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 2)), alpha = FALSE)
```

```{r}
pt2 <- table(all_niches_f$named_niches_final, all_niches_f$orig.ident) %>% as.data.frame()
pt2$Var1 <- as.character(pt2$Var1)

cols <- as.data.frame(ArchR::paletteDiscrete(all_niches_f@meta.data[, "named_niches_final"]))
colnames(cols)<-"colors"

pt2 <- pt2 %>% filter(Freq > 0)


ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1))  + 
        geom_bar(stat = 'identity', position = position_fill())+
        theme() +
        coord_flip() +  scale_fill_manual(values = cols$colors)+
        theme(legend.position = 'bottom', panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank()) + 
        labs(fill = 'Tissue-Defined Cluster', y = 'Cluster Frequency in Sample') + 
        guides(fill = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 2)), alpha = FALSE)
```






```{r}

#do same for visium

pt3 <- table(merging_vis$named_clusters, merging_vis$diseasestate) %>% as.data.frame()
pt3$Var1 <- as.character(pt3$Var1)

cols <- as.data.frame(ArchR::paletteDiscrete(merging_vis@meta.data[, "named_clusters"]))
colnames(cols)<-"colors"

pt3 <- pt3 %>% filter(Freq > 0)


ggplot(pt3, aes(x = Var2, y = Freq, fill = Var1))  + 
        geom_bar(stat = 'identity', position = position_fill())+
        theme() +
        coord_flip() +  scale_fill_manual(values = cols$colors)+
        theme(legend.position = 'bottom', panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank()) + 
        labs(fill = 'Tissue-Defined Cluster', y = 'Cluster Frequency in Sample') + 
        guides(fill = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 2)), alpha = FALSE)

test2 <- sc_utils(merging_vis)
Idents(merging_vis) <- 'diseasestate'
levels2 <- (levels(merging_vis)[-2])
prop.test_list2 <- list()
for (i in 1:length(levels2)){
prop.test_list2[[i]] <- permutation_test(test2, cluster_identity = "named_clusters", sample_1="OA", sample_2=levels2[[i]], sample_identity="diseasestate", n_permutations=10000)
}

library(scProportionTest)
for (i in 1:length(prop.test_list2)){
  print(permutation_plot(prop.test_list2[[i]], FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = F))}


```













```{r}

x_y_test <- data_x_y$eRA_BX020.txt


x_y_test_f <- x_y_test[rownames(x_y_test) %in% colnames(nn_obj),]

x_y_test_f$niches_spatial <- nn_obj@meta.data$RNA_snn_res.0.1
#cols <- ArchR::paletteDiscrete(s_obj_meta[[i]][, "named_niches"]);
ggplot(x_y_test_f, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=niches_spatial)) +
    geom_point(size=0.1)+theme_classic()+ggtitle("niches_spatial")





```
```{r}


Idents(nn_obj) <- 'RNA_snn_res.0.1'
DotPlot(nn_obj, features=rownames(nn_obj))+RotatedAxis()

#add RNA_snn_res.0.1 to s_obj and look at expression


```

```{r}


library(devtools)
library(spatula)
library(furrr)

names(xy_list) <- names(data_x_y)

xy_list_coloc_disease <- xy_list

coloc_res_coarse_list <- list()

for (i in 1:length(xy_list_coloc_disease)){
coloc_res_coarse<-tryCatch({
    coloc_res_coarse_list[[i]] = coloc_all_types(
        index_type = unique(xy_list_coloc_disease[[i]]$named),
        coords = xy_list_coloc_disease[[i]][, c("Centroid.X.µm", "Centroid.Y.µm")],
        y = xy_list_coloc_disease[[i]]$named,
        compartments = NULL,
        max_dist = 40,
        nperm = 1000,
        parallel = TRUE
    )}, error = function(e) {
    message('Failed run on: ', xy_list_coloc_disease[[i]])
    return(NULL)
    }
)
}

names(coloc_res_coarse_list) <- names(data_x_y)

for (i in 1:length(names(data_x_y))){
  coloc_res_coarse_list[[i]]$sample <- names(data_x_y)[[i]]}

library(data.table)
library(splitstackshape)

library(ggpubr)
library(rstatix)
library(DescTools)


coloc_res_coarse_list_bind <- rbindlist(coloc_res_coarse_list) %>% cSplit(splitCols="sample", sep="_")
coloc_res_coarse_list_bind %>% filter(index_type == "Lymphocytic1" & type != "low_quality" & type != "RBC_contam") %>% as.data.frame() %>%  ggplot(aes(x=factor(sample_1, levels=c("OA", "eRA", "EstRA", "Res")), y=freq, fill=sample_1)) +
  geom_boxplot() + geom_jitter(width=0.1,alpha=0.2)+facet_wrap(vars(type))+RotatedAxis()


new_fre_df <- coloc_res_coarse_list_bind %>% filter(index_type == "Lymphocytic1" & type != "low_quality" & type != "RBC_contam") %>% as.data.frame() %>% select(sample_2, zscore, type) %>%
  pivot_wider(names_from = sample_2, values_from = zscore) %>% 
  as.data.frame()  

library(dplyr)
library(tidyverse)

new_fre_df <- column_to_rownames(new_fre_df,'type')  
new_fre_df$mean_early <- rowMeans(new_fre_df[,1:3])
new_fre_df$mean_Est <- rowMeans(new_fre_df[,4:6])
new_fre_df$mean_Oa <- rowMeans(new_fre_df[,7:9])
new_fre_df$mean_Res <- rowMeans(new_fre_df[,10:12])

new_fre_df[-4,] %>% dplyr::select(mean_early, mean_Est, mean_Oa, mean_Res) %>% as.matrix() %>% Heatmap(col = circlize::colorRamp2(c(-20, 0, 20), c('blue', 'white', 'red')), border = T, cluster_columns = F, cluster_rows = F)


```

```{r}

# of lymphatic cells in vascular niche

all_niches_f$niches_named <- paste(all_niches_f$spatial_niches, all_niches_f$named_niches_final, sep="_")

nums <- table(all_niches_f2$spatial_niches, all_niches_f2$condition) %>% as.data.frame()
#nums <- filter(nums, grepl("1_",Var1, ignore.case = TRUE) & Var1 != "1_low_quality")
nums <- nums %>% filter(Freq >0)

nums$Var1 <- as.character(nums$Var1)

cols <- as.data.frame(ArchR::paletteDiscrete(all_niches_f2@meta.data[, "spatial_niches"]))
colnames(cols)<-"colors"

ggplot(nums, aes(x = Var2, y = Freq, fill = Var1))  + 
        geom_bar(stat = 'identity', position = position_fill()) +
        coord_flip() +  scale_fill_manual(values = cols$colors)+
        theme(legend.position = 'bottom', panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank()) + 
        labs(fill = 'Tissue-Defined Cluster', y = 'Cluster Frequency in Sample') + 
        guides(fill = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 2)), alpha = FALSE)

```

```{r}


all_niches_f$sample_condition<-paste(all_niches_f$orig.ident, all_niches_f$condition, sep=".") #condition would be resting/peak/resolved etc

nums <- table(all_niches_f$sample_condition, all_niches_f$niches_named) %>% as.data.frame()
nums$Var2 <- paste("n", nums$Var2,sep="" )

nums <- filter(nums, grepl("4_",Var2, ignore.case = TRUE) & Var2 != "4_low_quality")

nums<-nums %>%
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  as.data.frame() 

library(tidyverse)
nums <- nums %>% remove_rownames %>% column_to_rownames(var="Var1")
nums <- nums/rowSums(nums)

nums$condition<-rownames(nums)

library(splitstackshape)
nums<-cSplit(nums, splitCols="condition", sep="_")

cluster1 <- ggplot(nums, aes(x=condition_1, y=n4_Lymphocytic1)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("1_Lymphocytic1")+theme_classic()+RotatedAxis()
cluster1
res.aov_cluster1 <- aov(Lymphocytic1 ~ condition_2, data = pt)
summary(res.aov_cluster1)
stats_cluster1 <- TukeyHSD(res.aov_cluster1)
stats_cluster1
```
```{r}

#expresison of collagens etc in peri aggr/perivasc fibs

table(merging_vis$diseasestate)
table(merging_vis$Annotations_new_new)

merging_vis$prox_disease <- paste(merging_vis$Annotations_new_new, merging_vis$diseasestate, sep="_")


Idents(merging_vis) <- 'prox_disease'
levels <- levels(merging_vis)[grep("periVasc_fibro|periAggr_fibro", levels(merging_vis))]
DotPlot(merging_vis, features = c("COL1A1", "COL4A1", "COL6A1", "MMP3", "MMP1"), idents = levels)+RotatedAxis()

SpatialDimPlot(merging_vis, group.by = "Annotations_new_new",images="slice1")
SpatialDimPlot(merging_vis, group.by = "Annotations_new_new",images="slice1.13")



SpatialFeaturePlot(merging_vis, images="slice1", features =c("COL1A1", "COL4A1", "COL6A1") )
SpatialFeaturePlot(merging_vis, images="slice1.17", features =c("COL1A1", "COL4A1", "COL6A1") )


table(paste(merging_vis$orig.ident, merging_vis$diseasestate))[19]


```
```{r}

all_niches_f2$niche_celltype <- paste(all_niches_f2$spatial_niches, all_niches_f2$named_niches_final, sep="_")
Idents(all_niches_f2) <- "niche_celltype"
levels(all_niches_f2)

table(all_niches_f2$niche_celltype)

DotPlot(all_niches_f2, idents=levels(all_niches_f2)[grep("Fibrotic", levels(all_niches_f2))], features=c("Cell..COL6A1.mean", "Cell..COMP.mean"))+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))


Idents(all_niches_f2) <- 'spatial_niches'

DefaultAssay(all_niches_f2) <- "SCT"
DotPlot(all_niches_f2, features = c( "Cell..COL6A1.mean", "Cell..COMP.mean",  "Cell..CD31.mean" ,  "Cell..CD3.mean"), idents=levels(all_niches_f2))+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))



all_niches_f2$celltype_group <- paste(all_niches_f2$named_niches_final, all_niches_f2$group, sep="_")
Idents(all_niches_f2) <- "celltype_group"
levels(all_niches_f2)


DotPlot(all_niches_f2, features = c( "Cell..COL6A1.mean", "Cell..COMP.mean",  "Cell..COL1.mean" ,  "Cell..POSTN.mean", "Cell..SPARC.mean"), idents=levels(all_niches_f2)[grep("Fibrotic_", levels(all_niches_f2))])+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))



```

```{r}


#in niche 1 colagen expression across timpoints


Idents(all_niches_f2) <- 'spatial_niches'
levels(all_niches_f2)

all_niches_f2$niche_condition <- paste(all_niches_f2$spatial_niches, all_niches_f2$group, sep="_")

Idents(all_niches_f2) <- 'niche_condition'



DefaultAssay(all_niches_f2) <- "SCT"
DotPlot(all_niches_f2, features = c( "Cell..COL6A1.mean", "Cell..COMP.mean",  "Cell..COL1.mean" ,  "Cell..POSTN.mean", "Cell..SPARC.mean"), idents=levels(all_niches_f2)[grep("1_", levels(all_niches_f2))])+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))



```






```{r}
table(all_niches_f$orig.ident)

all_niches_f$pathotype <- all_niches_f@meta.data[["orig.ident"]]
Idents(all_niches_f) <- 'pathotype'
levels(all_niches_f)
current.sample.ids <- c("eRA_BX020.txt" ,  "eRA_BX031.txt"  , "eRA_BX115.txt" ,  "EstRA_BX127.txt", "EstRA_BX195.txt", "EstRA_BX240.txt", "OA_JRP127.txt"  , "OA_JRP146.txt" ,  "OA_JRP149.txt"  , "Res_BX026.txt" ,  "Res_BX028.txt" ,  "Res_BX202.txt"  )
new.sample.ids <- c("PI" ,  "PI"  , "Diffuse" ,  "PI", "Diffuse", "Lymphoid", "JRP"  , "JRP" ,  "JRP"  , "PI" ,  "PI" ,  "Lymphoid")

all_niches_f@meta.data[["pathotype"]] <- plyr::mapvalues(x = all_niches_f@meta.data[["pathotype"]], from = current.sample.ids, to = new.sample.ids)


pt2 <- table(all_niches_f2$spatial_niches, all_niches_f2$condition) %>% as.data.frame()
pt2$Var1 <- as.character(pt2$Var1)

cols <- as.data.frame(ArchR::paletteDiscrete(all_niches_f@meta.data[, "spatial_niches"]))
colnames(cols)<-"colors"

pt2 <- pt2 %>% filter(Freq > 0)

ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1))  + 
        geom_bar(stat = 'identity', position = position_fill())+
        theme() +
        coord_flip() +  scale_fill_manual(values = cols$colors)+
        theme(legend.position = 'bottom', panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank()) + 
        labs(fill = 'Tissue-Defined Cluster', y = 'Cluster Frequency in Sample') + 
        guides(fill = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 2)), alpha = FALSE)


test <- sc_utils(all_niches_f2)
Idents(all_niches_f2) <- 'pathotype'

prop.test <- permutation_test(test, cluster_identity = "spatial_niches", sample_1="Diffuse", sample_2="Lymphoid", sample_identity="pathotype", n_permutations=10000)
p1 <- permutation_plot(prop.test, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = F)

prop.test <- permutation_test(test, cluster_identity = "spatial_niches", sample_1="PI", sample_2="Lymphoid", sample_identity="pathotype", n_permutations=10000)
p2 <-permutation_plot(prop.test, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = F)

prop.test <- permutation_test(test, cluster_identity = "spatial_niches", sample_1="Lymphoid", sample_2="Diffuse", sample_identity="pathotype", n_permutations=10000)
p3 <-permutation_plot(prop.test, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = F)

prop.test <- permutation_test(test, cluster_identity = "spatial_niches", sample_1="PI", sample_2="Diffuse", sample_identity="pathotype", n_permutations=10000)
p4 <-permutation_plot(prop.test, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = F)
```


```{r}


#named_niches_final
Idents(all_niches_f2) <- "named_niches_final"
levels(all_niches_f2)[c(1,2,4,5,6)]

all_niches_f2_f <- subset(all_niches_f2, idents=levels(all_niches_f2)[c(1,2,4,5,6)])
test <- sc_utils(all_niches_f2_f)


library(scProportionTest)
prop.test1 <- permutation_test(test, cluster_identity = "named_niches_final", sample_1="Res", sample_2="eRA", sample_identity="condition", n_permutations=10000)
p5 <- permutation_plot(prop.test1, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = F)

prop.test2 <- permutation_test(test, cluster_identity = "named_niches_final", sample_1="eRA", sample_2="EstRA", sample_identity="condition", n_permutations=10000)
p6 <-permutation_plot(prop.test2, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = F)

prop.test3 <- permutation_test(test, cluster_identity = "named_niches_final", sample_1="OA", sample_2="EstRA", sample_identity="condition", n_permutations=10000)
p7 <-permutation_plot(prop.test3, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = F)
#

p5+geom_point(size=5,color = c( "grey", "red", "red", "red", "grey"))

p6+geom_point(size=5,color = c( "red", "red", "red", "red", "grey"))

p7+geom_point(size=5,color = c( "grey", "red", "grey", "grey", "grey"))


```

```{r}


all_niches_f2$niche_pathotype <- paste(all_niches_f2$spatial_niches, all_niches_f2$pathotype, sep="_")
Idents(all_niches_f2) <- "niche_pathotype"
to_plot <- levels(all_niches_f2)[grep("PI|Lymphoid", levels(all_niches_f2))]
to_plot <- to_plot[grep("1_|0_|2_|4_", to_plot)]
rownames(all_niches_f)

dotplot<-DotPlot(all_niches_f2, features =c("Cell..COL1.mean","Cell..COL4A1.mean"   ,  "Cell..COL6A1.mean","Cell..COMP.mean", "Cell..SPARC.mean") , idents = to_plot)+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

dotplot

dotplot<-dotplot$data

dotplot<-dotplot %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()


to_plot_df <- to_plot %>% as.data.frame()
library(splitstackshape)
to_plot_df <- cSplit(to_plot_df, splitCols = ".", sep="_")
to_plot_df$._1 <- paste("S",to_plot_df$._1, sep="")

col_ann <- HeatmapAnnotation(df = to_plot_df)


library(ComplexHeatmap)

Heatmap(dotplot, top_annotation = col_ann)

```
```{r}

to_plot_s <- levels(all_niches_f)[grep("PI|Lymphoid", levels(all_niches_f))]
to_plot_s <- to_plot_s[grep("4_", to_plot_s)]

DotPlot(all_niches_f, features =c("Cell..COL1.mean", "Cell..COL4A1.mean", "Cell..COL6A1.mean") , idents =to_plot_s)+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))


all$pathotype <- all@meta.data[["orig.ident"]]
Idents(all) <- 'pathotype'
levels(all)
current.sample.ids <- c("eRA_BX020.txt" ,  "eRA_BX031.txt"  , "eRA_BX115.txt" ,  "EstRA_BX127.txt", "EstRA_BX195.txt", "EstRA_BX240.txt", "OA_JRP127.txt"  , "OA_JRP146.txt" ,  "OA_JRP149.txt"  , "Res_BX026.txt" ,  "Res_BX028.txt" ,  "Res_BX202.txt"  )
new.sample.ids <- c("PI" ,  "PI"  , "Diffuse" ,  "PI", "Diffuse", "Lymphoid", "JRP"  , "JRP" ,  "JRP"  , "PI" ,  "PI" ,  "Lymphoid")

all@meta.data[["pathotype"]] <- plyr::mapvalues(x = all@meta.data[["pathotype"]], from = current.sample.ids, to = new.sample.ids)


all$niche_pathotype <- paste(all$named_niches_final, all$pathotype, sep="_")

Idents(all) <- 'niche_pathotype'
idents_use <- levels(all)[grep("_PI|Lymphoid", levels(all))]
idents_use[grep("Fibrotic", idents_use)]

DefaultAssay(all) <- 'SCT'
DotPlot(all, features =c("Cell..COL6A1.mean", "Cell..COMP.mean") , idents =c("Fibrotic_PI"    ,   "Fibrotic_Lymphoid"))+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

```



```{r}


all_niches_f2$niche_pathotype <- paste(all_niches_f2$spatial_niches, all_niches_f2$pathotype, sep="_")

Idents(all_niches_f2) <- 'spatial_niches'

DotPlot(all, features =c("Cell..COL6A1.mean", "Cell..COMP.mean") , idents =c("Fibrotic_PI"    ,   "Fibrotic_Lymphoid"))+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))




```



