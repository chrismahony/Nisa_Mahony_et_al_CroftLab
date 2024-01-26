
hg19.v27.samplefiltered <- read.delim("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/Bulk_human_data/hum0207_RNAseq/hum0207.v1.RNA.v1/RNAseq/hg19.v27.samplefiltered.htseq")

hg19.v27.samplefiltered<-hg19.v27.samplefiltered  %>% distinct(gname, .keep_all = TRUE)

meta.data <- read.delim("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/Bulk_human_data/hum0207_RNAseq/hum0207.v1.RNA.v1/RNAseq/meta.data.txt", header=FALSE)


hg19.v27.samplefiltered_all<-hg19.v27.samplefiltered[,colnames(hg19.v27.samplefiltered) %in% meta.data$V1,]

rownames(hg19.v27.samplefiltered_all)=hg19.v27.samplefiltered$gname

hg19.v27.samplefiltered_all$ENSGID<-NULL
hg19.v27.samplefiltered_all$gname<-NULL

library(splitstackshape)
meta_data_dds=colnames(hg19.v27.samplefiltered_all)
meta_data_dds<-as.data.frame(meta_data_dds)
meta_data_dds<-cSplit(meta_data_dds, "meta_data_dds", sep="_", type.convert=FALSE)
meta_data_dds$sample<-colnames(hg19.v27.samplefiltered_all)
meta_data_dds<- meta_data_dds %>% remove_rownames %>% column_to_rownames(var="sample")
meta_data_dds$meta_data_dds_1<-NULL
colnames(meta_data_dds)<-('group')



dds <- DESeqDataSetFromMatrix(countData=hg19.v27.samplefiltered_all, 
                              colData=meta_data_dds, 
                              design=~group, tidy = F)

dds$group <- relevel(dds$group, ref = "NS")

#mingenecount <- quantile(rowSums(counts(dds)), 0.5)
mingenecount <- 200
maxgenecount <- quantile(rowSums(counts(dds)), 0.999)
dim(counts(dds))
# Subset low-expressed genes
keep <- rowSums(counts(dds)) > mingenecount & rowSums(counts(dds)) < maxgenecount
dds <- dds[keep, ]
print(quantile(rowSums(counts(dds))))
dim(dds)

dds@colData[['group']] <- as.factor(dds@colData[['group']])

design(dds) <- formula(~ group)
print(design(dds))
dds <- DESeq(dds, test = "Wald")


targetvar <- "group"

#comps1 <- data.frame(t(combn(unique(as.character(meta_data_dds[[targetvar]])), 2)))
     # head(comps1)
      
      comps1 <- data.frame(X1=unique(meta.data$V2)[-c(6,11)] ,X2= rep("NS",9))      

      
      ress <- apply(comps1, 1, function(cp) {
        print(cp)
        res <- data.frame(results(dds, contrast=c(targetvar, cp[1], cp[2])))
        res[["gene"]] <- rownames(res)
        res[["comparison"]] <- paste0(cp[1], "_vs_", cp[2])
        res
      })

     
            
      res1 <- Reduce(rbind, ress)

res1 %>% 
      filter(padj < 0.05) %>%
      mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
      arrange(desc(abs(score))) -> subres

for (i in 1:length(unique(res1$comparison))){res1 %>% filter(comparison == unique(res1$comparison)[[i]] & padj < 0.05 & log2FoldChange > 2) %>% write.csv(paste0("/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/paper_figures/", unique(res1$comparison)[[i]], "_new.csv"))}



