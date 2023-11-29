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

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


dds$group <- relevel(dds$group, ref = "NS")

dds <- DESeq(dds)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)


deseq2Results <- results(dds)
deseq2ResDF <- as.data.frame(deseq2Results)
summary(deseq2Results)
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA)
#strict
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .01,])


deseq2VST <- vst(dds)
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

library(ComplexHeatmap)
Heatmap(deseq2VST,cluster_columns =F)

library(DESeq2)
resultsNames(dds)
IL6_vs_NS <- lfcShrink(dds, coef="group_IL6_vs_NS", type="apeglm")
IL6_vs_NS<-as.data.frame(IL6_vs_NS)
IL6_vs_NS$gene<-rownames(IL6_vs_NS)
IL6_vs_NS<-IL6_vs_NS[order(IL6_vs_NS$padj),]
#IL6_vs_NS<-IL6_vs_NS[IL6_vs_NS$log2FoldChange > 0.2,]
#write.csv(IL6_vs_NS, "/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/IL6_vs_NS.csv")

#write.csv(IL6_vs_NS, "/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/zone_annotationIL6_vs_NS_sorted.csv")

TNFa_vs_NS <- lfcShrink(dds, coef="group_TNFa_vs_NS", type="apeglm")
TNFa_vs_NS<-as.data.frame(TNFa_vs_NS)
TNFa_vs_NS$gene<-rownames(TNFa_vs_NS)
TNFa_vs_NS<-TNFa_vs_NS[order(TNFa_vs_NS$padj),]
#TNFa_vs_NS<-TNFa_vs_NS[TNFa_vs_NS$log2FoldChange > 0.2,]
#write.csv(TNFa_vs_NS, "/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/TNFa_vs_NS.csv")

TGFb_vs_NS <- lfcShrink(dds, coef="group_TGFb_vs_NS", type="apeglm")
TGFb_vs_NS<-as.data.frame(TGFb_vs_NS)
TGFb_vs_NS$gene<-rownames(TGFb_vs_NS)
TGFb_vs_NS<-TGFb_vs_NS[order(TGFb_vs_NS$padj),]
TGFb_vs_NS<-TGFb_vs_NS[TGFb_vs_NS$log2FoldChange > 0.2,]
write.csv(TGFb_vs_NS, "/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/TGFb_vs_NS.csv")

IFNg_vs_NS <- lfcShrink(dds, coef="group_IFNg_vs_NS", type="apeglm")
IFNg_vs_NS<-as.data.frame(IFNg_vs_NS)
IFNg_vs_NS$gene<-rownames(IFNg_vs_NS)
IFNg_vs_NS<-IFNg_vs_NS[order(IFNg_vs_NS$padj),]
IFNg_vs_NS<-IFNg_vs_NS[IFNg_vs_NS$log2FoldChange > 0.2,]
write.csv(IFNg_vs_NS, "/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/IFNg_vs_NS.csv")

IL1b_vs_NS <- lfcShrink(dds, coef="group_IL1b_vs_NS", type="apeglm")
IL1b_vs_NS<-as.data.frame(IL1b_vs_NS)
IL1b_vs_NS$gene<-rownames(IL1b_vs_NS)
IL1b_vs_NS<-IL1b_vs_NS[order(IL1b_vs_NS$padj),]
IL1b_vs_NS<-IL1b_vs_NS[IL1b_vs_NS$log2FoldChange > 0.2,]
write.csv(IL1b_vs_NS, "/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/IL1b_vs_NS.csv")

