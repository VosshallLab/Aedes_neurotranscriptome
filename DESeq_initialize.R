library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(gplots)

source('ntx_deseq_functions.R')

setwd('~/Dropbox/ntx_deseq/')

PE <- read.csv('PE_featurecounts.txt',sep = "\t",comment.char = "#")
SE <- read.csv('SE_featurecounts.txt',sep = "\t",comment.char = "#")

merge(PE,SE,by=c("Geneid","Chr","Start","End","Strand","Length")) -> test

columns <- read.csv('test.csv')
colnames(test) <- columns$y
row.names(test) <- test$gene
raw_FC <- test[,7:140]

##### read in library information from CSV file

libprop <- read.csv('libraries.csv')
row.names(libprop) <- libprop$library
raw_FC_reordered <- raw_FC[row.names(libprop)]


##### define libraries for SF and male for comparisons

libprop_SF <- libprop[libprop$condition %in% c("SF") & libprop$keep == 1,]
libprop_male <- libprop[libprop$condition %in% c("male") & libprop$keep == 1,]
libprop_SF_and_male <- libprop[libprop$condition %in% c("male","SF") & libprop$keep == 1,]

raw_FC_reordered_SF <- subset(raw_FC_reordered,select=row.names(libprop_SF))
raw_FC_reordered_male <- subset(raw_FC_reordered,select=row.names(libprop_male))
raw_FC_reordered_SF_and_male <- subset(raw_FC_reordered,select=row.names(libprop_SF_and_male))

############################# example gonotrophic analysis (brain)

# # single end libraries only, with keep flag == 1, and appropriate tissue
# libprop_brain_gono <- libprop[libprop$tissue %in% c("brain") & libprop$keep == 1 & libprop$type %in% c("single"),]   
# 
# # count matrix from selected libraries
# raw_FC_reordered_brain_gono <- subset(raw_FC_reordered,select=row.names(libprop_brain_gono))
# 
# # model design - condition is the variable
# dds_brain_gono <- DESeqDataSetFromMatrix(countData = raw_FC_reordered_brain_gono,
#                                          colData = libprop_brain_gono,
#                                          design = ~ condition
# )
# 
# # perform DESeq2
# dds_brain_gono <- DESeq(dds_brain_gono)
# res_brain_gono <- results(dds_brain_gono)
# rld_brain_gono <- rlog(dds_brain_gono)
# vsd_brain_gono <- varianceStabilizingTransformation(dds_brain_gono)
# rlogMat_brain <- assay(rld_brain_gono)
# vstMat_brain <- assay(vsd_brain_gono)
# 
# ## PCA plot of brain samples
# 
# data_brain_gono <- plotPCA(vsd_brain_gono, intgroup=c("condition", "type"), returnData=TRUE)
# percentVar <- round(100 * attr(data_brain_gono, "percentVar"))
# ggplot(data_brain_gono, aes(PC1, PC2, color=condition)) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance"))
# 
# ## euclidian distance heatmap plot
# 
# hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
# distsRL <- dist(t(assay(vsd_brain_gono)))
# mat <- as.matrix(distsRL)
# rownames(mat) <- colnames(mat) <- with(colData(dds_brain_gono), paste(condition, type, sep="-"))
# hc <- hclust(distsRL)
# heatmap.2(mat,Rowv=as.dendrogram(hc),symm=TRUE,trace="none",col=rev(hmcol),margin=c(13,13))
# 
# 


############################ end of brain example


########### SF, Male, SF_Male

dds_SF <- DESeqDataSetFromMatrix(countData = raw_FC_reordered_SF,
                              colData = libprop_SF,
                              design = ~ type + tissue
                              )
dds_male <- DESeqDataSetFromMatrix(countData = raw_FC_reordered_male,
                              colData = libprop_male,
                              design = ~ type + tissue
                                  )

dds_SF_male <- DESeqDataSetFromMatrix(countData = raw_FC_reordered_SF_and_male,
                                                colData = libprop_SF_and_male,
                                                design = ~ tissue + condition + type
                                      )

dds_male <- DESeqDataSetFromMatrix(countData = raw_FC_reordered_male,
                                   colData = libprop_male,
                                   design = ~ type + tissue
)



dds_SF <- DESeq(dds_SF)
dds_male <- DESeq(dds_male)
dds_SF_male <- DESeq(dds_SF_male)


res_SF <- results(dds_SF)
res_male <- results(dds_male)
res_SF_male <- results(dds_SF_male)

rld_male <- rlog(dds_male)
vsd_male <- varianceStabilizingTransformation(dds_male)
rlogMat_male <- assay(rld_male)
vstMat_male <- assay(vsd_male)

rld_SF <- rlog(dds_SF)
vsd_SF <- varianceStabilizingTransformation(dds_SF)
rlogMat_SF <- assay(rld_SF)
vstMat_SF <- assay(vsd_SF)

rld_SF_male <- rlog(dds_SF_male)
vsd_SF_male <- varianceStabilizingTransformation(dds_SF_male)
rlogMat_SF_male <- assay(rld_SF_male)
vstMat_SF_male <- assay(vsd_SF_male)

## PCA plots of rld and vsd transformed data

data <- plotPCA(rld_SF_male, intgroup=c("tissue", "type","condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=tissue, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

data_vsd <- plotPCA(vsd_SF_male, intgroup=c("tissue", "type","condition"), returnData=TRUE)
percentVar <- round(100 * attr(data_vsd, "percentVar"))
ggplot(data_vsd, aes(PC1, PC2, color=tissue, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))


## heatmap of sample-sample euclidian distances

hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
distsRL <- dist(t(assay(vsd_SF)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, type, sep="-"))
hc <- hclust(distsRL)
heatmap.2(mat,Rowv=as.dendrogram(hc),symm=TRUE,trace="none",col=rev(hmcol),margin=c(13,13))


