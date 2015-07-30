#### set up to look at gene expression in 3 tissues we have three time-points for: brain, antenna, hindlegs

## DESeq variables are:
# dds_hindlegs_gono
# dds_brain_gono
# dds_antenna_gono

####################
# brain as example
resovaries_sfo <- results(dds_ovaries_gono,contrast=c("condition","SF","O"))

ovaries_modulated <- vstMat_brain[resovaries_sfo$padj < 0.01,]
row.names(ovaries_modulated[!is.na(row.names(ovaries_modulated)),]) -> modGenesovaries

ovaries_FCs <- data.frame("row.names" = row.names(resovaries_sfo), "sfo" = resovaries_sfo$log2FoldChange)


####################
# brain as example
resbrain_sfbf <- results(dds_brain_gono,contrast=c("condition","SF","BF"))
resbrain_sfo <- results(dds_brain_gono,contrast=c("condition","SF","O"))

brain_modulated <- vstMat_brain[resbrain_sfbf$padj < 0.01 | resbrain_sfo$padj < 0.01,]
row.names(brain_modulated[!is.na(row.names(brain_modulated)),]) -> modGenesbrain

brain_FCs <- data.frame("row.names" = row.names(resbrain_sfbf), "sfbf" = resbrain_sfbf$log2FoldChange, "sfo" = resbrain_sfo$log2FoldChange)

#heatmap.2(as.matrix(brain_FCs[row.names(brain_FCs) %in% modGenesbrain,]),trace="none",col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

#################
## antenna
resantenna_sfbf <- results(dds_antenna_gono,contrast=c("condition","SF","BF"))
resantenna_sfo <- results(dds_antenna_gono,contrast=c("condition","SF","O"))


antenna_modulated <- vstMat_antenna[resantenna_sfbf$padj < 0.01 | resantenna_sfo$padj < 0.01,]
row.names(antenna_modulated[!is.na(row.names(antenna_modulated)),]) -> modGenesantenna

antenna_FCs <- data.frame("row.names" = row.names(resantenna_sfbf), "sfbf" = resantenna_sfbf$log2FoldChange, "sfo" = resantenna_sfo$log2FoldChange)

#heatmap.2(as.matrix(antenna_FCs[row.names(antenna_FCs) %in% modGenesantenna,]),trace="none",breaks=seq(-2,2,length.out=256),col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

#################
## hindlegs
reshindlegs_sfbf <- results(dds_hindlegs_gono,contrast=c("condition","SF","BF"))
reshindlegs_sfo <- results(dds_hindlegs_gono,contrast=c("condition","SF","O"))


hindlegs_modulated <- vstMat_hindlegs[reshindlegs_sfbf$padj < 0.01 | reshindlegs_sfo$padj < 0.01,]
row.names(hindlegs_modulated[!is.na(row.names(hindlegs_modulated)),]) -> modGeneshindlegs

hindlegs_FCs <- data.frame("row.names" = row.names(reshindlegs_sfbf), "sfbf" = reshindlegs_sfbf$log2FoldChange, "sfo" = reshindlegs_sfo$log2FoldChange)

#heatmap.2(as.matrix(hindlegs_FCs[row.names(hindlegs_FCs) %in% modGeneshindlegs,]),trace="none",breaks=seq(-2,2,length.out=256),col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

#################################
# z-score clustering and heatmaps
########## brain

expMat <- data.frame("sf" = resbrain_sfbf$baseMean,"bf"=resbrain_sfbf$baseMean * (2^resbrain_sfbf$log2FoldChange),"o"=resbrain_sfo$baseMean * (2^resbrain_sfo$log2FoldChange))
row.names(expMat) = row.names(reshindlegs_sfbf)

## z-score heatmap
z.expMat <- as.data.frame(t(apply((expMat[row.names(expMat) %in% modGenesbrain,]), MARGIN = 1, FUN = scale )))

heatmap.2(as.matrix(z.expMat),
          trace="none",symm=TRUE,Colv="none",
          col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),dendrogram="row",
          labCol = c("non-blood-fed","blood-fed","gravid"),main="brain")

#################################
# z-score clustering and heatmaps
########## antenna

expMat <- data.frame("sf" = resantenna_sfbf$baseMean,"bf"=resantenna_sfbf$baseMean * (2^resantenna_sfbf$log2FoldChange),"o"=resantenna_sfo$baseMean * (2^resantenna_sfo$log2FoldChange))
row.names(expMat) = row.names(reshindlegs_sfbf)

## z-score heatmap
z.expMat <- as.data.frame(t(apply((expMat[row.names(expMat) %in% modGenesantenna,]), MARGIN = 1, FUN = scale )))

heatmap.2(as.matrix(z.expMat),
          trace="none",symm=TRUE,Colv="none",
          col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),dendrogram="row",
          labCol = c("non-blood-fed","blood-fed","gravid"),main="antenna")

#################################
# z-score clustering and heatmaps
########## hindlegs

expMat <- data.frame("sf" = reshindlegs_sfbf$baseMean,"bf"=reshindlegs_sfbf$baseMean * (2^reshindlegs_sfbf$log2FoldChange),"o"=reshindlegs_sfo$baseMean * (2^reshindlegs_sfo$log2FoldChange))
row.names(expMat) = row.names(reshindlegs_sfbf)

## z-score heatmap
z.expMat <- as.data.frame(t(apply((expMat[row.names(expMat) %in% modGeneshindlegs,]), MARGIN = 1, FUN = scale )))

heatmap.2(as.matrix(z.expMat),
          trace="none",symm=TRUE,Colv="none",
          col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),dendrogram="row",
          labCol = c("non-blood-fed","blood-fed","gravid"),main="hindlegs")

#################################
# z-score clustering and heatmaps
########## ovaries

# expMat <- data.frame("sf" = resovaries_sfo$baseMean,"o"=resovaries_sfo$baseMean * (2^resovaries_sfo$log2FoldChange))
# row.names(expMat) = row.names(resovaries_sfo)
# 
# ## z-score heatmap
# z.expMat <- as.data.frame(t(apply((expMat[row.names(expMat) %in% modGenesovaries,]), MARGIN = 1, FUN = scale )))
# 
# heatmap.2(as.matrix(z.expMat),
#           trace="none",symm=TRUE,Colv="none",
#           col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),dendrogram="none",
#           labCol = c("non-blood-fed","gravid"),main="ovaries")



#### trying to cluster as well  - http://stackoverflow.com/questions/22278508/how-to-add-colsidecolors-on-heatmap-2-after-performing-bi-clustering-row-and-co

hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

mydata = z.expMat

cl.row <- hclustfunc(distfunc(mydata))

# extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
gr.row <- cutree(cl.row, 7)

# require(RColorBrewer)
col1 <- brewer.pal(7, "Set1")

# require(gplots)    
heatmap.2(as.matrix(mydata), hclustfun=hclustfunc, distfun=distfunc,   
          RowSideColors=col1[gr.row],trace="none",dendrogram="row",
          Rowv = reorder(as.dendrogram(cl.row),c(1,2,3,4,5,6,7)),
          col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          Colv = "none",labCol = c("non-blood-fed","blood-fed","gravid")
          )

names(gr.row[gr.row == 1]) -> clust1
names(gr.row[gr.row == 2]) -> clust2
names(gr.row[gr.row == 3]) -> clust3
names(gr.row[gr.row == 4]) -> clust4
names(gr.row[gr.row == 5]) -> clust5
names(gr.row[gr.row == 6]) -> clust6
names(gr.row[gr.row == 7]) -> clust7

vbase <- read.csv('~/Downloads/mart_export-2.txt',sep="\t")
description7 <- vbase[vbase$Gene.stable.ID %in% gene_annotations[gene_annotations$internal.gene_id %in% clust7,]$vectorbase.RU,]
description6 <- vbase[vbase$Gene.stable.ID %in% gene_annotations[gene_annotations$internal.gene_id %in% clust6,]$vectorbase.RU,]
description5 <- vbase[vbase$Gene.stable.ID %in% gene_annotations[gene_annotations$internal.gene_id %in% clust5,]$vectorbase.RU,]
description4 <- vbase[vbase$Gene.stable.ID %in% gene_annotations[gene_annotations$internal.gene_id %in% clust4,]$vectorbase.RU,]
description3 <- vbase[vbase$Gene.stable.ID %in% gene_annotations[gene_annotations$internal.gene_id %in% clust3,]$vectorbase.RU,]
description2 <- vbase[vbase$Gene.stable.ID %in% gene_annotations[gene_annotations$internal.gene_id %in% clust2,]$vectorbase.RU,]
description1 <- vbase[vbase$Gene.stable.ID %in% gene_annotations[gene_annotations$internal.gene_id %in% clust1,]$vectorbase.RU,]

