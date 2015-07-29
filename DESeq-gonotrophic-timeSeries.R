#### set up to look at gene expression in 3 tissues we have three time-points for: brain, antenna, hindlegs

## DESeq variables are:
# dds_hindlegs_gono
# dds_brain_gono
# dds_antenna_gono

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
          col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),dendrogram="none",
          labCol = c("non-blood-fed","blood-fed","gravid"),main="brain",trace="none")

#################################
# z-score clustering and heatmaps
########## antenna

expMat <- data.frame("sf" = resantenna_sfbf$baseMean,"bf"=resantenna_sfbf$baseMean * (2^resantenna_sfbf$log2FoldChange),"o"=resantenna_sfo$baseMean * (2^resantenna_sfo$log2FoldChange))
row.names(expMat) = row.names(reshindlegs_sfbf)

## z-score heatmap
z.expMat <- as.data.frame(t(apply((expMat[row.names(expMat) %in% modGenesantenna,]), MARGIN = 1, FUN = scale )))

heatmap.2(as.matrix(z.expMat),
          trace="none",symm=TRUE,Colv="none",
          col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),dendrogram="none",
          labCol = c("non-blood-fed","blood-fed","gravid"),main="antenna",trace="none")

#################################
# z-score clustering and heatmaps
########## hindlegs

expMat <- data.frame("sf" = reshindlegs_sfbf$baseMean,"bf"=reshindlegs_sfbf$baseMean * (2^reshindlegs_sfbf$log2FoldChange),"o"=reshindlegs_sfo$baseMean * (2^reshindlegs_sfo$log2FoldChange))
row.names(expMat) = row.names(reshindlegs_sfbf)

## z-score heatmap
z.expMat <- as.data.frame(t(apply((expMat[row.names(expMat) %in% modGeneshindlegs,]), MARGIN = 1, FUN = scale )))

heatmap.2(as.matrix(z.expMat),
          trace="none",symm=TRUE,Colv="none",
          col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),dendrogram="none",
          labCol = c("non-blood-fed","blood-fed","gravid"),main="hindlegs",trace="none")
