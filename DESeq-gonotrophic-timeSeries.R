#### set up to look at gene expression in 3 tissues we have three time-points for: brain, antenna, hindlegs

source('~/bioinfo/github/ntx_deseq/ntx_deseq_functions.R')

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

##### re-do Z-score based on TPM

### antenna
cond_list <- c("SF","BF","O")
antenna.tpm <- apply(counts(dds_antenna_gono),2,countToTpm,gene_annotations$Length)
antenna_gono_mean <- gene_annotations
### cycle through tissues and add means/medians to appropriate vectors

library(plyr)

for (cond in cond_list)
{
  cond_select <- subset_TPM(cond,antenna.tpm)
  
  antenna_gono_mean$temp_mean <- cond_select$mean

  antenna_gono_mean <- rename(antenna_gono_mean,replace = c("temp_mean" = cond))
}

tpm.Z <- as.data.frame(t(apply((antenna_gono_mean[,c("SF","BF","O")][row.names(antenna_gono_mean[,c("SF","BF","O")]) %in% modGenesantenna,]), MARGIN = 1, FUN = scale )))
colnames(tpm.Z) <- c("SF","BF","O")

heatmap_with_cluster(tpm.Z,'antenna',8)

### brain
cond_list <- c("SF","BF","O")
tpm_brain_gono <- apply(counts(dds_brain_gono),2,countToTpm,gene_annotations$Length)
brain_gono_mean <- gene_annotations
### cycle through tissues and add means/medians to appropriate vectors

library(plyr)

for (cond in cond_list)
{
  cond_select <- subset_TPM(cond,tpm_brain_gono)
  
  brain_gono_mean$temp_mean <- cond_select$mean
  
  brain_gono_mean <- rename(brain_gono_mean,replace = c("temp_mean" = cond))
}

tpm.Z <- as.data.frame(t(apply((brain_gono_mean[,c("SF","BF","O")][row.names(brain_gono_mean[,c("SF","BF","O")]) %in% modGenesbrain,]), MARGIN = 1, FUN = scale )))
colnames(tpm.Z) <- c("SF","BF","O")

heatmap_with_cluster(tpm.Z,'brain',8)

### hindlegs
cond_list <- c("SF","BF","O")
tpm_hindlegs_gono <- apply(counts(dds_hindlegs_gono),2,countToTpm,gene_annotations$Length)
hindlegs_gono_mean <- gene_annotations
### cycle through tissues and add means/medians to appropriate vectors

library(plyr)

for (cond in cond_list)
{
  cond_select <- subset_TPM(cond,tpm_hindlegs_gono)
  
  hindlegs_gono_mean$temp_mean <- cond_select$mean
  
  hindlegs_gono_mean <- rename(hindlegs_gono_mean,replace = c("temp_mean" = cond))
}



tpm.Z <- as.data.frame(t(apply((hindlegs_gono_mean[,c("SF","BF","O")][row.names(hindlegs_gono_mean[,c("SF","BF","O")]) %in% modGeneshindlegs,]), MARGIN = 1, FUN = scale )))
colnames(tpm.Z) <- c("SF","BF","O")

heatmap_with_cluster(tpm.Z,'hindlegs',8)
