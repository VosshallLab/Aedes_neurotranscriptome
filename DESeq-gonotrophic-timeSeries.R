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


#################################
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

antenna_reg_genes <- heatmap_with_cluster(tpm.Z,'antenna',8)

antenna_up_early <- subset(antenna_reg_genes,gr.row==1)
antenna_up_late <- subset(antenna_reg_genes,gr.row==3)
antenna_down_early <- subset(antenna_reg_genes,gr.row==4)
antenna_down_late <- subset(antenna_reg_genes,gr.row==7)

antenna_up_transient <- subset(antenna_reg_genes,gr.row %in% c(2,5))
antenna_down_transient <- subset(antenna_reg_genes,gr.row %in% c(6,8))

all_interesting_antenna <- subset(antenna_reg_genes,gr.row %in% c(1,3,4,7))

#heatmap_with_cluster(log10(1+subset(antenna_gono_mean[row.names(antenna_gono_mean) %in% row.names(all_interesting_antenna),],select=c("SF","BF","O"))),'antenna',4)
tpm.Z.interesting <- as.data.frame(t(apply((antenna_gono_mean[,c("SF","BF","O")][row.names(antenna_gono_mean[,c("SF","BF","O")]) %in% all_interesting_antenna$internal.gene_id,]), MARGIN = 1, FUN = scale )))
colnames(tpm.Z.interesting) <- c("SF","BF","O")
heatmap_with_cluster(tpm.Z.interesting,'antenna',4)


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

brain_reg_genes <- heatmap_with_cluster(tpm.Z,'brain',8)

brain_up_early <- subset(brain_reg_genes,gr.row %in% c(8,7))
brain_up_late <- subset(brain_reg_genes,gr.row==6)
brain_down_early <- subset(brain_reg_genes,gr.row==5)
brain_down_late <- subset(brain_reg_genes,gr.row==3)

brain_up_transient <- subset(brain_reg_genes,gr.row %in% c(1))
brain_down_transient <- subset(brain_reg_genes,gr.row %in% c(2,4))

all_interesting_brain <- subset(brain_reg_genes,gr.row %in% c(3,8,5,6,7))

#heatmap_with_cluster(log10(1+subset(brain_gono_mean[row.names(brain_gono_mean) %in% row.names(all_interesting_brain),],select=c("SF","BF","O"))),'brain',4)
tpm.Z.interesting <- as.data.frame(t(apply((brain_gono_mean[,c("SF","BF","O")][row.names(brain_gono_mean[,c("SF","BF","O")]) %in% all_interesting_brain$internal.gene_id,]), MARGIN = 1, FUN = scale )))
colnames(tpm.Z.interesting) <- c("SF","BF","O")
heatmap_with_cluster(tpm.Z.interesting,'brain',4)





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

hindlegs_reg_genes <- heatmap_with_cluster(tpm.Z,'hindlegs',8)

hindlegs_up_early <- subset(hindlegs_reg_genes,gr.row ==8)
hindlegs_up_late <- subset(hindlegs_reg_genes,gr.row==6)
hindlegs_down_early <- subset(hindlegs_reg_genes,gr.row==3)
hindlegs_down_late <- subset(hindlegs_reg_genes,gr.row==7)

hindlegs_up_transient <- subset(hindlegs_reg_genes,gr.row %in% c(1,4,5))
hindlegs_down_transient <- subset(hindlegs_reg_genes,gr.row==2)

all_interesting_hindlegs <- subset(hindlegs_reg_genes,gr.row %in% c(3,6,7,8))

#heatmap_with_cluster(log10(1+subset(hindlegs_gono_mean[row.names(hindlegs_gono_mean) %in% row.names(all_interesting_hindlegs),],select=c("SF","BF","O"))),'hindlegs',4)
tpm.Z.interesting <- as.data.frame(t(apply((hindlegs_gono_mean[,c("SF","BF","O")][row.names(hindlegs_gono_mean[,c("SF","BF","O")]) %in% all_interesting_hindlegs$internal.gene_id,]), MARGIN = 1, FUN = scale )))
colnames(tpm.Z.interesting) <- c("SF","BF","O")
heatmap_with_cluster(tpm.Z.interesting,'hindlegs',4)

##############
antenna_gene4252.tpm <- antenna.tpm["gene4252",]



## OR TPM plots
orlist <- c("gene13247","gene14944","gene15193","gene15909","gene2482","gene2976","gene3236","gene4613","gene5735")
or_tpm <- melt(antenna.tpm[which(row.names(antenna.tpm) %in% orlist),] )
or_tpm$condition <- c(rep("BF",36),rep("O",27),rep("SF",36))
or_tpm$condition_ordered <- factor(or_tpm$condition, levels = c("SF","BF","O"))
or_tpm$Var1 <- factor(or_tpm$Var1)
or_tpm.plot <- ggplot(data=or_tpm,aes(x=condition_ordered,y=value,colour=Var1)) + geom_boxplot() + facet_wrap(~Var1, scales="free")


#### OR by theme
# up early

orlist <- c("gene13247","gene14944","gene15193","gene15909","gene2482","gene2976","gene3236","gene4613","gene5735")
or_list_up <- c("gene13247","gene15193","gene2482","gene2976","gene3236")
or_list_uplate <- c("gene14944","gene4613","gene5735")
or_list_down <- c("gene15909")
or_tpm <- melt(antenna.tpm[which(row.names(antenna.tpm) %in% orlist),] )
or_tpm$condition <- c(rep("BF",36),rep("O",27),rep("SF",36))
or_tpm$condition_ordered <- factor(or_tpm$condition, levels = c("SF","BF","O"))
or_tpm$Var1 <- factor(or_tpm$Var1,levels = c(or_list_down,or_list_uplate,or_list_up))
or_tpm.plot <- ggplot(data=or_tpm,aes(x=condition_ordered,y=value)) + geom_boxplot() + geom_point() + facet_wrap(~Var1, scales="free")



## OR4 and OR5 outliers
orlist <- c("gene7459","gene9065")
or_tpm <- melt(antenna.tpm[which(row.names(antenna.tpm) %in% orlist),] )
or_tpm$condition <- c(rep("BF",8),rep("O",6),rep("SF",8))
or_tpm$condition_ordered <- factor(or_tpm$condition, levels = c("SF","BF","O"))
or_tpm.plot <- ggplot(data=or_tpm,aes(x=condition_ordered,y=value)) +geom_boxplot() + geom_point() + facet_wrap(~Var1,scales = "free")

