#### set up to look at gene expression in 3 tissues we have three time-points for: brain, antenna, hindlegs

## DESeq variables are:
# dds_hindlegs_gono
# dds_brain_gono
# dds_antenna_gono


# brain as example
resBrain_sfbf <- results(dds_brain_gono,contrast=c("condition","SF","BF"))
resBrain_sfo <- results(dds_brain_gono,contrast=c("condition","SF","O"))

brain_modulated <- vstMat_brain[resBrain_sfbf$padj < 0.01 | resBrain_sfo$padj < 0.01,]
brain_modulated <- brain_modulated[,c("Fe_Br_SF_5","Fe_Br_SF_6","Fe_Br_SF_7","Fe_Br_SF_8","Fe_Br_BF_1","Fe_Br_BF_2","Fe_Br_BF_3","Fe_Br_BF_4","Fe_Br_O_1","Fe_Br_O_2","Fe_Br_O_3","Fe_Br_O_4","Fe_Br_O_5")]

row.names(brain_modulated[!is.na(row.names(brain_modulated)),]) -> modGenesBrain

heatmap.2(log10(1+vstMat_brain[row.names(vstMat_brain) %in% modGenesBrain,]))
