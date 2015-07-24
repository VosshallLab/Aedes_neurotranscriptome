##### define libraries for SF and male for comparisons

# libprop_SF <- libprop[libprop$condition %in% c("SF") & libprop$keep == 1,]
# libprop_male <- libprop[libprop$condition %in% c("male") & libprop$keep == 1,]
libprop_SF_and_male <- libprop[libprop$condition %in% c("male","SF") & libprop$keep == 1,]

# raw_FC_reordered_SF <- subset(raw_FC_reordered,select=row.names(libprop_SF))
# raw_FC_reordered_male <- subset(raw_FC_reordered,select=row.names(libprop_male))
raw_FC_reordered_SF_and_male <- subset(raw_FC_reordered,select=row.names(libprop_SF_and_male))

########### SF, Male, SF_Male

# dds_SF <- DESeqDataSetFromMatrix(countData = raw_FC_reordered_SF,
#                               colData = libprop_SF,
#                               design = ~ type + tissue
#                               )
# 
# dds_male <- DESeqDataSetFromMatrix(countData = raw_FC_reordered_male,
#                               colData = libprop_male,
#                               design = ~ type + tissue
#                                   )

dds_SF_male <- DESeqDataSetFromMatrix(countData = raw_FC_reordered_SF_and_male,
                                      colData = libprop_SF_and_male,
                                      design = ~ tissue + condition + type
)



# dds_SF <- DESeq(dds_SF)
# dds_male <- DESeq(dds_male)
dds_SF_male <- DESeq(dds_SF_male)


# res_SF <- results(dds_SF)
# res_male <- results(dds_male)
res_SF_male <- results(dds_SF_male)
# 
# rld_male <- rlog(dds_male)
# vsd_male <- varianceStabilizingTransformation(dds_male)
# rlogMat_male <- assay(rld_male)
# vstMat_male <- assay(vsd_male)
# 
# rld_SF <- rlog(dds_SF)
# vsd_SF <- varianceStabilizingTransformation(dds_SF)
# rlogMat_SF <- assay(rld_SF)
# vstMat_SF <- assay(vsd_SF)

rld_SF_male <- rlog(dds_SF_male,blind=FALSE)
vsd_SF_male <- varianceStabilizingTransformation(dds_SF_male,blind=FALSE)
rlogMat_SF_male <- assay(rld_SF_male)
vstMat_SF_male <- assay(vsd_SF_male)

## PCA plots of rld and vsd transformed data

SF_male_PCA_vsd <- PCA_tissue_condition(vsd_SF_male)
SF_male_PCA_rld <- PCA_tissue_condition(rld_SF_male)

## euclidean heatmap plot of rld and vsd transformed data

SF_male_heatmap_vsd <- euclidean_heatmap(vsd_SF_male,dds_SF_male)
SF_male_heatmap_rld <- euclidean_heatmap(rld_SF_male,dds_SF_male)