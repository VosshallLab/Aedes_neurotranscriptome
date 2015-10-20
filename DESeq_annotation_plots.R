##### define libraries for SF and male for comparisons

# only keep libraries that have a keep == 1 in the library key
libprop_SF_and_male <- libprop[libprop$condition %in% c("male","SF") & libprop$keep == 1,]

# add a column for tissue and condition to the matrix
libprop_SF_and_male$tissue_condition <- paste(libprop_SF_and_male$condition,libprop_SF_and_male$tissue,sep="-")

# get the raw counts so we can process in DESeq2
rawcounts_reordered_SF_and_male <- subset(rawcounts_with_annotation_reordered,select=row.names(libprop_SF_and_male))

########### take SF, Male libraries through DESeq2

dds_SF_male_combined <- DESeqDataSetFromMatrix(countData = rawcounts_reordered_SF_and_male,
                                      colData = libprop_SF_and_male,
                                      design = ~ type + tissue_condition )

### perform DESeq

dds_SF_male_combined <- DESeq(dds_SF_male_combined)

### perform variance stabilizing tranformation
vst_SF_male <- varianceStabilizingTransformation(dds_SF_male_combined,blind=FALSE)
vstMat_SF_male <- assay(vst_SF_male)

## PCA plots of vst transformed data

SF_male_PCA_vst <- PCA_tissue_condition(vst_SF_male); SF_male_PCA_vst + theme_bw()

## euclidean heatmap plot of vst transformed data

SF_male_heatmap_vst <- euclidean_heatmap(vst_SF_male,dds_SF_male_combined)

