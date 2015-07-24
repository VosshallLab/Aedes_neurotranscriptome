###### DESeq_gonotrophic.R

##### hindlegs gonotrophic
libprop_hindlegs_gono <- libprop[libprop$tissue %in% c("hindlegs") & libprop$keep == 1 & libprop$type %in% c("single"),]   
# 
# # count matrix from selected libraries
raw_FC_reordered_hindlegs_gono <- subset(raw_FC_reordered,select=row.names(libprop_hindlegs_gono))
# 
# # model design - condition is the variable
dds_hindlegs_gono <- DESeqDataSetFromMatrix(countData = raw_FC_reordered_hindlegs_gono, colData = libprop_hindlegs_gono, design = ~ condition)
# 
# # perform DESeq2
dds_hindlegs_gono <- DESeq(dds_hindlegs_gono)
res_hindlegs_gono <- results(dds_hindlegs_gono)
rld_hindlegs_gono <- rlog(dds_hindlegs_gono)
vsd_hindlegs_gono <- varianceStabilizingTransformation(dds_hindlegs_gono)
rlogMat_hindlegs <- assay(rld_hindlegs_gono)
vstMat_hindlegs <- assay(vsd_hindlegs_gono)

gono_hindlegs_PCA_vsd <- PCA_tissue_condition(vsd_hindlegs_gono)
gono_hindlegs_PCA_rld <- PCA_tissue_condition(rld_hindlegs_gono)
############


############################# example gonotrophic analysis (brain)

# # single end libraries only, with keep flag == 1, and appropriate tissue
libprop_brain_gono <- libprop[libprop$tissue %in% c("brain") & libprop$keep == 1 & libprop$type %in% c("single"),]   
# 
# # count matrix from selected libraries
raw_FC_reordered_brain_gono <- subset(raw_FC_reordered,select=row.names(libprop_brain_gono))
# 
# # model design - condition is the variable
dds_brain_gono <- DESeqDataSetFromMatrix(countData = raw_FC_reordered_brain_gono, colData = libprop_brain_gono, design = ~ condition)
# 
# # perform DESeq2
dds_brain_gono <- DESeq(dds_brain_gono)
res_brain_gono <- results(dds_brain_gono)
rld_brain_gono <- rlog(dds_brain_gono)
vsd_brain_gono <- varianceStabilizingTransformation(dds_brain_gono)
rlogMat_brain <- assay(rld_brain_gono)
vstMat_brain <- assay(vsd_brain_gono)

gono_brain_PCA_vsd <- PCA_tissue_condition(vsd_brain_gono)
gono_brain_PCA_rld <- PCA_tissue_condition(rld_brain_gono)

############################ end of brain example

#### antenna gono
libprop_antenna_gono <- libprop[libprop$tissue %in% c("antenna") & libprop$keep == 1 & libprop$type %in% c("single"),]   
# 
# # count matrix from selected libraries
raw_FC_reordered_antenna_gono <- subset(raw_FC_reordered,select=row.names(libprop_antenna_gono))
# 
# # model design - condition is the variable
dds_antenna_gono <- DESeqDataSetFromMatrix(countData = raw_FC_reordered_antenna_gono, colData = libprop_antenna_gono, design = ~ condition)
# 
# # perform DESeq2
dds_antenna_gono <- DESeq(dds_antenna_gono)
res_antenna_gono <- results(dds_antenna_gono)
rld_antenna_gono <- rlog(dds_antenna_gono)
vsd_antenna_gono <- varianceStabilizingTransformation(dds_antenna_gono)
rlogMat_antenna <- assay(rld_antenna_gono)
vstMat_antenna <- assay(vsd_antenna_gono)

gono_antenna_PCA_vsd <- PCA_tissue_condition(vsd_antenna_gono)
gono_antenna_PCA_rld <- PCA_tissue_condition(rld_antenna_gono)
##################


###### forelegs

libprop_forelegs_gono <- libprop[libprop$tissue %in% c("forelegs") & libprop$keep == 1 & libprop$type %in% c("single"),]   
# 
# # count matrix from selected libraries
raw_FC_reordered_forelegs_gono <- subset(raw_FC_reordered,select=row.names(libprop_forelegs_gono))
# 
# # model design - condition is the variable
dds_forelegs_gono <- DESeqDataSetFromMatrix(countData = raw_FC_reordered_forelegs_gono, colData = libprop_forelegs_gono, design = ~ condition)
# 
# # perform DESeq2
dds_forelegs_gono <- DESeq(dds_forelegs_gono)
res_forelegs_gono <- results(dds_forelegs_gono)
rld_forelegs_gono <- rlog(dds_forelegs_gono)
vsd_forelegs_gono <- varianceStabilizingTransformation(dds_forelegs_gono)
rlogMat_forelegs <- assay(rld_forelegs_gono)
vstMat_forelegs <- assay(vsd_forelegs_gono)

gono_forelegs_PCA_vsd <- PCA_tissue_condition(vsd_forelegs_gono)
gono_forelegs_PCA_rld <- PCA_tissue_condition(rld_forelegs_gono)

#######################


#### midlegs

libprop_midlegs_gono <- libprop[libprop$tissue %in% c("midlegs") & libprop$keep == 1 & libprop$type %in% c("single"),]   
# 
# # count matrix from selected libraries
raw_FC_reordered_midlegs_gono <- subset(raw_FC_reordered,select=row.names(libprop_midlegs_gono))
# 
# # model design - condition is the variable
dds_midlegs_gono <- DESeqDataSetFromMatrix(countData = raw_FC_reordered_midlegs_gono, colData = libprop_midlegs_gono, design = ~ condition)
# 
# # perform DESeq2
dds_midlegs_gono <- DESeq(dds_midlegs_gono)
res_midlegs_gono <- results(dds_midlegs_gono)
rld_midlegs_gono <- rlog(dds_midlegs_gono)
vsd_midlegs_gono <- varianceStabilizingTransformation(dds_midlegs_gono)
rlogMat_midlegs <- assay(rld_midlegs_gono)
vstMat_midlegs <- assay(vsd_midlegs_gono)

gono_midlegs_PCA_vsd <- PCA_tissue_condition(vsd_midlegs_gono)
gono_midlegs_PCA_rld <- PCA_tissue_condition(rld_midlegs_gono)
#########