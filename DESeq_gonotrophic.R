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
vst_hindlegs_gono <- varianceStabilizingTransformation(dds_hindlegs_gono)
rlogMat_hindlegs <- assay(rld_hindlegs_gono)
vstMat_hindlegs <- assay(vst_hindlegs_gono)

gono_hindlegs_PCA_vst <- PCA_tissue_condition(vst_hindlegs_gono)
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
vst_brain_gono <- varianceStabilizingTransformation(dds_brain_gono)
rlogMat_brain <- assay(rld_brain_gono)
vstMat_brain <- assay(vst_brain_gono)

gono_brain_PCA_vst <- PCA_tissue_condition(vst_brain_gono)
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
vst_antenna_gono <- varianceStabilizingTransformation(dds_antenna_gono)
rlogMat_antenna <- assay(rld_antenna_gono)
vstMat_antenna <- assay(vst_antenna_gono)

gono_antenna_PCA_vst <- PCA_tissue_condition(vst_antenna_gono)
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
vst_forelegs_gono <- varianceStabilizingTransformation(dds_forelegs_gono)
rlogMat_forelegs <- assay(rld_forelegs_gono)
vstMat_forelegs <- assay(vst_forelegs_gono)

gono_forelegs_PCA_vst <- PCA_tissue_condition(vst_forelegs_gono)
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
vst_midlegs_gono <- varianceStabilizingTransformation(dds_midlegs_gono)
rlogMat_midlegs <- assay(rld_midlegs_gono)
vstMat_midlegs <- assay(vst_midlegs_gono)

gono_midlegs_PCA_vst <- PCA_tissue_condition(vst_midlegs_gono)
gono_midlegs_PCA_rld <- PCA_tissue_condition(rld_midlegs_gono)
#########


## rostrum
libprop_rostrum_gono <- libprop[libprop$tissue %in% c("rostrum") & libprop$keep == 1 & libprop$type %in% c("single"),]   
# 
# # count matrix from selected libraries
raw_FC_reordered_rostrum_gono <- subset(raw_FC_reordered,select=row.names(libprop_rostrum_gono))
# 
# # model design - condition is the variable
dds_rostrum_gono <- DESeqDataSetFromMatrix(countData = raw_FC_reordered_rostrum_gono, colData = libprop_rostrum_gono, design = ~ condition)
# 
# # perform DESeq2
dds_rostrum_gono <- DESeq(dds_rostrum_gono)
res_rostrum_gono <- results(dds_rostrum_gono)
rld_rostrum_gono <- rlog(dds_rostrum_gono)
vst_rostrum_gono <- varianceStabilizingTransformation(dds_rostrum_gono)
rlogMat_rostrum <- assay(rld_rostrum_gono)
vstMat_rostrum <- assay(vst_rostrum_gono)

gono_rostrum_PCA_vst <- PCA_tissue_condition(vst_rostrum_gono)
gono_rostrum_PCA_rld <- PCA_tissue_condition(rld_rostrum_gono)

########

###### ovaries

libprop_ovaries_gono <- libprop[libprop$tissue %in% c("ovaries") & libprop$keep == 1 & libprop$type %in% c("paired"),]   
# 
# # count matrix from selected libraries
raw_FC_reordered_ovaries_gono <- subset(raw_FC_reordered,select=row.names(libprop_ovaries_gono))
# 
# # model design - condition is the variable
dds_ovaries_gono <- DESeqDataSetFromMatrix(countData = raw_FC_reordered_ovaries_gono, colData = libprop_ovaries_gono, design = ~ condition)
# 
# # perform DESeq2
dds_ovaries_gono <- DESeq(dds_ovaries_gono)
res_ovaries_gono <- results(dds_ovaries_gono)
rld_ovaries_gono <- rlog(dds_ovaries_gono)
vst_ovaries_gono <- varianceStabilizingTransformation(dds_ovaries_gono)
rlogMat_ovaries <- assay(rld_ovaries_gono)
vstMat_ovaries <- assay(vst_ovaries_gono)

gono_ovaries_PCA_vst <- PCA_tissue_condition(vst_ovaries_gono)
gono_ovaries_PCA_rld <- PCA_tissue_condition(rld_ovaries_gono)

##########

##### abd.tip

libprop_abdominaltip_gono <- libprop[libprop$tissue %in% c("abdominaltip") & libprop$keep == 1 & libprop$type %in% c("single"),]   
# 
# # count matrix from selected libraries
raw_FC_reordered_abdominaltip_gono <- subset(raw_FC_reordered,select=row.names(libprop_abdominaltip_gono))
# 
# # model design - condition is the variable
dds_abdominaltip_gono <- DESeqDataSetFromMatrix(countData = raw_FC_reordered_abdominaltip_gono, colData = libprop_abdominaltip_gono, design = ~ condition)
# 
# # perform DESeq2
dds_abdominaltip_gono <- DESeq(dds_abdominaltip_gono)
res_abdominaltip_gono <- results(dds_abdominaltip_gono)
rld_abdominaltip_gono <- rlog(dds_abdominaltip_gono)
vst_abdominaltip_gono <- varianceStabilizingTransformation(dds_abdominaltip_gono)
rlogMat_abdominaltip <- assay(rld_abdominaltip_gono)
vstMat_abdominaltip <- assay(vst_abdominaltip_gono)

gono_abdominaltip_PCA_vst <- PCA_tissue_condition(vst_abdominaltip_gono)
gono_abdominaltip_PCA_rld <- PCA_tissue_condition(rld_abdominaltip_gono)