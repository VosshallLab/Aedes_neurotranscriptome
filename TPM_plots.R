libprop_SF_and_male <- libprop[libprop$condition %in% c("male","SF") & libprop$keep == 1,]

tpm_male_SF <- subset(tpm_all_with_annotation,select=c(row.names(libprop_SF_and_male),colnames(gene_annotations)))


### list of tissues of interest for the annotation plots

tissue_list <- c("Fe_Ov_SF","Fe_Pa_SF","Fe_Os_SF","Fe_Br_SF","Fe_An_SF","Fe_Rs_SF","Fe_FL_SF","Fe_ML_SF","FE_HL_SF","Fe_At_SF","Ma_Br","Ma_An","Ma_Rs","Ma_FL","Ma_ML","Ma_HL","Ma_At")


### function to find a specific tissue from the whole TPM matrix

subset_TPM <- function(tissue_to_find,tpm){
  
  selection <- grepl(tissue_to_find,colnames(tpm))
  tpm_sub <- tpm[,selection]
  temp_mean <- apply(tpm_sub,1,mean)
  temp_median <- apply(tpm_sub,1,median)
  tpm_sub$median <- temp_median
  tpm_sub$mean <- temp_mean
  tpm_sub
}

### create appropriate empty vectors for tissue means and medians

compiled_means <- gene_annotations
compiled_medians <- gene_annotations

### cycle through tissues and add means/medians to appropriate vectors

library(plyr)

for (tissue in tissue_list)
{
  tissue_select <- subset_TPM(tissue,tpm_male_SF)
  
  compiled_means$temp_mean <- tissue_select$mean
  compiled_medians$temp_median <- tissue_select$median
  
  compiled_means <- rename(compiled_means,replace = c("temp_mean" = tissue))
  compiled_medians <- rename(compiled_medians,replace = c("temp_median" = tissue))
  }

library(RColorBrewer)

genefamily_heatmap <- function(genefam,tpmMat) {
 
  family_subset <- tpmMat[tpmMat$gene.family %in% genefam,]
#  row.names(family_subset) <- family_subset$display.name
  disp_name <- family_subset$display.name
  family_subset <- subset(family_subset,select=tissue_list)
  column_labels <- c("Ovaries","Palp","Proboscis","Brain","Antenna","Rostrum","Forelegs","Midlegs","Hindlegs","Abdominal tip","Brain","Antenna","Rostrum","Forelegs","Midlegs","Hindlegs","Abdominal tip")
 
  
  
  heatmap.2(log10(as.matrix(family_subset)+1),dendrogram="none",Colv=FALSE,colsep=c(3,10),sepwidth=c(.5,.5),
            trace="none",labCol=column_labels,srtCol=45,density.info="none",labRow=disp_name,
            key.xlab="Log10(TPM)",key.title=NA,
           # col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
            col = colorRampPalette(rev(brewer.pal(9,"Blues")))(255) 
           )
}

genefamily_heatmap("OR",compiled_means)
