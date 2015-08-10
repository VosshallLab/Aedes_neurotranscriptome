# male/female fold change calculations
library(plyr)

#write.csv(compiled_means,'compiled_means.csv')

compiled_means <- read.csv('compiled_means.csv',header=TRUE,row.names = 1)

rm(tpm_dimorphic_tpm)
rm(tpm_dimorphic_tpm.Z)

tpm_dimorphic_tpm <- as.matrix(compiled_means[,c(16:27,29:32)])  
tpm_dimorphic_tpm.Z <- tpm_dimorphic_tpm

#tpm_dimorphic.Z <- apply(tpm_dimorphic_tpm,1,sum)  # no rostrum

 for (i in 1:length(row.names(tpm_dimorphic_tpm))) {
   
  tpm_dimorphic_tpm.Z[i,] <- scale(as.numeric(tpm_dimorphic_tpm[i,]))[,1]
   
 }

 euclid <- function(testRow,refRow) {
   
   sum(testRow - refRow)
   
 }

rm(tpm_dimorphic_tpm.Z.metric)
 
tpm_dimorphic_tpm.Z.metric[1:length(tpm_dimorphic_tpm.Z)] <- tpm_dimorphic_tpm.Z

tpm_male.sum.Z <- apply(tpm_dimorphic_tpm.Z[,10:15],1,sum)
tpm_female.sum.Z <- apply(tpm_dimorphic_tpm.Z[,1:9],1,sum)

tpm_female.sum <- apply(tpm_dimorphic_tpm[,1:9],1,sum)
tpm_male.sum <- apply(tpm_dimorphic_tpm[,10:15],1,sum)

temp.sum <- data.frame("male" = tpm_male.sum,"female" = tpm_female.sum)
temp.sum.Z <- data.frame("male" = tpm_male.sum.Z,"female" = tpm_female.sum.Z)

tpm_mean.Z <- data.frame("male" = tpm_male.Z,"female" = tpm_female.Z)
temp <- na.omit(tpm_mean.Z)

genes_temp <- row.names(temp[abs(temp$male - temp$female) > 1.5,])

dimorphism_heatmap(genes_temp,compiled_means,3)

#tpm_dimorphic.Z <- as.data.frame(tpm_dimorphic)
row.names(tpm_dimorphic) <- compiled_means$internal.gene_id

tpm_dimorphic.Z$FM <- apply(tpm_dimorphic.Z[,1:10],2,mean) /  apply(tpm_dimorphic.Z[,c(11:16)],2,mean) 
tpm_dimorphic.Z$diff <- apply(tpm_dimorphic.Z[,1:10],2,sum) / apply(tpm_dimorphic.Z[,c(11:16)],2,sum)

fold_change <- subset(tpm_dimorphic.Z,select=c("FM","diff"))
fold_change <- fold_change[order(fold_change$FM),]

male_genes <- head(row.names(fold_change), 25)
female_genes <- tail(row.names(fold_change),25)

tpm_male <- tpm_male_SF[tpm_male_SF$internal.gene_id %in% male_genes,]
tpm_female <- tpm_male_SF[tpm_male_SF$internal.gene_id %in% female_genes,]

mean_tpm_male <- compiled_means[row.names(compiled_means) %in% male_genes,]
mean_tpm_female <- compiled_means[row.names(compiled_means) %in% female_genes,]

dimorphism_heatmap(female_genes,compiled_means,5)



dds_SF_male_sex <- DESeqDataSetFromMatrix(countData = rawcounts_reordered_SF_and_male,
                                          colData = libprop_SF_and_male,
                                          design = ~ type + condition )

dds_SF_male_sex <- DESeq(dds_SF_male_sex)
dimorphism_DESeq <- results(dds_SF_male_sex,contrast=c("condition","SF","male"))

dimorphism_sig <- dimorphism_DESeq[!is.na(dimorphism_DESeq$padj) & dimorphism_DESeq$padj < 0.01,]
dimorphism_sig[order(dimorphism_sig$log2FoldChange),] -> dimorphism_sig

male_specific <- head(dimorphism_sig,30)
female_specific <- tail(dimorphism_sig,30)

dimorphism_heatmap(row.names(female_specific),compiled_means,5)
dimorphism_heatmap(row.names(male_specific),compiled_means,5)

