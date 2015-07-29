# male/female fold change


tpm_dimorphic <- apply(compiled_means[,15:31],2,scale)
tpm_dimorphic <- as.data.frame(tpm_dimorphic)
tpm_dimorphic$internal.gene_id <- compiled_means$internal.gene_id
row.names(tpm_dimorphic) <- row.names(compiled_means)
tpm_dimorphic$FM <- apply(tpm_dimorphic[,1:10],1,mean) - apply(tpm_dimorphic[,c(11:12,14:17)],1,mean) 
tpm_dimorphic$diff <- apply(tpm_dimorphic[,1:10],1,sum) - apply(tpm_dimorphic[,c(11:12,14:17)],1,sum)
fold_change <- subset(tpm_dimorphic,select=c("internal.gene_id","diff"))
fold_change <- fold_change[order(fold_change$diff),]

male_genes <- head(fold_change$internal.gene_id, 25)
female_genes <- tail(fold_change$internal.gene_id,25)

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
dimorphism_sig[order(dimorphism_sig$log2FoldChange),] -> dimorphism_DESeq

male_specific <- head(dimorphism_DESeq,30)
female_specific <- tail(dimorphism_DESeq,30)
female_specific
