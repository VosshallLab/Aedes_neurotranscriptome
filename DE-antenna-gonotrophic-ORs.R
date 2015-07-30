#### antenna OR gonotrophic

reg_ORs <- c("gene15193","gene2976","gene11112","gene13247","gene2482","gene3236","gene8278","gene7459","gene9065","gene14816","gene15909","gene14944","gene4613","gene5735","gene15967","gene14494")

antenna_gono_TPM <- countToTpm(counts(dds_antenna_gono),gene_annotations$Length)

heatmap.2(log10(1+as.matrix(antenna_gono_TPM[row.names(antenna_gono_TPM) %in% reg_ORs,])),Colv="none")




