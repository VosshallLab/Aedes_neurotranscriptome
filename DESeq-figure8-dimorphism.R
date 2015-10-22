### sexual dimorphism plots, October 21, 2015

"gene16037" -> nix
"gene16140" -> myosex

FC_cutoff <- 3
padj_cutoff <- 0.01

antenna <- res_SF_male_antenna[!is.na(res_SF_male_antenna$padj) & res_SF_male_antenna$padj < padj_cutoff,]
brain <- res_SF_male_brain[!is.na(res_SF_male_brain$padj) & res_SF_male_brain$padj < padj_cutoff,]
hl <- res_SF_male_hindlegs[!is.na(res_SF_male_hindlegs$padj) & res_SF_male_hindlegs$padj < padj_cutoff,]
fl <- res_SF_male_forelegs[!is.na(res_SF_male_forelegs$padj) & res_SF_male_forelegs$padj < padj_cutoff,]
ml <- res_SF_male_midlegs[!is.na(res_SF_male_midlegs$padj) & res_SF_male_midlegs$padj < padj_cutoff,]
abdtip <- res_SF_male_abdominal_tip[!is.na(res_SF_male_abdominal_tip$padj) & res_SF_male_abdominal_tip$padj < padj_cutoff,]

dimorphic_fem <- list("antenna" = row.names(antenna[antenna$log2FoldChange >= FC_cutoff,]), 
                            "brain" = row.names(brain[brain$log2FoldChange >= FC_cutoff,]), 
                            "abdominal_tip" = row.names(abdtip[abdtip$log2FoldChange >= FC_cutoff,]),
                            "hindlegs" = row.names(hl[hl$log2FoldChange >= FC_cutoff,]),
                      "midlegs" = row.names(ml[ml$log2FoldChange >= FC_cutoff,]),
                      "forelegs" = row.names(fl[fl$log2FoldChange >= FC_cutoff,]))

dimorphic_male <- list("antenna" = row.names(antenna[antenna$log2FoldChange <= -1*FC_cutoff,]), 
                      "brain" = row.names(brain[brain$log2FoldChange <= -1*FC_cutoff,]), 
                      "abdominal_tip" = row.names(abdtip[abdtip$log2FoldChange <= -1*FC_cutoff,]),
                      "hindlegs" = row.names(hl[hl$log2FoldChange <= -1*FC_cutoff,]),
                      "forelegs" = row.names(fl[fl$log2FoldChange <= -1*FC_cutoff,]),
                      "midlegs" = row.names(ml[ml$log2FoldChange <= -1*FC_cutoff,])
                      )

melt.male <- melt(dimorphic_male)
melt.fem <- melt(dimorphic_fem)

melt.male.legs <- melt.male
melt.male.legs[melt.male.legs$L1 %in% c("forelegs","midlegs","hindlegs"),]$L1 <- "legs"

melt.fem.legs <- melt.fem
melt.fem.legs[melt.fem.legs$L1 %in% c("forelegs","midlegs","hindlegs"),]$L1 <- "legs"

melt.fem.legs <- unique(melt.fem.legs)
melt.male.legs <- unique(melt.male.legs)

rownames(melt.fem.legs) <- NULL
rownames(melt.male.legs) <- NULL


count(melt.male.legs$value)[count(melt.male.legs$value)$freq %in% c(3,4),] -> male_all
count(melt.fem.legs$value)[count(melt.fem.legs$value)$freq %in% c(3,4),] -> female_all

dimorphism_heatmap_fig8_M(gene_list = male_all$x,tpmMat = compiled_means,maxTPM = 2)
dimorphism_heatmap_fig8_F(gene_list = female_all$x,tpmMat = compiled_means,maxTPM = 4)

# generate venn diagrams (as TIFF files)

venn.diagram(unstack(melt.fem.legs),euler.d=TRUE,scaled=TRUE,'femalevenn.tiff',imagetype='tiff')
venn.diagram(unstack(melt.male.legs),euler.d=TRUE,scaled=TRUE,'malevenn.tiff',imagetype='tiff')


##### to export as supplementals

annotation_subset <- subset(gene_annotations,select=c("vectorbase.RU","display.name"))

write.csv(merge(annotation_subset,as.data.frame(antenna),by="row.names"),'sex_dimorphism_antenna.csv')
write.csv(merge(annotation_subset,as.data.frame(brain),by="row.names"),'sex_dimorphism_brain.csv')
write.csv(merge(annotation_subset,as.data.frame(hl),by="row.names"),'sex_dimorphism_hindlegs.csv')
write.csv(merge(annotation_subset,as.data.frame(fl),by="row.names"),'sex_dimorphism_forelegs.csv')
write.csv(merge(annotation_subset,as.data.frame(ml),by="row.names"),'sex_dimorphism_midlegs.csv')
write.csv(merge(annotation_subset,as.data.frame(abdtip),by="row.names"),'sex_dimorphism_abdominal_tip.csv')