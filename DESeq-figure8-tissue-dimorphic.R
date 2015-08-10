# for individual tissues
###########################################################

"gene16037" -> nix
"gene16140" -> myosex

FC_cutoff <- 4
padj_cutoff <- 0.1
TPM_cutoff <- 10
maxTPM <- 5

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

tissues <- c("antenna","brain","abdominal_tip","hindlegs","forelegs","midlegs")

melt.male <- melt(dimorphic_male)
melt.fem <- melt(dimorphic_fem)

antenna <- subset(merge(as.data.frame(antenna),compiled_means, by="row.names"),select=c("Fe_An_SF","Ma_An","log2FoldChange","internal.gene_id","vectorbase.RU","display.name"))
brain <- subset(merge(compiled_means,as.data.frame(brain), by="row.names"),select=c("Fe_Br_SF","Ma_Br","log2FoldChange","internal.gene_id","vectorbase.RU","display.name"))
hindlegs <- subset(merge(compiled_means,as.data.frame(hl), by="row.names"),select=c("FE_HL_SF","Ma_HL","log2FoldChange","internal.gene_id","vectorbase.RU","display.name"))
midlegs <- subset(merge(compiled_means,as.data.frame(ml), by="row.names"),select=c("Fe_ML_SF","Ma_ML","log2FoldChange","internal.gene_id","vectorbase.RU","display.name"))
forelegs <- subset(merge(compiled_means,as.data.frame(fl), by="row.names"),select=c("Fe_FL_SF","Ma_FL","log2FoldChange","internal.gene_id","vectorbase.RU","display.name"))
abdominal_tip <- subset(merge(compiled_means,as.data.frame(abdtip), by="row.names"),select=c("Fe_At_SF","Ma_At","log2FoldChange","internal.gene_id","vectorbase.RU","display.name"))



#### antenna
antenna.male <- antenna[(antenna$Ma_An > TPM_cutoff | antenna$Fe_An_SF > TPM_cutoff) & (antenna$log2FoldChange <= -1 * FC_cutoff | antenna$log2FoldChange >= FC_cutoff),]
antenna.male <- antenna.male[order(antenna.male$log2FoldChange),]
antenna.male <- subset(antenna.male,select=c("Fe_An_SF","Ma_An","vectorbase.RU","display.name"))
row.names(antenna.male) <- paste(antenna.male$vectorbase.RU,antenna.male$display.name)
antenna.male <- antenna.male[,1:2]
pdf('antenna_dimorphic.pdf',width=6)
heatmap.2(as.matrix(log10(1+antenna.male),trace="none",dendrogram="none",main="antenna"),Colv="none", breaks= seq(0,5,length.out=256),
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),trace="none")
dev.off()

##### brain
brain.male <- brain[(brain$Ma_Br > TPM_cutoff | brain$Fe_Br_SF > TPM_cutoff) & (brain$log2FoldChange <= -1 * FC_cutoff | brain$log2FoldChange >= FC_cutoff),]
brain.male <- brain.male[order(brain.male$log2FoldChange),]
brain.male <- subset(brain.male,select=c("Fe_Br_SF","Ma_Br","vectorbase.RU","display.name"))
row.names(brain.male) <- paste(brain.male$vectorbase.RU,brain.male$display.name)
brain.male <- brain.male[,1:2]
pdf('brain_dimorphic.pdf',width=6)
heatmap.2(as.matrix(log10(1+brain.male),trace="none",dendrogram="none",main="brain"),Colv="none", breaks= seq(0,3,length.out=256),
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),trace="none")
dev.off()

# ##### hindlegs
hindlegs.male <- hindlegs[(hindlegs$Ma_HL > TPM_cutoff | hindlegs$FE_HL_SF > TPM_cutoff) & (hindlegs$log2FoldChange <= -1 * FC_cutoff | hindlegs$log2FoldChange >= FC_cutoff),]
hindlegs.male <- hindlegs.male[order(hindlegs.male$log2FoldChange),]
hindlegs.male <- subset(hindlegs.male,select=c("FE_HL_SF","Ma_HL","vectorbase.RU","display.name"))
row.names(hindlegs.male) <- paste(hindlegs.male$vectorbase.RU,hindlegs.male$display.name)
hindlegs.male <- hindlegs.male[,1:2]
pdf('hindlegs_dimorphic.pdf',width=6)
heatmap.2(as.matrix(log10(1+hindlegs.male),trace="none",dendrogram="none",main="hindlegs"),Colv="none", breaks= seq(0,3,length.out=256),
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),trace="none")
dev.off()

# ##### forelegs
forelegs.male <- forelegs[(forelegs$Ma_FL > TPM_cutoff | forelegs$Fe_FL_SF > TPM_cutoff) & (forelegs$log2FoldChange <= -1 * FC_cutoff | forelegs$log2FoldChange >= FC_cutoff),]
forelegs.male <- forelegs.male[order(forelegs.male$log2FoldChange),]
forelegs.male <- subset(forelegs.male,select=c("Fe_FL_SF","Ma_FL","vectorbase.RU","display.name"))
row.names(forelegs.male) <- paste(forelegs.male$vectorbase.RU,forelegs.male$display.name)
forelegs.male <- forelegs.male[,1:2]
pdf('forelegs_dimorphic.pdf',width=6)
heatmap.2(as.matrix(log10(1+forelegs.male),trace="none",dendrogram="none",main="forelegs"),Colv="none", breaks= seq(0,3,length.out=256),
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),trace="none")
dev.off()

##### midlegs
midlegs.male <- midlegs[(midlegs$Ma_ML > TPM_cutoff | midlegs$Fe_ML_SF > TPM_cutoff) & (midlegs$log2FoldChange <= -1 * FC_cutoff | midlegs$log2FoldChange >= FC_cutoff),]
midlegs.male <- midlegs.male[order(midlegs.male$log2FoldChange),]
midlegs.male <- subset(midlegs.male,select=c("Fe_ML_SF","Ma_ML","vectorbase.RU","display.name"))
row.names(midlegs.male) <- paste(midlegs.male$vectorbase.RU,midlegs.male$display.name)
midlegs.male <- midlegs.male[,1:2]
pdf('midlegs_dimorphic.pdf',width=6)
heatmap.2(as.matrix(log10(1+midlegs.male),trace="none",dendrogram="none",main="midlegs"),Colv="none", breaks= seq(0,3,length.out=256),
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),trace="none")
dev.off()

##### abdominal_tip
abdominal_tip.male <- abdominal_tip[(abdominal_tip$Fe_At_SF > TPM_cutoff | abdominal_tip$Ma_At > TPM_cutoff) & (abdominal_tip$log2FoldChange <= -1 * FC_cutoff | abdominal_tip$log2FoldChange >= FC_cutoff),]
abdominal_tip.male <- abdominal_tip.male[order(abdominal_tip.male$log2FoldChange),]
abdominal_tip.male <- subset(abdominal_tip.male,select=c("Fe_At_SF","Ma_At","vectorbase.RU","display.name"))
row.names(abdominal_tip.male) <- paste(abdominal_tip.male$vectorbase.RU,abdominal_tip.male$display.name)
abdominal_tip.male <- abdominal_tip.male[,1:2]
pdf('abdominal_tip_dimorphic.pdf',width=6)
heatmap.2(as.matrix(log10(1+abdominal_tip.male),trace="none",dendrogram="none",main="abdominal_tip"),Colv="none", breaks= seq(0,4,length.out=256),
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),trace="none")
dev.off()


