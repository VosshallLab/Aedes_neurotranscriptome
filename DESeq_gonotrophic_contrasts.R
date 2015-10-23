#### contrasts, SF vs. BF

ALPHA = 0.01

## HL

BF_hindlegs <- results(dds_hindlegs_gono,alpha=ALPHA,contrast=c("condition","SF","BF"))
BF_hindlegs[!is.na(BF_hindlegs$padj) & BF_hindlegs$padj < ALPHA,] -> sigBF_hindlegs


## brain

BF_brain <- results(dds_brain_gono,alpha=ALPHA,contrast=c("condition","SF","BF"))
BF_brain[!is.na(BF_brain$padj) & BF_brain$padj < ALPHA,] -> sigBF_brain


## rostrum

BF_rostrum <- results(dds_rostrum_gono,alpha=ALPHA,contrast=c("condition","SF","BF"))
BF_rostrum[!is.na(BF_rostrum$padj) & BF_rostrum$padj < ALPHA,] -> sigBF_rostrum

## antenna

BF_antenna <- results(dds_antenna_gono,alpha=ALPHA,contrast=c("condition","SF","BF"))
BF_antenna[!is.na(BF_antenna$padj) & BF_antenna$padj < ALPHA,] -> sigBF_antenna



#### contrasts, SF vs. O

## brain

O_brain <- results(dds_brain_gono,alpha=ALPHA,contrast=c("condition","SF","O"))
O_brain[!is.na(O_brain$padj) & O_brain$padj < ALPHA,] -> sigO_brain


## antenna

O_antenna <- results(dds_antenna_gono,alpha=ALPHA,contrast=c("condition","SF","O"))
O_antenna[!is.na(O_antenna$padj) & O_antenna$padj < ALPHA,] -> sigO_antenna

## foreleg

O_forelegs <- results(dds_forelegs_gono,alpha=ALPHA,contrast=c("condition","SF","O"))
O_forelegs[!is.na(O_forelegs$padj) & O_forelegs$padj < ALPHA,] -> sigO_foreleg


## midleg

O_midlegs <- results(dds_midlegs_gono,alpha=ALPHA,contrast=c("condition","SF","O"))
O_midlegs[!is.na(O_midlegs$padj) & O_midlegs$padj < ALPHA,] -> sigO_midleg

## hindleg

O_hindlegs <- results(dds_hindlegs_gono,alpha=ALPHA,contrast=c("condition","SF","O"))
O_hindlegs[!is.na(O_hindlegs$padj) & O_hindlegs$padj < ALPHA,] -> sigO_hindlegs

## abdominal tip

O_at <- results(dds_abdominaltip_gono,alpha=ALPHA,contrast=c("condition","SF","O"))
O_at[!is.na(O_at$padj) & O_at$padj < ALPHA,] -> sigO_at

## ovaries

O_ovaries <- results(dds_ovaries_gono,alpha=ALPHA,contrast=c("condition","SF","O"))
O_ovaries[!is.na(O_ovaries$padj) & O_ovaries$padj < ALPHA,] -> sigO_ovaries

#################
# print MA plots and volcano plots for all


#pdf("plots/gono/MA_BF_brain.pdf",width=8,height=6)
plotMA(BF_brain,ylim=c(-5,5),alpha=ALPHA)
#dev.off()

#pdf("plots/gono/MA_BF_antenna.pdf",width=8,height=6)
plotMA(BF_antenna,ylim=c(-5,5),alpha=ALPHA)
#dev.off()

#pdf("plots/gono/MA_BF_hindlegs.pdf",width=8,height=6)
plotMA(BF_hindlegs,ylim=c(-5,5),alpha=ALPHA)
#dev.off()

#pdf("plots/gono/MA_BF_rostrum.pdf",width=8,height=6)
plotMA(BF_rostrum,ylim=c(-5,5),alpha=ALPHA)
#dev.off()

##### O

#pdf("plots/gono/MA_O_brain.pdf",width=8,height=6)
plotMA(O_brain,ylim=c(-5,5),alpha=ALPHA)
#dev.off()

#pdf("plots/gono/MA_O_antenna.pdf",width=8,height=6)
plotMA(O_antenna,ylim=c(-5,5),alpha=ALPHA)
#dev.off()

#pdf("plots/gono/MA_O_hindlegs.pdf",width=8,height=6)
plotMA(O_hindlegs,ylim=c(-5,5),alpha=ALPHA)
#dev.off()

#pdf("plots/gono/MA_O_midlegs.pdf",width=8,height=6)
plotMA(O_midlegs,ylim=c(-5,5),alpha=ALPHA)
#dev.off()

#pdf("plots/gono/MA_O_forelegs.pdf",width=8,height=6)
plotMA(O_forelegs,ylim=c(-5,5),alpha=ALPHA)
#dev.off()

#pdf("plots/gono/MA_O_ovaries.pdf",width=8,height=6)
plotMA(O_ovaries,ylim=c(-5,5),alpha=ALPHA)
#dev.off()

#pdf("plots/gono/MA_O_abdominal_tip.pdf",width=8,height=6)
plotMA(O_at,ylim=c(-5,5),alpha=ALPHA)
#dev.off()



#### write csv files for supplemental files

annotation_subset <- subset(gene_annotations,select=c("vectorbase.RU","display.name"))
padj_cutoff <- 0.1

write.csv(merge(annotation_subset,as.data.frame(O_ovaries[!is.na(O_ovaries$padj) & O_ovaries$padj < padj_cutoff,]),by="row.names"),'O_ovaries.csv')
write.csv(merge(annotation_subset,as.data.frame(O_at[!is.na(O_at$padj) & O_at$padj < padj_cutoff,]),by="row.names"),'O_at.csv')
write.csv(merge(annotation_subset,as.data.frame(O_forelegs[!is.na(O_forelegs$padj) & O_forelegs$padj < padj_cutoff,]),by="row.names"),'O_forelegs.csv')
write.csv(merge(annotation_subset,as.data.frame(O_midlegs[!is.na(O_midlegs$padj) & O_midlegs$padj < padj_cutoff,]),by="row.names"),'O_midlegs.csv')
write.csv(merge(annotation_subset,as.data.frame(O_hindlegs[!is.na(O_hindlegs$padj) & O_hindlegs$padj < padj_cutoff,]),by="row.names"),'O_hindlegs.csv')
write.csv(merge(annotation_subset,as.data.frame(O_antenna[!is.na(O_antenna$padj) & O_antenna$padj < padj_cutoff,]),by="row.names"),'O_antenna.csv')
write.csv(merge(annotation_subset,as.data.frame(O_brain[!is.na(O_brain$padj) & O_brain$padj < padj_cutoff,]),by="row.names"),'O_brain.csv')

write.csv(merge(annotation_subset,as.data.frame(BF_brain[!is.na(BF_brain$padj) & BF_brain$padj < padj_cutoff,]),by="row.names"),'BF_brain.csv')
write.csv(merge(annotation_subset,as.data.frame(BF_antenna[!is.na(BF_antenna$padj) & BF_antenna$padj < padj_cutoff,]),by="row.names"),'BF_antenna.csv')
write.csv(merge(annotation_subset,as.data.frame(BF_rostrum[!is.na(BF_rostrum$padj) & BF_rostrum$padj < padj_cutoff,]),by="row.names"),'BF_rostrum.csv')
write.csv(merge(annotation_subset,as.data.frame(BF_hindlegs[!is.na(BF_hindlegs$padj) & BF_hindlegs$padj < padj_cutoff,]),by="row.names"),'BF_hindlegs.csv')



