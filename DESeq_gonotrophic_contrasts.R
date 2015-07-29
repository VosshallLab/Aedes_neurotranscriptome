#### contrasts, SF vs. BF

## HL

BF_hindlegs <- results(dds_hindlegs_gono,alpha=0.05,contrast=c("condition","SF","BF"))
BF_hindlegs[!is.na(BF_hindlegs$padj) & BF_hindlegs$padj < 0.05,] -> sigBF_hindlegs


## brain

BF_brain <- results(dds_brain_gono,alpha=0.05,contrast=c("condition","SF","BF"))
BF_brain[!is.na(BF_brain$padj) & BF_brain$padj < 0.05,] -> sigBF_brain


## rostrum

BF_rostrum <- results(dds_rostrum_gono,alpha=0.05,contrast=c("condition","SF","BF"))
BF_rostrum[!is.na(BF_rostrum$padj) & BF_rostrum$padj < 0.05,] -> sigBF_rostrum

## antenna

BF_antenna <- results(dds_antenna_gono,alpha=0.05,contrast=c("condition","SF","BF"))
BF_antenna[!is.na(BF_antenna$padj) & BF_antenna$padj < 0.05,] -> sigBF_antenna



#### contrasts, SF vs. O

## brain

O_brain <- results(dds_brain_gono,alpha=0.05,contrast=c("condition","SF","O"))
O_brain[!is.na(O_brain$padj) & O_brain$padj < 0.05,] -> sigO_brain


## antenna

O_antenna <- results(dds_antenna_gono,alpha=0.05,contrast=c("condition","SF","O"))
O_antenna[!is.na(O_antenna$padj) & O_antenna$padj < 0.05,] -> sigO_antenna

## foreleg

O_forelegs <- results(dds_forelegs_gono,alpha=0.05,contrast=c("condition","SF","O"))
O_forelegs[!is.na(O_forelegs$padj) & O_forelegs$padj < 0.05,] -> sigO_foreleg


## midleg

O_midlegs <- results(dds_midlegs_gono,alpha=0.05,contrast=c("condition","SF","O"))
O_midlegs[!is.na(O_midlegs$padj) & O_midlegs$padj < 0.05,] -> sigO_midleg

## hindleg

O_hindlegs <- results(dds_hindlegs_gono,alpha=0.05,contrast=c("condition","SF","O"))
O_hindlegs[!is.na(O_hindlegs$padj) & O_hindlegs$padj < 0.05,] -> sigO_hindlegs

## abdominal tip

O_at <- results(dds_abdominaltip_gono,alpha=0.05,contrast=c("condition","SF","O"))
O_at[!is.na(O_at$padj) & O_at$padj < 0.05,] -> sigO_at

## ovaries

O_ovaries <- results(dds_ovaries_gono,alpha=0.05,contrast=c("condition","SF","O"))
O_ovaries[!is.na(O_ovaries$padj) & O_ovaries$padj < 0.05,] -> sigO_ovaries

#################
# print MA plots and volcano plots for all


pdf("plots/gono/MA_BF_brain.pdf",width=6)
plotMA(BF_brain,ylim=c(-5,5))
dev.off()

pdf("plots/gono/MA_BF_antenna.pdf",width=6)
plotMA(BF_antenna,ylim=c(-5,5))
dev.off()

pdf("plots/gono/MA_BF_hindlegs.pdf",width=6)
plotMA(BF_hindlegs,ylim=c(-5,5))
dev.off()

pdf("plots/gono/MA_BF_rostrum.pdf",width=6)
plotMA(BF_rostrum,ylim=c(-5,5))
dev.off()

#### O

pdf("plots/gono/MA_O_brain.pdf",width=6)
plotMA(O_brain,ylim=c(-5,5))
dev.off()

pdf("plots/gono/MA_O_antenna.pdf",width=6)
plotMA(O_antenna,ylim=c(-5,5))
dev.off()

pdf("plots/gono/MA_O_hindlegs.pdf",width=6)
plotMA(O_hindlegs,ylim=c(-5,5))
dev.off()

pdf("plots/gono/MA_O_midlegs.pdf",width=6)
plotMA(O_midlegs,ylim=c(-5,5))
dev.off()

pdf("plots/gono/MA_O_forelegs.pdf",width=6)
plotMA(O_forelegs,ylim=c(-5,5))
dev.off()

#pdf("plots/gono/MA_O_ovaries.pdf",width=6)
plotMA(O_ovaries,ylim=c(-5,5))
#dev.off()

#pdf("plots/gono/MA_O_abdominal_tip.pdf",width=6)
plotMA(O_at,ylim=c(-5,5))
#dev.off()

