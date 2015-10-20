# code to merge FlyBase OrthoDB results
##source('~/bioinfo/github/ntx_deseq/annotation.R')

gene_annotations <- rawcounts

# set the row names as the gene ID and create a column called internal.gene_id
row.names(gene_annotations) <- gene_annotations$Geneid
gene_annotations$internal.gene_id <- gene_annotations$Geneid

# delete extraneous columns to generate an annotation template
gene_annotations <- subset(gene_annotations, select=c("internal.gene_id","Chr","Start","End","Strand","Length"))

# read in manual annotations - currently a copy of latest draft from Dropbox (2015.10.20)
manual_annotations <- read.csv('annotations/Additional file 3 - AaegLRU annotation.csv')
row.names(manual_annotations) <- manual_annotations$internal.gene_id

# merge manual and count annotations
gene_annotations <- merge(gene_annotations, manual_annotations,by="internal.gene_id")
row.names(gene_annotations) <- gene_annotations$internal.gene_id

# write out finalized annotation file as .csv
write.csv(gene_annotations,'annotations/AaegL2.1RUv2_annotations.csv')