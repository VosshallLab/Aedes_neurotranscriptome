setwd('~/bioinfo/github/ntx_deseq/')

PE <- read.csv('raw_gene_counts/PE_featurecounts.txt',sep = "\t",comment.char = "#")
SE <- read.csv('raw_gene_counts/SE_featurecounts.txt',sep = "\t",comment.char = "#")

gene_annotations <- merge(PE,SE,by=c("Geneid","Chr","Start","End","Strand","Length"))
row.names(gene_annotations) <- gene_annotations$Geneid
gene_annotations$internal.gene_id <- gene_annotations$Geneid
gene_annotations <- subset(gene_annotations, select=c("internal.gene_id","Chr","Start","End","Strand","Length"))

manual_annotations <- read.csv('~/Dropbox/writing/Aedes neurotranscriptome paper/CURRENT DRAFT/Supplemental data/Additional File 3 - AaegL2.1RU gene description DONE/Additional file 3 - AaegL2.1RU annotation.csv')
row.names(manual_annotations) <- manual_annotations$internal.gene_id

gene_annotations <- merge(gene_annotations, manual_annotations,by="internal.gene_id")
row.names(gene_annotations) <- gene_annotations$internal.gene_id

write.csv(gene_annotations,'annotations/AaegL2.1RUv2_annotations.csv')