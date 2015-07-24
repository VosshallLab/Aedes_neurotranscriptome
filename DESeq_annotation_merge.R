setwd('~/bioinfo/github/ntx_deseq/')

PE <- read.csv('raw_gene_counts/PE_featurecounts.txt',sep = "\t",comment.char = "#")
SE <- read.csv('raw_gene_counts/SE_featurecounts.txt',sep = "\t",comment.char = "#")

gene_annotations <- merge(PE,SE,by=c("Geneid","Chr","Start","End","Strand","Length"))
gene_annotations <- subset(gene_annotations, select=c("Geneid","Chr","Start","End","Strand","Length"))

manual_annotations <- read.csv('~/Dropbox/writing/Aedes neurotranscriptome paper/CURRENT DRAFT/Supplemental data/Additional File 4 - AaegL2.1RU gene description DONE/Additional file 4 - AaegL2.1RU.csv')
row.names(manual_annotations) <- manual_annotations$internal.gene_id

gene_annotations <- merge(gene_annotations, manual_annotations)

write.csv(gene_annotations,'annotations/AaegL2.1RUv2_annotations.csv')