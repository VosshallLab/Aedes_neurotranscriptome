updated_annotations <- read.csv('~/Dropbox/writing/Aedes neurotranscriptome paper/CURRENT DRAFT/Supplemental data/Additional File 3 - AaegLRU gene description DONE/Additional file 3 - AaegLRU annotation.csv')


setwd('/Users/ben/Dropbox/writing/Aedes neurotranscriptome paper/Data/Figure 3/Genome annotation/cuffcompare/cuffcompare_08192015/')

loci <- read.csv('3_ref_RU.loci',header=F,sep="\t")
loci[loci$V3 %in% "-",] -> new_loci

loci_2 <- read.csv('2_ref_RU.loci',header=F,sep="\t")
loci_2[loci_2$V3 %in% "-",] -> new_loci_2

write.csv(new_loci,'3_ref_new_loci.csv',quote = F)
write.csv(new_loci_2, '2_ref_new_loci.csv',quote = F)

genes_only <- read.csv('new_loci_genes_3',header=T)
genes_only_2 <- read.csv('new_loci_genes_2',header=T)

novel_loci_genes <- updated_annotations[updated_annotations$internal.gene_id %in% as.character(genes_only$V4),]
novel_loci_genes$vectorbase.RU

novel_loci_genes_2 <- updated_annotations[updated_annotations$internal.gene_id %in% as.character(genes_only_2$V4),]
novel_loci_genes_2$vectorbase.RU
