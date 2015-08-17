blastx <- read.csv('~/Dropbox/blast_stupid/L2.1CDS_blastX',sep="\t",header=FALSE)
colnames(blastx) <- c("transcript","fbpp","pident","length","mismatch","gapopen","qstart","qend","sstart","ssend","e_value","bitscore")

flybase1 <- read.csv('~/Dropbox/blast_stupid/fbgn_annotation_ID_fb_2015_03.tsv.gz.tsv',sep="\t")
flybase2 <- read.csv('~/Dropbox/blast_stupid/fbgn_fbtr_fbpp_fb_2015_03.tsv.gz.tsv',sep="\t")

flybase <- merge(flybase1,flybase2,by=c("fbgn"),all=TRUE)
rm(flybase1); rm(flybase2)

blastx_flybase <- merge(blastx,flybase, by="fbpp")


library(plyr)

test.gene <- strsplit(as.character(blastx_flybase$transcript), ".",fixed=TRUE)
df <- ldply(test.gene)
colnames(df) <- c("gene", "isoform")
blastx_flybase <- cbind(df,blastx_flybase)

blastx_flybase <- blastx_flybase[order(blastx_flybase$"e_value"),]

blastx_reduced <- blastx_flybase[!duplicated(blastx_flybase[c(1,18)]),]

blastx.combined <- ddply(blastx_reduced,'gene',summarize,
                   flybase.fbpp.collapsed = paste((fbpp), sep="-", collapse = "-"),
                   flybase.fbgn.collapsed = paste((fbgn), sep="-", collapse = "-"),
                   flybase.annot.collapsed = paste((annotation_ID), sep="-", collapse = "-"),
                   flybase.gene.collapsed = paste((gene_symbol), sep="-", collapse = "-"),
                   blastx.evalue = paste((e_value), sep="-", collapse = "-"),
                   .progress = "text")

colnames(blastx.combined) <- c("internal.gene_id","dmel.blastx.fbpp","dmel.blastx.fbgn","dmel.blastx.annot","dmel.blastx.gene","dmel.blastx.evalues")

updated_annotations_w_blastx <- merge(updated_annotations,blastx.combined,by="internal.gene_id",all.x=TRUE)


## below, write .CSV file
# write.csv(updated_annotations_w_blastx,'~/Dropbox/writing/Aedes neurotranscriptome paper/CURRENT DRAFT/Supplemental data/Additional File 3 - AaegLRU gene description DONE/Additional file 3 - AaegLRU annotation_withBLASTX.csv')
# blastx.combined.selected <- subset(blastx.combined,select=c("gene","flybase.annot.collapsed",flybase.fbpp.collapsed))