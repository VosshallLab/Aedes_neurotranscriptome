# annotation


updated_annotations <- read.csv('~/Dropbox/writing/Aedes neurotranscriptome paper/CURRENT DRAFT/Supplemental data/Additional File 4 - AaegL2.1RU gene description DONE/Additional file 4 - AaegL2.1RU annotation.csv')

###### everything below has been run already and consolidated into the annotation file
###### run on 2015-07-30 and saved in the Supplemental Data File 4 as the following
###### '~/Dropbox/writing/Aedes aegypti neurotranscriptome/CURRENT DRAFT/Supplemental data/Additional file 4 - AaegL2.1RU gene description/Additional file 4 - AaegL2.1RU - with OrthoDB.csv'



flybase1 <- read.csv('~/Dropbox/blast_stupid/fbgn_annotation_ID_fb_2015_03.tsv.gz.tsv',sep="\t")
flybase2 <- read.csv('~/Dropbox/blast_stupid/fbgn_fbtr_fbpp_fb_2015_03.tsv.gz.tsv',sep="\t")
orthodb <- read.csv('~/Dropbox/blast_stupid/ODB_with_simple.txt',sep="\t")

flybase <- merge(flybase1,flybase2,by=c("fbgn"),all=TRUE)
rm(flybase1); rm(flybase2)

fly_ortho <- subset(orthodb[which(orthodb$organism %in% "Drosophila melanogaster"),],select=c("vectorbase.RU","odb8_og_id"))
colnames(fly_ortho) <- c("fbpp","odb8_og_id")
aedes_orth <- subset(orthodb[which(orthodb$organism %in% "Aedes aegypti"),],select=c("vectorbase.RU","odb8_og_id"))

# pick out essential annotation fields
minimal_set <- subset(gene_annotations,select=c("internal.gene_id","display.name","gene.family","gene.family.specific","vectorbase.RU"))

# merge OrthoDB data
combined <- merge(minimal_set,aedes_orth,by="vectorbase.RU",all=TRUE)
combined <- merge(combined,fly_ortho,by="odb8_og_id",all=TRUE)
combined <- combined[which(!is.na(combined$fbpp)),]

flybase <- subset(flybase,select=c("fbpp","fbgn","annotation_ID","gene_symbol"))
flybase <- flybase[which(!is.na(flybase$fbpp)),]

final <- merge(combined,flybase,by="fbpp",all=TRUE)
final_filtered <- final[which(!is.na(final$vectorbase.RU)),]
annotation_combo <- merge(minimal_set,final_filtered,by="vectorbase.RU",all=TRUE)

write.csv(annotation_combo,'annotation_with_orthodb.csv')
###################### 



#### merge fields into human-readable, one line per gene
library(plyr)

read.csv('annotation_with_orthodb.csv') -> raw

raw.added <- ddply(raw,'internal.gene_id.x',summarize,
                   flybase.fbpp.collapsed = paste(flybase.fbpp, sep="-", collapse = "-"),
                   flybase.fbgn.collapsed = paste(flybase.fbgn, sep="-", collapse = "-"),
                   flybase.annot.collapsed = paste(flybase.annot, sep="-", collapse = "-"),
                   flybase.gene.collapsed = paste(flybase.gene, sep="-", collapse = "-"),
                   orthodb.collapsed = paste(orthodb.group, sep = "-", collapse = "-")
)

raw.original <- unique(subset(raw,select=c("vectorbase.RU","internal.gene_id.x",
                           "display.name","gene.family","gene.family.specific")))

raw.added <- merge(raw.added,
                   raw.original,
                                    by="internal.gene_id.x")
