# annotation

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

