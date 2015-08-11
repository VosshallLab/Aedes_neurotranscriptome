#### TPM_vectorbase calculations

vec_tpm_PE <- read.csv('~/Dropbox/writing/Aedes neurotranscriptome paper/Data/Additional_file_vectorbase_TPM/PE-vectorbase-renamed.txt',sep = "\t",comment.char = "#",header=TRUE)
vec_tpm_SE <- read.csv('~/Dropbox/writing/Aedes neurotranscriptome paper/Data/Additional_file_vectorbase_TPM/SE-vectorbase-renamed.txt',sep = "\t",comment.char = "#",header=TRUE)
merge(vec_tpm_PE,vec_tpm_SE,by=c("Geneid","Chr","Start","End","Strand","Length")) -> vec_rawcounts

columns <- read.csv('annotations/library_key.csv')
colnames(vec_rawcounts) <- columns$y
row.names(vec_rawcounts) <- vec_rawcounts$gene

write.csv(vec_rawcounts,'vectorbase_counts_all.csv')

vec_tpm <- apply(vec_rawcounts[,7:length(colnames(vec_rawcounts))],2,countToTpm,vec_rawcounts$len)
row.names(vec_tpm) <- vec_rawcounts$gene

libprop_all <- libprop[ libprop$keep == 1 ,]

vec_tpm_keep <- subset(vec_tpm,select=row.names(libprop_all))
write.csv(vec_tpm_keep,'~/Dropbox/writing/Aedes neurotranscriptome paper/Data/Additional_file_vectorbase_TPM/vectorbase_tpm.csv')

