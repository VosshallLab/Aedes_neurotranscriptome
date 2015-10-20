### load a set of core functions
source('ntx_deseq_functions.R')

### load the feature counts data and generate count matrix (rawcounts)

PE <- read.csv('raw_gene_counts/PE_featurecounts.txt',sep = "\t",comment.char = "#")
SE <- read.csv('raw_gene_counts/SE_featurecounts.txt',sep = "\t",comment.char = "#")
rawcounts <- merge(PE,SE,by=c("Geneid","Chr","Start","End","Strand","Length"))

### generate or load gene annotations
# generate
#source('DESeq_annotation_merge.R')

# or load
gene_annotations <- read.csv('annotations/AaegL2.1RUv2_annotations.csv')
row.names(gene_annotations) <- gene_annotations$internal.gene_id

# read in library information
columns <- read.csv('annotations/library_key.csv')
colnames(rawcounts) <- columns$y
row.names(rawcounts) <- rawcounts$gene

# generate merged file with annotation information and raw counts
rawcounts_only <- rawcounts[,7:length(colnames(rawcounts))]
rawcounts_with_annotation <- merge(rawcounts_only,gene_annotations,by="row.names")
row.names(rawcounts_with_annotation) <- rawcounts_with_annotation$internal.gene_id

# convert counts to TPM
tpm_all <- apply(rawcounts_only,2,countToTpm,rawcounts$len)
tpm_all_with_annotation <- merge(tpm_all,gene_annotations,by="row.names")


##### read in library information from CSV file and add info

libprop <- read.csv('annotations/library_info.csv')
row.names(libprop) <- libprop$library

# pull out counts and TPM for each library
rawcounts_with_annotation_reordered <- rawcounts_with_annotation[row.names(libprop)]
tpm_all_with_annotation_reordered <- tpm_all_with_annotation[row.names(libprop)]
