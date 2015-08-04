library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(reshape2)

setwd('~/bioinfo/github/ntx_deseq/')


### load a set of core functions
source('ntx_deseq_functions.R')

### load gene annotations
source('DESeq_annotation_merge.R')

### load the feature counts data and generate count matrices

PE <- read.csv('raw_gene_counts/PE_featurecounts.txt',sep = "\t",comment.char = "#")
SE <- read.csv('raw_gene_counts/SE_featurecounts.txt',sep = "\t",comment.char = "#")

merge(PE,SE,by=c("Geneid","Chr","Start","End","Strand","Length")) -> rawcounts

columns <- read.csv('annotations/library_key.csv')
colnames(rawcounts) <- columns$y
row.names(rawcounts) <- rawcounts$gene

rawcounts_only <- rawcounts[,7:length(colnames(rawcounts))]
rawcounts_with_annotation <- merge(rawcounts_only,gene_annotations,by="row.names")
row.names(rawcounts_with_annotation) <- rawcounts_with_annotation$internal.gene_id

tpm_all <- apply(rawcounts_only,2,countToTpm,rawcounts$len)
tpm_all_with_annotation <- merge(tpm_all,gene_annotations,by="row.names")

##### read in library information from CSV file

libprop <- read.csv('annotations/library_info.csv')
row.names(libprop) <- libprop$library

rawcounts_with_annotation_reordered <- rawcounts_with_annotation[row.names(libprop)]
tpm_all_with_annotation_reordered <- tpm_all_with_annotation[row.names(libprop)]

### load gene annotations
