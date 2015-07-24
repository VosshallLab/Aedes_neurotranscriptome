library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(gplots)

source('ntx_deseq_functions.R')

setwd('~/bioinfo/github/ntx_deseq/')

PE <- read.csv('raw_gene_counts/PE_featurecounts.txt',sep = "\t",comment.char = "#")
SE <- read.csv('raw_gene_counts/SE_featurecounts.txt',sep = "\t",comment.char = "#")

merge(PE,SE,by=c("Geneid","Chr","Start","End","Strand","Length")) -> test

subset(test,select=c("Geneid","Chr","Start","End","Strand","Length")) -> gene_annotations

columns <- read.csv('annotations/library_key.csv')
colnames(test) <- columns$y
row.names(test) <- test$gene
raw_FC <- test[,7:140]

##### read in library information from CSV file

libprop <- read.csv('annotations/library_info.csv')
row.names(libprop) <- libprop$library
raw_FC_reordered <- raw_FC[row.names(libprop)]

