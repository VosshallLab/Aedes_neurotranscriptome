##### flowchart to produce figures for Matthews et al., BMC Genomics, 2015

# load required libraries
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(reshape2)
library(plyr)
library(VennDiagram)

# set working directory to github repository
setwd('~/bioinfo/github/ntx_deseq/')

# initialize environment by loading/generating counts and annotations
source('~/bioinfo/github/ntx_deseq/DESeq_initialize.R')

# manipulate matrices to generate annotation plots
source('~/bioinfo/github/ntx_deseq/DESeq_annotation_plots.R')

## TPM heatmaps (figures 4-7)
source('~/bioinfo/github/ntx_deseq/TPM_plots.R')

## figure 8 - male vs. sugar-fed female
source('~/bioinfo/github/ntx_deseq/DESeq_annotation_results.R')
source('~/bioinfo/github/ntx_deseq/DESeq-figure8-dimorphism.R')

## figure 9 - gonotrophic
source('~/bioinfo/github/ntx_deseq/DESeq_gonotrophic.R')
source('~/bioinfo/github/ntx_deseq/DESeq_gonotrophic_contrasts.R')
source('~/bioinfo/github/ntx_deseq/DESeq-gonotrophic-timeSeries.R')

## orco analysis (figure 6)
source('~/bioinfo/github/ntx_deseq/DESeq-orcoFigure.R')