##### flowchart to produce figures

source('~/bioinfo/github/ntx_deseq/DESeq_initialize.R')

source('~/bioinfo/github/ntx_deseq/DESeq_annotation_merge.R')

#source('~/bioinfo/github/ntx_deseq/annotation.R')
updated_annotations <- read.csv('~/Dropbox/writing/Aedes neurotranscriptome paper/CURRENT DRAFT/Supplemental data/Additional File 4 - AaegL2.1RU gene description DONE/Additional file 4 - AaegL2.1RU annotation.csv')

source('~/bioinfo/github/ntx_deseq/DESeq_annotation_plots.R')

source('~/bioinfo/github/ntx_deseq/DESeq_annotation_results.R')

source('~/bioinfo/github/ntx_deseq/TPM_plots.R')

source('~/bioinfo/github/ntx_deseq/DESeq_gonotrophic.R')

source('~/bioinfo/github/ntx_deseq/DESeq_gonotrophic_contrasts.R')

source('~/bioinfo/github/ntx_deseq/DESeq-orcoFigure.R')

source('~/bioinfo/github/ntx_deseq/DESeq-gonotrophic-timeSeries.R')
