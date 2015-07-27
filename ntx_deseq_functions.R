### PCA plot by tissue and condition


euclidean_heatmap <- function(deseq,cols)
{
  hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
  distsRL <- dist(t(assay(deseq)))
  mat <- as.matrix(distsRL)
  rownames(mat) <- colnames(mat) <- with(colData(cols), paste(tissue,condition, type, sep="-"))
  hc <- hclust(distsRL)
  heatmap.2(mat,Rowv=as.dendrogram(hc),symm=TRUE,trace="none",col=rev(hmcol),margin=c(13,13))

}


PCA_tissue_condition <- function(deseq)
{
data <- plotPCA(deseq, intgroup=c("tissue", "type","condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=tissue, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_shape_manual(values=c(21,16))
}

#### TPM/RPKM functions from 
#### https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/

### convert raw counts for a single library to TPM - takes gene counts and effective length as vectors

countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}



### count to effective count

countToEffCounts <- function(counts, len, effLen)
{
  counts * (len / effLen)
}

##################################################
## plotting functions from: https://gist.github.com/stephenturner/f60c1934405c127f09a6

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## I like mine better:

maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=display.name, cex=textcx, col=2))
  }
}


## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=display.name, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

