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
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
}

### convert raw counts for a single library to TPM - takes gene counts and effective length as vectors

countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}


### convert raw counts for a single library to TPM - takes gene counts and effective length as vectors

countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}


### convert FPKM to TPM

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}


### count to effective count

countToEffCounts <- function(counts, len, effLen)
{
  counts * (len / effLen)
}

