### PCA plot by tissue and condition

dimorphism_heatmap <- function(gene_list,tpmMat,maxTPM) {
  
  family_subset <- tpmMat[tpmMat$internal.gene_id %in% gene_list,]
  #  row.names(family_subset) <- family_subset$display.name
  disp_name <- paste(family_subset$vectorbase.RU,family_subset$display.name,sep="-")
  family_subset <- subset(family_subset,select=tissue_list)
  column_labels <- c("Ovaries","Palp","Proboscis","Brain","Antenna","Rostrum","Forelegs","Midlegs","Hindlegs","Abdominal tip","Brain","Antenna","Rostrum","Forelegs","Midlegs","Hindlegs","Abdominal tip")
  heatmap.2(log10(as.matrix(family_subset)+1),dendrogram="none",Colv=FALSE,sepwidth=c(.5,.5),
            trace="none",labCol=column_labels,srtCol=45,density.info="none",labRow=disp_name,
            key.xlab="Log10(TPM+1)",key.title=NA,
            #scale = "row",
            breaks= seq(0,maxTPM,length.out=256),
            col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
            # col = colorRampPalette(rev(brewer.pal(9,"Blues")))(255) 
            # col = rev(viridis(100))
  )
}

dimorphism_heatmap_fig8_F <- function(gene_list,tpmMat,maxTPM) {
  
  family_subset <- tpmMat[tpmMat$internal.gene_id %in% gene_list,]
  family_subset$disp_name <- paste(family_subset$vectorbase.RU,family_subset$display.name,sep="-")
  row.names(family_subset) <- family_subset$disp_name
  family_subset <- subset(family_subset,select=c("Fe_Br_SF","Fe_An_SF","Fe_FL_SF","Fe_ML_SF","FE_HL_SF","Fe_At_SF","Ma_Br","Ma_An","Ma_FL","Ma_ML","Ma_HL","Ma_At"))
  
  fem_sum <- apply(family_subset[,1:6],1,sum)
  male_sum <- apply(family_subset[,7:12],1,sum)
  family_subset <- family_subset[order(-fem_sum),]
  
  #  row.names(family_subset) <- family_subset$display.name
 
  column_labels <- c("Brain","Antenna","Forelegs","Midlegs","Hindlegs","Abdominal tip","Brain","Antenna","Forelegs","Midlegs","Hindlegs","Abdominal tip")
  heatmap.2(log10(as.matrix(family_subset)+1),dendrogram="none",Colv=FALSE,Rowv = FALSE,sepwidth=c(.5,.5),
            trace="none",labCol=column_labels,srtCol=45,density.info="none",
            key.xlab="Log10(TPM+1)",key.title=NA,
            #scale = "row",
            breaks= seq(0,maxTPM,length.out=256),
            col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
            
            # col = colorRampPalette(rev(brewer.pal(9,"Blues")))(255) 
            # col = rev(viridis(100))
  )
}

dimorphism_heatmap_fig8_M <- function(gene_list,tpmMat,maxTPM) {
  
  family_subset <- tpmMat[tpmMat$internal.gene_id %in% gene_list,]
  family_subset$disp_name <- paste(family_subset$vectorbase.RU,family_subset$display.name,sep="-")
  row.names(family_subset) <- family_subset$disp_name
  family_subset <- subset(family_subset,select=c("Fe_Br_SF","Fe_An_SF","Fe_FL_SF","Fe_ML_SF","FE_HL_SF","Fe_At_SF","Ma_Br","Ma_An","Ma_FL","Ma_ML","Ma_HL","Ma_At"))
  
  fem_sum <- apply(family_subset[,1:6],1,sum)
  male_sum <- apply(family_subset[,7:12],1,sum)
  family_subset <- family_subset[order(male_sum),]

  
  column_labels <- c("Brain","Antenna","Forelegs","Midlegs","Hindlegs","Abdominal tip","Brain","Antenna","Forelegs","Midlegs","Hindlegs","Abdominal tip")
  heatmap.2(log10(as.matrix(family_subset)+1),dendrogram="none",Colv=FALSE,Rowv = FALSE,sepwidth=c(.5,.5),
            trace="none",labCol=column_labels,srtCol=45,density.info="none",
            key.xlab="Log10(TPM+1)",key.title=NA,
            #scale = "row",
            breaks= seq(0,maxTPM,length.out=256),
            col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
            
            # col = colorRampPalette(rev(brewer.pal(9,"Blues")))(255) 
            # col = rev(viridis(100))
  )
}
### function to find a specific tissue from the whole TPM matrix

subset_TPM <- function(tissue_to_find,tpm){
  
  selection <- grepl(tissue_to_find,colnames(tpm))
  tpm_sub <- tpm[,selection]
  temp_mean <- apply(tpm_sub,1,mean)
  temp_median <- apply(tpm_sub,1,median)
  tpm_sub$median <- temp_median
  tpm_sub$mean <- temp_mean
  tpm_sub
}


euclidean_heatmap <- function(deseq,cols)
{
  hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
  distsRL <- dist(t(assay(deseq)))
  mat <- as.matrix(distsRL)
  rownames(mat) <- colnames(mat) <- with(colData(cols), paste(tissue,condition, type, sep="-"))
  hc <- hclust(distsRL)
  heatmap.2(mat,Rowv=as.dendrogram(hc),symm=TRUE,trace="none",col=rev(hmcol),margin=c(13,13))

}

TPM_one_gene <- function(dds,gene,group){
  
  d <- plotCounts(dds, gene=gene, intgroup=group, returnData=TRUE)
  ggplot(d, aes(x=condition, y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + scale_y_log10(breaks=c(0,100,1000))

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

### cluster heatmap plots

hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")


# inspired by http://stackoverflow.com/questions/22278508/how-to-add-colsidecolors-on-heatmap-2-after-performing-bi-clustering-row-and-co

heatmap_with_cluster <- function(mydata,tissue,clustNum) 
{
  
  
  
  rowLabels <- updated_annotations[updated_annotations$internal.gene_id %in% row.names(mydata),]
  
  cl.row <- hclustfunc(distfunc(mydata))
  
  # extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
  gr.row <- cutree(cl.row, clustNum)
  
  # require(RColorBrewer)
  col1 <- brewer.pal(clustNum, "Set1")
  
  # require(gplots)
  
  pdf(file = paste(tissue,'_cluster_plot_zscore.pdf'),width=12,height=24)
  
  heatmap.2(as.matrix(mydata), hclustfun=hclustfunc, distfun=distfunc,   
            RowSideColors=col1[gr.row],trace="none",dendrogram="row",
            Rowv = reorder(as.dendrogram(cl.row),1:clustNum),
            col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),main=tissue,
            Colv = "none")#,labCol = c("non-blood-fed","blood-fed","gravid"))
  
  dev.off()
  
  cluster_data <- as.data.frame(gr.row)
  cluster_data$internal.gene_id = row.names(cluster_data)
  merge(cluster_data,updated_annotations,by="internal.gene_id") -> cluster_data
  
  write.csv(cluster_data,paste(tissue,'_cluster_data.csv',sep=''))
  
  return(cluster_data)
  
}