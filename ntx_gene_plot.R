library(ggplot2)
library(reshape2)

ntx_gene_plot = function(gene_list){
  
  meanSubset = compiled_means[compiled_means$internal.gene_id %in% gene_list,9:25]
  
  colnames(meanSubset) <- c("Female Ovaries","Female Palps","Female Proboscis","Female Brain","Female Antenna","Female Rostrum","Female Forelegs","Female Midlegs","Female Hindlegs","Female Abdominal Tip","Male Brain","Male Antenna","Male Rostrum","Male Forelegs","Male Midlegs","Male Hindlegs","Male Abdominal Tip")
  
  meanPlot <- ggplot(data = melt(as.matrix(meanSubset)), aes(x=Var2,y=value,group=Var1,colour=Var1)) + geom_point(size=5) + facet_grid(Var1~., scales="free_y") +theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  print(meanPlot)
  return(meanPlot)

}

gene_lookup_AAEL = function(gene_list){
  
  gene_num <- gene_annotations[gene_annotations$vectorbase.RU %in% gene_list,]$internal.gene_id
  
  return(gene_num)
}


ammonium_plot <- ntx_gene_plot(c("gene13896","gene9231","gene9235"))

rh50_plot <- ntx_gene_plot(c("gene13896"))