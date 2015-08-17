grab_genefamily <- function(genefam,tpmMat) {
  
  family_subset <- tpmMat[tpmMat$gene.family %in% genefam,]
  #  row.names(family_subset) <- family_subset$display.name
}

OR <- grab_genefamily("OR",compiled_means)

OR_plot_An <- ggplot(data=OR,aes(x=log10(Fe_An_SF+1),y=log10(Ma_An+1))) + geom_point() + theme_bw()
OR_plot_Rs <- ggplot(data=OR,aes(x=log10(Fe_Rs_SF+1),y=log10(Ma_Rs+1))) + geom_point() + theme_bw()

OR$status <- rep(0,length(OR$vectorbase.RU))

OR[OR$vectorbase.RU %in% c("AAEL017043","AAEL011409","AAEL018095","AAEL016966","AAEL017347","AAEL017123","AAEL017505","AAEL017000","AAEL018070","AAEL017201","AAEL014197","AAEL017014"),]$status <- 1


OR_plot_An_subset <- ggplot(data=OR,aes(x=log10(Fe_An_SF+1),y=log10(Ma_An+1),colour=status)) + geom_point() + theme_bw()


OR$FM <- OR$Fe_An_SF / OR$Ma_An

FMFC <- OR[order(-OR$FM),]