grab_genefamily <- function(genefam,tpmMat) {
  
  family_subset <- tpmMat[tpmMat$gene.family %in% genefam,]
  #  row.names(family_subset) <- family_subset$display.name
}

OR <- grab_genefamily("OR",compiled_means)

OR_plot_An <- ggplot(data=OR,aes(x=log10(Fe_An_SF+1),y=log10(Ma_An+1))) + geom_point() + theme_bw()
  OR_plot_Rs <- ggplot(data=OR,aes(x=log10(Fe_Rs_SF+1),y=log10(Ma_Rs+1))) + geom_point() + theme_bw()
