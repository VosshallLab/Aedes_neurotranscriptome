high_brain_male <- subset(compiled_means,select=c("Ma_Br","Ma_Rs","Ma_HL"))
high_brain_female <- subset(compiled_means,select=c("Fe_Br_SF","Fe_Rs_SF","FE_HL_SF"))
high_brain_male_rev <- subset(compiled_means,select=c("Ma_Br","Ma_Rs","Ma_HL"))
high_brain_female_rev <- subset(compiled_means,select=c("Fe_Br_SF","Fe_Rs_SF","FE_HL_SF"))


# pairs.panels(high_brain_male)
# pairs.panels(high_brain_male_rev)
# pairs.panels(high_brain_female)
# pairs.panels(high_brain_female_rev)

brain_rs <- grab_genefamily("neuropeptide",compiled_means)

brain_rs <- subset(brain_rs,select=c("Ma_Br","Fe_Br_SF","Ma_Rs","Fe_Rs_SF"))

heatmap.2(log10(1+as.matrix(brain_rs)),scale="row",dendrogram="none")
