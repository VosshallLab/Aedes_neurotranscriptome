or_mean = compiled_means[compiled_means$gene.family %in% c("OR"),]
ir_mean = compiled_means[compiled_means$gene.family %in% c("IR"),]

combo_mean = c(or_mean,ir_mean)

combo_mean_df = data.frame(gene=combo_mean$vectorbase.RU,family = combo_mean$gene.family,
                           fe_an = combo_mean$Fe_An_SF,ma_an = combo_mean$Ma_An,
                           fe_rs = combo_mean$Fe_Rs_SF,ma_rs = combo_mean$Ma_Rs)
