### generate results

## antenna
res_SF_male_antenna <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-antenna","male-antenna"))
res_SF_male_antenna$internal.gene_id <- gene_annotations$internal.gene_id
res_SF_male_antenna$display <- gene_annotations$display.name

# res_SF_male_antenna$Gene <- res_SF_male_antenna$internal.gene_id
res_SF_male_antenna$Gene <- res_SF_male_antenna$display

# volcanoplot(res_SF_male_antenna)


## rostrum
res_SF_male_rostrum <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-rostrum","male-rostrum"))
res_SF_male_rostrum$internal.gene_id <- gene_annotations$internal.gene_id
res_SF_male_rostrum$display <- gene_annotations$display.name

# res_SF_male_rostrum$Gene <- res_SF_male_rostrum$internal.gene_id
res_SF_male_rostrum$Gene <- res_SF_male_rostrum$display

# volcanoplot(res_SF_male_rostrum)

## hindlegs
res_SF_male_hindlegs <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-hindlegs","male-hindlegs"))
res_SF_male_hindlegs$internal.gene_id <- gene_annotations$internal.gene_id
res_SF_male_hindlegs$display <- gene_annotations$display.name

# res_SF_male_hindlegs$Gene <- res_SF_male_hindlegs$internal.gene_id
res_SF_male_hindlegs$Gene <- res_SF_male_hindlegs$display

# volcanoplot(res_SF_male_hindlegs)

## midlegs
res_SF_male_midlegs <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-midlegs","male-midlegs"))
res_SF_male_midlegs$internal.gene_id <- gene_annotations$internal.gene_id
res_SF_male_midlegs$display <- gene_annotations$display.name

# res_SF_male_midlegs$Gene <- res_SF_male_midlegs$internal.gene_id
res_SF_male_midlegs$Gene <- res_SF_male_midlegs$display

# volcanoplot(res_SF_male_midlegs)

## forelegs
res_SF_male_forelegs <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-forelegs","male-forelegs"))
res_SF_male_forelegs$internal.gene_id <- gene_annotations$internal.gene_id
res_SF_male_forelegs$display <- gene_annotations$display.name

# res_SF_male_forelegs$Gene <- res_SF_male_forelegs$internal.gene_id
res_SF_male_forelegs$Gene <- res_SF_male_forelegs$display

# volcanoplot(res_SF_male_forelegs)


res_SF_male_antenna <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-antenna","male-antenna"))
plotMA(res_SF_male_antenna,ylim=c(-5,5))

res_SF_male_brain <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-brain","male-brain"))
plotMA(res_SF_male_brain,ylim=c(-5,5))

res_SF_male_hindlegs <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-hindlegs","male-hindlegs"))
plotMA(res_SF_male_hindlegs,ylim=c(-5,5))

res_SF_male_midlegs <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-midlegs","male-midlegs"))
plotMA(res_SF_male_midlegs,ylim=c(-5,5))

res_SF_male_forelegs <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-forelegs","male-forelegs"))
plotMA(res_SF_male_forelegs,ylim=c(-5,5))

res_SF_male_rostrum <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-rostrum","male-rostrum"))
plotMA(res_SF_male_rostrum,ylim=c(-5,5))

res_SF_male_abdominal_tip <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-abdominaltip","male-abdominaltip"))
plotMA(res_SF_male_abdominal_tip,ylim=c(-5,5))



#### nix, myo-sex
res_SF_male_brain["gene16037",] -> nix
res_SF_male_brain["gene16140",] -> myosex




