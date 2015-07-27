### generate results


## antenna
res_SF_male_antenna <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-antenna","male-antenna"))
res_SF_male_antenna$internal.gene_id <- gene_annotations$internal.gene_id
res_SF_male_antenna$display <- gene_annotations$display.name

# res_SF_male_antenna$Gene <- res_SF_male_antenna$internal.gene_id
res_SF_male_antenna$Gene <- res_SF_male_antenna$display

volcanoplot(res_SF_male_antenna)


## rostrum
res_SF_male_rostrum <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-rostrum","male-rostrum"))
res_SF_male_rostrum$internal.gene_id <- gene_annotations$internal.gene_id
res_SF_male_rostrum$display <- gene_annotations$display.name

# res_SF_male_rostrum$Gene <- res_SF_male_rostrum$internal.gene_id
res_SF_male_rostrum$Gene <- res_SF_male_rostrum$display

volcanoplot(res_SF_male_rostrum)

## hindlegs
res_SF_male_hindlegs <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-hindlegs","male-hindlegs"))
res_SF_male_hindlegs$internal.gene_id <- gene_annotations$internal.gene_id
res_SF_male_hindlegs$display <- gene_annotations$display.name

# res_SF_male_hindlegs$Gene <- res_SF_male_hindlegs$internal.gene_id
res_SF_male_hindlegs$Gene <- res_SF_male_hindlegs$display

volcanoplot(res_SF_male_hindlegs)

## midlegs
res_SF_male_midlegs <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-midlegs","male-midlegs"))
res_SF_male_midlegs$internal.gene_id <- gene_annotations$internal.gene_id
res_SF_male_midlegs$display <- gene_annotations$display.name

# res_SF_male_midlegs$Gene <- res_SF_male_midlegs$internal.gene_id
res_SF_male_midlegs$Gene <- res_SF_male_midlegs$display

volcanoplot(res_SF_male_midlegs)

## forelegs
res_SF_male_forelegs <- results(dds_SF_male_combined,contrast=c("tissue_condition","SF-forelegs","male-forelegs"))
res_SF_male_forelegs$internal.gene_id <- gene_annotations$internal.gene_id
res_SF_male_forelegs$display <- gene_annotations$display.name

# res_SF_male_forelegs$Gene <- res_SF_male_forelegs$internal.gene_id
res_SF_male_forelegs$Gene <- res_SF_male_forelegs$display

volcanoplot(res_SF_male_forelegs)


