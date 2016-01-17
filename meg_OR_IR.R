require(gridExtra)
require(dplyr)
require(reshape2)
require(ggplot2)

or_mean = compiled_means[compiled_means$gene.family %in% c("OR"),]
ir_mean = compiled_means[compiled_means$gene.family %in% c("IR"),]

#combo_mean = compiled_means[compiled_means$gene.family %in% c("IR","OR"),]
combo_mean = compiled_means[compiled_means$gene.family %in% c("IR","OR","GR"),]
combo_mean_sub = combo_mean[,c(3:4,9:25)]
combo_mean_sub_melt = melt(combo_mean_sub)

total_thresh = data.frame(thresh = 0,ir=74,or=129)
temp_thresh = data.frame(thresh = 0,ir=0,or=0)

# ECDF plots
THRESH <- 0
temp <- as.data.frame(combo_mean_sub_melt[combo_mean_sub_melt$value > THRESH,])

temp_sub <- temp[temp$variable %in% c("Fe_An_SF","Ma_An","Fe_Pa_SF"),]
temp_sub <- transform(temp_sub,y=rep(1,nrow(temp_sub)))
temp_sub$combo = paste(temp_sub$gene.family,temp_sub$variable)
temp_sub <- temp_sub[order(temp_sub$value),]

require(plyr)

DF.t <- ddply(temp_sub, .(combo), transform, cy = cumsum(y))
DF.t <- ddply(DF.t, .(combo), transform, cyn = max(cy)-cy)

ggplot(DF.t, aes(x=log10(value), y=cy, colour=combo, group=combo)) + geom_point() + facet_grid(variable ~ gene.family) + theme_bw()

final_meg_plot <- ggplot(DF.t, aes(x=log10(value), y=cyn, colour=variable, group=combo)) + geom_point(aes(shape=gene.family)) + facet_grid(. ~ variable,scales="free_y") + theme_bw()

temp_facet <- ggplot(temp_sub, aes(log10(value), fill = variable)) + geom_point(y=) + facet_grid(gene.family * variable ~ . )
dat <-  ggplot_build(temp_facet)$data[[1]]
library(data.table)
ggplot(setDT(dat)[,y:=cumsum(y),"PANEL"],aes(x=x)) +
  geom_point(aes(y=y,fill=PANEL),stat="identity")+facet_grid(PANEL~.) +
  guides(title="Modul")


THRESH <- 0
temp <- as.data.frame(combo_mean_sub_melt[combo_mean_sub_melt$value > THRESH,])
t0 <- ggplot(temp,aes(x=value,group=variable,fill=gene.family)) + geom_histogram() + ylim(0,200) + theme_classic()
temp_thresh$thresh = 0.0000000001
temp_thresh$or = length(unique(temp[temp$gene.family %in% c("OR"),]$display.name))
temp_thresh$ir = length(unique(temp[temp$gene.family %in% c("IR"),]$display.name))
total_thresh <- rbind(total_thresh,temp_thresh)

THRESH <- 0.5
temp <- as.data.frame(combo_mean_sub_melt[combo_mean_sub_melt$value > THRESH,])
t0.5 <- ggplot(temp,aes(x=variable,group=gene.family,fill=gene.family)) + geom_histogram() + ylim(0,200) + theme_classic()
temp_thresh$thresh = THRESH
temp_thresh$or = length(unique(temp[temp$gene.family %in% c("OR"),]$display.name))
temp_thresh$ir = length(unique(temp[temp$gene.family %in% c("IR"),]$display.name))
total_thresh <- rbind(total_thresh,temp_thresh)
                      
THRESH <- 1
temp <- as.data.frame(combo_mean_sub_melt[combo_mean_sub_melt$value > THRESH,])
t1 <- ggplot(temp,aes(x=variable,group=gene.family,fill=gene.family)) + geom_histogram() + ylim(0,200) + theme_classic()
temp_thresh$thresh = THRESH
temp_thresh$or = length(unique(temp[temp$gene.family %in% c("OR"),]$display.name))
temp_thresh$ir = length(unique(temp[temp$gene.family %in% c("IR"),]$display.name))
total_thresh <- rbind(total_thresh,temp_thresh)

THRESH <- 2
temp <- as.data.frame(combo_mean_sub_melt[combo_mean_sub_melt$value > THRESH,])
t2 <- ggplot(temp,aes(x=variable,group=gene.family,fill=gene.family)) + geom_histogram() + ylim(0,200) + theme_classic()
temp_thresh$thresh = THRESH
temp_thresh$or = length(unique(temp[temp$gene.family %in% c("OR"),]$display.name))
temp_thresh$ir = length(unique(temp[temp$gene.family %in% c("IR"),]$display.name))
total_thresh <- rbind(total_thresh,temp_thresh)

THRESH <- 5
temp <- as.data.frame(combo_mean_sub_melt[combo_mean_sub_melt$value > THRESH,])
t5 <- ggplot(temp,aes(x=variable,group=gene.family,fill=gene.family)) + geom_histogram() + ylim(0,200) + theme_classic()
temp_thresh$thresh = THRESH
temp_thresh$or = length(unique(temp[temp$gene.family %in% c("OR"),]$display.name))
temp_thresh$ir = length(unique(temp[temp$gene.family %in% c("IR"),]$display.name))
total_thresh <- rbind(total_thresh,temp_thresh)

THRESH <- 10
temp <- as.data.frame(combo_mean_sub_melt[combo_mean_sub_melt$value > THRESH,])
t10 <- ggplot(temp,aes(x=variable,group=gene.family,fill=gene.family)) + geom_histogram(stat="count") + ylim(0,200) + theme_classic()
temp_thresh$thresh = THRESH
temp_thresh$or = length(unique(temp[temp$gene.family %in% c("OR"),]$display.name))
temp_thresh$ir = length(unique(temp[temp$gene.family %in% c("IR"),]$display.name))
total_thresh <- rbind(total_thresh,temp_thresh)

grid.arrange(t0,t0.5,t1,t2,t5,t10)

total_num <- ggplot(data=melt(total_thresh,id.vars="thresh"),aes(x=thresh,y=value,group=variable,fill=variable)) + geom_point() + facet_grid(.~variable) + ylim(0,150)
