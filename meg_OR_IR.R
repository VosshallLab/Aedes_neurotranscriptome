require(gridExtra)

or_mean = compiled_means[compiled_means$gene.family %in% c("OR"),]
ir_mean = compiled_means[compiled_means$gene.family %in% c("IR"),]

combo_mean = compiled_means[compiled_means$gene.family %in% c("IR","OR"),]
combo_mean_sub = combo_mean[,c(3:4,9:25)]
combo_mean_sub_melt = melt(combo_mean_sub)

total_thresh = data.frame(thresh = 0,ir=74,or=129)
temp_thresh = data.frame(thresh = 0,ir=0,or=0)


THRESH <- 0
temp <- as.data.frame(combo_mean_sub_melt[combo_mean_sub_melt$value > THRESH,])
t0 <- ggplot(temp,aes(x=variable,group=gene.family,fill=gene.family)) + geom_histogram() + ylim(0,200) + theme_classic()
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
t10 <- ggplot(temp,aes(x=variable,group=gene.family,fill=gene.family)) + geom_histogram() + ylim(0,200) + theme_classic()
temp_thresh$thresh = THRESH
temp_thresh$or = length(unique(temp[temp$gene.family %in% c("OR"),]$display.name))
temp_thresh$ir = length(unique(temp[temp$gene.family %in% c("IR"),]$display.name))
total_thresh <- rbind(total_thresh,temp_thresh)

grid.arrange(t0,t0.5,t1,t2,t5,t10)

total_num <- ggplot(data=melt(total_thresh,id.vars="thresh"),aes(x=thresh,y=value,group=variable,fill=variable)) + geom_point() + facet_grid(.~variable) + ylim(0,150)
