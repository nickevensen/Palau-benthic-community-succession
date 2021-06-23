######################################################
#### Micro-benthic community composition heatmaps #### 
######################################################

library(pheatmap)
library(RColorBrewer)

## R script to create Figure S9

#input data

exp2_HM <- read.csv('HM_micro_benthos_exp2_order.csv', row.names=1, header=TRUE) # contains only orders that are present at >1% in at least one of the samples
exp2_rowmeta <- read.csv('HM_micro_benthos_exp2_rowmeta.csv', row.names=1, header=TRUE)
exp2_colmeta <- read.csv('HM_micro_benthos_exp2_colmeta.csv', row.names=1, header=TRUE)


#create heatmap

breaks=c(0.1,1,2,4,6,10,15,20)

exp2_heatmap <- pheatmap(exp2_HM,color = brewer.pal(n = 8, name = "YlGnBu"),
                         cluster_cols=FALSE,cluster_rows=FALSE,scale="none",breaks=breaks,
                         annotation_row=exp2_rowmeta,annotation_col=exp2_colmeta,
                         fontsize_row=8,fontsize_col=6,angle_col=45,cellwidth=8,border_color="grey87",
                         main="Micro-benthic community composition")


#the heatmap was further modified and colours updated in Illustrator to create Figure S9 