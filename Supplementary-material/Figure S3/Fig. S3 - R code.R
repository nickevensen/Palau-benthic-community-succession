######################################################
#### Macro-benthic community composition heatmaps #### 
######################################################

library(pheatmap)
library(RColorBrewer)

## R script to create Figure S3

#input data

benthic_combined2 <- read.csv('HM_macro_benthos_combined2.csv', row.names=1, header=TRUE)
benthic_rows2 <- read.csv('HM_macro_benthos_rownames2.csv', row.names=1, header=TRUE)
benthic_cols2 <- read.csv('HM_macro_benthos_colnames2.csv', row.names=1, header=TRUE)


#create heatmap

breaks=c(0,1,20,40,60,80,100)
HM_combined <- pheatmap(benthic_combined2,color = brewer.pal(n = 6, name = "YlGnBu"),
                        cluster_cols=FALSE,cluster_rows=FALSE,treeheight_row=100,scale="none",
                        breaks=breaks,annotation_row=benthic_rows2,annotation_col=benthic_cols2,
                        fontsize_row=8,fontsize_col=6,angle_col=45,cellwidth=5,cellheight=9,border_color="grey87",
                        main="Macro-benthic community composition")


#the heatmap was further modified and colours updated in Illustrator to create Figure S3 






