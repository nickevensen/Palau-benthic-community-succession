#############################
#### Ordination analyses #### 
#############################

library(vegan)
library(ggplot2)

## R script to create Figure S5

#input data

exp1_asv <- read.csv('micro_ASVs_exp1.csv', header=TRUE, row.names=1)
exp1_asv <- t(sqrt(exp1_asv))
exp1_metadata <- read.csv('micro_metadata_exp1.csv', header=TRUE, row.names=1)


#obtain pca scores

exp1_rda <- (rda(exp1_asv)) #PC1 = 14.1%, PC2 = 7.2%
exp1.scores <- as.data.frame(scores(exp1_rda)$sites)
exp1.scores[,"Treatment"]<-as.character(NA)
exp1.data.scores <- cbind(exp1.scores, exp1_metadata)
exp1.data.scores$Treatment[Name=(exp1.data.scores$treatment == "AO")]<-"Ambient-Open"
exp1.data.scores$Treatment[Name=(exp1.data.scores$treatment == "AC")]<-"Ambient-Caged"
exp1.data.scores$Treatment[Name=(exp1.data.scores$treatment == "NO")]<-"Enriched-Open"
exp1.data.scores$Treatment[Name=(exp1.data.scores$treatment == "NC")]<-"Enriched-Caged"


# plot figure S3a 

cols_treatment <- c("Ambient-Caged" = "navy", "Ambient-Open" = "steelblue1", "Enriched-Caged" = "darkgreen", "Enriched-Open" = "palegreen3")
cent_treatment <- aggregate(cbind(PC1, PC2) ~ Treatment, data = exp1.data.scores, FUN = mean)
ggplot(exp1.data.scores,aes(x = PC1, y = PC2, col=Treatment)) + scale_colour_manual(values=cols_treatment) +
  geom_point(size = 2.5, stroke=1, aes(shape=treatment)) +  scale_shape_manual(values=c(15, 16, 17, 18)) + 
  stat_ellipse(geom = "polygon", alpha = 0.25, fill="grey", size= 1, level = 0.50) + coord_fixed() + theme_bw() + 
  geom_point(data = cent_treatment, shape=(values=c(22, 21, 24, 23)), size = 4, stroke = 2, colour=(values=c("navy","steelblue1","darkgreen","palegreen3")), fill="black") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme_classic() + 
  theme(text = element_text(size=15, face="bold")) + theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background = element_rect(colour = "black", size=1)) + theme(axis.ticks.length=unit(.2,"cm")) + 
  scale_x_continuous(limits=c(-2,2)) + scale_y_continuous(limits=c(-2,2))

# plot figure S3b

cols_benthos <- c("Bare" = "#FDE0EF", "CCA" = "#C51B7D", "CCA_EAM" = "#DE77AE", "EAM" = "#F1B6DA", "Turf" = "#4D9221", "Thick_turf" = "#7FBC41", "Long_turf" = "#E6F5D0", "Macroalgae" = "#B8E186")
cent_dominant <- aggregate(cbind(PC1, PC2) ~ dominant, data = exp1.data.scores, FUN = mean)
ggplot(exp1.data.scores,aes(x = PC1, y = PC2, col=dominant)) + scale_colour_manual(values=cols_benthos) +
  geom_point(size = 2.5, stroke=1, aes(shape=Treatment)) +  scale_shape_manual(values=c(15, 16, 17, 18)) + 
  stat_ellipse(geom = "polygon", alpha = 0.25, fill="grey", size= 1, level = 0.50) + coord_fixed() + theme_bw() + 
  geom_point(data = cent_dominant, size = 4, stroke = 2, colour=(values=c("#FDE0EF","#C51B7D","#DE77AE","#F1B6DA","#E6F5D0","#B8E186","#7FBC41","#4D9221")), fill="black") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme_classic() + 
  theme(text = element_text(size=15, face="bold")) + theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background = element_rect(colour = "black", size=1)) + theme(axis.ticks.length=unit(.2,"cm")) + 
  scale_x_continuous(limits=c(-2,2)) + scale_y_continuous(limits=c(-2,2))

# plot figure S3c 

cols_time <- c("T1" = "#f3ebdf", "T2" = "#ddd4c5", "T3" = "#c7bdab", "T4" = "#aa9e8a", "T5" = "#6a5e4c", "T6" = "#443928")
sizes <- c("Ambient-Caged" = 15, "Ambient-Open" = 16, "Enriched-Caged" = 17, "Enriched-Open" = 18)
cent_time <- aggregate(cbind(PC1, PC2) ~ time, data = exp1.data.scores, FUN = mean)
ggplot(exp1.data.scores,aes(x = PC1, y = PC2, col=time)) + scale_colour_manual(values=cols_time) +
  geom_point(size = 2.5, stroke=1, aes(shape=Treatment)) +  scale_shape_manual(values=c(15, 16, 17, 18)) + 
  stat_ellipse(geom = "polygon", alpha = 0.25, fill="grey", size= 1, level = 0.50) + coord_fixed() + theme_bw() + 
  geom_point(data = cent_time, size = 4, stroke = 2, colour=(values=c("#f3ebdf","#ddd4c5","#c7bdab","#aa9e8a","#6a5e4c","#443928")), fill="black") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme_classic() + 
  theme(text = element_text(size=15, face="bold")) + theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background = element_rect(colour = "black", size=1)) + theme(axis.ticks.length=unit(.2,"cm")) + 
  scale_x_continuous(limits=c(-2,2)) + scale_y_continuous(limits=c(-2,2))


#PCAs were combined and modified in Illustrator to create Figure S5

