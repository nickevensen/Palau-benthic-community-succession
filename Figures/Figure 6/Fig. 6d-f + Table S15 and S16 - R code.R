#############################################
#### Ordination and statistical analyses #### 
#############################################

library(vegan)
library(ggplot2)

## R script to create Figures 6d, 6e and 6f

#input data

exp2_asv <- read.csv('micro_ASVs_exp2.csv',header=TRUE,row.names=1)
exp2_asv <- t(sqrt(exp2_asv))
exp2_metadata <- read.csv('micro_metadata_exp2.csv',header=TRUE,row.names=1)
exp2_metadata$dominant <- factor(exp2_metadata$dominant)
exp2_benthos <- read.csv('micro_macro_correlation_exp2.csv',header=TRUE,row.names=1)
exp2_settlement <- read.csv('micro_settlement_correlation_exp2.csv',header=TRUE,row.names=1)
exp2_asv_orders <- read.csv('micro_ASVs_exp2_orders.csv',header=TRUE,row.names=1) #includes all orders with relative abundance >= 1% in at least one of the samples
exp2_asv_orders <- t(sqrt(exp2_asv_orders))


#permanova

adonis2(exp2_asv~factor(nutrients)*factor(herbivore),data=exp2_metadata,method="bray",perm=999) #results are shown in Table S15


#ordination analysis

exp2_rda <- rda(exp2_asv) #PC1 = 20.0%, PC2 = 11.0%

exp2_rda_treatment <- rda(exp2_asv ~ factor(treatment), data=exp2_metadata) # 34.5% of community variation explained by treatment
exp2_rda_herbivore <- rda(exp2_asv ~ factor(herbivore), data=exp2_metadata) # 18.9% of community variation explained by caging
exp2_rda_nutrient <- rda(exp2_asv ~ factor(nutrients), data=exp2_metadata) # 8.6% of community variation explained by nutrients


#environmental parameter fitting of macro- and micro-benthic community and settlement

exp2_rda_benthos <- envfit(exp2_rda, exp2_benthos, perm=999, na.rm = TRUE) #results are shown in Table S16

exp2_rda_orders <- envfit(exp2_rda, exp2_asv_orders, perm=999, na.rm = TRUE) #top 10 orders with p < 0.05 and highest R2 values can be found in file micro_ASVs_exp2_orders_top10.csv

exp2_rda_settlement <- envfit(exp2_rda, exp2_settlement, perm=999, na.rm = TRUE) #p > 0.05


#obtain pca scores

exp2.scores <- as.data.frame(scores(exp2_rda)$sites)
exp2.scores[,"Treatment"]<-as.character(NA)
exp2.data.scores <- cbind(exp2.scores,exp2_metadata)
exp2.data.scores$Treatment[Name=(exp2.data.scores$treatment == "AO")]<-"Ambient-Open"
exp2.data.scores$Treatment[Name=(exp2.data.scores$treatment == "AC")]<-"Ambient-Caged"
exp2.data.scores$Treatment[Name=(exp2.data.scores$treatment == "NO")]<-"Enriched-Open"
exp2.data.scores$Treatment[Name=(exp2.data.scores$treatment == "NC")]<-"Enriched-Caged"


# plot figure 6d

cols_treatment <- c("Ambient-Caged" = "navy", "Ambient-Open" = "steelblue1", "Enriched-Caged" = "darkgreen", "Enriched-Open" = "palegreen3")
sizes <- c("Ambient-Caged" = 15, "Ambient-Open" = 16, "Enriched-Caged" = 17, "Enriched-Open" = 18)
cent_treatment <- aggregate(cbind(PC1, PC2) ~ Treatment, data = exp2.data.scores, FUN = mean)
ggplot(exp2.data.scores,aes(x = PC1, y = PC2, col=Treatment)) + scale_colour_manual(values=cols_treatment) +
  geom_point(size = 2.5, stroke=1, aes(shape=Treatment)) +  scale_shape_manual(values=c(15, 16, 17, 18)) + 
  stat_ellipse(geom = "polygon", alpha = 0.25, fill="grey", size= 1, level = 0.50) + coord_fixed() + theme_bw() + 
  geom_point(data = cent_treatment, shape=(values=c(22, 21, 24, 23)), size = 4, stroke = 2, colour=(values=c("navy","steelblue1","darkgreen","palegreen3")), fill="black") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme_classic() + 
  theme(text = element_text(size=15, face="bold")) + theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background = element_rect(colour = "black", size=1)) + theme(axis.ticks.length=unit(.2,"cm")) + 
  scale_x_continuous(limits=c(-3,3)) + scale_y_continuous(limits=c(-3,3))

# plot figure 6e

cols_benthos <- c("EAM" = "#FDE0EF", "NCCA" = "#C51B7D", "CCA" = "#DE77AE", "Lobophora" = "#F1B6DA", "Turf" = "#4D9221", "unknown" = "#F2EADE")
cent_dominant <- aggregate(cbind(PC1, PC2) ~ dominant, data = exp2.data.scores, FUN = mean)
ggplot(exp2.data.scores,aes(x = PC1, y = PC2, col=dominant)) + scale_colour_manual(values=cols_benthos) +
  geom_point(size = 2.5, stroke=1, aes(shape=Treatment)) +  scale_shape_manual(values=c(15, 16, 17, 18)) + 
  stat_ellipse(geom = "polygon", alpha = 0.25, fill="grey", size= 1, level = 0.50) + coord_fixed() + theme_bw() + 
  geom_point(data = cent_dominant, size = 4, stroke = 2, colour=(values=c("#DE77AE","#FDE0EF","#F1B6DA","#C51B7D","#4D9221","#F2EADE")), fill="black") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme_classic() + 
  theme(text = element_text(size=15, face="bold")) + theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background = element_rect(colour = "black", size=1)) + theme(axis.ticks.length=unit(.2,"cm")) + 
  scale_x_continuous(limits=c(-3,3)) + scale_y_continuous(limits=c(-3,3))


#create and plot benthic, settlement and microbial vectors from environmental parameter fitting

#benthic vectors

plot(exp2_rda, xlim=c(-2,2), ylim=c(-2,2))
plot(exp2_rda_benthos, col='black') + abline(h = 0, v = 0, col = "black", lty = 2) 

#settlement vector

plot(exp2_rda, xlim=c(-2,2), ylim=c(-2,2))
plot(exp2_rda_settlement, col='black') + abline(h = 0, v = 0, col = "black", lty = 2) 

#otu vectors

exp2_order_signfit <- read.csv('micro_ASVs_exp2_orders_top10.csv', header=TRUE, row.names=1)
exp2_order_signfit <- t(sqrt(exp2_order_signfit))
exp2_rda_orders_signfit <- envfit(exp2_rda, exp2_order_signfit, perm=999, na.rm = TRUE) 

plot(exp2_rda, xlim=c(-2,2), ylim=c(-2,2))
plot(exp2_rda_orders_signfit, col='black') + abline(h = 0, v = 0, col = 'black', lty = 2) 


#PCAs and vector plots were merged and modified in Illustrator to create Figures 6d, 6e and 6f, and combined with Figures 6a, 6b and 6c from the macro-benthic community analysis



