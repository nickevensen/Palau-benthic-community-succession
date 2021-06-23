############################################# 
#### Ordination and statistical analyses #### 
#############################################

library(vegan)
library(ggplot2)

## R script to create Figures 3b and 3d

#input data

exp1_asv <- read.csv('micro_ASVs_exp1.csv', header=TRUE, row.names=1)
exp1_asv <- t(sqrt(exp1_asv))
exp1_metadata <- read.csv('micro_metadata_exp1.csv', header=TRUE, row.names=1)
exp1_metadata$dominant <- factor(exp1_metadata$dominant)
exp1_benthos <- read.csv('micro_macro_correlation_exp1.csv', header=TRUE, row.names=1)
exp1_asv_orders <- read.csv('micro_ASVs_exp1_orders.csv', header=TRUE, row.names=1) #includes all orders with relative abundance >= 1% in at least one of the samples
exp1_asv_orders <- t(sqrt(exp1_asv_orders))


#permanova

adonis2(exp1_asv ~ factor(nutrients) * factor(herbivore) * factor(time), data=exp1_metadata, method="bray", perm=999) #results are shown in Table S3
adonis2(exp1_asv ~ factor(dominant), data=exp1_metadata, method="bray", perm=999) #results are shown in Table S4

#ordination analysis

exp1_rda <- (rda(exp1_asv)) #PC1 = 14.1%, PC2 = 7.2%

exp1_rda_treatment <- rda(exp1_asv ~ factor(treatment), data=exp1_metadata) # 10.6% of community variation explained by treatment
exp1_rda_herbivore <- rda(exp1_asv ~ factor(herbivore), data=exp1_metadata) # 7.6% of community variation explained by caging
exp1_rda_nutrients <- rda(exp1_asv ~ factor(nutrients), data=exp1_metadata) # 1.7% of community variation explained by nutrients
exp1_rda_time <- rda(exp1_asv ~ factor(time), data=exp1_metadata) # 18.6% of community variation explained by time point
exp1_rda_dominant <- rda(exp1_asv ~ factor(dominant), data=exp1_metadata) # 19.2% of community variation explained by dominant benthic group


#environmental parameter fitting of macro- and micro-benthic community

exp1_rda_benthos <- envfit(exp1_rda, exp1_benthos, perm=999, na.rm = TRUE) #results are shown in Table S5; significant correlation for CCA p = 0.001; EAM p = 0.014; Inverts p = 0.022; Macroalgae p = 0.001; Turf p = 0.001

exp1_rda_orders <- envfit(exp1_rda, exp1_asv_orders, perm=999, na.rm = TRUE) #top 10 orders with p < 0.05 and highest R2 values can be found in file micro_ASVs_exp1_orders_top10.csv


#Tukey's HSD - results are shown in Table S6

exp1_rda_score <- scores(exp1_rda)$sites

exp1_rda_tuk1_herbivore <- TukeyHSD(aov(exp1_rda_score[,1] ~ factor(herbivore),data=exp1_metadata))
exp1_rda_tuk2_herbivore <- TukeyHSD(aov(exp1_rda_score[,2] ~ factor(herbivore),data=exp1_metadata))

exp1_rda_tuk1_nutrient <- TukeyHSD(aov(exp1_rda_score[,1] ~ factor(nutrients),data=exp1_metadata))
exp1_rda_tuk2_nutrient <- TukeyHSD(aov(exp1_rda_score[,2] ~ factor(nutrients),data=exp1_metadata))

exp1_rda_tuk <- TukeyHSD(aov(exp1_rda_score[,1] ~ factor(treatment), data=exp1_metadata))
exp1_rda_tuk <- TukeyHSD(aov(exp1_rda_score[,2] ~ factor(treatment), data=exp1_metadata))


#plot PCA

exp1.scores <- as.data.frame(scores(exp1_rda)$sites)
exp1.scores[,"Treatment"]<-as.character(NA)
exp1.data.scores <- cbind(exp1.scores, exp1_metadata)
exp1.data.scores$Treatment[Name=(exp1.data.scores$treatment == "AO")]<-"Ambient-Open"
exp1.data.scores$Treatment[Name=(exp1.data.scores$treatment == "AC")]<-"Ambient-Caged"
exp1.data.scores$Treatment[Name=(exp1.data.scores$treatment == "NO")]<-"Enriched-Open"
exp1.data.scores$Treatment[Name=(exp1.data.scores$treatment == "NC")]<-"Enriched-Caged"

cols <- c("Ambient-Caged" = "navy", "Ambient-Open" = "steelblue1", "Enriched-Caged" = "darkgreen", "Enriched-Open" = "palegreen3")
sizes <- c("Ambient-Caged" = 15, "Ambient-Open" = 16, "Enriched-Caged" = 17, "Enriched-Open" = 18)
cent <- aggregate(cbind(PC1, PC2) ~ time:Treatment, data = exp1.data.scores, FUN = mean)
ggplot(cent, aes(x=PC1, y=PC2, col=Treatment, shape=Treatment)) + 
  scale_colour_manual(values=cols) + geom_point(size = 3) + 
  scale_shape_manual(values=c(15,17,15,17)) + coord_fixed() + 
  geom_text(aes(label=factor(time)),hjust=0,vjust=0) + 
  geom_segment(aes(xend=c(tail(PC1, n=-1), NA), yend=c(tail(PC2, n=-1), NA)), arrow=arrow(length=unit(0.3,"cm")), size=0.7) +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme_classic() + theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) + scale_x_continuous(limits=c(-3,3)) + scale_y_continuous(limits=c(-2,2))


#create and plot benthic and microbial vectors from environmental parameter fitting

#benthic vectors

plot(exp1_rda, xlim=c(-2,2), ylim=c(-2,2))
plot(exp1_rda_benthos, col='black') + abline(h = 0, v = 0, col = "black", lty = 2) 

#otu vectors

exp1_order_signfit <- read.csv('micro_ASVs_exp1_orders_top10.csv', header=TRUE, row.names=1)
exp1_order_signfit <- t(sqrt(exp1_order_signfit))
exp1_rda_orders_signfit <- envfit(exp1_rda, exp1_order_signfit, perm=999, na.rm = TRUE) 

plot(exp1_rda, xlim=c(-2,2), ylim=c(-2,2))
plot(exp1_rda_orders_signfit, col='black') + abline(h = 0, v = 0, col = 'black', lty = 2) 


#PCAs and vector plots were merged and modified in Illustrator to create Figures 3b and 3d, and combined with Figures 3a and 3c from the macro-benthic community analysis
