#Code for Figure 6a-c of "Benthic micro- and macro-community succession and coral recruitment under overfishing and nutrient enrichment"
#Also includes statistical output from Table S13 and S14

#setwd to source file location (data file and R script are in the same folder)
Settlement_cover_full<-read.csv("Settlement_cover_full.csv")

dataIR<-transform(Settlement_cover_full, 
                  Nutrient=as.factor(Nutrient), 
                  Caging=as.factor(Caging),
                  Microhabitat=as.factor(Microhabitat))
str(dataIR)
names(dataIR)
attach(dataIR)

##create separate data sets
library(dplyr)
## sometimes the select function is acting weird, if errors message, close R and re opens and start again from here after reading dataset, should work
dataIR_com<-dplyr::select(dataIR,CCA,Lobophora,Sand,Turf,Thick_turf,
                          NCCA,Macroalgae,EAM,Invert)
dataIR_meta<-dplyr::select(dataIR,Nutrient,Caging,Microhabitat)

#make distance matrix
library(vegan)
DISTdataIR_com<- vegdist(dataIR_com, method="bray")

# ADONIS function
IR<-adonis(DISTdataIR_com~Nutrient*Caging*Microhabitat, data=dataIR_meta)
IR

library(rtf)
tableS13<-as.data.frame(IR$aov.tab, keep.rownames = TRUE)
tableS13$Source <- rownames(tableS13)
tableS13<-format(tableS13, digits = 3)

rtffile13 <- RTF("rtf13.doc")
addParagraph(rtffile13, "This is the output for Table S13:\n")
addTable(rtffile13, tableS13)
done(rtffile13)

#pairwise comparisons
library(pairwiseAdonis)
str(dataIR)

pairwise.adonis(dataIR_com,factors=dataIR$Caging,sim.function='vegdist',sim.method='bray',p.adjust.m='bonferroni') 

#Prepare data for plotting
Sol <- metaMDS(DISTdataIR_com)
Sol
stressplot(Sol)

data.scores <- as.data.frame(scores(Sol))
new.data.scores <- cbind(data.scores,dataIR_meta)

#Create new column to pair treatment combinations to plot all treatments separately
new.data.scores[,"treatment"]<-as.character(NA)
#crowns
new.data.scores$treatment[Name=(new.data.scores$Nutrient == "Enriched") & (new.data.scores$Caging=="Open") & (new.data.scores$Microhabitat=="Crown")]<-"Enriched-Open-Crown"
new.data.scores$treatment[Name=(new.data.scores$Nutrient == "Enriched") & (new.data.scores$Caging=="Caged") & (new.data.scores$Microhabitat=="Crown")]<-"Enriched-Caged-Crown"

new.data.scores$treatment[Name=(new.data.scores$Nutrient == "Ambient") & (new.data.scores$Caging=="Open") & (new.data.scores$Microhabitat=="Crown")]<-"Ambient-Open-Crown"
new.data.scores$treatment[Name=(new.data.scores$Nutrient == "Ambient") & (new.data.scores$Caging=="Caged") & (new.data.scores$Microhabitat=="Crown")]<-"Ambient-Caged-Crown"

#crevices
new.data.scores$treatment[Name=(new.data.scores$Nutrient == "Enriched") & (new.data.scores$Caging=="Open") & (new.data.scores$Microhabitat=="Crevice")]<-"Enriched-Open-Crevice"
new.data.scores$treatment[Name=(new.data.scores$Nutrient == "Enriched") & (new.data.scores$Caging=="Caged") & (new.data.scores$Microhabitat=="Crevice")]<-"Enriched-Caged-Crevice"

new.data.scores$treatment[Name=(new.data.scores$Nutrient == "Ambient") & (new.data.scores$Caging=="Open") & (new.data.scores$Microhabitat=="Crevice")]<-"Ambient-Open-Crevice"
new.data.scores$treatment[Name=(new.data.scores$Nutrient == "Ambient") & (new.data.scores$Caging=="Caged") & (new.data.scores$Microhabitat=="Crevice")]<-"Ambient-Caged-Crevice"

#subset for plotting
crown.data.scores<-subset(new.data.scores, Microhabitat=="Crown")
crevice.data.scores<-subset(new.data.scores, Microhabitat=="Crevice")

#Plot for crown communities
cent <- aggregate(cbind(NMDS1, NMDS2) ~ treatment, data = crown.data.scores, FUN = mean)
segs <- merge(crown.data.scores, setNames(cent, c('treatment','oNMDS1','oNMDS2')),
              by = 'treatment', sort = FALSE)

cols <- c("Ambient-Caged-Crown" = "navy","Ambient-Open-Crown" = "steelblue1", 
          "Enriched-Caged-Crown" = "darkgreen", "Enriched-Open-Crown" = "palegreen3")
library(ggplot2)
ggplot(crown.data.scores, aes(x = NMDS1, y = NMDS2, col=treatment)) + scale_colour_manual(values=cols) +
  stat_ellipse(geom = "polygon", alpha = 0.25, fill="grey", size= 1, level = 0.50) + # ellipses encompassing 50% of replicates
  geom_point(data = cent, shape=(values=c(22, 24, 22, 24)), size = 5, stroke = 2, colour=(values=c("navy","steelblue1","darkgreen","palegreen3")), fill="black") +       # centroids
  geom_point(size = 2.5, stroke=1, aes(shape=treatment)) +  scale_shape_manual(values=c(15, 17, 15, 17)) +                                     # sample scores
  coord_fixed() +                                              # same axis scaling
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + geom_vline(xintercept=0, linetype="dashed") + geom_hline(yintercept=0, linetype="dashed") +
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) + scale_y_continuous(limits=c(-0.4,0.5), breaks = c(-0.4,-0.2,0,0.2,0.4)) + scale_x_continuous(limits=c(-0.6,0.6), breaks = c(-0.4,-0.2,0,0.2,0.4,0.6))

#Plot for crevice communities
cent <- aggregate(cbind(NMDS1, NMDS2) ~ treatment, data = crevice.data.scores, FUN = mean)
segs <- merge(crevice.data.scores, setNames(cent, c('treatment','oNMDS1','oNMDS2')),
              by = 'treatment', sort = FALSE)

cols <- c("Ambient-Caged-Crevice" = "navy","Ambient-Open-Crevice" = "steelblue1", 
          "Enriched-Caged-Crevice" = "darkgreen", "Enriched-Open-Crevice" = "palegreen3")

library(ggplot2)
ggplot(crevice.data.scores, aes(x = NMDS1, y = NMDS2, col=treatment)) + scale_colour_manual(values=cols) +
  stat_ellipse(geom = "polygon", alpha = 0.25, fill="grey", size= 1, level = 0.50) + # ellipses encompassing 50% of replicates
  geom_point(data = cent, shape=(values=c(22, 24, 22, 24)), size = 5, stroke = 2, colour=(values=c("navy","steelblue1","darkgreen","palegreen3")), fill="black") +       # centroids
  geom_point(size = 2.5, stroke=1, aes(shape=treatment)) +  scale_shape_manual(values=c(15, 17, 15, 17)) +                                     # sample scores
  coord_fixed() +                                              # same axis scaling
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + geom_vline(xintercept=0, linetype="dashed") + geom_hline(yintercept=0, linetype="dashed") +
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) + scale_y_continuous(limits=c(-0.4,0.5), breaks = c(-0.4,-0.2,0,0.2,0.4)) + scale_x_continuous(limits=c(-0.5,0.6), breaks = c(-0.4,-0.2,0,0.2,0.4,0.6))

# That's the main plots done. Now to plot the vectors to show what's driving difference in ordinations (i.e., coordinates of the point)
library(MASS)
library(vegan)
ord_IR<-metaMDS(DISTdataIR_com,distance = "bray",
                k=2,trymax=1000,autotransform=TRUE,expand=FALSE, plot=FALSE)
ord_IR
fit<-envfit(ord_IR, dataIR_com, perm=999)
fit
scores(fit, "vectors")
plot(ord_IR,col="white",axes=FALSE, ylim=c(-0.4, 0.5), xlim=c(-0.5, 0.6))
plot(fit, p.max = 0.95, col = "black") + abline(h = 0, v = 0, col = "black", lty = 2)
axis(side=1, at=seq(-0.5, 0.6, by=0.1))
axis(side=2, at=seq(-0.4, 0.5, by=0.1))
box()

tableS14 <- data.frame((fit$vectors)$arrows, (fit$vectors)$r, (fit$vectors)$pvals)
tableS14$Source <- rownames(tableS14)
tableS14<-format(tableS14, digits = 3)

rtffile14 <- RTF("rtf14.doc")
addParagraph(rtffile14, "This is the output for Table S14:\n")
addTable(rtffile14, tableS14)
done(rtffile14)

#END