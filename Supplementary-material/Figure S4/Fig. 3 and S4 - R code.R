#Code for Figure 3 of "Benthic micro- and macro-community succession and coral recruitment under overfishing and nutrient enrichment"

#setwd to source file Microhabitat & read in data
Macro_communities<-read.csv("Macro_communities.csv")

dataIR<-transform(Macro_communities, 
                  Nutrients=as.factor(Nutrients), 
                  Caging=as.factor(Caging),
                  Microhabitat=as.factor(Microhabitat),
                  Timepoint=as.factor(Timepoint))

# prepare data sets for PERMANOVA
library(dplyr)
## sometimes the select function is acting weird, if errors message, close R and re opens and start again from here after reading dataset, should work
dataIR_com<-dplyr::select(dataIR, Bare, CCA, Cyanobacteria, EAM, Inverts, Lobophora, 
                          Long_turf, Macroalgae, NCCA, Sand, Turf, Thick_turf)
dataIR_meta<-dplyr::select(dataIR,Nutrients,Caging,Timepoint,Microhabitat)

#Bray–Curtis dissimilarity
DISTdataIR_com<- vegdist(dataIR_com, method="bray")

#PERMANOVA - adonis function
Perm <- adonis(DISTdataIR_com ~ Caging*Nutrients*Microhabitat*Timepoint, data=dataIR_meta)
Perm

Perm2<- adonis(DISTdataIR_com ~ Caging*Nutrients*Microhabitat + Caging*Nutrients*Timepoint +
                 Caging*Microhabitat*Timepoint + Nutrients*Microhabitat*Timepoint, data=dataIR_meta)
Perm2

Perm3<- adonis(DISTdataIR_com ~ Caging*Nutrients + Caging*Microhabitat + Nutrients*Microhabitat + Caging*Timepoint +
                 Nutrients*Timepoint + Microhabitat*Timepoint, data=dataIR_meta)
Perm3 #this is the simplest model when doing a simple drop of non-significant terms in R

#save model output to table - Table S1
library(rtf)
tableS1<-as.data.frame(Perm3$aov.tab, keep.rownames = TRUE)
tableS1$Source <- rownames(tableS1)
tableS1<-format(tableS1, digits = 2)

rtffile <- RTF("rtf.doc")
addParagraph(rtffile, "This is the output for Table S1:\n")
addTable(rtffile, tableS1)
done(rtffile)

#### Plotting the tile community data through time  ####
## average data set for ordination
library(plyr)
dataIR_AV<-ddply(dataIR, .(Timepoint,Nutrients,Caging,Microhabitat),summarize,
                 Bare=mean(Bare),
                 CCA=mean(CCA),
                 Cyanobacteria=mean(Cyanobacteria),
                 EAM=mean(EAM),
                 Inverts=mean(Inverts),
                 Lobophora=mean(Lobophora),
                 Long_turf=mean(Long_turf),
                 Macroalgae=mean(Macroalgae),
                 NCCA=mean(NCCA),
                 Sand=mean(Sand),
                 Thick_turf=mean(Thick_turf),
                 Turf=mean(Turf))

library(dplyr)
dataIRAV_com<-dplyr::select(dataIR_AV, Bare, CCA, Cyanobacteria, EAM, Inverts, Lobophora, 
                            Long_turf, Macroalgae, NCCA, Sand, Turf, Thick_turf)
dataIRAV_meta<-dplyr::select(dataIR_AV,Timepoint,Nutrients,Caging,Microhabitat)

#Analyse all together, but present separating crowns and crevices - presentation matches up with KM plots.
Sol <- metaMDS(dataIRAV_com)
Sol
stressplot(Sol)

data.scores <- as.data.frame(scores(Sol))
new.data.scores <- cbind(data.scores,dataIRAV_meta)

#merge data to plot the combination of the plot x tile treatments
new.data.scores[,"treatment"]<-as.character(NA)
new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Normal") & (new.data.scores$Caging=="Open")]<-"Normal-Open"
new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Normal") & (new.data.scores$Caging=="Caged")]<-"Normal-Caged"
new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Normal") & (new.data.scores$Caging=="Partial")]<-"Normal-Partial"
new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Nutrient") & (new.data.scores$Caging=="Open")]<-"Nutrient-Open"
new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Nutrient") & (new.data.scores$Caging=="Caged")]<-"Nutrient-Caged"
new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Nutrient") & (new.data.scores$Caging=="Partial")]<-"Nutrient-Partial"

#Here - subset and plot crowns and crevices separately.
Crown.data.scores<-subset(new.data.scores, Microhabitat=="Crown")
attach(Crown.data.scores)

Crevice.data.scores<-subset(new.data.scores, Microhabitat=="Crevice")
attach(Crevice.data.scores)

# plot with all treatments with lines joining the timepoints within a treatment and labels showing what timepoint each point is (colours represent treatments)
cols <- c("Normal-Caged" = "navy", "Normal-Partial" = "dodgerblue4", "Normal-Open" = "steelblue1",
          "Nutrient-Caged" = "darkgreen", "Nutrient-Partial" = "palegreen4", "Nutrient-Open" = "palegreen3")
# colour list: 'navy','steelblue1','dodgerblue4','darkgreen','palegreen3','palegreen4'
library(ggplot2)

#crowns - Fig. 3a
ggplot(Crown.data.scores, aes(x = NMDS1, y = NMDS2, col=treatment, shape=treatment)) + scale_colour_manual(values=cols) + 
  geom_point(size=3) + scale_shape_manual(values=c(15, 17, 16, 15, 17, 16)) +
  coord_fixed() + geom_text(aes(label=Timepoint),hjust=0.5, vjust=-1) + 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme_classic() + theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  geom_vline(xintercept=c(0,0), linetype="dotted") + geom_hline(yintercept=c(0,0), linetype="dotted") +
  theme(axis.ticks.length=unit(.2,"cm")) + scale_x_continuous(limits=c(-1.2,1)) + scale_y_continuous(limits=c(-1,1))

#crevices - Fig. S4a
ggplot(Crevice.data.scores, aes(x = NMDS1, y = NMDS2, col=treatment, shape=treatment)) + scale_colour_manual(values=cols) + 
  geom_point(size=3) + scale_shape_manual(values=c(15, 17, 16, 15, 17, 16)) +
  coord_fixed() + geom_text(aes(label=Timepoint),hjust=0.5, vjust=-1) + 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme_classic() + theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  geom_vline(xintercept=c(0,0), linetype="dotted") + geom_hline(yintercept=c(0,0), linetype="dotted") +
  theme(axis.ticks.length=unit(.2,"cm")) + scale_x_continuous(limits=c(-1.2,1)) + scale_y_continuous(limits=c(-1,1))

### NOTE: Vectors connecting the timepoints for each treatment were later added in using Adobe Illustrator.

#plot vectors driving the MDS - Figs. 3c and S4b
library(MASS)
library(vegan)
ord_IR<-metaMDS(dataIRAV_com,distance = "bray",
                k=2,trymax=1000,autotransform=TRUE,expand=FALSE, plot=FALSE)
ord_IR
fit<-envfit(ord_IR, dataIRAV_com, perm=999)
fit
scores(fit, "vectors")
plot(ord_IR,col="white")
plot(fit, p.max = 0.9, col = "black") + abline(h = 0, v = 0, col = "black", lty = 2) + scale_x_continuous(limits=c(-1.2,1)) + scale_y_continuous(limits=c(-1,1))

#save model output to table - Table S2
tableS2 <- data.frame((fit$vectors)$arrows, (fit$vectors)$r, (fit$vectors)$pvals)
tableS2$Source <- rownames(tableS2)
tableS2<-format(tableS2, digits = 2)

rtffile2 <- RTF("rtf2.doc")
addParagraph(rtffile2, "This is the output for Table S2:\n")
addTable(rtffile2, tableS2)
done(rtffile2)

#END