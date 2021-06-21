#Code for Figure 3a,c and Fig. S4 of "Benthic micro- and macro-community succession and coral recruitment under overfishing and nutrient enrichment"
#Also includes statistical output from Table S1, S2, and S12

#setwd to source file location & read in data
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

#### Multiple regressions analyses for Exp. 1 - summary output included in Table S12
#Contains two sections - one looking at all tiles, but not including grazing effects, and one including only open tiles (w/ grazing)

library(corrplot)
library(mctest)
library(nlme)
library(lme4)
library(Hmisc)
library(dplyr)
library(MASS)

Survival<-read.csv("Regression_survival.csv")
Macro<-read.csv("Regression_macro.csv")
Micro<-read.csv("Regression_micro.csv")

number.of.failures <- Survival$dead
number.of.successes <- Survival$alive
binomial.denominator <- Survival$total

number.of.failures <- binomial.denominator - number.of.successes
y <- cbind(number.of.successes, number.of.failures)

Survival_MRA<- cbind(y,Macro,Micro)
Survival_MRA<- Survival_MRA[-c(19:20)]

Survival_MRA$Tile <- as.factor(Survival_MRA$Tile)
Survival_MRA$Timepoint <- as.factor(Survival_MRA$Timepoint)

Survival_MRA$number.of.successes<-as.character(Survival_MRA$number.of.successes)
Survival_MRA$number.of.failures<-as.character(Survival_MRA$number.of.failures)

str(Survival_MRA)

Survival_MRA1<-Survival_MRA %>% mutate_if(is.numeric, scale)
Survival_MRA1$number.of.successes<-as.numeric(Survival_MRA1$number.of.successes)
Survival_MRA1$number.of.failures<-as.numeric(Survival_MRA1$number.of.failures)
str(Survival_MRA1)

#### First Mult. Reg. section - All tiles ####

#Macro
# Run multicollinearity test
colnames(Survival_MRA1)
Macros<-dplyr::select(Survival_MRA, Bare, CCA, Cyanobacteria, EAM, Inverts, Lobophora, Long_turf, Macroalgae, NCCA,Sand, Thick_turf, Turf)

res2 <- rcorr(as.matrix(Macros))
-corrplot(res2$r, type="upper", order="hclust", p.mat = res2$P, sig.level = 0.01, insig = "blank") 
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(res2$r, method = "color",
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res2$P, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)

#Drop Macroalgae and cyanobacteria
Macros<-dplyr::select(Survival_MRA, Bare, CCA, EAM, Inverts, Lobophora, Long_turf, Turf, NCCA, Sand, Thick_turf)

Macro.mod<- glmer(y ~ Bare+ CCA+ EAM+ Inverts+ Lobophora+ Long_turf+ Turf+ NCCA+ Sand+ Thick_turf + (1|Timepoint/Tile), data=Survival_MRA, family="binomial")
drop1(Macro.mod, test="Chi") # drop Thick_turf

Macro.mod.1<-update(Macro.mod, . ~ . -(Thick_turf))
drop1(Macro.mod.1, test="Chi") #drop Sand

Macro.mod.2<-update(Macro.mod.1, . ~ . -(Sand))
drop1(Macro.mod.2, test="Chi") #drop Long_turf

Macro.mod.3<-update(Macro.mod.2, . ~ . -(Long_turf))
drop1(Macro.mod.3, test="Chi") #drop Turf

Macro.mod.4<-update(Macro.mod.3, . ~ . -(Turf))
drop1(Macro.mod.4, test="Chi") #drop Inverts

Macro.mod.5<-update(Macro.mod.4, . ~ . -(Inverts))
drop1(Macro.mod.5, test="Chi") #drop EAM

Macro.mod.6<-update(Macro.mod.5, . ~ . -(EAM))
drop1(Macro.mod.6, test="Chi") #drop NCCA

Macro.mod.7<-update(Macro.mod.6, . ~ . -(NCCA))
drop1(Macro.mod.7, test="Chi") #Final - Bare, CCA, Lobophora

#microbial
Micros<-dplyr::select(Survival_MRA, Rhodobacterales, Chitinophagales, Flavobacteriales, Cytophagales, Alteromonadales, Nostocales, Cellvibrionales, Vibrionales, Caulobacterales,Phormidesmiales)

res2 <- rcorr(as.matrix(Micros))
-corrplot(res2$r, type="upper", order="hclust", p.mat = res2$P, sig.level = 0.01, insig = "blank") 
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(res2$r, method = "color",
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res2$P, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)

#Drop Cellvibrionales

Micros<-dplyr::select(Survival_MRA1, Rhodobacterales, Chitinophagales, Flavobacteriales, Cytophagales, Alteromonadales, Nostocales, Phormidesmiales, Vibrionales, Caulobacterales)

Micro.mod<- glmer(y ~ Rhodobacterales+ Chitinophagales+ Flavobacteriales+ Cytophagales+ Alteromonadales+ Nostocales+ Phormidesmiales+ Vibrionales+ Caulobacterales + (1|Timepoint/Tile), data=Survival_MRA, family="binomial")
drop1(Micro.mod, test="Chi") # drop Vibrionales

Micro.mod.1<-update(Micro.mod, . ~ . -(Vibrionales))
drop1(Micro.mod.1, test="Chi") #drop Caulobacterales

Micro.mod.2<-update(Micro.mod.1, . ~ . -(Caulobacterales))
drop1(Micro.mod.2, test="Chi") #drop Rhodobacterales

Micro.mod.3<-update(Micro.mod.2, . ~ . -(Rhodobacterales))
drop1(Micro.mod.3, test="Chi") #drop Phormidesmiales

Micro.mod.4<-update(Micro.mod.3, . ~ . -(Phormidesmiales))
drop1(Micro.mod.4, test="Chi") #drop Flavobacteriales

Micro.mod.5<-update(Micro.mod.4, . ~ . -(Flavobacteriales))
drop1(Micro.mod.5, test="Chi") #drop Nostocales

Micro.mod.6<-update(Micro.mod.5, . ~ . -(Nostocales))
drop1(Micro.mod.6, test="Chi") #drop Chitinophagales

Micro.mod.7<-update(Micro.mod.6, . ~ . -(Chitinophagales))
drop1(Micro.mod.7, test="Chi") #Final - Alteromonadales, Cytophagales

#combined mod
Combined.mod<-glmer(y ~ Alteromonadales + Cytophagales + Bare + CCA + Lobophora + (1|Timepoint/Tile), data=Survival_MRA1, family="binomial")

drop1(Combined.mod, test="Chi") # drop Cytophagales

Combined.mod.1<-update(Combined.mod, . ~ . -(Cytophagales))
drop1(Combined.mod.1, test="Chi") #drop Alteromonadales

Combined.mod.2<-update(Combined.mod.1, . ~ . -(Alteromonadales))
drop1(Combined.mod.2, test="Chi") #all macro groups significant

#### Model output for "All tiles" section of table S12 ####
sjPlot::tab_model(Combined.mod.2); summary(Combined.mod.2)

#### Second Mult. Reg. section - Open tiles ####

#data file containing survival rates, standardised % cover of macro- and micro-communities, and grazing rates on open tiles
Survival_MRA_open<-read.csv("Regression_open_tiles.csv")
colnames(Survival_MRA_open)

y <- cbind(Survival_MRA_open$number.of.successes, Survival_MRA_open$number.of.failures)

#grazers
Grazing<-dplyr::select(Survival_MRA_open,Carnivore,Detritivore_Grazer,Excavator,Grazer,Omnivore,Scraper,Total,Herbivores)

# Run collinearity test
res2 <- rcorr(as.matrix(Grazing))
-corrplot(res2$r, type="upper", order="hclust", p.mat = res2$P, sig.level = 0.01, insig = "blank") 
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(res2$r, method = "color",
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res2$P, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)

#remove Total and Herbivores

Grazing<-dplyr::select(Survival_MRA_open,Carnivore,Detritivore_Grazer,Excavator,Grazer,Omnivore,Scraper)

Grazing.mod<- glmer(y ~ Carnivore+Detritivore_Grazer+Excavator+Grazer+Omnivore+Scraper + (1|Timepoint/Tile), data=Survival_MRA_open, family="binomial")
drop1(Grazing.mod, test="Chi") # drop Grazer

Grazing.mod.1<-update(Grazing.mod, . ~ . -(Grazer))
drop1(Grazing.mod.1, test="Chi") #drop Omnivore

Grazing.mod.2<-update(Grazing.mod.1, . ~ . -(Omnivore))
drop1(Grazing.mod.2, test="Chi") #drop Carnivore

Grazing.mod.3<-update(Grazing.mod.2, . ~ . -(Carnivore))
drop1(Grazing.mod.3, test="Chi") #drop Scraper

Grazing.mod.4<-update(Grazing.mod.3, . ~ . -(Scraper))
drop1(Grazing.mod.4, test="Chi") # Detritivore_Grazer & Excavator

#Open macros

#Note: Macroalgae, Cyanobacteria, Long and thick turf removed as these were all zeros on open tiles.
Macros_O<-dplyr::select(Survival_MRA_open, Bare, Inverts, CCA, EAM, Inverts, Lobophora, Turf, NCCA, Sand)

res2 <- rcorr(as.matrix(Macros_O))
-corrplot(res2$r, type="upper", order="hclust", p.mat = res2$P, sig.level = 0.01, insig = "blank") 
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(res2$r, method = "color",
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res2$P, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)

#remove EAM and Sand

Open.macro.mod<- glmer(y ~ CCA + Lobophora + Turf + NCCA + Bare + Inverts + (1|Timepoint/Tile), data=Survival_MRA_open, family="binomial")
drop1(Open.macro.mod, test="Chi") # drop Inverts

Open.macro.mod.1<-update(Open.macro.mod, . ~ . -(Inverts))
drop1(Open.macro.mod.1, test="Chi") #drop Lobophora

Open.macro.mod.2<-update(Open.macro.mod.1, . ~ . -(Lobophora))
drop1(Open.macro.mod.2, test="Chi") #drop NCCA

Open.macro.mod.3<-update(Open.macro.mod.2, . ~ . -(NCCA))
drop1(Open.macro.mod.3, test="Chi") #Keep CCA, Turf, and Bare for main model with Micro and grazing predictors

#Open micros

Micros_O<-dplyr::select(Survival_MRA_open, Rhodobacterales, Chitinophagales, Flavobacteriales, Cytophagales, Alteromonadales, Nostocales, Cellvibrionales, Vibrionales, Caulobacterales)

res2 <- rcorr(as.matrix(Micros_O))
-corrplot(res2$r, type="upper", order="hclust", p.mat = res2$P, sig.level = 0.01, insig = "blank") 
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(res2$r, method = "color",
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res2$P, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)

#remove Nostocales

Micros_O<-dplyr::select(Survival_MRA_open, Rhodobacterales, Chitinophagales, Flavobacteriales, Cytophagales, Alteromonadales, Cellvibrionales, Vibrionales, Caulobacterales)

Open.Micro.mod<- glmer(y ~ Rhodobacterales+ Chitinophagales+ Flavobacteriales+ Cytophagales+ Alteromonadales+ Cellvibrionales+ Vibrionales+ Caulobacterales + (1|Timepoint/Tile), data=Survival_MRA_open, family="binomial")
drop1(Open.Micro.mod, test="Chi") # drop Flavobacteriales

Open.Micro.mod.1<-update(Open.Micro.mod, . ~ . -(Flavobacteriales))
drop1(Open.Micro.mod.1, test="Chi") #drop Chitinophagales

Open.Micro.mod.2<-update(Open.Micro.mod.1, . ~ . -(Chitinophagales))
drop1(Open.Micro.mod.2, test="Chi") #all sig. - Rhodobacterales, Cytophagales, Alteromonadales, Cellvibrionales, Vibrionales, Caulobacterales


sjPlot::tab_model(Open.Micro.mod.2); summary(Open.Micro.mod.2)
sjPlot::plot_model(Open.Micro.mod.2, show.values = TRUE, value.offset = .3)

# FULL MOD with macro, micro and grazing predictors
Full_O<-dplyr::select(Survival_MRA_open, CCA, Turf, Bare, Rhodobacterales, Cytophagales, Alteromonadales, Cellvibrionales, Vibrionales, Caulobacterales, Detritivore_Grazer, Excavator)

Full.mod<-glmer(y ~ Lobophora + Turf + Sand + Rhodobacterales+ Cytophagales+ Alteromonadales+ Cellvibrionales+ Vibrionales+ Caulobacterales + Detritivore_Grazer + Excavator + (1|Timepoint/Tile), data=Survival_MRA_open, family="binomial")

drop1(Full.mod, test="Chi") # drop Turf

Full.mod.1<-update(Full.mod, . ~ . -(Turf))
drop1(Full.mod.1, test="Chi") #drop Sand

Full.mod.2<-update(Full.mod.1, . ~ . -(Sand))
drop1(Full.mod.2, test="Chi") #drop Lobophora

Full.mod.3<-update(Full.mod.2, . ~ . -(Lobophora))
drop1(Full.mod.3, test="Chi") #drop Excavator

Full.mod.4<-update(Full.mod.3, . ~ . -(Excavator))
drop1(Full.mod.4, test="Chi") #drop Caulobacterales 

Full.mod.5<-update(Full.mod.4, . ~ . -(Caulobacterales))
drop1(Full.mod.5, test="Chi") #drop Cellvibrionales

Full.mod.6<-update(Full.mod.5, . ~ . -(Cellvibrionales))
drop1(Full.mod.6, test="Chi") #drop Cytophagales

Full.mod.7<-update(Full.mod.6, . ~ . -(Cytophagales))
drop1(Full.mod.7, test="Chi") #drop Alteromonadales

Full.mod.8<-update(Full.mod.7, . ~ . -(Alteromonadales))
drop1(Full.mod.8, test="Chi") #final - Rhodobacterales, Vibrionales, Detritivore_Grazer

#### Model output for "Open tiles" section of table S12 ####
sjPlot::tab_model(Full.mod.8); summary(Full.mod.8)

#END