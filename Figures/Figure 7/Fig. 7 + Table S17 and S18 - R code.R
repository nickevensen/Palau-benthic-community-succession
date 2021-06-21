#Code for Figure 7 of "Benthic micro- and macro-community succession and coral recruitment under overfishing and nutrient enrichment"
#Also includes statistical output from Tables S17 and S18

#setwd to source file location (data file and R script are in the same folder)

#Data and analysis used to plot Fig. 7a - the final figure was plotted using GraphPad Prism Software, 
#thus no R code to plot the figure is provided.
Settlement_per_cm2<-read.csv("Settlement_rates_cm2.csv")

library(car)

Settlement_mod<-glm(Settlers ~ Nutrients*Caging*Microhabitat, data=Settlement_per_cm2)
summary(Settlement_mod)
Anova(Settlement_mod)

tableS17<-Anova(Settlement_mod)
tableS17$Source <- rownames(tableS17)
tableS17<-format(tableS17, digits = 3)

rtffile17 <- RTF("rtf17.doc")
addParagraph(rtffile17, "This is the output for Table S17:\n")
addTable(rtffile17, tableS17)
done(rtffile17)

#data from post hoc comparisons below also included in Table s17
library(emmeans)
emmeans(Settlement_mod, pairwise ~ Caging | Microhabitat)
emmeans(Settlement_mod, pairwise ~ Microhabitat | Caging)

#### Code for Multiple Regression analysis (Table S18) and Figure 7c analysis and plot
library(corrplot)
library(mctest)
library(nlme)
library(lme4)
library(Hmisc)
library(dplyr)
library(MASS)

#### MICROBIAL ONLY ####
Microbial_settlement<-read.csv("Microbes_vs_settlement.csv")

MicroAll<-dplyr::select(Microbial_settlement,Rhodobacterales,Flavobacteriales,Cytophagales,Cellvibrionales,
                 Chitinophagales,Nostocales,Alteromonadales,Caulobacterales, Gammaproteobacteria_unclassified, Vibrionales)

res2 <- rcorr(as.matrix(MicroAll))

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(res2$r, method = "color",
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res2$P, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)

-corrplot(res2$r, type="upper", p.mat = res2$P, sig.level = 0.01, insig = "blank")

#remove Gammaproteobacteria_unclassified & Vibrionales due to high correlation with other microbial groups

Microbial_select<-dplyr::select(Microbial_settlement,Settlement,Rhodobacterales,Flavobacteriales,Cytophagales,Cellvibrionales,
                                Chitinophagales,Nostocales,Alteromonadales,Caulobacterales)

str(Microbial_select)

#Model selection
Micro.mod<-lm(Settlement~Rhodobacterales+Flavobacteriales+Cytophagales+Cellvibrionales+
                Chitinophagales+Nostocales+Alteromonadales+Caulobacterales,data=Microbial_select)
drop1(Micro.mod, test="Chi") #drop Alteromonadales

Micro.mod.1<-update(Micro.mod, . ~ . -(Alteromonadales))
drop1(Micro.mod.1, test="Chi") #drop Nostocales

Micro.mod.2<-update(Micro.mod.1, . ~ . -(Nostocales))
drop1(Micro.mod.2, test="Chi") #drop Caulobacterales

Micro.mod.3<-update(Micro.mod.2, . ~ . -(Caulobacterales))
drop1(Micro.mod.3, test="Chi") #drop Cellvibrionales

Micro.mod.4<-update(Micro.mod.3, . ~ . -(Cellvibrionales))
drop1(Micro.mod.4, test="Chi") #drop Chitinophagales

Micro.mod.5<-update(Micro.mod.4, . ~ . -(Chitinophagales))
drop1(Micro.mod.5, test="Chi") #drop Cytophagales

Micro.mod.6<-update(Micro.mod.5, . ~ . -(Cytophagales))
drop1(Micro.mod.6, test="Chi") #final Rhodobacterales & Flavobacteriales

#For Table S18
sjPlot::tab_model(Micro.mod.6); summary(Micro.mod.6)

#individual microbial groups
Rhodobacterales.mod<-lm(Settlement~Rhodobacterales,data=Microbial_select)
summary(Rhodobacterales.mod) #N.S.

Flavobacteriales.mod<-lm(Settlement~Flavobacteriales,data=Microbial_select)
summary(Flavobacteriales.mod) #N.S.

#### MACRO ONLY ####
Benthos_vs_settlement<-read.csv("Benthos_vs_settlement.csv")

MacroAll<-dplyr::select(Benthos_vs_settlement,CCA,EAM,Lobophora,NCCA,Sand,Thick_turf,Turf,Macroalgae,Invert)

res2 <- rcorr(as.matrix(MacroAll))

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(res2$r, method = "color",
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res2$P, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)

-corrplot(res2$r, type="upper", p.mat = res2$P, sig.level = 0.01, insig = "blank") 
#drop turf

###
Macro_select<-dplyr::select(Benthos_vs_settlement,Settlement,CCA,EAM,Lobophora,NCCA,Sand,Thick_turf,Macroalgae,Invert)


Macro.mod<-lm(Settlement~CCA+EAM+Lobophora+NCCA+Sand+Thick_turf+Macroalgae+Invert,data=Macro_select)
drop1(Macro.mod, test="Chi") #drop Thick_turf

Macro.mod.1<-update(Macro.mod, . ~ . -(Thick_turf))
drop1(Macro.mod.1, test="Chi") #drop Invert

Macro.mod.2<-update(Macro.mod.1, . ~ . -(Invert))
drop1(Macro.mod.2, test="Chi") #drop Macroalgae

Macro.mod.3<-update(Macro.mod.2, . ~ . -(Macroalgae))
drop1(Macro.mod.3, test="Chi") #drop Lobophora

Macro.mod.4<-update(Macro.mod.3, . ~ . -(Lobophora))
drop1(Macro.mod.4, test="Chi") #drop NCCA

Macro.mod.5<-update(Macro.mod.4, . ~ . -(NCCA))
drop1(Macro.mod.5, test="Chi") #all sig

#final - CCA, EAM, Sand

#For Table S18
sjPlot::tab_model(Macro.mod.5); summary(Macro.mod.5)

#individual
CCA.mod<-lm(Settlement~CCA,data=Macro_select)
summary(CCA.mod) #Sig. - p < 0.01 / Adj. R2 = 0.42

EAM.mod<-lm(Settlement~EAM,data=Macro_select)
summary(EAM.mod) #N.S.

Sand.mod<-lm(Settlement~Sand,data=Macro_select)
summary(Sand.mod) #N.S.

#### Plot Fig. 7c - CCA vs. Settlement ####
ggplot(Macro_select, aes(x = CCA, y = Settlement)) + 
  stat_smooth(method = "glm.nb", col = "grey40", fill="coral", alpha=0.2) + 
  geom_point(shape=21, col="black", fill="black", alpha=0.8, stroke = 0.8) +
  theme_classic() + theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2, "cm")) + 
  scale_x_continuous(name="CCA cover (%)", breaks = round(seq(min(Macro_select$CCA), max(Macro_select$CCA), by = 10),60)) + 
  scale_y_continuous(name="No. of settlers") + scale_y_continuous(limits=c(0,50))

#END