#Code for Figure S2 of "Benthic micro- and macro-community succession and coral recruitment under overfishing and nutrient enrichment"

#setwd to source file location (data file and R script are in the same folder)

#Proportion of coral survivors of vertical surfaces at final timepoint - for Fig. S2a (plotted using GraphPad Prism)
Final_survival<-read.csv("Final_survival_vertical.csv")

#organise data
number.of.failures <- Final_survival$dead
number.of.successes <- Final_survival$alive
binomial.denominator <- Final_survival$total

number.of.failures <- binomial.denominator - number.of.successes
y <- cbind(number.of.successes, number.of.failures)

prop_alive <- Final_survival$proportion
Caging <- as.factor(Final_survival$Caging)
Nutrients <- as.factor(Final_survival$Nutrients)
Microhabitat <- as.factor(Final_survival$Microhabitat)
Tile <- as.factor(Final_survival$Tile)

#Analysis
library(MASS)
library(car)

Final.surv.mod <- glm(y~Caging*Nutrients, family=quasibinomial)
summary(Final.surv.mod)
Anova(Final.surv.mod) # Caging is significant

library(emmeans)
emmeans(Final.surv.mod, pairwise ~ Caging)

library(plyr)
summary_data <- ddply(Final_survival, c("Caging"), summarise, N= length(proportion), mean = mean(proportion),
                      sd   = sd(proportion),
                      se   = sd / sqrt(N))
summary_data


##### PERMANOVA to compare final cover on vertical surfaces #####
Vertical_cover<-read.csv("Vertical_cover_for_Fig.S2a.csv")

dataIR<-transform(Vertical_cover, 
                  Nutrients=as.factor(Nutrients), 
                  Caging=as.factor(Caging))
str(dataIR)
names(dataIR)
attach(dataIR)


library(dplyr)
## Create a community data frame and a predictor variable data frame (sometimes the select function is acting weird, if errors message, close R and re opens and start again from here after reading dataset)
dataIR_com<-dplyr::select(dataIR,CCA,EAM,Inverts,Lobophora,Macroalgae,Turf,Thick_turf,
                   NCCA)
dataIR_meta<-dplyr::select(dataIR,Nutrients,Caging)

#Make distance matrix
library(vegan)
DISTdataIR_com<- vegdist(dataIR_com, method="bray")

#Run adonis (i.e., PERMANOVA)
IR<-adonis(dataIR_com~Nutrients*Caging, data=dataIR_meta)
IR #only caging treatment is significant

library(pairwiseAdonis)
str(dataIR)

pairwise.adonis(dataIR_com,factors=dataIR$Caging,sim.function='vegdist',sim.method='bray',p.adjust.m='bonferroni')

#This produce the coordinates for each point
Sol <- metaMDS(dataIR_com)
Sol
stressplot(Sol)

data.scores <- as.data.frame(scores(Sol))
new.data.scores <- cbind(data.scores,dataIR_meta)

new.data.scores[,"treatment"]<-as.character(NA)
new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Ambient") & (new.data.scores$Caging=="Open")]<-"Ambient-Open"
new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Ambient") & (new.data.scores$Caging=="Caged")]<-"Ambient-Caged"
new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Ambient") & (new.data.scores$Caging=="Partial")]<-"Ambient-Partial"
new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Nutrient") & (new.data.scores$Caging=="Open")]<-"Nutrient-Open"
new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Nutrient") & (new.data.scores$Caging=="Caged")]<-"Nutrient-Caged"
new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Nutrient") & (new.data.scores$Caging=="Partial")]<-"Nutrient-Partial"

cent <- aggregate(cbind(NMDS1, NMDS2) ~ treatment, data = new.data.scores, FUN = mean)
segs <- merge(new.data.scores, setNames(cent, c('treatment','oNMDS1','oNMDS2')),
              by = 'treatment', sort = FALSE)

# plot with all treatments with lines joining the timepoints within a treatment and labels showing what timepoint each point is (colours represent treatments)
cols <- c("Ambient-Caged" = "navy", "Ambient-Partial" = "dodgerblue4", "Ambient-Open" = "steelblue1",
          "Nutrient-Caged" = "darkgreen", "Nutrient-Partial" = "palegreen4", "Nutrient-Open" = "palegreen3")
# 'navy','steelblue1','dodgerblue4','darkgreen','palegreen3','palegreen4'

#Plot nMDS of final community structure on vertical surfaces - Fig. S2b-c
library(ggplot2)
ggplot(new.data.scores, aes(x = NMDS1, y = NMDS2, col=treatment)) + scale_colour_manual(values=cols) + 
  stat_ellipse(geom = "polygon", alpha = 0.25, fill="grey", size= 1, level = 0.66) + # ellipses
  geom_point(data = cent, shape=(values=c(22, 24, 21, 22, 24, 21)), size = 3, stroke = 1, colour=(values=c("navy","steelblue1","dodgerblue4","darkgreen","palegreen3","palegreen4")), fill="black") +       # centroids
  geom_point(aes(shape=treatment)) +  scale_shape_manual(values=c(15, 17, 16, 15, 17, 16)) +                                               # sample scores
  coord_fixed() +                                              # same axis scaling
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + geom_vline(xintercept=0, linetype="dashed") + geom_hline(yintercept=0, linetype="dashed") +
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) + 
  scale_y_continuous(limits=c(-1.1,1.2), breaks = c(-1.0,-0.5,0,0.5,1)) + scale_x_continuous(limits=c(-1.3,1.5), breaks = c(-1,0,1.0))

#plot vectors driving the MDS
library(MASS)
library(vegan)
ord_IR<-metaMDS(dataIR_com,distance = "bray",
                k=2,trymax=1000,autotransform=TRUE,expand=FALSE, plot=FALSE)
ord_IR
fit<-envfit(ord_IR, dataIR_com, perm=999)
fit
scores(fit, "vectors")
plot(ord_IR,col="white", ylim=c(-1.1, 1.2), xlim=c(-1.3, 1.5))
plot(fit, p.max = 0.9, col = "black") + abline(h = 0, v = 0, col = "black", lty = 2)

#END
