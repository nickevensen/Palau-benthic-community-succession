#Code for Figure 4 of "Benthic micro- and macro-community succession and coral recruitment under overfishing and nutrient enrichment"
#Also includes statistical output from Tables S9 and S10

#setwd to source file location (data file and R script are in the same folder)

library(vegan)
library(dplyr)
Fish_bites<-read.csv("Fish_bites.csv")


dataIR<-transform(Fish_bites, 
                  Nutrients=as.factor(Nutrients), 
                  Caging=as.factor(Caging),
                  Microhabitat=as.factor(Microhabitat))
str(dataIR)
names(dataIR)
attach(dataIR)

#create separate dataframes
## sometimes the select function is acting weird, if errors message, close R and re opens and start again from here after reading dataset, should work
dataIR_com<-dplyr::select(dataIR,Detritivore_Grazer,Excavator,Grazer,Omnivore,Scraper)
dataIR_meta<-dplyr::select(dataIR,Nutrients,Caging,Microhabitat)

#make distance matrix
DISTdataIR_com<- vegdist(dataIR_com, method="bray")

# ADONIS function
IR<-adonis(DISTdataIR_com ~ Nutrients*Caging*Microhabitat, data=dataIR_meta)
IR

library(rtf)
tables9<-as.data.frame(IR$aov.tab, keep.rownames = TRUE)
tables9$Source <- rownames(tables9)
tables9<-format(tables9, digits = 2)

rtffile9 <- RTF("rtf9.doc")
addParagraph(rtffile9, "This is the output for Table s9:\n")
addTable(rtffile9, tables9)
done(rtffile9)

#pairwise comparisons
library(pairwiseAdonis)
str(dataIR)

Nutrients<-dataIR$Nutrients
pairwise.adonis2(dataIR_com ~ Microhabitat, data = dataIR_meta,sim.function = "vegdist",sim.method = "bray",p.adjust.m = "bonferroni", strata = 'Nutrients') 

#Create new column to pair treatment combinations to Nutrients all your treatments in the same figure - remember that R is case-sensitive, below 'treatment' is lower case
new.data.scores[,"treatment"]<-as.character(NA)

#Crowns
new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Nutrient") & (new.data.scores$Microhabitat=="Crown")]<-"Nutrient-Crown"

new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Ambient") & (new.data.scores$Microhabitat=="Crown")]<-"Ambient-Crown"

#Crevices
new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Nutrient") & (new.data.scores$Microhabitat=="Crevice")]<-"Nutrient-Crevice"

new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Ambient") & (new.data.scores$Microhabitat=="Crevice")]<-"Ambient-Crevice"

#Sides
new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Nutrient") & (new.data.scores$Microhabitat=="Side")]<-"Nutrient-Side"

new.data.scores$treatment[Name=(new.data.scores$Nutrients == "Ambient") & (new.data.scores$Microhabitat=="Side")]<-"Ambient-Side"


#Nutrients for all, pooling partial and open as they didn't differ
cent <- aggregate(cbind(NMDS1, NMDS2) ~ treatment, data = new.data.scores, FUN = mean)
segs <- merge(new.data.scores, setNames(cent, c('treatment','oNMDS1','oNMDS2')),
              by = 'treatment', sort = FALSE)

new.data.scores$treatment<-as.factor(new.data.scores$treatment)
print(levels(new.data.scores$treatment))

#"Ambient-Crevice"   "Ambient-Side"     "Ambient-Crown"      "Nutrient-Crevice" "Nutrient-Side"   "Nutrient-Crown"
# 'navy','steelblue1','dodgerblue4','darkgreen','palegreen3','palegreen4'
cols <- c("Ambient-Crevice" = "navy","Ambient-Side" = "dodgerblue4", "Ambient-Crown" = "steelblue1",
          "Nutrient-Crevice" = "darkgreen", "Nutrient-Side" = "palegreen4", "Nutrient-Crown" = "palegreen3")

library(ggplot2)
ggplot(new.data.scores, aes(x = NMDS1, y = NMDS2, col=treatment)) + scale_colour_manual(values=cols) +
  stat_ellipse(geom = "polygon", alpha = 0.25, fill="grey", size= 1, level = 0.50) + # ellipses encompassing 50% of replicates
  #geom_segment(data = segs, mapping = aes(xend = oNMDS1, yend = oNMDS2), colour="black") +  # spiders
  geom_point(data = cent, shape=(values=c(24, 21, 22, 24, 21, 22)), size = 5, stroke = 2, colour=(values=c("navy","steelblue1","dodgerblue4","darkgreen","palegreen3","palegreen4")), fill="black") +       # centroids
  geom_point(size = 2.5, stroke=1, aes(shape=treatment)) +  scale_shape_manual(values=c(17, 15, 16, 17, 15, 16)) +                                     # sample scores
  coord_fixed() +                                              # same axis scaling
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + geom_vline(xintercept=0, linetype="dashed") + geom_hline(yintercept=0, linetype="dashed") +
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) + scale_y_continuous(limits=c(-0.65,0.2), breaks = c(-0.6,-0.4,-0.2,0,0.2)) + scale_x_continuous(limits=c(-0.51,0.4), breaks = c(-0.4,-0.2,0,0.2,0.4))

# That's the main plot done. Now to plot the vectors to show what's driving difference in ordinations (i.e., coordinates of the point)
library(MASS)
library(vegan)
ord_IR<-metaMDS(DISTdataIR_com,distance = "bray",
                k=2,trymax=1000,autotransform=TRUE,expand=FALSE, plot=FALSE)
ord_IR
fit<-envfit(ord_IR, dataIR_com, perm=999)
fit
scores(fit, "vectors")
plot(ord_IR,col="white",axes=FALSE, ylim=c(-0.65,0.2), xlim=c(-0.51,0.4))
plot(fit, p.max = 0.95, col = "black") + abline(h = 0, v = 0, col = "black", lty = 2)
axis(side=1, at=seq(-0.5, 0.6, by=0.1))
axis(side=2, at=seq(-0.4, 0.5, by=0.1))
box()

tables10 <- data.frame((fit$vectors)$arrows, (fit$vectors)$r, (fit$vectors)$pvals)
tables10$Source <- rownames(tables10)
tables10<-format(tables10, digits = 3)

rtffile10 <- RTF("rtf10.doc")
addParagraph(rtffile10, "This is the output for Table s10:\n")
addTable(rtffile10, tables10)
done(rtffile10)

#END