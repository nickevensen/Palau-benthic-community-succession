#Code for analysis represented in Figure S8 of "Benthic micro- and macro-community succession and coral recruitment under overfishing and nutrient enrichment"
#Also includes statistical output from Table S7 and S8
#Note Figure S8 was created using GraphPad Prism software, and therefore no code to plot the data is included in this script


#setwd to source file location (data file and R script are in the same folder)
Fish_bite_means<-read.csv("Fish_bite_means.csv")


Fish_bite_means$N<-NA
cdata<-aggregate(Fish_bite_means["N"], by=Fish_bite_means[c("Nutrients","Caging","Microhabitat")], FUN=length)
cdata.means<-aggregate(Fish_bite_means[c("Total")], by=Fish_bite_means[c("Nutrients","Caging","Microhabitat")], FUN=mean)
cdata<-merge(cdata, cdata.means)
cdata
cdata.sd<-aggregate(Fish_bite_means[c("Total")], by=Fish_bite_means[c("Nutrients","Caging","Microhabitat")], FUN=sd)
cdata.sd<-aggregate(Fish_bite_means[c("Total")], by=Fish_bite_means[c("Nutrients","Caging","Microhabitat")], FUN=sd)
cdata.se<-cdata.sd$Total/15

cdata<-cbind(cdata,cdata.se)

cdata[,"treatment"]<-as.character(NA)

cdata$treatment[Name=(cdata$Nutrients == "Nutrient") & (cdata$Caging=="Open")]<-"Nutrient-Open"
cdata$treatment[Name=(cdata$Nutrients == "Nutrient") & (cdata$Caging=="Partial")]<-"Nutrient-Partial"

cdata$treatment[Name=(cdata$Nutrients == "Ambient") & (cdata$Caging=="Open")]<-"Ambient-Open"
cdata$treatment[Name=(cdata$Nutrients == "Ambient") & (cdata$Caging=="Partial")]<-"Ambient-Partial"

cdata$Microhabitat<-as.factor(cdata$Microhabitat)
levels(cdata$Microhabitat)
cdata$Microhabitat = factor(cdata$Microhabitat,levels(cdata$Microhabitat)[c(2,1,3)]) 


library(ggplot2)
ggplot(cdata, aes(fill=treatment, y=Total, x=Microhabitat)) + 
  geom_bar(position="dodge", stat="identity") + scale_fill_manual(values=c("steelblue1","navy","palegreen3","darkgreen"))  + ylim(0,0.8) + theme_bw() +
  geom_errorbar(aes(ymin=Total-cdata.se, ymax=Total+cdata.se), width=.2,
                position=position_dodge(.9)) 

###Fish bite data is in bites per min per tile for each functional group

library(lme4)
library(car)
library(emmeans)

#All fish
mod1<-aov(Total ~ Caging*Nutrients*Microhabitat + Error(Week), data=Fish_bite_means)
summary(mod1)

emmeans(mod1, pairwise ~ Microhabitat)

tableS7<-summary(mod1)[["Error: Within"]][[1]]
tableS7$Source <- rownames(tableS7)
tableS7<-format(tableS7, digits = 2)

rtffile7 <- RTF("rtf7.doc")
addParagraph(rtffile7, "This is the output for Table S7:\n")
addTable(rtffile7, tableS7)
done(rtffile7)

#All herbivores
mod1<-aov(Herbivores ~ Caging*Nutrients*Microhabitat + Error(Week), data=Fish_bite_means)
summary(mod1)

emmeans(mod1, pairwise ~ Caging | Microhabitat)
emmeans(mod1, pairwise ~ Microhabitat | Caging)

#Detritivore_grazer category
mod1<-aov(Detritivore_Grazer ~ Caging*Nutrients*Microhabitat + Error(Week), data=Fish_bite_means)
summary(mod1)

emmeans(mod1, pairwise ~ Microhabitat)
emmeans(mod1, pairwise ~ Nutrients)

tableS7<-summary(mod1)[["Error: Within"]][[1]]
tableS7$Source <- rownames(tableS7)
tableS7<-format(tableS7, digits = 2)

rtffile7 <- RTF("rtf7.doc")
addParagraph(rtffile7, "This is the output for Table S7:\n")
addTable(rtffile7, tableS7)
done(rtffile7)

#Grazer category
mod1<-aov(Grazer ~ Caging*Nutrients*Microhabitat + Error(Week), data=Fish_bite_means)
summary(mod1)

emmeans(mod1, list(pairwise ~ Nutrients))
emmeans(mod1, pairwise ~ Caging | Microhabitat)
emmeans(mod1, pairwise ~ Microhabitat | Caging)

tableS7<-summary(mod1)[["Error: Within"]][[1]]
tableS7$Source <- rownames(tableS7)
tableS7<-format(tableS7, digits = 2)

rtffile7 <- RTF("rtf7.doc")
addParagraph(rtffile7, "This is the output for Table S7:\n")
addTable(rtffile7, tableS7)
done(rtffile7)

#Omnivore category
mod1<-aov(Omnivore ~ Caging*Nutrients*Microhabitat + Error(Week), data=Fish_bite_means)
summary(mod1)

emmeans(mod1, list(pairwise ~ Nutrients))

tableS7<-summary(mod1)[["Error: Within"]][[1]]
tableS7$Source <- rownames(tableS7)
tableS7<-format(tableS7, digits = 2)

rtffile7 <- RTF("rtf7.doc")
addParagraph(rtffile7, "This is the output for Table S7:\n")
addTable(rtffile7, tableS7)
done(rtffile7)

#Scraper category
mod1<-aov(Scraper ~ Caging*Nutrients*Microhabitat + Error(Week), data=Fish_bite_means)
summary(mod1)

emmeans(mod1, list(pairwise ~ Nutrients))

tableS7<-summary(mod1)[["Error: Within"]][[1]]
tableS7$Source <- rownames(tableS7)
tableS7<-format(tableS7, digits = 2)

rtffile7 <- RTF("rtf7.doc")
addParagraph(rtffile7, "This is the output for Table S7:\n")
addTable(rtffile7, tableS7)
done(rtffile7)

#Excavator category
mod1<-aov(Excavator ~ Caging*Nutrients*Microhabitat + Error(Week), data=Fish_bite_means)
summary(mod1)

emmeans(mod1, list(pairwise ~ Microhabitat))
emmeans(mod1, list(pairwise ~ Nutrients))

tableS7<-summary(mod1)[["Error: Within"]][[1]]
tableS7$Source <- rownames(tableS7)
tableS7<-format(tableS7, digits = 2)

rtffile7 <- RTF("rtf7.doc")
addParagraph(rtffile7, "This is the output for Table S7:\n")
addTable(rtffile7, tableS7)
done(rtffile7)

#Corallivore category
mod1<-aov(Corallivore ~ Caging*Nutrients*Microhabitat + Error(Week), data=Fish_bite_means)
summary(mod1)

emmeans(mod1, list(pairwise ~ Nutrients))

tableS7<-summary(mod1)[["Error: Within"]][[1]]
tableS7$Source <- rownames(tableS7)
tableS7<-format(tableS7, digits = 2)

rtffile7 <- RTF("rtf7.doc")
addParagraph(rtffile7, "This is the output for Table S7:\n")
addTable(rtffile7, tableS7)
done(rtffile7)

#Carnivore category
mod1<-aov(Carnivore ~ Caging*Nutrients*Microhabitat + Error(Week), data=Fish_bite_means)
summary(mod1)

tableS7<-summary(mod1)[["Error: Within"]][[1]]
tableS7$Source <- rownames(tableS7)
tableS7<-format(tableS7, digits = 3)

rtffile7 <- RTF("rtf7.doc")
addParagraph(rtffile7, "This is the output for Table S7:\n")
addTable(rtffile7, tableS7)
done(rtffile7)

#END