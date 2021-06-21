#Code for Figure 5 of "Benthic micro- and macro-community succession and coral recruitment under overfishing and nutrient enrichment"
#Also includes statistical output from Table S11

#setwd to source file location (data file and R script are in the same folder)

Survival_data<-read.csv("Survival_data.csv")
# Organize data for Coxme
Long_format<-data.frame(Survival_data[rep(seq_len(dim(Survival_data)[1]), Survival_data$Count), , drop = FALSE], row.names=NULL)
Long_format[,"binomial"]<-as.character(NA)
Long_format$binomial[grep("Alive", Long_format$Status)]<-"0"
Long_format$binomial[grep("Dead", Long_format$Status)]<-"1"

#Now running coxme
library(survival)
library(coxme)

Long_format$binomial<-as.numeric(Long_format$binomial)

library(AICcmodavg)
Can.mod<-list()
Can.mod[[1]]<-coxme(Surv(Time, binomial) ~ Caging*Nutrients*Microhabitat + (1 | Plot), data = Long_format)
Can.mod[[2]]<-coxme(Surv(Time, binomial) ~ Caging*Nutrients + Caging*Microhabitat + Nutrients*Microhabitat + (1 | Plot), data = Long_format)
Can.mod[[3]]<-coxme(Surv(Time, binomial) ~ Caging*Nutrients + Caging*Microhabitat + (1 | Plot), data = Long_format)
Can.mod[[4]]<-coxme(Surv(Time, binomial) ~ Caging*Nutrients + Nutrients*Microhabitat + (1 | Plot), data = Long_format)
Can.mod[[5]]<-coxme(Surv(Time, binomial) ~ Nutrients*Microhabitat + Caging*Microhabitat + (1 | Plot), data = Long_format)
Can.mod[[6]]<-coxme(Surv(Time, binomial) ~ Caging + Nutrients + Microhabitat + (1 | Plot), data = Long_format)
Can.mod[[7]]<-coxme(Surv(Time, binomial) ~ Caging + Nutrients + (1 | Plot), data = Long_format)
Can.mod[[8]]<-coxme(Surv(Time, binomial) ~ Caging + Microhabitat + (1 | Plot), data = Long_format)
Can.mod[[9]]<-coxme(Surv(Time, binomial) ~ Nutrients + Microhabitat + (1 | Plot), data = Long_format)
Can.mod[[10]]<-coxme(Surv(Time, binomial) ~ Nutrients + (1 | Plot), data = Long_format)
Can.mod[[11]]<-coxme(Surv(Time, binomial) ~ Caging + (1 | Plot), data = Long_format)
Can.mod[[12]]<-coxme(Surv(Time, binomial) ~ Microhabitat + (1 | Plot), data = Long_format)
Can.mod[[13]]<-coxme(Surv(Time, binomial) ~ (1 | Plot), data = Long_format)

aic_out<-aictab(cand.set = Can.mod)
aic_out
evidence(aic_out) #best model fit (per AIC) includes interaction of Nutrients x Caging x Microhabitat

anova(Can.mod[[1]]) 

tableS11<-anova(Can.mod[[1]]) 
tableS11$Source <- rownames(tableS11)
tableS11<-format(tableS11, digits = 3)

rtffile11 <- RTF("rtf11.doc")
addParagraph(rtffile11, "This is the output for Table S11:\n")
addTable(rtffile11, tableS11)
done(rtffile11)

#KM plots
library(survminer)
library(ggpubr)

#plot by nutrient + caging treatment for crowns
crown_surv <- Long_format[which(Long_format$Microhabitat=='Crown'), ]
crown_fit <- survfit(Surv(Time, binomial) ~ Nutrients+Caging, data = crown_surv)
ggsurvplot(crown_fit, conf.int = TRUE, conf.int.style = c("ribbon"), conf.int.alpha = 0.2, palette = c('navy','steelblue1','dodgerblue4','darkgreen','palegreen3','palegreen4'), 
           censor = FALSE, linetype = c("solid", "dotted","dashed","solid", "dotted","dashed")) 

#plot by nutrient + caging treatment for crevices
crevice_surv <- Long_format[which(Long_format$Microhabitat=='Crevice'), ]
crevice_fit <- survfit(Surv(Time, binomial) ~ Nutrients+Caging, data = crevice_surv)
ggsurvplot(crevice_fit, conf.int = TRUE, conf.int.style = c("ribbon"), conf.int.alpha = 0.2, palette = c('navy','steelblue1','dodgerblue4','darkgreen','palegreen3','palegreen4'), 
           censor = FALSE, linetype = c("solid", "dotted","dashed","solid", "dotted","dashed"))

#END