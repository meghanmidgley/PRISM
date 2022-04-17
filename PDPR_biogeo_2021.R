#you can always start over again by clearing the memory
rm(list = ls(all=TRUE))
#try

# -----------------------------------------------------------
# load packages and relevant file paths
# -----------------------------------------------------------
library(tidyverse)
library(car)
library(multcomp)
library(lmerTest)
library(phytools)
library(vegan)
library(ggtree)

#Load data
#file.choose()
PDPR2021.df<-read.csv("C:\\Users\\mortonarb\\Dropbox\\NSF-ROL-MOR\\DATA\\MORTON\\Soils\\PDPR2021_MasterData.csv")
PDPRtree2016<-read.tree("C:\\Users\\mortonarb\\Dropbox\\NSF-ROL-MOR\\DATA\\MORTON\\Soils\\phylogeny.analyzed.2016-01-05b.tre")

# -----------------------------------------------------------
# transform data and create new data frames
# -----------------------------------------------------------
PDPR2021.df$sqrtOM<-sqrt(PDPR2021.df$X..OM)
PDPR2021.df$logTOC<-log(PDPR2021.df$TOC..ug.C.g.soil.)
PDPR2021.df$logTN<-log(PDPR2021.df$TN..ug.N.g.soil.)
PDPR2021.df$logNH4<-log(PDPR2021.df$Initial.Ammonium..ug.g.)
PDPR2021.df$logNO3<-log(PDPR2021.df$Initial.Nitrate..ug.g.+1)
PDPR2021.df$loginorgN<-log(PDPR2021.df$Inorganic.N..ug.N.g.soil.)
PDPR2021.df$sqrtON<-sqrt(PDPR2021.df$Organic.N..ug.N..g.soil.)
PDPR2021.df$sqrttOCTN<-sqrt(PDPR2021.df$TOC.TN)
PDPR2021.df$logTOCTN<-log(PDPR2021.df$TOC.TN)
PDPR2021.df$logNmin<-log(PDPR2021.df$N.min..ug.N.g.soil.1.day.1.+1)
PDPR2021.df$logNitrification<-log(PDPR2021.df$Nitrification..ug.N.g.soil.1.day.1.+1)

#subset polycultures and monocultures
PDPR2021poly.df<-PDPR2021.df%>%
  filter(Type=="poly")
PDPR2021mono.df<-PDPR2021.df%>%
  filter(Type=="mono")
  
PDPR2021mono.df<-PDPR2021mono.df%>%
  mutate(Clade= case_when (PDPR2021mono.df$Group=="grammoid"|PDPR2021mono.df$Family=="liliaceae"|PDPR2021mono.df$Family=="commelinaceae"|PDPR2021mono.df$Family=="iridaceae"~ "monocot",  
                           PDPR2021mono.df$Group!="grammoid"|PDPR2021mono.df$Family!="liliaceae"|PDPR2021mono.df$Family!="commelinaceae"|PDPR2021mono.df$Family!="iridaceae"~"eudicot")) 

#Soil properties to map onto phylogenetic tree
MonoSoils.df<-PDPR2021mono.df %>% 
  group_by(Species)%>%                       
  summarise(mean(TOC..ug.C.g.soil.), mean(Inorganic.N..ug.N.g.soil.), mean(N.min..ug.N.g.soil.1.day.1.), mean(Nitrification..ug.N.g.soil.1.day.1.))%>%
  remove_rownames %>% 
  column_to_rownames(var="Species")%>%
  rename(TOC="mean(TOC..ug.C.g.soil.)", InorgN="mean(Inorganic.N..ug.N.g.soil.)", Nmin="mean(N.min..ug.N.g.soil.1.day.1.)", Nitrification="mean(Nitrification..ug.N.g.soil.1.day.1.)")

# Soil properties for phylogenetic signal detection
TOC<-MonoSoils.df[,1]
names(TOC)<-rownames(MonoSoils.df)
InorgN<-MonoSoils.df[,2]
names(InorgN)<-rownames(MonoSoils.df)
Nmin<-MonoSoils.df[,3]
names(Nmin)<-rownames(MonoSoils.df)
Nitrification<-MonoSoils.df[,4]
names(Nitrification)<-rownames(MonoSoils.df)

#Fix and prune tree
PDPRtree2016$tip.label[which(PDPRtree2016$tip.label == "Symphyotrichum_novae-angliae")] <-
  "Symphyotrichum_novae.angliae"
PDPRtree2016 <- drop.tip(PDPRtree2016, which(!PDPRtree2016$tip.label %in% PDPR2021.df$Species))

# -----------------------------------------------------------
# Analyses
# -----------------------------------------------------------

#BASIC ASSESSMENTS AND VISUALIZATIONS TO CHECK DATA

#WATER CONTENT
#evaluate normality
hist(PDPR2021.df$X..water)
qqnorm(PDPR2021.df$X..water)
qqline(PDPR2021.df$X..water)
#identify outliers
boxplot(PDPR2021.df$X..water)
boxplot(PDPR2021.df$X..water~ Block, PDPR2021.df)
#run ANOVA
water.aov<-aov(X..water~Block, PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(water.aov))
qqline(residuals(water.aov))
shapiro.test(residuals(water.aov))
#homogenity of variances
plot(residuals(water.aov)~fitted(water.aov))
bartlett.test(X..water~ Block, data = PDPR2021.df)
leveneTest(X..water ~ Block, data = PDPR2021.df)
#model results
summary(water.aov)
# get (adjusted) weight means per group
water_means <- emmeans(object = water.aov,
                       specs = "Block")
# add letters to each mean
water_means_cld <- cld(object = water_means,
                       adjust = "Tukey",
                       Letters = letters,
                       alpha = 0.05)
# show output
water_means_cld

#PH
#evaluate normality
hist(PDPR2021.df$pH)
qqnorm(PDPR2021.df$pH)
qqline(PDPR2021.df$pH)
#identify outliers
boxplot(PDPR2021.df$pH)
boxplot(PDPR2021.df$pH~ Block, PDPR2021.df)
#run ANOVA
pH.aov<-aov(pH~Block, PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(pH.aov))
qqline(residuals(pH.aov))
shapiro.test(residuals(pH.aov))
#homogenity of variances
plot(residuals(pH.aov)~fitted(pH.aov))
bartlett.test(pH~ Block, data = PDPR2021.df)
leveneTest(pH ~ Block, data = PDPR2021.df)
#model results
summary(pH.aov)
# get (adjusted) weight means per group
pH_means <- emmeans(object = pH.aov,
                       specs = "Block")
# add letters to each mean
pH_means_cld <- cld(object = pH_means,
                       adjust = "Tukey",
                       Letters = letters,
                       alpha = 0.05)
# show output
pH_means_cld

#ORGANIC MATTER
#evaluate normality
hist(PDPR2021.df$sqrtOM)
qqnorm(PDPR2021.df$sqrtOM)
qqline(PDPR2021.df$sqrtOM)
#identify outliers
boxplot(PDPR2021.df$sqrtOM)
boxplot(PDPR2021.df$sqrtOM~ Block, PDPR2021.df)
#run ANOVA
OM.aov<-aov(sqrtOM~Block, PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(OM.aov))
qqline(residuals(OM.aov))
shapiro.test(residuals(OM.aov))
#homogenity of variances
plot(residuals(OM.aov)~fitted(OM.aov))
bartlett.test(sqrtOM~ Block, data = PDPR2021.df)
leveneTest(sqrtOM ~ Block, data = PDPR2021.df)
#model results
summary(OM.aov)
# get (adjusted) weight means per group
OM_means <- emmeans(object = OM.aov,
                    specs = "Block")
# add letters to each mean
OM_means_cld <- cld(object = OM_means,
                    adjust = "Tukey",
                    Letters = letters,
                    alpha = 0.05)
# show output
OM_means_cld

#TOC
#evaluate normality
hist(PDPR2021.df$logTOC)
qqnorm(PDPR2021.df$logTOC)
qqline(PDPR2021.df$logTOC)
#identify outliers
boxplot(PDPR2021.df$logTOC)
boxplot(PDPR2021.df$logTOC~ Block, PDPR2021.df)
#run ANOVA
TOC.aov<-aov(logTOC~Block, PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(TOC.aov))
qqline(residuals(TOC.aov))
shapiro.test(residuals(TOC.aov))
#homogenity of variances
plot(residuals(TOC.aov)~fitted(TOC.aov))
bartlett.test(logTOC~ Block, data = PDPR2021.df)
leveneTest(logTOC ~ Block, data = PDPR2021.df)
#model results
summary(TOC.aov)
# get (adjusted) weight means per group
TOC_means <- emmeans(object = TOC.aov,
                    specs = "Block")
# add letters to each mean
TOC_means_cld <- cld(object = TOC_means,
                    adjust = "Tukey",
                    Letters = letters,
                    alpha = 0.05)
# show output
TOC_means_cld

#TN
#evaluate normality
hist(PDPR2021.df$logTN)
qqnorm(PDPR2021.df$logTN)
qqline(PDPR2021.df$logTN)
#identify outliers
boxplot(PDPR2021.df$logTN)
boxplot(PDPR2021.df$logTN~ Block, PDPR2021.df)
#run ANOVA
TN.aov<-aov(logTN~Block, PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(TN.aov))
qqline(residuals(TN.aov))
shapiro.test(residuals(TN.aov))
#homogenity of variances
plot(residuals(TN.aov)~fitted(TN.aov))
bartlett.test(logTN~ Block, data = PDPR2021.df)
leveneTest(logTN ~ Block, data = PDPR2021.df)
#model results
summary(TN.aov)
# get (adjusted) weight means per group
TN_means <- emmeans(object = TN.aov,
                     specs = "Block")
# add letters to each mean
TN_means_cld <- cld(object = TN_means,
                     adjust = "Tukey",
                     Letters = letters,
                     alpha = 0.05)
# show output
TN_means_cld

#NH4
#evaluate normality
hist(PDPR2021.df$logNH4)
qqnorm(PDPR2021.df$logNH4)
qqline(PDPR2021.df$logNH4)
#identify outliers
boxplot(PDPR2021.df$logNH4)
boxplot(PDPR2021.df$logNH4~ Block, PDPR2021.df)
#run ANOVA
NH4.aov<-aov(logNH4~Block, PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(NH4.aov))
qqline(residuals(NH4.aov))
shapiro.test(residuals(NH4.aov))
#homogenity of variances
plot(residuals(NH4.aov)~fitted(NH4.aov))
bartlett.test(logNH4~ Block, data = PDPR2021.df)
leveneTest(logNH4 ~ Block, data = PDPR2021.df)
#model results
summary(NH4.aov)
# get (adjusted) weight means per group
NH4_means <- emmeans(object = NH4.aov,
                    specs = "Block")
# add letters to each mean
NH4_means_cld <- cld(object = NH4_means,
                    adjust = "Tukey",
                    Letters = letters,
                    alpha = 0.05)
# show output
NH4_means_cld

#NO3
#evaluate normality
hist(PDPR2021.df$logNO3)
qqnorm(PDPR2021.df$logNO3)
qqline(PDPR2021.df$logNO3)
#identify outliers
boxplot(PDPR2021.df$logNO3)
boxplot(PDPR2021.df$logNO3~ Block, PDPR2021.df)
#run ANOVA
NO3.aov<-aov(logNO3~Block, PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(NO3.aov))
qqline(residuals(NO3.aov))
shapiro.test(residuals(NO3.aov))
#homogenity of variances
plot(residuals(NO3.aov)~fitted(NO3.aov))
bartlett.test(logNO3~ Block, data = PDPR2021.df)
leveneTest(logNO3 ~ Block, data = PDPR2021.df)
#model results
summary(NO3.aov)
# get (adjusted) weight means per group
NO3_means <- emmeans(object = NO3.aov,
                     specs = "Block")
# add letters to each mean
NO3_means_cld <- cld(object = NO3_means,
                     adjust = "Tukey",
                     Letters = letters,
                     alpha = 0.05)
# show output
NO3_means_cld

#inorgN
#evaluate normality
hist(PDPR2021.df$loginorgN)
qqnorm(PDPR2021.df$loginorgN)
qqline(PDPR2021.df$loginorgN)
#identify outliers
boxplot(PDPR2021.df$loginorgN)
boxplot(PDPR2021.df$loginorgN~ Block, PDPR2021.df)
#run ANOVA
inorgN.aov<-aov(loginorgN~Block, PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(inorgN.aov))
qqline(residuals(inorgN.aov))
shapiro.test(residuals(inorgN.aov))
#homogenity of variances
plot(residuals(inorgN.aov)~fitted(inorgN.aov))
bartlett.test(loginorgN~ Block, data = PDPR2021.df)
leveneTest(loginorgN ~ Block, data = PDPR2021.df)
#model results
summary(inorgN.aov)
# get (adjusted) weight means per group
inorgN_means <- emmeans(object = inorgN.aov,
                     specs = "Block")
# add letters to each mean
inorgN_means_cld <- cld(object = inorgN_means,
                     adjust = "Tukey",
                     Letters = letters,
                     alpha = 0.05)
# show output
inorgN_means_cld

#orgN
#evaluate normality
hist(PDPR2021.df$sqrtON)
qqnorm(PDPR2021.df$sqrtON)
qqline(PDPR2021.df$sqrtON)
#identify outliers
boxplot(PDPR2021.df$sqrtON)
boxplot(PDPR2021.df$sqrtON~ Block, PDPR2021.df)
#run ANOVA
orgN.aov<-aov(sqrtON~Block, PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(orgN.aov))
qqline(residuals(orgN.aov))
shapiro.test(residuals(orgN.aov))
#homogenity of variances
plot(residuals(orgN.aov)~fitted(orgN.aov))
bartlett.test(sqrtON~ Block, data = PDPR2021.df)
leveneTest(sqrtON ~ Block, data = PDPR2021.df)
#model results
summary(orgN.aov)
# get (adjusted) weight means per group
orgN_means <- emmeans(object = orgN.aov,
                        specs = "Block")
# add letters to each mean
orgN_means_cld <- cld(object = orgN_means,
                        adjust = "Tukey",
                        Letters = letters,
                        alpha = 0.05)
# show output
orgN_means_cld

#%inorg
#evaluate normality
hist(PDPR2021.df$X.inorg.N)
qqnorm(PDPR2021.df$X.inorg.N)
qqline(PDPR2021.df$X.inorg.N)
#identify outliers
boxplot(PDPR2021.df$X.inorg.N)
boxplot(PDPR2021.df$X.inorg.N~ Block, PDPR2021.df)
#run ANOVA
X.inorg.aov<-aov(X.inorg.N~Block, PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(X.inorg.aov))
qqline(residuals(X.inorg.aov))
shapiro.test(residuals(X.inorg.aov))
#homogenity of variances
plot(residuals(X.inorg.aov)~fitted(X.inorg.aov))
bartlett.test(X.inorg.N~ Block, data = PDPR2021.df)
leveneTest(X.inorg.N ~ Block, data = PDPR2021.df)
#model results
summary(X.inorg.aov)
# get (adjusted) weight means per group
Xinorg_means <- emmeans(object = X.inorg.aov,
                      specs = "Block")
# add letters to each mean
Xinorg_means_cld <- cld(object = Xinorg_means,
                      adjust = "Tukey",
                      Letters = letters,
                      alpha = 0.05)
# show output
Xinorg_means_cld

#avaliable C:N
#evaluate normality
hist(PDPR2021.df$sqrttOCTN)
qqnorm(PDPR2021.df$sqrttOCTN)
qqline(PDPR2021.df$sqrttOCTN)
#identify outliers
boxplot(PDPR2021.df$sqrttOCTN)
boxplot(PDPR2021.df$sqrttOCTN~ Block, PDPR2021.df)
#run ANOVA
TOCTN.aov<-aov(sqrttOCTN~Block, PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(TOCTN.aov))
qqline(residuals(TOCTN.aov))
shapiro.test(residuals(TOCTN.aov))
#homogenity of variances
plot(residuals(TOCTN.aov)~fitted(TOCTN.aov))
bartlett.test(sqrttOCTN~ Block, data = PDPR2021.df)
leveneTest(sqrttOCTN ~ Block, data = PDPR2021.df)
#model results
summary(TOCTN.aov)

#Nmin
#evaluate normality
hist(PDPR2021.df$logNmin)
qqnorm(PDPR2021.df$logNmin)
qqline(PDPR2021.df$logNmin)
#identify outliers
boxplot(PDPR2021.df$logNmin)
boxplot(PDPR2021.df$logNmin~ Block, PDPR2021.df)
#run ANOVA
Nmin.aov<-aov(logNmin~Block, PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(Nmin.aov))
qqline(residuals(Nmin.aov))
shapiro.test(residuals(Nmin.aov))
#homogenity of variances
plot(residuals(Nmin.aov)~fitted(Nmin.aov))
bartlett.test(logNmin~ Block, data = PDPR2021.df)
leveneTest(logNmin ~ Block, data = PDPR2021.df)
#model results
summary(Nmin.aov)

#Nitrification
#evaluate normality
hist(PDPR2021.df$logNitrification)
qqnorm(PDPR2021.df$logNitrification)
qqline(PDPR2021.df$logNitrification)
#identify outliers
boxplot(PDPR2021.df$logNitrification)
boxplot(PDPR2021.df$logNitrification~ Block, PDPR2021.df)
#run ANOVA
Nit.aov<-aov(logNitrification~Block, PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(Nit.aov))
qqline(residuals(Nit.aov))
shapiro.test(residuals(Nit.aov))
#homogenity of variances
plot(residuals(Nit.aov)~fitted(Nit.aov))
bartlett.test(logNitrification~ Block, data = PDPR2021.df)
leveneTest(logNitrification ~ Block, data = PDPR2021.df)
#model results
summary(Nit.aov)

#DO MONOCULTURES AND POLYCULTURES DIFFER IN THEIR SOIL PROPERTIES?

water.type.lmer<-lmer(X..water~Type+ (1|Block) + (1|Mix), data = PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(water.type.lmer))
qqline(residuals(water.type.lmer))
shapiro.test(residuals(water.type.lmer))
#homogenity of variances
plot(residuals(water.type.lmer)~fitted(water.type.lmer))
#model results
summary(water.type.lmer)
anova(water.type.lmer)

TOC.type.lmer<-lmer(logTOC~Type+ (1|Block) + (1|Mix), data = PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(TOC.type.lmer))
qqline(residuals(TOC.type.lmer))
shapiro.test(residuals(TOC.type.lmer))
#homogenity of variances
plot(residuals(TOC.type.lmer)~fitted(TOC.type.lmer))
#model results
summary(TOC.type.lmer)
anova(TOC.type.lmer)

NO3.type.lmer<-lmer(logNO3~Type+ (1|Block) + (1|Mix), data = PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(NO3.type.lmer))
qqline(residuals(NO3.type.lmer))
shapiro.test(residuals(NO3.type.lmer))
#homogenity of variances
plot(residuals(NO3.type.lmer)~fitted(NO3.type.lmer))
#model results
summary(NO3.type.lmer)
anova(NO3.type.lmer)

TOCTN.type.lmer<-lmer(logTOCTN~Type+ (1|Block) + (1|Mix), data = PDPR2021.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(TOCTN.type.lmer))
qqline(residuals(TOCTN.type.lmer))
shapiro.test(residuals(TOCTN.type.lmer))
#homogenity of variances
plot(residuals(TOCTN.type.lmer)~fitted(TOCTN.type.lmer))
#model results
summary(TOCTN.type.lmer)
anova(TOCTN.type.lmer)

#PERMANOVA?
# adonis2, data~Type+ random?

#DO MONOCULTURE SOIL PROPERTIES EXHIBIT A PHYLOGENETIC SIGNAL?

phylosig(PDPRtree2016, TOC, method="lambda", test=TRUE, nsim=999)
phylosig(PDPRtree2016, TOC, method="K", test=TRUE, nsim=999)
phylosig(PDPRtree2016, InorgN, method="lambda", test=TRUE, nsim=999)
phylosig(PDPRtree2016, InorgN, method="K", test=TRUE, nsim=999)
phylosig(PDPRtree2016, Nmin, method="lambda", test=TRUE, nsim=999)
phylosig(PDPRtree2016, Nmin, method="K", test=TRUE, nsim=999)

#DO MONOCULTURE CLADES DIFFER?

TOC.clade.lmer<-lmer(logTOC~Clade+ (1|Block)+ (1|Mix), data = PDPR2021mono.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(TOC.clade.lmer))
qqline(residuals(TOC.clade.lmer))
shapiro.test(residuals(TOC.clade.lmer))
#homogenity of variances
plot(residuals(TOC.clade.lmer)~fitted(TOC.clade.lmer))
#model results
summary(TOC.clade.lmer)
anova(TOC.clade.lmer)

NO3.clade.lmer<-lmer(logNO3~Clade+ (1|Block)+ (1|Mix), data = PDPR2021mono.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(NO3.clade.lmer))
qqline(residuals(NO3.clade.lmer))
shapiro.test(residuals(NO3.clade.lmer))
#homogenity of variances
plot(residuals(NO3.clade.lmer)~fitted(NO3.clade.lmer))
#model results
summary(NO3.clade.lmer)
anova(NO3.clade.lmer)

Nmin.clade.lmer<-lmer(logNmin~Clade+ (1|Block)+ (1|Mix), data = PDPR2021mono.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(Nmin.clade.lmer))
qqline(residuals(Nmin.clade.lmer))
shapiro.test(residuals(Nmin.clade.lmer))
#homogenity of variances
plot(residuals(Nmin.clade.lmer)~fitted(Nmin.clade.lmer))
#model results
summary(Nmin.clade.lmer)
anova(Nmin.clade.lmer)

#DO MONOCULTURE FUNCTIONAL GROUPS DIFFER?

TOC.group.lmer<-lmer(logTOC~Group+ (1|Block)+ (1|Mix), data = PDPR2021mono.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(TOC.group.lmer))
qqline(residuals(TOC.group.lmer))
shapiro.test(residuals(TOC.group.lmer))
#homogenity of variances
plot(residuals(TOC.group.lmer)~fitted(TOC.group.lmer))
#model results
summary(TOC.group.lmer)
anova(TOC.group.lmer)

NO3.group.lmer<-lmer(logNO3~Group+ (1|Block)+ (1|Mix), data = PDPR2021mono.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(NO3.group.lmer))
qqline(residuals(NO3.group.lmer))
shapiro.test(residuals(NO3.group.lmer))
#homogenity of variances
plot(residuals(NO3.group.lmer)~fitted(NO3.group.lmer))
#model results
summary(NO3.group.lmer)
anova(NO3.group.lmer)

Nmin.group.lmer<-lmer(logNmin~Group+ (1|Block)+ (1|Mix), data = PDPR2021mono.df)
#check model assumptions
#normality of residuals
qqnorm(residuals(Nmin.group.lmer))
qqline(residuals(Nmin.group.lmer))
shapiro.test(residuals(Nmin.group.lmer))
#homogenity of variances
plot(residuals(Nmin.group.lmer)~fitted(Nmin.group.lmer))
#model results
summary(Nmin.group.lmer)
anova(Nmin.group.lmer)

#DO C3 and C4 grasses differ?

#DO POLYCULTURE LEVELS DIFFER?  
#lmers and PERMANOVA

#variability differ? hypothesis: lower diversity will have greater dispersion than high diversity, but centroids will be the same


# -----------------------------------------------------------
# figures
# -----------------------------------------------------------

#TOC across the experiment (map), ideal for both 2016 and 2021

#TOC for mono and poly
TOCmp.plot<- ggplot(PDPR2021.df, aes(x = Type, y = TOC..ug.C.g.soil.)) +
  theme_linedraw() +
  geom_boxplot(notch=F, 
               lwd=.5, 
               colour='black', 
               stat="boxplot", 
               outlier.shape=NA)+
  guides(fill=FALSE)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill=NA, colour = "black", size=.7), 
        axis.title.x = element_text(margin = margin(t = 10, b=5), size=16),
        axis.title.y = element_text(margin = margin(l = 5, r=5), size=16), 
        axis.text.x= element_text(margin = margin(t = 10), size=14),
        axis.text.y=element_text(margin = margin(r = 10), size=12), 
        axis.ticks.length=unit(-0.25, "cm"),
        axis.ticks.margin=unit(0.5, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.4), 
        axis.ticks.x = element_blank())+
  labs(x='Plot type', 
       y= expression(atop("Total organic C",paste(~(mu*g~C~g~soil^{-1})))))+
  scale_x_discrete(limits = c("mono","poly"), 
                   labels = c("Mono", "Poly"))+
  geom_jitter(position=position_jitter(0.07)) 

TOCmp.plot

#NO3 for mono and poly
NO3mp.plot<- ggplot(PDPR2021.df, aes(x = Type, y = Initial.Nitrate..ug.g.)) +
  theme_linedraw() +
  geom_boxplot(notch=F, 
               lwd=.5, 
               colour='black', 
               stat="boxplot", 
               outlier.shape=NA)+
  guides(fill=FALSE)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill=NA, colour = "black", size=.7), 
        axis.title.x = element_text(margin = margin(t = 10, b=5), size=16),
        axis.title.y = element_text(margin = margin(l = 5, r=5), size=16), 
        axis.text.x= element_text(margin = margin(t = 10), size=14),
        axis.text.y=element_text(margin = margin(r = 10), size=12), 
        axis.ticks.length=unit(-0.25, "cm"),
        axis.ticks.margin=unit(0.5, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.4), 
        axis.ticks.x = element_blank())+
  labs(x='Plot type', 
       y= expression(atop("Nitrate",paste(~(mu*g~N-NO[3]~g~soil^{-1})))))+
  scale_x_discrete(limits = c("mono","poly"), 
                   labels = c("Mono", "Poly"))+
  geom_jitter(position=position_jitter(0.07)) 

NO3mp.plot

#TOC:TN for mono and poly
TOCTNmp.plot<- ggplot(PDPR2021.df, aes(x = Type, y = TOC.TN)) +
  theme_linedraw() +
  geom_boxplot(notch=F, 
               lwd=.5, 
               colour='black', 
               stat="boxplot", 
               outlier.shape=NA)+
  guides(fill=FALSE)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill=NA, colour = "black", size=.7), 
        axis.title.x = element_text(margin = margin(t = 10, b=5), size=16),
        axis.title.y = element_text(margin = margin(l = 5, r=5), size=16), 
        axis.text.x= element_text(margin = margin(t = 10), size=14),
        axis.text.y=element_text(margin = margin(r = 10), size=12), 
        axis.ticks.length=unit(-0.25, "cm"),
        axis.ticks.margin=unit(0.5, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.4), 
        axis.ticks.x = element_blank())+
  labs(x='Plot type', 
       y= 'TOC:TN')+
  scale_x_discrete(limits = c("mono","poly"), 
                   labels = c("Mono", "Poly"))+
  geom_jitter(position=position_jitter(0.07)) 

TOCTNmp.plot

#TOC for monocots and dicots
TOCmd.plot<- ggplot(PDPR2021mono.df, aes(x = Clade, y = TOC..ug.C.g.soil.)) +
  theme_linedraw() +
  geom_boxplot(notch=F, 
               lwd=.5, 
               colour='black', 
               stat="boxplot", 
               outlier.shape=NA)+
  guides(fill=FALSE)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill=NA, colour = "black", size=.7), 
        axis.title.x = element_text(margin = margin(t = 10, b=5), size=16),
        axis.title.y = element_text(margin = margin(l = 5, r=5), size=16), 
        axis.text.x= element_text(margin = margin(t = 10), size=14),
        axis.text.y=element_text(margin = margin(r = 10), size=12), 
        axis.ticks.length=unit(-0.25, "cm"),
        axis.ticks.margin=unit(0.5, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.4), 
        axis.ticks.x = element_blank())+
  labs(x='Clade', 
       y= expression(atop("Total organic C",paste(~(mu*g~C~g~soil^{-1})))))+
  scale_x_discrete(limits = c("monocot","eudicot"), 
                   labels = c("Monocot", "Eudicot"))+
  geom_jitter(position=position_jitter(0.07)) 

TOCmd.plot

#phylogeny with heatmap of soil properties (TOC, inorg N, N min?)
p1 <- ggtree(PDPRtree2016, layout="circular") +
  geom_cladelabel(node=192, label="Asteraceae", offset = 70) +
  geom_cladelabel(node=233, label="Lamiaceae", offset = 70, hjust = 1.1) +
  geom_cladelabel(node=173, "Fabaceae", offset = 70, hjust = 1) +
  geom_cladelabel(node=169, "Rosaceae", offset = 70, hjust =1.1) +
  geom_cladelabel(node=147, "Poaceae", offset = 70, hjust = -.1)

legend_title1 <- c("TOC")
fig1 <- p1 + new_scale_fill()
fig12 <- gheatmap(fig1, MonoSoils.df[,"TOC", drop=FALSE], offset=0, width=.175,
                  colnames = FALSE) +
  scale_fill_gradient(na.value = "white", low = "white", high = "#ae3918", legend_title1) + #blue
  theme(legend.position="none", legend.title=element_text() )  

legend_title2 <- c("inorg N")
p2 <- fig12 + new_scale_fill()
fig13 <- gheatmap(p2, MonoSoils.df[,"InorgN", drop=FALSE], offset=0.8, width=.1,
                  colnames = FALSE) +
  scale_fill_gradient(na.value = "white", low = "white", high = "#027ab0", legend_title2) + #blue
  theme(legend.position="none", legend.title=element_text() ) 

legend_title3 <- c("N min")
p3 <- fig13 + new_scale_fill()
fig14 <- gheatmap(p3, MonoSoils.df[,"Nmin", drop=FALSE], offset=40, width=.1,
                  colnames = FALSE) +
  scale_fill_gradient(na.value = "white", low = "white", high = "#49997c", legend_title2) + #blue
  theme(legend.position="none", legend.title=element_text() ) 

fig14

#NMDS ordinations of monocultures vs polycultures, monoculture FGs, polyculture levels
#using vegan function 'metaMDS' with Bray-Curtis dissimilarity 




