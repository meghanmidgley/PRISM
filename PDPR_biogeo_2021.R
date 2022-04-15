#you can always start over again by clearing the memory
rm(list = ls(all=TRUE))
#try

# -----------------------------------------------------------
# load packages and relevant file paths
# -----------------------------------------------------------
library(tidyverse)
library(car)
library(multcomp)

#Load data
#file.choose()
PDPR2021.df<-read.csv("C:\\Users\\mortonarb\\Dropbox\\NSF-ROL-MOR\\DATA\\MORTON\\Soils\\PDPR2021_MasterData.csv")

# -----------------------------------------------------------
# Transform data
# -----------------------------------------------------------
PDPR2021.df$sqrtOM<-sqrt(PDPR2021.df$X..OM)
PDPR2021.df$logTOC<-log(PDPR2021.df$TOC..ug.C.g.soil.)
PDPR2021.df$logTN<-log(PDPR2021.df$TN..ug.N.g.soil.)
PDPR2021.df$logNH4<-log(PDPR2021.df$Initial.Ammonium..ug.g.)
PDPR2021.df$logNO3<-log(PDPR2021.df$Initial.Nitrate..ug.g.+1)
PDPR2021.df$loginorgN<-log(PDPR2021.df$Inorganic.N..ug.N.g.soil.)
PDPR2021.df$sqrtON<-sqrt(PDPR2021.df$Organic.N..ug.N..g.soil.)
PDPR2021.df$sqrttOCTN<-sqrt(PDPR2021.df$TOC.TN)
PDPR2021.df$logNmin<-log(PDPR2021.df$N.min..ug.N.g.soil.1.day.1.+1)
PDPR2021.df$logNitrification<-log(PDPR2021.df$Nitrification..ug.N.g.soil.1.day.1.+1)

PDPR2021poly.df<-PDPR2021.df%>%
  filter(Type=="poly")
PDPR2021mono.df<-PDPR2021.df%>%
  filter(Type=="mono")
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


