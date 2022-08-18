#

#     BrainScaleANOVAREV.R
#     /DS/RhDec20/BrainScale4/newHippoScale
#
#  #  November 17, 2021 
#
# ANOVA for Brain Scale values 
#  TrtShort only covariate
#  including year would not be workable
#
#    REVised having switched to Science paper DEGs


# Following:
# http://www.sthda.com/english/wiki/one-way-anova-test-in-r

library(ggplot2)
library(car)
# install.packages("multcomp")
library(multcomp)

#   set location of the directory containing the datasets
in1DataPath="/DS/RhDec20/BrainScale4/newHippoScale/"
outDataPath="/DS/RhDec20/BrainScale4/newHippoScale/"
setwd(in1DataPath)

#  = = = = LOAD RMM INFO FOR LUNG = = = = = 

f=paste0(in1DataPath,"BrainScalePCAstuffREV.RData")
load(f)   #  data for 46 samples

# Compute the analysis of variance
res.aov <- aov(PC1 ~ TrtShort, data = Brain_ScoreAllSamples)
# Summary of the analysis
summary(res.aov)
#  overall test p<0.001

# 1. Homogeneity of variances
plot(res.aov, 1)  # looks ~OK
leveneTest(PC1 ~ TrtShort, data = Brain_ScoreAllSamples) # OK
#  Normality of residuals
plot(res.aov, 2) # not bad

#TukeyHSD(res.aov)  #  all pairwise

#  just compare to PTC
Brain_ScoreAllSamples$TrtShortF <- relevel(Brain_ScoreAllSamples$TrtShortF, ref = "PTC")

# Check that PTC is the reference category (appears 1st):
levels(Brain_ScoreAllSamples$TrtShortF)

# Compute the analysis of variance
res.aov <- aov(PC1 ~ TrtShortF, data = Brain_ScoreAllSamples)
# Summary of the analysis
summary(res.aov)
#  overall test p ~ 10 e-7 --checks

# Dunnett's test:
post_test <- glht(res.aov,linfct = mcp(TrtShortF = "Dunnett"))

summry=summary(post_test)
summry

#   TO PASTE INTO SPREADSHEET

res.aov <- aov(PC1 ~ TrtShortF, data = Brain_ScoreAllSamples)

summary(res.aov)

post_test <- glht(res.aov,linfct = mcp(TrtShortF = "Dunnett"))

summary(post_test)




1
2
3
quit()
