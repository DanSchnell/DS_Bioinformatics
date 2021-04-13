#

#     LungScaleANOVA.R
#     /DS/RhDec20/LungScale4/
#
#  #  March 18, 2021 
#
# ANOVA for Lung Scale values 
#  TrtShort only covariate
#  including year would not be workable

# Following:
# http://www.sthda.com/english/wiki/one-way-anova-test-in-r

library(ggplot2)
library(car)
# install.packages("multcomp")
library(multcomp)

#   set location of the directory containing the datasets
in1DataPath="/DS/RhDec20/LungScale4/Results/"
outDataPath="/DS/RhDec20/LungScale4/Results/"
setwd(in1DataPath)

#  = = = = LOAD RMM INFO FOR LUNG = = = = = 

f=paste0(in1DataPath,"LungScalePCAstuff.RData")
load(f)   #  data for 35 samples

# Compute the analysis of variance
res.aov <- aov(PC1 ~ TrtShort, data = ScoreAllSamples)
# Summary of the analysis
summary(res.aov)
#  overall test p<0.001

# 1. Homogeneity of variances
plot(res.aov, 1)  # looks ~OK
leveneTest(PC1 ~ TrtShort, data = ScoreAllSamples) # OK
#  Normality of residuals
plot(res.aov, 2) # not bad

TukeyHSD(res.aov)  #  all pairwise

#  just compare to PTC
ScoreAllSamples$TrtShortF <- relevel(ScoreAllSamples$TrtShortF, ref = "PTC")

# Check that PTC is the reference category (appears 1st):
levels(ScoreAllSamples$TrtShortF)

# Compute the analysis of variance
res.aov <- aov(PC1 ~ TrtShortF, data = ScoreAllSamples)
# Summary of the analysis
summary(res.aov)
#  overall test p<0.001--checks

# Dunnett's test:
post_test <- glht(res.aov,linfct = mcp(TrtShortF = "Dunnett"))

summry=summary(post_test)
summry

1
2
3
quit()
