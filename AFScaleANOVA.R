#

#     AFScaleANOVA.R
#     /DS/RhDec20/AFScale4/
#
#  #  March 19, 2021 
#
# ANOVA for AF Scale values 
#  TrtShort only covariate
#  including year would not be workable

# Following:
# http://www.sthda.com/english/wiki/one-way-anova-test-in-r

library(ggplot2)
library(car)
# install.packages("multcomp")
library(multcomp)

#   set location of the directory containing the datasets
in1DataPath="/DS/RhDec20/AFScale4/Results/"
outDataPath="/DS/RhDec20/AFScale4/Results/"
setwd(in1DataPath)

#  = = = = LOAD RMM INFO FOR AF = = = = = 

f=paste0(in1DataPath,"AFScalePCAstuff.RData")
load(f)   #  data for 35 samples

# Compute the analysis of variance
res.aov <- aov(PC1 ~ TrtShort, data = AF_ScoreAllSamples)
# Summary of the analysis
summary(res.aov)
#  overall test p<0.001

# 1. Homogeneity of variances
plot(res.aov, 1)  # TC outliers
leveneTest(PC1 ~ TrtShort, data = AF_ScoreAllSamples) # signif
#  Normality of residuals
plot(res.aov, 2) # TC outliers again

# GO AHEAD WITH THIS ANALYSIS METHOD
# POOLED SIGMA2 is probably an overestimate due to TC
# outliers, but this should just make analysis more conservative


#TukeyHSD(res.aov)  #  all pairwise

#  just compare to PTC
AF_ScoreAllSamples$TrtShortF <- relevel(AF_ScoreAllSamples$TrtShortF, ref = "PTC")

# Check that PTC is the reference category (appears 1st):
levels(AF_ScoreAllSamples$TrtShortF)

# Compute the analysis of variance
res.aov <- aov(PC1 ~ TrtShortF, data = AF_ScoreAllSamples)
# Summary of the analysis
summary(res.aov)
#  overall test p<0.001--checks

# Dunnett's test:
post_test <- glht(res.aov,linfct = mcp(TrtShortF = "Dunnett"))

summry=summary(post_test)
summry

#   TO PASTE INTO SPREADSHEET

res.aov <- aov(PC1 ~ TrtShortF, data = AF_ScoreAllSamples)

summary(res.aov)

post_test <- glht(res.aov,linfct = mcp(TrtShortF = "Dunnett"))

summary(post_test)

1
2
3
quit()
