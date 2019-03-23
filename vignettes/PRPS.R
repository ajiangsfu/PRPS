## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)
library(dplyr)
library(ggplot2)
set.seed(1014)

## ------------------------------------------------------------------------
library(PRPS)
data("rosenwald.cli")
str(rosenwald.cli)

## ------------------------------------------------------------------------
data("rosenwald.expr")
str(rosenwald.expr)

## ------------------------------------------------------------------------
data("WrightCOO")
str(WrightCOO)
head(WrightCOO)

## ------------------------------------------------------------------------
data("DHITsigGenes")
str(DHITsigGenes)
head(DHITsigGenes)

## ------------------------------------------------------------------------
group = rosenwald.cli$group
dat = rosenwald.expr

## ------------------------------------------------------------------------
nin = 100 
trainLPS = LPStraining (trainDat = dat, groupInfo = group, refGroup = "GCB", topN = nin,
                      weightMethod = "ttest")
str(trainLPS)


## ------------------------------------------------------------------------
trainLPS$classTable

## ------------------------------------------------------------------------
?(LPStraining)
help(LPStraining)

## ------------------------------------------------------------------------
set.seed(1027562)
N = 200  ### sample size for testing
testdat = sample(1:dim(dat)[2], size = N, replace = TRUE)
testdat = dat[, testdat]
str(testdat)


## ------------------------------------------------------------------------
testLPS = LPStesting(LPStrainObj = trainLPS, newdat = testdat)
str(testLPS)
table(testLPS$LPS_class)

## ------------------------------------------------------------------------
fullLPS = LPStesting(LPStrainObj = trainLPS, newdat = cbind(testdat, dat))
str(fullLPS)
table(fullLPS$LPS_class)

## ------------------------------------------------------------------------
table(testLPS$LPS_class, fullLPS$LPS_class[1:N])

## ------------------------------------------------------------------------
?(LPStesting)
help(LPStesting)

## ------------------------------------------------------------------------
nin = 100 
trainPRPS = PRPStraining (trainDat = dat, groupInfo = group, refGroup = "GCB", topN = nin,
                      weightMethod = "ttest")
str(trainPRPS)


## ------------------------------------------------------------------------
trainPRPS$classTable

## ------------------------------------------------------------------------
?(PRPStraining)
help(PRPStraining)

## ------------------------------------------------------------------------
testPRPS = PRPStesting(PRPStrainObj = trainPRPS, newdat = testdat)
str(testPRPS)
table(testPRPS$PRPS_class)

## ------------------------------------------------------------------------
fullPRPS = PRPStesting(PRPStrainObj = trainPRPS, newdat = cbind(testdat, dat))
str(fullPRPS)
table(fullPRPS$PRPS_class)

## ------------------------------------------------------------------------
table(testPRPS$PRPS_class, fullPRPS$PRPS_class[1:N])

## ------------------------------------------------------------------------
nin = 100 
trainPS = PStraining (trainDat = dat, groupInfo = group, refGroup = "GCB", topN = nin,
                         weightMethod = "ttest")
str(trainPS)


## ------------------------------------------------------------------------
trainPS$classTable

## ------------------------------------------------------------------------
?(PStraining)
help(PStraining)

## ------------------------------------------------------------------------
testPS = PStesting(PStrainObj = trainPS, newdat = testdat)
str(testPS)
table(testPS$PS_class)

## ------------------------------------------------------------------------
fullPS = PStesting(PStrainObj = trainPS, newdat = cbind(testdat, dat))
str(fullPS)
table(fullPS$PS_class)

## ------------------------------------------------------------------------
table(testPS$PS_class, fullPS$PS_class[1:N])

## ------------------------------------------------------------------------
lpstrain = LPStraining (trainDat = dat,standardization=TRUE, groupInfo = group, refGroup = "GCB", topN = nin,weightMethod = "ttest")
str(trainLPS)

lpstest = LPStesting(lpstrain, newdat = testdat, standardization = TRUE)
str(lpstest)

## ------------------------------------------------------------------------
lpswts = trainLPS$LPS_pars$weights
str(lpswts)

## ------------------------------------------------------------------------
genes = rownames(lpswts)
lpswts = lpswts[,1]
names(lpswts) = genes
head(lpswts)


## ------------------------------------------------------------------------
lpsscores = getClassScores(testdat, classMethod = "LPS", weights = lpswts)
head(lpsscores)

## ------------------------------------------------------------------------
plotNames = c("lpsPlots","prpsPlots","psPlots")
trainNames = paste("train",plotNames, sep="_")
plotTraining(trainObj = trainLPS, plotName = trainNames[1])
plotTraining(trainObj = trainPRPS, plotName = trainNames[2])
plotTraining(trainObj = trainPS, plotName = trainNames[3])

## ------------------------------------------------------------------------
testNames = paste("test",plotNames, sep="_")
plotTesting(testObj = testLPS, plotName = testNames[1])
plotTesting(testObj = testPRPS, plotName = testNames[2])
plotTesting(testObj = testPS, plotName = testNames[3])

