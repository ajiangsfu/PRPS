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
data("DHITsig")
str(DHITsig)
head(DHITsig)

## ------------------------------------------------------------------------
dat = rosenwald.expr
trainset = subset(rosenwald.cli, rosenwald.cli$set == "Training")
testset = subset(rosenwald.cli, rosenwald.cli$set == "Validation")

## ------------------------------------------------------------------------
nin = 100 
trainLPS = LPStraining (trainDat = dat[,rownames(trainset)], groupInfo = trainset$group, refGroup = "GCB", topN = nin,
                      weightMethod = "ttest")
str(trainLPS)


## ------------------------------------------------------------------------
trainLPS$classTable

## ------------------------------------------------------------------------
?(LPStraining)
help(LPStraining)

## ------------------------------------------------------------------------
testLPS = LPStesting(LPStrainObj = trainLPS, newdat = dat[,rownames(testset)])
str(testLPS)
table(testLPS$LPS_class)
table(testLPS$LPS_class, testset$group)

## ------------------------------------------------------------------------
fullLPS = LPStesting(LPStrainObj = trainLPS, newdat = dat)
str(fullLPS)
table(fullLPS$LPS_class)

## ------------------------------------------------------------------------
table(testLPS$LPS_class, fullLPS[rownames(testset), "LPS_class"])


## ------------------------------------------------------------------------
?(LPStesting)
help(LPStesting)

## ------------------------------------------------------------------------
nin = 100 
trainPRPS = PRPStraining (trainDat = dat[,rownames(trainset)], groupInfo = trainset$group, refGroup = "GCB", topN = nin,
                      weightMethod = "ttest")
str(trainPRPS)


## ------------------------------------------------------------------------
trainPRPS$classTable

## ------------------------------------------------------------------------
?(PRPStraining)
help(PRPStraining)

## ------------------------------------------------------------------------
testPRPS = PRPStesting(PRPStrainObj = trainPRPS, newdat = dat[,rownames(testset)])
str(testPRPS)
table(testPRPS$PRPS_class)
table(testPRPS$PRPS_class, testset$group)



## ------------------------------------------------------------------------
fullPRPS = PRPStesting(PRPStrainObj = trainPRPS, newdat = dat)
str(fullPRPS)
table(fullPRPS$PRPS_class)

## ------------------------------------------------------------------------
table(testPRPS$PRPS_class, fullPRPS[rownames(testset), "PRPS_class"])

## ------------------------------------------------------------------------
nin = 100 
trainPS = PStraining (trainDat = dat[,rownames(trainset)], groupInfo = trainset$group, refGroup = "GCB", topN = nin,
                         weightMethod = "ttest")
str(trainPS)


## ------------------------------------------------------------------------
trainPS$classTable

## ------------------------------------------------------------------------
?(PStraining)
help(PStraining)

## ------------------------------------------------------------------------
testPS = PStesting(PStrainObj = trainPS, newdat = dat[,rownames(testset)])
str(testPS)
table(testPS$PS_class)
table(testPS$PS_class, testset$group)

## ------------------------------------------------------------------------
fullPS = PStesting(PStrainObj = trainPS, newdat = dat)
str(fullPS)
table(fullPS$PS_class)

## ------------------------------------------------------------------------
table(testPS$PS_class, fullPS[rownames(testset), "PS_class"])

## ------------------------------------------------------------------------
lpstrain = LPStraining (trainDat = dat[,rownames(trainset)], standardization=TRUE, groupInfo = trainset$group, refGroup = "GCB", topN = nin,weightMethod = "ttest")
str(trainLPS)

lpstest = LPStesting(lpstrain, newdat = dat[,rownames(testset)], standardization = TRUE)
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
lpsscores = getClassScores(testdat = dat[,rownames(testset)], classMethod = "LPS", weights = lpswts)
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

