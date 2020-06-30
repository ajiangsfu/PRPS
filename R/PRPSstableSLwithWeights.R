
#' PRPS stable self-training 
#' @description This function is to calculate PRPS (Probability ratio based classification predication score) scores and 
#' make binary classification calls for a testing data set without PRPS training object. It involves a self-training 
#' process with given features and their weights.
#'
#' @details This function is trying to get reasonable PRPS based classification without training data set, but with 
#' selected features and their weights. The actual steps are as following:
#' 1) assume that we have a pool for group ratio priors such as seq(0.05, 0.95, by = 0.05) for default ratioRange = c(0.05, 0.95)
#' 2) With given features and their weights
#'    a) for each prior in 1), call PSSLwithWeightsPrior with given features and weights to achieve PRPS scores
#'       apply EM on PRPS scores with Mclust, get 2 group classification
#'    b) define the samples that are always in the same classes across searching range as stable classes
#' 3) repeat step 2) but this time with opposite signs in the given weights, result in another set of stable classes
#' 4) get final stable classes that are common in 2) and 3)
#' 5) use final stable classes to get group means and sds for each feature and for each group
#' 5) calculate PRPS scores
#' 6) Once we have PRPS scores, we could use the theoretic natual cutoff 0 to make classification calls, which may or may not appropriate. 
#' Alternatively, with two groups based on stable classes assuming that PRPS score is a mixture of two normal distributions, 
#' we can get Empirical Bayesian probability and make calls
#' @param newdat a input data matrix or data frame, columns for samples and rows for features
#' @param weights a numeric vector with selected features (as names of the vector) and their weights
#' @param plotName a pdf file name with full path and is ended with ".pdf", which is used to save multiple pages 
#'  of PRPS histgrams with distribution densities. Default value us NULL, no plot is saved.
#' @param ratioRange a numeric vector with two numbers, which indicates ratio search range. The default is
#'  c(0.05, 0.95), which should NOT be changed in most of situations. However, if your classification is very
#'  unbalanced such as one group is much smaller than the other, and/or sample variation is quite big,
#'  and/or classification results are far away from what you expect, you might want to change the default values.
#'  c(0.15, 0.85) is recommended as an alternative setting other than default. In an extreme rare situation, c(0.4, 0,6) could a good try.
#' @param stepby a numeric parameter for distance between percentage searching step, it should be within (0,1), default value is 0.05, 
#'  but a user can change it to other values such as 0.01
#' @param standardization a logic variable to indicate if standardization is needed before classification 
#'  score calculation
#' @param classProbCut a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
#'  often used values are 0.8 and 0.9, default value is 0.9. Only one value is used for both groups, 
#'  the samples that are not included in either group will be assigned as UNCLASS
#' @param PRPShighGroup a string to indicate group name with high PRPS score
#' @param PRPSlowGroup a string to indicate group name with low PRPS score
#' @param breaks a integer to indicate number of bins in histogram, default is 50
#' @param imputeNA a logic variable to indicate if NA imputation is needed, if it is TRUE, NA imputation is 
#'  processed before any other steps, the default is FALSE
#' @param byrow a logic variable to indicate direction for imputation, default is TRUE, 
#'  which will use the row data for imputation
#' @param imputeValue a character variable to indicate which value to be used to replace NA, default is "median", 
#'  the median value of the chose direction with "byrow" data to be used
#' @return A list with two items is returned: PRPS parameters for selected features, PRPS scores and classifications for the given samples.
#' \item{PRPS_pars}{a list of 3 items, the 1st item is a data frame with weights of each selected features for PRPS
#'  calculation, the 2nd item is a numeric vector containing PRPS mean and sd for two groups，the 3rd item is a data frame contains mean and sd
#'   for each group and for each selected feature based on stable classes}
#' \item{PRPS_test}{a data frame of PRPS score, classification and two groups' Empirical Bayesian probabilites based on stable classes, 
#' and classification with natural 0 cutoff}
#' @keywords PRPS EM self-learning
#' @author Aixiang Jiang
#' @references Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P, 
#'  Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, Kridel R, 
#'  Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. Double-Hit Gene Expression Signature Defines
#'  a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell Lymphoma. J Clin Oncol. 
#'  2018 Dec 3:JCO1801583. doi: 10.1200/JCO.18.01583.
#' 
#' Ultsch, A., Thrun, M.C., Hansen-Goos, O., Loetsch, J.: Identification of Molecular Fingerprints
#' in Human Heat Pain Thresholds by Use of an Interactive Mixture Model R Toolbox(AdaptGauss),
#' International Journal of Molecular Sciences, doi:10.3390/ijms161025897, 2015.
#' 
#' Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and 
#' density estimation using Gaussian finite mixture models, The R Journal, 8/1, pp. 205-233.
#' 
#' @import mclust
#' @export

PRPSstableSLwithWeights = function(newdat, weights, plotName = NULL, ratioRange = c(0.05, 0.95), stepby = 0.05, standardization=FALSE, classProbCut = 0.9, 
                    PRPShighGroup = "PRPShigh", PRPSlowGroup = "PRPSlow", breaks = 50, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  imputeValue = imputeValue[1]
  ## imputee NA if imputeNA is true
  if(imputeNA){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for PRPS approach, it does not require standardization, however, if standardization = TRUE, do the standardization
  if(standardization){newdat = standardize(newdat)}
  
  ### in case that some features in weight list are actually not in data, do the following
  tmp = intersect(names(weights), rownames(newdat))
  weights = weights[tmp]
  newdat = newdat[tmp,]
  
  rps = seq(ratioRange[1], ratioRange[2], by = stepby)
  
  if(!is.null(plotName)){
    pdf(plotName)
  }
  
  rpsres = sapply(rps, FUN = function(xx){
    tmp = PRPSSLwithWeightsPrior(newdat=newdat, weights=weights, ratioPrior = xx, PRPShighGroup = PRPShighGroup, PRPSlowGroup = PRPSlowGroup)
    ### note on 20200306, since I only keep the PRPS score from PRPSSLwithWeightsPrior and ignore its classification, I am fine
    ###       otherwise I might have to consider if theoretical score cut at 0 is good or not
    
    ######  add plot EM step back for this local function on 20191206 ######################
    emsearch = plotHistEM(tmp$PRPS_test$PRPS_score, G = 2, breaks = breaks, 
                          scoreName = paste("PRPS_score with rho = ", xx, sep=""))
    #########################################################################
    
    mcls = mclust::Mclust(tmp$PRPS_test$PRPS_score, G=2, warn = TRUE)
    ### note on 20200306, since G = 2 is set in plotHistEM, the above line is the same as shown in plot
    ###  this means that I actually did the same Mclust for two times, which I do not need to bather to export the Mclust project
    ###   -> which is reasonable, since plot is usually the end-function without any return
    ### alternatively, I can run Mclust first, and then inout Mclust into plot function
    ###   -> but the problem with this is: sometimes I do not have Mclust before I call the plot function
    ### putting together, the current approach is acceptable, so no change is needed
    
    return(mcls$classification)
  })
  
  rownames(rpsres) = colnames(newdat)
  ### the classification is coded with 1 and 2
  ### in order to find samples that are always 1 and always 2, what should I do?
  ### for always 1, if I extract 1 for each cell, then row sum for the sample should be 0
  ### for always 2, if I extract 2 for each cell, then row sum for the sample should be 0
  
  res1 = rpsres - 1
  res2 = rpsres - 2
  
  grp1 = rownames(res1[which(rowSums(res1) == 0),])
  grp2 = rownames(res1[which(rowSums(res2) == 0),])
  
  ########### changes on Oct 31, 2019 #######################
  orpsres = sapply(rps, FUN = function(xx){
    tmp = PRPSSLwithWeightsPrior(newdat=newdat, weights=-weights, ratioPrior = xx, PRPShighGroup = PRPSlowGroup, PRPSlowGroup = PRPShighGroup)
    
    ######  add plot EM step back for this local function on 20191206 ######################
    emsearch = plotHistEM(tmp$PRPS_test$PRPS_score, G = 2, breaks = breaks, 
                          scoreName = paste("Reverse weight PRPS_score with rho = ", xx, sep=""))
    #########################################################################
    
    mcls = mclust::Mclust(tmp$PRPS_test$PRPS_score, G=2, warn = TRUE)
    return(mcls$classification)
  })
  
  rownames(orpsres) = colnames(newdat)
  ### the classification is coded with 1 and 2
  ### in order to find samples that are always 1 and always 2, what should I do?
  ### for always 1, if I extract 1 for each cell, then row sum for the sample should be 0
  ### for always 2, if I extract 2 for each cell, then row sum for the sample should be 0
  
  ores1 = orpsres - 1
  ores2 = orpsres - 2
  
  ogrp1 = rownames(ores1[which(rowSums(ores1) == 0),])
  ogrp2 = rownames(ores1[which(rowSums(ores2) == 0),])
  
  grp1 = intersect(grp1, ogrp2)
  grp2 = intersect(grp2, ogrp1)
  
  #########################################################
  
  ### now, I need to work on group mean and sd for each feature
  datgrp1 = newdat[,grp1]
  datgrp2 = newdat[,grp2]
  
  means1 = rowMeans(datgrp1) 
  means2 = rowMeans(datgrp2)
  
  sds1 = apply(datgrp1,1,sd)
  sds2 = apply(datgrp2,1,sd)
  
  #### I need to determine the direction for allmeans, and allsds
  #### maybe I can use the middle prior to make the decision
  xx=median(rps)
  tmp = PRPSSLwithWeightsPrior(newdat=newdat, weights=weights, ratioPrior = xx, PRPShighGroup = PRPShighGroup, PRPSlowGroup = PRPSlowGroup)
  mcls = mclust::Mclust(tmp$PRPS_test$PRPS_score, G=2, warn = TRUE)
  
  mmean = mcls$parameters$mean
  if(mmean[2] > mmean[1]){
    allmeans = cbind(means2, means1)
    allsds = cbind(sds2, sds1)
    tt = grp1
    grp1 = grp2
    grp2 = tt
  }else{
    allmeans = cbind(means1, means2)
    allsds = cbind(sds1, sds2)
  }
  
  df1 = length(grp1) - 1
  df2 = length(grp2) - 1
  dfs = c(df1, df2)
  
  PRPS_score = weightedLogProbClass(newdat = newdat, topTraits=names(weights), weights=weights,
                                    classMeans = allmeans, classSds = allsds, dfs = dfs)
  
  ###################################################################
  # #### 20190503, call plotHistEM , make changes of scoreName on 20191206 for this local function
  emsearch = plotHistEM(PRPS_score, G = 2, breaks = breaks, scoreName = "Final PRPS_score")
  ####################################################################
  
  bestG = emsearch$bestG
  emcut = emsearch$emcut
  #
  # ### no matter how many bestG, only keep the 1st and last one for the following
  gmeans = c(emcut$Means[1], emcut$Means[bestG])
  gsds = c(emcut$SDs[1], emcut$SDs[bestG])
  #
  if(gmeans[1]> gmeans[2]){
    name1 = PRPShighGroup
    name2 = PRPSlowGroup
  }else{
    name1 = PRPSlowGroup
    name2 = PRPShighGroup
  }
  
  grp1score = PRPS_score[grp1]
  grp2score = PRPS_score[grp2]
  #
  scoreMeans = c(mean(grp1score), mean(grp2score))
  scoreSds = c(sd(grp1score), sd(grp2score))
  #
  PRPS_prob1 = getProb(PRPS_score, groupMeans = scoreMeans, groupSds = scoreSds)
  PRPS_prob2 = getProb(PRPS_score, groupMeans = rev(scoreMeans), groupSds = rev(scoreSds))
  
  PRPS_class = rep("UNCLASS",length(PRPS_score))
  PRPS_class[which(PRPS_prob1 >= classProbCut)] = PRPShighGroup
  PRPS_class[which(PRPS_prob2 >= classProbCut)] = PRPSlowGroup
  
  PRPS_class0 = ifelse(PRPS_score > 0,  PRPShighGroup, PRPSlowGroup)
  PRPS_score = data.frame(PRPS_score)
  PRPS_test = cbind(PRPS_score, PRPS_class, PRPS_prob1, PRPS_prob2, PRPS_class0,stringsAsFactors =F)
  
  
  ### 20190503: add the stable classification 
  PRPS_test$stable_class = "UNCLASS"
  PRPS_test[grp1,"stable_class"] = PRPShighGroup
  PRPS_test[grp2,"stable_class"] = PRPSlowGroup
  
  
  weights = data.frame(weights)
  scoreMeanSds = c(scoreMeans, scoreSds)
  
  PRPS_pars =  list(weights, meansds = scoreMeanSds, traitsmeansds = cbind(allmeans, allsds), dfs)
  names(PRPS_pars) = c("weights","meansds","traitsmeansds", "dfs")
  
  if(!is.null(plotName)){
    dev.off()
  }
  
  outs = list(PRPS_pars, PRPS_test)
  names(outs) = c("PRPS_pars","PRPS_test")
  return(outs)
  
}
