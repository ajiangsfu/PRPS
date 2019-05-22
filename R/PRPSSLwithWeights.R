
#' Self learning binary classification with selected features and their weights
#' @description This function is to calculate PRPS (Probability ratio based classification predication score) scores and 
#' make binary classificationcalls for a testing data set without PRPS training object. This function involves a self learning 
#' process with weights + priors + EM + Bayes but no need to input a group ratio prior. However, we do need selected feature list
#' with the feature weights, and we use self learning method to estimate mean and sd for two groups for each selected feature 
#' in order to calculate PRPS scores. Our original idea in JCO 2018 does requrire a group ratio prior to do so,
#' which is written as PRPSSLwithWeightsPrior function in this package. We also have another function: PRPSSLwithWeightsEM, 
#' which uses EM to estimate group mean and sd for each selected feature. In this fucnction, however, we combine the ideas in both
#' of the above two functions to achieve our new version of self learning algorithm.
#' @details This function is trying to get reasonable PRPS based classification without training data set, but with 
#' selected features and their weights. It combines ideas from PRPSSLwithWeightsPrior and
#' PRPSSLwithWeightsEM, which we call it as self learning algorithm, the actual steps are as following:
#' 1) assume that we have a pool for group ratio priors such as: seq(0.05, 0.95, by = 0.05), this will give us 19 ratio priors
#' 2) for each prior in 1), call PRPSSLwithWeightsPrior to achieve PRPS scores
#' 3) apply EM on PRPS scores from 2) with Mclust, which includes 2 group classification
#' 4) use the samples that are always in the same groups to get group means and sds for each feature and for each group
#' 5) calculate PRPS scores
#' 6) Once we have PRPS scores, we could use the theoretic natual cutoff 0 to make classification calls, which may or may not appropriate. 
#' Alternatively, we can also apply EM to calcualate mean and sd for the
#' two groups assuming that PRPS score is a mixture of two normal distributions, followed by Empirical Bayes' probability 
#' calculation and final binary classification calls.
#' PRPS calculation is based on Ennishi 2018, its formula is:
#' \eqn{PRPS(X_i) = \sum (|a_j| log(P1(x_ij)/P0(x_ij)))}
#' Here, a_j represents the jth selected feature weights, and x_ij is the corresponding feature value
#'  for the ith sample, 
#' P1 and P0 are the probabilities that the ith sample belongs to two different group.
#' However, in order to calculate P1 and P0, we need to have two group mean and sd for each selected feature. Although there are multiple way
#' to obtain these values, in this function, we design to use EM algorithm to achieve group mean and sd assuming that each selected feature
#' is a mixture of two normal distributions. 
#' After we have PRPS scores, we provide classification calls based on theoretic 0 cutoff and based on probability, which may or may not appropriate
#' for this function since it is based on self learning. To calculate a Empirical Bayes' 
#' probability to make classification calls, we also need to apply EM to get PRPS score mean and sd for two groups. 
#' After that, we can calcualte probability that a sample belongs to either group,
#' and then use the following formula to get Empirical Bayes' probability:
#' \eqn{prob(x) = p_test(x)/(p_test(x) + p_ref(x))}
#' Here prob(x) is the Empirical Bayes' probability of a given sample, p_test(x) is the probability that a given sample
#' belongs to the test group, p_ref(x) is the probability that a given sample belongs to the reference group.
#' Notice that the test and reference group is just the relative grouping, in fact, for this step, we often need
#'  to calculate Empirical Bayes' probabilities for a given sample from two different standing points.
#' @param newdat a input data matrix or data frame, columns for samples and rows for features
#' @param weights a numeric vector with selected features (as names of the vector) and their weights
#' @param standardization a logic variable to indicate if standardization is needed before classification 
#'  score calculation
#' @param classProbCut a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
#'  often used values are 0.8 and 0.9, default value is 0.9. Only one value is used for both groups, 
#'  the samples that are not included in either group will be assigned as UNCLASS
#' @param PRPShighGroup a string to indicate group name with high PRPS score
#' @param PRPSlowGroup a string to indicate group name with low PRPS score
#' @param breaks a integer to indicate number of bins in histogram, default is 50
#' @param EMmaxRuns number of Iterations for EM searching; default=50
#' @param imputeNA a logic variable to indicate if NA imputation is needed, if it is TRUE, NA imputation is 
#'  processed before any other steps, the default is FALSE
#' @param byrow a logic variable to indicate direction for imputation, default is TRUE, 
#'  which will use the row data for imputation
#' @param imputeValue a character variable to indicate which value to be used to replace NA, default is "median", 
#'  the median value of the chose direction with "byrow" data to be used
#' @return A list with two items is returned: PRPS parameters for selected features, PRPS scores and classifications for the given samples.
#' \item{PRPS_pars}{a list of 3 items, the 1st item is a data frame with weights of each selected features for PRPS
#'  calculation, the 2nd item is a numeric vector containing PRPS mean and sd for two groups，and the 3rd item is a data frame contains mean and sd
#'   for each group and for each selected feature}
#' \item{PRPS_test}{a data frame of PRPS score,  classification and two groups' Empirical Bayesian probabilites based on EM, and
#' classification with natural 0 cutoff}
#' @keywords PRPS EM 
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

#' @export

PRPSSLwithWeights = function(newdat, weights, standardization=FALSE, classProbCut = 0.9, PRPShighGroup = "PRPShigh", 
                    PRPSlowGroup = "PRPSlow", breaks = 50, EMmaxRuns = 50, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  require(mclust)
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
  
  rps = seq(0.05, 0.95, by = 0.05)
  
  rpsres = sapply(rps, FUN = function(xx){
    tmp = PRPSSLwithWeightsPrior(newdat=newdat, weights=weights, ratioPrior = xx, PRPShighGroup = PRPShighGroup, PRPSlowGroup = PRPSlowGroup)
    mcls = mclust::Mclust(tmp$PRPS_test$PRPS_score, G=2)
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
  
  ### now, I need to work on group mean and sd for each feature
  datgrp1 = newdat[,grp1]
  datgrp2 = newdat[,grp2]
  
  means1 = rowMeans(datgrp1) 
  means2 = rowMeans(datgrp2)
  
  sds1 = apply(datgrp1,1,sd)
  sds2 = apply(datgrp2,1,sd)
  
  # if(length(grp1) >= length(grp2)){
  #   allmeans = cbind(means2, means1)
  #   allsds = cbind(sds2, sds1)
  # }else{
  #   allmeans = cbind(means1, means2)
  #   allsds = cbind(sds1, sds2)
  # }
  
  #### I need to determine the direction for allmeans, and allsds
  #### maybe I can use the middle prior to make the decision
  xx=median(rps)
  tmp = PRPSSLwithWeightsPrior(newdat=newdat, weights=weights, ratioPrior = xx, PRPShighGroup = PRPShighGroup, PRPSlowGroup = PRPSlowGroup)
  mcls = mclust::Mclust(tmp$PRPS_test$PRPS_score, G=2)
  #### in fact, there are group labels in mcls$mean
  #### get the two groups score means and decide which group is which
  
  mmean = mcls$parameters$mean
  if(mmean[2] > mmean[1]){
    allmeans = cbind(means2, means1)
    allsds = cbind(sds2, sds1)
  }else{
    allmeans = cbind(means1, means2)
    allsds = cbind(sds1, sds2)
  }
  
  PRPS_score = weightedLogProbClass(newdat = newdat, topTraits=names(weights), weights=weights,
                                    classMeans = allmeans, classSds = allsds)
  
  #### 20190503, call plotHistEM 
  emsearch = plotHistEM(PRPS_score, G = 2:4, breaks = breaks, EMmaxRuns = EMmaxRuns, scoreName = "PRPS_score")
  bestG = emsearch$bestG
  emcut = emsearch$emcut
  
  ### no matter how many bestG, only keep the 1st and last one for the following
  gmeans = c(emcut$Means[1], emcut$Means[bestG])
  gsds = c(emcut$SDs[1], emcut$SDs[bestG])
  
  PRPS_prob1 = getProb(PRPS_score, groupMeans = gmeans, groupSds = gsds)
  PRPS_prob2 = getProb(PRPS_score, groupMeans = rev(gmeans), groupSds = rev(gsds))
  
  if(gmeans[1]> gmeans[2]){
    name1 = PRPShighGroup
    name2 = PRPSlowGroup
    scoreMeanSds = c(gmeans, gsds)
  }else{
    name1 = PRPSlowGroup
    name2 = PRPShighGroup
    scoreMeanSds = c(rev(gmeans), rev(gsds))
  }
  
  PRPS_class = rep("UNCLASS",length(PRPS_score))
  PRPS_class[which(PRPS_prob1 >= classProbCut)] = name1
  PRPS_class[which(PRPS_prob2 >= classProbCut)] = name2
  
  names(scoreMeanSds) = c("testPRPSmean","refPRPSmean","testPRPSsd","refPRPSsd")
  
  PRPS_class0 = ifelse(PRPS_score > 0,  PRPShighGroup, PRPSlowGroup)
  
  PRPS_score = data.frame(PRPS_score)
  PRPS_test = cbind(PRPS_score, PRPS_class, PRPS_prob1, PRPS_prob2, PRPS_class0,stringsAsFactors =F)

  weights = data.frame(weights)
  
  PRPS_pars =  list(weights, meansds = scoreMeanSds, traitsmeansds = cbind(allmeans, allsds))
  names(PRPS_pars) = c("weights","meansds","traitsmeansds")
  
  outs = list(PRPS_pars, PRPS_test)
  names(outs) = c("PRPS_pars","PRPS_test")
  return(outs)
  
}
