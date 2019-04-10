
#' PRPS score calculation and binary classification for a testing data set without PRPStraining output object, but with selected feature weights,
#' The current function uses weights + priors + EM + Bayes but no need to input a group ratio prior.
#' @description This function is to calculate PRPS (Probability ratio based classification predication score) scores for a testing data set 
#' without PRPS training object. However, we do need selected feature list with their weights, and we can estimate mean and sd for two groups for each
#' selected feature in order to calculate PRPS scores. Our original idea in JCO 2018 does requrire a group ratio prior to do so, which is 
#' written as PRPSClassWithWeightsPrior function in this package. We also have another function: PRPSClassWithWeightsEM, which uses EM to 
#' estimate group mean and sd for each selected feature. In this fucnction, however, we combine the ideas in both of the above two functions.
#' 
#' @details This function is trying to get reasonable PRPS based classification without training object and given group mean and sd for each selected feature,
#' but with selected features and their weights. It combined ideas from PRPSClassWithWeightsPrior and PRPSClassWithWeightsEM. 
#' The actual steps are as following:
#' 1) assume that we have a pool for group ratio priors such as: seq(0.05, 0.95, by = 0.05), this will give us 19 ratio priors
#' 2) for each prior in 1), call PRPSClassWithWeightsPrior, and we will get PRPS score and classification with natural cutoff 0
#' 3) apply EM on PRPS scores from 2) with Mclust
#' 4) use the samples that are always in the same group to get group means and sds for each feature
#' 5) calculate PRPS scores
#' 6) Once we have PRPS scores, we could use the theoretic natual cutoff 0 to make classification calls. 
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
#' After we have PRPS scores, we need to calculate a Empirical Bayes' probability to make classification calls. To do that, 
#' we also need to apply EM to get PRPS score mean and sd for two groups. 
#' After that, we can calcualte probability that a sample belongs to either group,
#' and then use the following formula to get Empirical Bayes' probability:
#' \eqn{prob(x) = p_test(x)/(p_test(x) + p_ref(x))}
#' Here prob(x) is the Empirical Bayes' probability of a given sample, p_test(x) is the probability that a given sample
#' belongs to the test group, p_ref(x) is the probability that a given sample belongs to the reference group.
#' Notice that the test and reference group is just the relative grouping, in fact, for this step, we often need
#'  to calculate Empirical Bayes' probabilities for a given sample from two different standing points.
#' @param newdat a new data matrix or data frame, columns for samples and rows for features
#' @param weights a numeric vector with selected features (as names of the vector) and their weights
#' @param standardization a logic variable to indicate if standardization is needed before classification 
#'  score calculation
#' @param classProbCut a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
#'  often used values are 0.8 and 0.9, default value is 0.8. Only one value is used for both groups, 
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
#' @return A data frame with PRPS score, Empirical Bayesian probabilites for two groups and classification
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
#' #' Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A trait expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.

#' @export

PRPSClassNoTraining = function(newdat, weights, standardization=FALSE, classProbCut = 0.8, PRPShighGroup = "PRPShigh", 
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
    tmp = PRPSClassWithWeightsPrior(newdat=newdat, weights=weights, ratioPrior = xx, PRPShighGroup = PRPShighGroup, PRPSlowGroup = PRPSlowGroup)
    mcls = mclust::Mclust(tmp$PRPS_score, G=2)
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
  
  if(length(grp1) >= length(grp2)){
    allmeans = cbind(means2, means1)
    allsds = cbind(sds2, sds1)
  }else{
    allmeans = cbind(means1, means2)
    allsds = cbind(sds1, sds2)
  }
  
  PRPS_score = weightedLogProbClass(newdat = newdat, topTraits=names(weights), weights=weights,
                                    classMeans = allmeans, classSds = allsds)
  
  emcut = AdaptGauss::EMGauss(PRPS_score, K = 2,fast=TRUE, MaxNumberofIterations = EMmaxRuns)
 
  ### add a plot, hist with two distribution lines, do not need to save, just plot it
  hist(PRPS_score, prob = TRUE, breaks = breaks)
  curve(emcut$Weights[1]*dnorm(x, mean=emcut$Means[1], sd=emcut$SDs[1]), 
        col="red", lwd=2, add=TRUE, yaxt="n")
  curve(emcut$Weights[2]*dnorm(x, mean=emcut$Means[2], sd=emcut$SDs[2]), 
        col="green", lwd=2, add=TRUE, yaxt="n")
  
  PRPS_prob1 = getProb(PRPS_score, groupMeans = emcut$Means, groupSds = emcut$SDs)
  PRPS_prob2 = getProb(PRPS_score, groupMeans = rev(emcut$Means), groupSds = rev(emcut$SDs))
  
  if(emcut$Means[1]> emcut$Means[2]){
    name1 = PRPShighGroup
    name2 = PRPSlowGroup
  }else{
    name1 = PRPSlowGroup
    name2 = PRPShighGroup
  }
  
  PRPS_class = rep("UNCLASS",length(PRPS_score))
  PRPS_class[which(PRPS_prob1 >= classProbCut)] = name1
  PRPS_class[which(PRPS_prob2 >= classProbCut)] = name2
  
  PRPS_score = data.frame(PRPS_score)
  PRPS_test = cbind(PRPS_score, PRPS_class, PRPS_prob1, PRPS_prob2, stringsAsFactors =F)
  
  return(PRPS_test)
  
}
