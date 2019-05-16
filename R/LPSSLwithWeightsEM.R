#' LPS score calculation and binary classification for a testing data set without training data set, but with selected feature weights
#' @description This is the function to calculate LPS (Linear Prediction Score) scores for a testing data set 
#' without training data set. However, we do need selected feature list with their weights. Once we have LPS scores, 
#' we then apply EM (Expectation-Maximization) to get LPS score mean and sd for the two groups assuming that LPS score
#' is a mixture of two normal distributions, followed by Empirical Bayes' probability calculation and final binary classification calls.
#' @details  This is the function to calculate LPS scores with given selected feature weights, use EM (Expectation-Maximization) to get LPS score
#' mean and sd for the two groups assuming that LPS score is a mixture of two normal distributions, and then make classification based on 
#' Empirical Bayesian probabilities for a testing new data set. 
#' LPS calculation is based on Wright 2003. The fomula is straightforward for a given testing data and selected feature weights:
#' \eqn{LPS(X) = \sum a_j x_ij}
#' Here a_j represents the jth selected feature weight, and x_ij is the corresponding feature value for the ith sample.
#' 
#' However, when there is no training data set, we do not have LPS mean and sd for the two groups. There are multiple ways to deal with this problem.
#' In this function, we design to use EM algorithm to get LPS score mean and sd calculation for the two groups assuming that LPS score
#' is a mixture of two normal distributions.
#' 
#' When calculate a Empirical Bayes' probability, we first calcualte probability that a sample belongs to either group,
#' and then use the following formula to get Empirical Bayes' probability:
#' \eqn{prob(x) = p_test(x)/(p_test(x) + p_ref(x))}
#' Here prob(x) is the Empirical Bayes' probability of a given sample, p_test(x) is the probability that a given sample
#' belongs to the test group, p_ref(x) is the probability that a given sample belongs to the reference group.
#' Notice that the test and reference group is just the relative grouping, in fact, for this step, we often need
#' to calculate Empirical Bayes' probabilities for a given sample from two different standing points.
#' 
#' Within this function, standardization step is included as an option to minimize the difference 
#' between training and testing data sets, be aware that standardization step should be consistent between the data used 
#' for weight calculation and the testing data set. 
#' Also, be aware that this standardization step is only done to make distributions of each selected
#' features comparable. Be aware that this feature-wise standardization cannot make the sample-wise distributions
#' comparable. For example, if the proportion of one group over all samples might be much lessin the testing data set 
#' compared to the data that were used for weight generation, or even worse, the testing data set might only contain one classification 
#' group only, in this case, feature-wise standardization can not solve the comparison problem. 
#' This is the common problem for classification and feature-wise standardization cannot solve the problem. 
#' In order to solve the problem, we should make data comparable as much as possbile before classification step.
#' 
#' @param newdat a new data matrix or data frame, columns for samples and rows for features
#' @param weights a numeric vector with selected features (as names of the vector) and their weights
#' @param standardization a logic variable to indicate if standardization is needed before classification 
#'  score calculation
#' @param classProbCut a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
#'  often used values are 0.8 and 0.9, default value is 0.8. Only one value is used for both groups, 
#'  the samples that are not included in either group will be assigned as UNCLASS
#' @param LPShighGroup a string to indicate group name with high LPS score
#' @param LPSlowGroup a string to indicate group name with low LPS score
#' @param breaks a integer to indicate number of bins in histogram, default is 50
#' @param EMmaxRuns number of Iterations for EM searching; default=50
#' @param imputeNA a logic variable to indicate if NA imputation is needed, if it is TRUE, NA imputation is 
#'  processed before any other steps, the default is FALSE
#' @param byrow a logic variable to indicate direction for imputation, default is TRUE, 
#'  which will use the row data for imputation
#' @param imputeValue a character variable to indicate which value to be used to replace NA, default is "median", 
#'  the median value of the chose direction with "byrow" data to be used
#' @return A data frame with LPS score, Empirical Bayesian probabilites for two groups and classification
#' 
#' #' @return A list with two items is returned: LPS parameters for selected features, LPS scores and classifications for the given samples.
#' \item{LPS_pars}{a list of 2 items, the 1st item is a data frame with weights of each selected features for LPS
#'  calculation, the 2nd item is a numeric vector containing LPS mean and sd for two groups based on EM}
#' \item{LPS_test}{a data frame of LPS score, Empirical Bayesian probabilites for both groups and classification}
#' @keywords LPS EM 
#' @author Aixiang Jiang
#' @references Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A gene expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.
#' 
#' Ultsch, A., Thrun, M.C., Hansen-Goos, O., Loetsch, J.: Identification of Molecular Fingerprints
#' in Human Heat Pain Thresholds by Use of an Interactive Mixture Model R Toolbox(AdaptGauss),
#' International Journal of Molecular Sciences, doi:10.3390/ijms161025897, 2015.
#' @export

LPSSLWithWeightsEM = function(newdat, weights, standardization=FALSE, classProbCut = 0.8, LPShighGroup = "LPShigh", 
                    LPSlowGroup = "LPSlow", breaks = 50, EMmaxRuns = 50, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  imputeValue = imputeValue[1]
  ## imputee NA if imputeNA is true
  if(imputeNA){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for LPS approach, it does not require standardization, however, if standardization = TRUE, do the standardization
  if(standardization){newdat = standardize(newdat)}
  
  ### in case that some features in weight list are actually not in data, do the following
  tmp = intersect(names(weights), rownames(newdat))
  weights = weights[tmp]

  # if NA is not imputed, remove it from LPS score calculation in the following getLPSscore function
  LPS_score = apply(data.matrix(newdat[names(weights),]), 2, getLPSscore, coefs = weights)
  
  ### now, use EM to define 1st draft of group, and then use it to calculate group mean and sd
  
  emcut = AdaptGauss::EMGauss(LPS_score, K = 2, fast=TRUE, MaxNumberofIterations = EMmaxRuns)
  
  # > emcut
  # $Means
  # [1] -106.36599  -39.04443
  # 
  # $SDs
  # [1] 23.62815 11.70491
  # 
  # $Weights
  # [1] 0.1139536 0.8860464
  #### this is exact what I need!
  
  ### add a plot, hist with two distribution lines, do not need to save, just plot it
  hist(LPS_score, prob = TRUE, breaks = breaks)
  curve(emcut$Weights[1]*dnorm(x, mean=emcut$Means[1], sd=emcut$SDs[1]), 
        col="red", lwd=2, add=TRUE, yaxt="n")
  curve(emcut$Weights[2]*dnorm(x, mean=emcut$Means[2], sd=emcut$SDs[2]), 
        col="green", lwd=2, add=TRUE, yaxt="n")
  
  LPS_prob1 = getProb(LPS_score, groupMeans = emcut$Means, groupSds = emcut$SDs)
  LPS_prob2 = getProb(LPS_score, groupMeans = rev(emcut$Means), groupSds = rev(emcut$SDs))
  
  if(emcut$Means[1]> emcut$Means[2]){
    name1 = LPShighGroup
    name2 = LPSlowGroup
  }else{
    name1 = LPSlowGroup
    name2 = LPShighGroup
  }
  
  LPS_class = rep("UNCLASS",length(LPS_score))
  LPS_class[which(LPS_prob1 >= classProbCut)] = name1
  LPS_class[which(LPS_prob2 >= classProbCut)] = name2
  
  LPS_score = data.frame(LPS_score)
  LPS_test = cbind(LPS_score, LPS_class, LPS_prob1, LPS_prob2, stringsAsFactors =F)
  
  #####################################################
  #### add more item to return
  
  weights = data.frame(weights)
  
  meansds = c(testLPSmean, refLPSmean, testLPSsd, refLPSsd)
  names(meansds) = c("testLPSmean","refLPSmean","testLPSsd","refLPSsd")
  
  LPS_pars =  list(weights,meansds)
  names(LPS_pars) = c("weights","meansds")
  
  outs = list(LPS_pars, LPS_test)
  names(outs) = c("LPS_pars","LPS_test")
  
  return(outs)
  
}
