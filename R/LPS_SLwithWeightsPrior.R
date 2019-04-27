
#' LPS score calculation and binary classification for a testing data set without training data set, but with selected features' weights and 
#' a ratio prior 
#' @description This is the function to calculate LPS (Linear Prediction Score) scores and make binary classification for a data set 
#' without LPS training set. However, we do need some prior information such as selected feature list with their weights, 
#' and ratio prior to calculate LPS score mean and sd for each group for Empirical Bayes' probability calculation.
#' @details  This is the function to calculate LPS scores and make binary classification based on given data set with selected features and 
#' their weights from elsewhere, and a ratio prior of one group over all samples. The given data and selected feature info are 
#' used for LPS score calculation, while the ratio prior is used to calculate Empirical Bayesian probabilities, 
#' followed by final classification calls.
#' 
#' Within this function, standardization step is included as an option to minimize the difference 
#' between training and testing data sets, be aware that standardization step should be consistent between the data used 
#' for weight calculation and the testing data set. 
#' Also, be aware that this step is only done to make distributions of each selected
#' features comparable. Be aware that this feature-wise standardization cannot make the sample-wise distributions
#' comparable. For example, if the proportion of one group over all samples might be much lessin the testing data set 
#' compared to the data that were used for weight generation, or even worse, the testing data set might only contain one classification 
#' group only, in this case, feature-wise standardization can not solve the comparison problem. 
#' This is the common problem for classification and feature-wise standardization cannot solve the problem. 
#' In order to solve the problem, we should make data comparable as much as possbile before classification step.

#' LPS calculation is based on Wright 2003. The fomula is straightforward:
#' \eqn{LPS(X) = \sum a_j x_ij}
#' Here a_j represents the jth selected feature weights, and x_ij is the corresponding feature value for the ith sample.

#' When calculate a Empirical Bayes' probability, the 1st group in the input mean and sd vectors is treated as the 
#' test group. When calculate the probabilities, we first calcualte probability that a sample belongs to either group,
#' and then use the following formula to get Empirical Bayes' probability:
#' \eqn{prob(x) = p_test(x)/(p_test(x) + p_ref(x))}
#' Here prob(x) is the Empirical Bayes' probability of a given sample, p_test(x) is the probability that a given sample
#' belongs to the test group, p_ref(x) is the probability that a given sample belongs to the reference group.
#' Notice that the test and reference group is just the relative grouping, in fact, for this step, we often need
#' to calculate Empirical Bayes' probabilities for a given sample from two different standing points.
#' @param newdat a new data matrix or data frame, which is comparable to training data set, 
#'  with columns for samples and rows for features
#' @param weights a numeric vector with selected features (as names of the vector) and their weights
#' @param ratioPrior A numeric vector with 2 items, ratio prior for test group, and ratio prior for reference group. The order is important
#' @param testGroup A string for test group name
#' @param refGroup A string for reference group name
#' @param isTestGroupHighLPS A logic variable to indicate if the test group has higher LPS score than the reference group
#' @param standardization a logic variable to indicate if standardization is needed before classification 
#'  score calculation
#' @param classProbCut a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
#'  often used values are 0.8 and 0.9, default value is 0.8. Only one value is used for both groups, 
#'  the samples that are not included in either group will be assigned as UNCLASS
#' @param imputeNA a logic variable to indicate if NA imputation is needed, if it is TRUE, NA imputation is 
#'  processed before any other steps, the default is FALSE
#' @param byrow a logic variable to indicate direction for imputation, default is TRUE, 
#'  which will use the row data for imputation
#' @param imputeValue a character variable to indicate which value to be used to replace NA, default is "median", 
#'  the median value of the chose direction with "byrow" data to be used
#' @return A list with two items is returned: LPS parameters for selected features, LPS scores and classifications for the given samples.
#' \item{LPS_pars}{a list of 2 items, the 1st item is a data frame with weights of each selected features for LPS
#'  calculation, the 2nd item is a numeric vector containing LPS mean and sd for two groups based on a given ratio prior}
#' \item{LPS_test}{a data frame of LPS score, and Empirical Bayesian probabilites for both groups and classification}
#' @keywords LPS
#' @author Aixiang Jiang
#' @references Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A gene expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.
#' @export

LPS_SLwithWeightsPrior = function(newdat, weights, ratioPrior = 1/3, testGroup, refGroup, isTestGroupHighLPS = TRUE,
                      standardization=FALSE, classProbCut = 0.8,  imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  imputeValue = imputeValue[1]
  ## imputee NA if imputeNA is true
  if(imputeNA){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for LPS approach, it does not require standardization, however, if standardization = TRUE, do the standardization
  if(standardization){newdat = standardize(newdat)}
  
  # if NA is not imputed, remove it from LPS score calculation
  LPS_score = apply(data.matrix(newdat[names(weights),]), 2, getLPSscore, coefs = weights)
  
  # for LPS, 0 is a NOT natural cutoff for two group classification
  # in order to get classification, need to get two groups' LPS mean and sd
  
  
  
  
  ##########################################################
  #### notice that I am using two ratio priors in this function
  #### should I change it to use one prior only? or should I change all other functions to be two priors as well?
  #### => think again, 1 ratio prior makes more sense, however, 
  #### it seems that it is more reasonable to use the two ends of the data to calculate mean and sd for feature level of LPS score level
  ####    => I should change all SL functions? but at least I should keep 1 cut instead of two cuts(to get two ends) to 
  ####        consistent to JCO paper for one PRPS function?
  
  #### in addition, here I use test group, reference group, and if test group score is high?
  ####  which is not consistent to PRPS SL functions as well, should I change all to setting as in here, or change all to be setting  as in PRPS SL?
  #### => 
  ###########################################################
  
  if(isTestGroupHighLPS == TRUE | isTestGroupHighLPS == T){
    ttmp = quantile(LPS_score, probs = 1-ratioPrior)
    rtmp = quantile(LPS_score, probs = ratioPrior)
    ttmp = LPS_score[which(LPS_score >= ttmp)]
    rtmp = LPS_score[which(LPS_score < rtmp)]
  } else{
    ttmp = quantile(LPS_score, probs = ratioPrior)
    rtmp = quantile(LPS_score, probs = 1- ratioPrior)
    ttmp = LPS_score[which(LPS_score <= ttmp)]
    rtmp = LPS_score[which(LPS_score > rtmp)]
  }
  ### notice that both ttmp and rtmp are the two ends, meaning that length(rtmp + ttmp) < length(all scores), it will be equal if ratioRrior is 0.5
  
  # ### make changes on 20190425 after comparisons on 20190424
  # ### the following 3 lines are actually consistent to the feature level setting for ratio prior usage
  # ###  more importantly, the following 3 lines are also consistent to JCO paper where the original PRPS was defined and used with ratio prior
  # tmpcut = quantile(PRPS_score, probs = 1-ratioPrior)
  # ttmp = PRPS_score[which(PRPS_score >= tmpcut)]
  # rtmp = PRPS_score[which(PRPS_score < tmpcut)]
  # 
  
  testLPSmean = mean(ttmp)
  refLPSmean = mean(rtmp)
  testLPSsd = sd(ttmp)
  refLPSsd = sd(rtmp)
  
  LPS_prob_test = getProb(LPS_score, groupMeans = c(testLPSmean, refLPSmean), groupSds = c(testLPSsd, refLPSsd))
  
  LPS_prob_ref = getProb(LPS_score, groupMeans = c(refLPSmean, testLPSmean), groupSds = c(refLPSsd, testLPSsd))
  
  LPS_class = rep("UNCLASS",length(LPS_score))
  LPS_class[which(LPS_prob_test >= classProbCut)] = testGroup
  LPS_class[which(LPS_prob_ref >= classProbCut)] = refGroup
  
  LPS_score = data.frame(LPS_score)
  LPS_test = cbind(LPS_score, LPS_class, LPS_prob_test, LPS_prob_ref, stringsAsFactors =F)
  
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
