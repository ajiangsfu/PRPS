
#' LPS score calculation for a testing data set without LPStraining output object, but with some priors
#' @description This is the function to calculate LPS (Linear Prediction Score) scores for a testing data set 
#' without LPS training object. However, we do need some prior information such as selected feature list with their weights, 
#' LPS score means and sds for two groups or priors for two groups' percentages.
#' @details  This is the function to calculate LPS scores and make classification and 
#' Empirical Bayesian probabilities for a testing new data set. However, this new data set should be 
#' comparable to the training data set that generates prior parameters as much as possible. When the testing data set is not comparable
#' to the training data set, calibration might be needed before calling this function. 
#' Within this LPStesting function, standardization step is included as an option to minimize the difference 
#' between training and testing data sets, be aware that standardization step should be consistent between training and testing data sets. 
#' Also, be aware that this step is only done to make distributions of each selected
#' features comparable. Be aware that this feature-wise standardization cannot make the sample-wise distributions
#' comparable. For example, the training data set must have two classification groups, however, 
#' the proportion of one group sample might be much less than the other group in the testing data set 
#' compared to the training data set, or even worse, the testing data set might only contain one classification 
#' group only. This is the common problem for classification and feature-wise standardization cannot solve the problem. 
#' In order to solve the problem, we should make data comparable as much as possbile before classification step, for example, calibration as mentioned above. 

#' LPS calculation is based on Wright 2003. The fomula is straightforward:
#'   \eqn{LPS(X) = \sum a_j x_ij}
#' Here a_j represents the jth selected feature weights, and x_ij is the corresponding feature value for the ith sample.
#' \eqn{PRPS(X_i) = \sum (|a_j| log(P1(x_ij)/P0(x_ij)))}
#' When calculate a Empirical Bayes' probability, the 1st group in the input mean and sd vectors is treated as the 
#' test group. When calculate the probabilities, we first calcualte probability that a sample belongs to either group,
#' and then use the following formula to get Empirical Bayes' probability:
#' \eqn{prob(x) = p_test(x)/(p_test(x) + p_ref(x))}
#' Here prob(x) is the Empirical Bayes' probability of a given sample, p_test(x) is the probability that a given sample
#' belongs to the test group, p_ref(x) is the probability that a given sample belongs to the reference group.
#' Notice that the test and reference group is just the relative grouping, in fact, for this step, we often need
#'  to calculate Empirical Bayes' probabilities for a given sample from two different standing points.
#' @param newdat a new data matrix or data frame, which is comparable to training data set, 
#'  with columns for samples and rows for features
#' @param weights a numeric vector with selected features (as names of the vector) and their weights
#' @param testGroup A string for test group name
#' @param refGroup A string for reference group name
#' @param LPSMeanSds A numeric vector with 4 items for LPS score distributions: test group mean, reference group mean, test group sd, reference group sd
#'   the names for these items are not necessary, but the order of these items are important
#' @param ratioPriors A numeric vector with 2 items, ratio prior for test group, and ratio prior for reference group. The order is important
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
#' @return A data frame with LPS score, Empirical Bayesian probabilites for two groups and classification
#' @keywords LPS
#' @author Aixiang Jiang
#' @references Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A gene expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.
#' @export

LPSClassWithPriors = function(newdat, weights, testGroup, refGroup, LPSMeanSds = NULL, ratioPriors = c(0.5, 0.5), isTestGroupHighLPS = TRUE,
                      standardization=FALSE, classProbCut = 0.8,  imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  imputeValue = imputeValue[1]
  ## imputee NA if imputeNA is true
  if(imputeNA == TRUE | imputeNA == T){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for LPS approach, it does not require standardization, however, if standardization = TRUE, do the standardization
  if(standardization){newdat = standardize(newdat)}
  
  # if NA is not imputed, remove it from LPS score calculation
  LPS_score = apply(data.matrix(newdat[names(weights),]), 2, getLPSscore, coefs = weights)
  
  # for LPS, 0 is a NOT natural cutoff for two group classification
  # in order to get classification, need to get two groups' LPS mean and sd
  
  if(is.null(LPSMeanSds)){
    if(isTestGroupHighLPS == TRUE | isTestGroupHighLPS == T){
      ttmp = quantile(LPS_score, probs = 1-ratioPriors[1])
      rtmp = quantile(LPS_score, probs = ratioPriors[2])
      ttmp = LPS_score[which(LPS_score >= ttmp)]
      rtmp = LPS_score[which(LPS_score < rtmp)]
    } else{
      ttmp = quantile(LPS_score, probs = ratioPriors[1])
      rtmp = quantile(LPS_score, probs = 1- ratioPriors[2])
      ttmp = LPS_score[which(LPS_score <= ttmp)]
      rtmp = LPS_score[which(LPS_score > rtmp)]
    }
    testLPSmean = mean(ttmp)
    refLPSmean = mean(rtmp)
    testLPSsd = sd(ttmp)
    refLPSsd = sd(rtmp)
  }else{
    testLPSmean = LPSMeanSds[1]
    refLPSmean = LPSMeanSds[2]
    testLPSsd = LPSMeanSds[3]
    refLPSsd = LPSMeanSds[4]
  }
  
  LPS_prob_test = getProb(LPS_score, groupMeans = c(testLPSmean, refLPSmean), groupSds = c(testLPSsd, refLPSsd))
  
  LPS_prob_ref = getProb(LPS_score, groupMeans = c(refLPSmean, testLPSmean), groupSds = c(refLPSsd, testLPSsd))
  
  LPS_class = rep("UNCLASS",length(LPS_score))
  LPS_class[which(LPS_prob_test >= classProbCut)] = testGroup
  LPS_class[which(LPS_prob_ref >= classProbCut)] = refGroup
  
  LPS_score = data.frame(LPS_score)
  LPS_test = cbind(LPS_score, LPS_class, LPS_prob_test, LPS_prob_ref, stringsAsFactors =F)
  
  return(LPS_test)
  
}
