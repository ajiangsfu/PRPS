
#' LPS score calculation and binary classification for a testing new data set based on a self learning object
#' @description This is the function to calculate LPS (Linear Prediction Score) scores for a testing data set 
#' with LPS self learning object. The selected feature list, these features' parameters are from the given LPS self learning object
#' that is an output from either LPSSLwithWeightsEM or LPSSLwithWeightsPrior.
#' @details  This is the function to calculate LPS scores and make classification and 
#' Empirical Bayesian probabilities for a testing new data set. However, this new data set should be 
#' comparable to the self learning data set as much as possible. Within LPS self learning and within this current 
#' functions, standardization step is included as an option to minimize the difference 
#' between self learning and testing data sets, but this step is only done to make distributions of each selected
#' features comparable. Be aware that this feature-wise standardization cannot make the sample-wise distributions
#' comparable. For example, the self learning data set must have two classification groups, however, 
#' the proportion of one group sample might be much less than the other group in the testing data set 
#' compared to the self learning data set, or even worse, the testing data set might only contain one classification 
#' group only. This is the common problem for classification and feature-wise standardization cannot solve the problem. 
#' In order to solve the problem, we should make data comparable as much as possbile before classification step. 
#' For example, use the same pre-processing settings and make suitable batch effect correction. 
#' 
#' LPS calculation is based on Wright 2003. The fomula is straightforward:
#'   \eqn{LPS(X) = \sum a_j x_ij}
#' Here a_j represents the jth selected feature weights, and x_ij is the corresponding feature value for the ith sample.
#' LPS mean and sd for two groups are needed for Empirical Bayes' probability calculation, which are from a LPS 
#' self learning object that is an output from either LPSSLwithWeightsEM or LPSSLwithWeightsPrior.
#' 
#' When calculate a Empirical Bayes' probability, the 1st group in the input mean and sd vectors is treated as the 
#' test group. When calculate the probabilities, we first calcualte probability that a sample belongs to either group,
#' and then use the following formula to get Empirical Bayes' probability:
#' \eqn{prob(x) = p_test(x)/(p_test(x) + p_ref(x))}
#' Here prob(x) is the Empirical Bayes' probability of a given sample, p_test(x) is the probability that a given sample
#' belongs to the test group, p_ref(x) is the probability that a given sample belongs to the reference group.
#' Notice that the test and reference group is just the relative grouping, in fact, for this step, we often need
#'  to calculate Empirical Bayes' probabilities for a given sample from two different standing points.
#' 
#' @param LPSSLObj a LPS self learning object, which is the output from either LPSSLwithWeightsEM or LPSSLwithWeightsPrior
#' @param newdat a new data matrix or data frame, which is comparable to self learning data set, 
#'  with columns for samples and rows for features
#' @param standardization a logic variable to indicate if standardization is needed before classification 
#'  score calculation
#' @param classProbCut a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
#'  often used values are 0.8 and 0.9, default value is 0.9. Only one value is used for both groups, 
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

LPSSLextension = function(LPSSLObj, newdat, standardization=FALSE, classProbCut = 0.9,  imputeNA = FALSE, byrow = TRUE, imputeValue =c("median","mean")){
  imputeValue = imputeValue[1]
  
  if(is.null(LPSSLObj)){print("Please input your LPS self learning object")}
  LPS_pars = LPSSLObj$LPS_pars
  weights = LPS_pars$weights
  
  ## imputee NA if imputeNA is true
  if(imputeNA){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for LPS approach, it does not require standardization, however, if standardization = TRUE, do the standardization
  if(standardization){newdat = standardize(newdat)}
  
  # if NA is not imputed, remove it from LPS score calculation
  LPS_score = apply(data.matrix(newdat[rownames(weights),]), 2, getLPSscore, coefs = weights[,1])
  
  # for LPS, 0 is a NOT natural cutoff for two group classification
  # in order to get classification, need to get two groups' LPS mean and sd
  
  ### get group info
  testres = LPSSLObj$LPS_test
  testres = testres[order(testres[,1], decreasing = T),]
  testGroup = testres[1,2]
  
  refGroup = setdiff(unique(LPSSLObj$LPS_test$PRPS_class),c(testGroup, "UNCLASS"))
  
  # for LPS, 0 is a NOT natural cutoff for two group classification
  # in order to get classification, need to get two groups' LPS mean and sd
  
  testLPSmean = LPS_pars$meansds[1]
  refLPSmean = LPS_pars$meansds[2]
  testLPSsd = LPS_pars$meansds[3]
  refLPSsd = LPS_pars$meansds[4]
 
  LPS_prob_test = getProb(LPS_score, groupMeans = c(testLPSmean, refLPSmean), groupSds = c(testLPSsd, refLPSsd))
  
  LPS_prob_ref = getProb(LPS_score, groupMeans = c(refLPSmean, testLPSmean), groupSds = c(refLPSsd, testLPSsd))
  
  LPS_class = rep("UNCLASS",length(LPS_score))
  LPS_class[which(LPS_prob_test >= classProbCut)] = testGroup
  LPS_class[which(LPS_prob_ref >= classProbCut)] = refGroup
  
  LPS_score = data.frame(LPS_score)
  LPS_test = cbind(LPS_score, LPS_class, LPS_prob_test, LPS_prob_ref, stringsAsFactors =F)
  
  return(LPS_test)
  
}
