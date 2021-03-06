
#' PRPS self-training (used to be called self-learning) extension 
#' @description This is a PRPS self-training extension function to calculate PRPS (Probability ratio based classification predication score)
#'  scores and make binary classification calls for a testing data set with a PRPS self-training object, e.g., output of PRPSstableSLwithWeights.
#'  The selected feature list, these features' parameters are extracted from the given PRPS self-training object.
#' @details  This is the function to calculate PRPS scores, Empirical Bayesian probabilities and make binary classification 
#' for a testing data set. This new testing data set should be comparable to the self-training data set as much as possible. 
#  
#' Within this current function, standardization step is included as an option to minimize the difference between self-training and testing data sets. 
#' Whether or not a user decides to do standardization, this should be consistent between self-training and testing data sets, 
#' otherwise this current testing function will not work.

#'  Notice that standardization step is only done to make distributions of each selected feature comparable within each data set. 
#'  Be aware that this feature-wise standardization cannot make the sample-wise distributions comparable. 
#'  For example, the self-training data set must have two classification groups, however, the proportion of one group 
#'  might be much smaller than the other group in the testing data set compared to the self-training data set, or even worse, 
#'  the testing data set might contain one classification group only. This is the common problem for classification 
#'  and feature-wise standardization cannot solve the problem. 
#'
#'  In order to solve the problem, we should make data comparable as much as possbile before classification step. 
#'  For example, use the same pre-processing settings and make suitable batch effect correction. 
#'  For classification with PRPS approach, we also suggest to combine self-training and testing data together as "newdat" 
#'  for this PRPSSLextension function, to avoid forcing samples into two groups while there is actual only one group
#'  in the testing data set.
#'  
#'  PRPS calculation is based on Ennishi 2018. The fomula is: 
#'  \eqn{PRPS(X_i) = \sum (|a_j| log10(P1(x_ij)/P0(x_ij)))}
#'  Here, a_j represents the jth selected feature weights, and x_ij is the corresponding feature value for the ith sample, 
#'  P1 and P0 are the probabilities that the ith sample belongs to two different groups. 
#'  The therotic cutoff is 0 to make classification calls based on PRPS score, alternatively, we can use empirical Bayesian approach to make calls.
#'  
#'  When a Empirical Bayesian probability is calculated, by default, the 1st group in the input mean and sd vectors is treated as the 
#' test group. When we calculate the probabilities, we first calcualte probability that a sample belongs to either group,
#' and then use the following formula to get Empirical Bayesian probability:
#' \eqn{prob(x) = d_test(x)/(d_test(x) + d_ref(x))}
#' Here prob(x) is the Empirical Bayesian probability of a given sample, d_test(x) is the density value assuming that a given sample
#' belongs to the test group, d_ref(x) is the density value assuming that a given sample belongs to the reference group.
#' In the current function, however, we calculate Empirical Bayesian probabilities for both directions.
#'  
#' @param PRPSSLObj a PRPS self-training object that is the output from function
#' PRPS_stableSLwithWeights or PRPSSLwithWeights or PRPSSLwithWeightsPrior or PRPSSLwithWeightsEM.
#' @param newdat a new data matrix or data frame, which is comparable to self-training data set, 
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
#' @return  A data frame with PRPS scores, Empirical Bayesian probabilites for two groups and classification, 
#'  and classification based on 0 natural cutoff on PRPS scores. 
#' @keywords PRPS
#' @author Aixiang Jiang
#' @references
#'  Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P, 
#'  Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, Kridel R,
#'  Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. 
#'  Double-Hit Gene Expression Signature Defines a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell
#'  Lymphoma. J Clin Oncol. 2019 Jan 20;37(3):190-201.
#'  
#' Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A gene expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.
#'  
#' @export


PRPSSLextension = function(PRPSSLObj, newdat, standardization=FALSE,  classProbCut = 0.9,
                           imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  imputeValue = imputeValue[1]
  
  if(is.null(PRPSSLObj)){print("Please input your PRPS self-training object")}
  PRPS_pars = PRPSSLObj$PRPS_pars
  weights = PRPS_pars$weights
  dfs = PRPS_pars$dfs
  
  ## impute NA if imputeNA is true
  if(imputeNA){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for PRPS approach, it does not require standardization, however, if standardization = TRUE, do the standardization
  if(standardization){newdat = standardize(newdat)}
  
  Traitsmeansds = PRPS_pars$traitsmeansds
  
  # if NA is not imputed, remove it from PRPS score calculation
  PRPS_score = weightedLogProbClass(newdat, topTraits=rownames(weights), weights=weights[,1],
                                    classMeans = Traitsmeansds[,1:2], classSds = Traitsmeansds[,3:4], dfs = dfs)
  
  ### get group info
  testres = PRPSSLObj$PRPS_test
  testres = testres[order(testres[,1], decreasing = T),]
  testGroup = testres[1,2]
  
  refGroup = setdiff(unique(PRPSSLObj$PRPS_test$PRPS_class),c(testGroup, "UNCLASS"))
  
  # for PRPS, 0 is a natural cutoff for two group classification
  # use 0 to get class0, refer the code in the training part
  PRPS_class0 = ifelse(PRPS_score>0, testGroup, refGroup)
  
  # alternatively, we can get classification based on prob, and we need to get two groups' PRPS mean and sd first from training
  # if testing and trainng data sets (more accurately: if the score are comparable) are comparable
  
  testPRPSmean = PRPS_pars$meansds[1]
  refPRPSmean = PRPS_pars$meansds[2]
  testPRPSsd = PRPS_pars$meansds[3]
  refPRPSsd = PRPS_pars$meansds[4]
  
  PRPS_prob_test = getProb(PRPS_score, groupMeans = c(testPRPSmean, refPRPSmean), groupSds = c(testPRPSsd, refPRPSsd))
  
  PRPS_prob_ref = getProb(PRPS_score, groupMeans = c(refPRPSmean, testPRPSmean), groupSds = c(refPRPSsd, testPRPSsd))
  
  PRPS_class = rep("UNCLASS",length(PRPS_score))
  PRPS_class[which(PRPS_prob_test >= classProbCut)] = testGroup
  PRPS_class[which(PRPS_prob_ref >= classProbCut)] = refGroup
  
  PRPS_score = data.frame(PRPS_score)
  PRPS_test = cbind(PRPS_score, PRPS_class, PRPS_prob_test, PRPS_prob_ref, PRPS_class0, stringsAsFactors =F)
  
  return(PRPS_test)
  
}

