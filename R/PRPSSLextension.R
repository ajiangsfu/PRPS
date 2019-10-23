
#' PRPS score calculation and binary classification for a testing new data set based on a self learning object
#' @description This is the function to calculate PRPS (Probability ratio based classification predication score)
#'  scores and make binary classification calls for a testing data set with PRPS_stableSLwithWeights or PRPSSLwithWeights or PRPSSLwithWeightsPrior or PRPSSLwithWeightsEM self learning object. 
#'  The selected feature list, these features' parameters are from the given PRPS_stableSLwithWeightsor PRPSSLwithWeights or PRPSSLwithWeightsPrior or PRPSSLwithWeightsEM object.
#' @details  This is the function to calculate PRPS scores, Empirical Bayesian probabilities and make classification 
#' for a testing new data set. However, this new data set should be comparable to the self learning data set used for PRPS_stableSLwithWeights or PRPSSLwithWeights
#' or PRPSSLwithWeightsPrior or PRPSSLwithWeightsEM as much as possible. Within PRPS_stableSLwithWeights/PRPSSLwithWeights/PRPSSLwithWeightsPrior/PRPSSLwithWeightsEM 
#' and within this current PRPSSLextension functions, 
#'  standardization step is included as an option to minimize the difference between self learning and testing data sets, 
#'  but this step is only done to make distributions of each selected feature comparable. 
#'  Be aware that this feature-wise standardization cannot make the sample-wise distributions comparable. 
#'  For example, the self learning data set must have two classification groups, however, the proportion of one group sample
#'  might be much less than the other group in the testing data set compared to the self learning data set, or even worse, 
#'  the testing data set might only contain one classification group only. This is the common problem for classification 
#'  and feature-wise standardization cannot solve the problem. 
#'  In order to solve the problem, we should make data comparable as much as possbile before classification step. 
#'  For example, use the same pre-processing settings and make suitable batch effect correction. 
#'  For classification with PRPS approach, we also suggest to combine self learning and testing data together as "newdat" 
#'  for this PRPSSLextension function, to avoid forcing two groups' classification while there is actual only one group
#'  in the testing data set.
#'  PRPS calculation is based on Ennishi 2018. The fomula is 
#'  \eqn{PRPS(X_i) = \sum (|a_j| log(P1(x_ij)/P0(x_ij)))}
#'  Here, a_j represents the jth selected feature weights, and x_ij is the corresponding feature value for the ith sample, 
#'  P1 and P0 are the probabilities that the ith sample belongs to two different group.
#'  The therotic cutoff is 0 to make classification calls based on PRPS score, alternatively, we can use Bayes approach to make calls.
#'  When we calculate a Empirical Bayes' probability, usually, the group with smaller proportion is treated as
#'  the test group while the other group with higher proportion is treated as reference group. 
#'  When calculate the probabilities, we first calcualte probability that a sample belongs to either group, 
#'  and then use the following formula to get Empirical Bayes' probability:
#' \eqn{prob(x) = d_test(x)/(d_test(x) + d_ref(x))}
#' Here prob(x) is the Empirical Bayes' probability of a given sample, d_test(x) is the density value
#'  that a given sample belongs to the test group, d_ref(x) is the density value that a given sample belongs
#'   to the reference group.
#'  Notice that the test and reference group is just the relative grouping, in fact, for this step, 
#'  we often need to calculate Empirical Bayes' probabilities for a given sample from two different standing points.
#' @param PRPSSLObj a PRPS self learning object that is the output from function
#' PRPS_stableSLwithWeights or PRPSSLwithWeights or PRPSSLwithWeightsPrior or PRPSSLwithWeightsEM.
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
#' @return  A data frame with PRPS scores, Empirical Bayesian probabilites for two groups and classification, and classification based on 0 natural cutoff on PRPS scores. 
#' @keywords PRPS
#' @author Aixiang Jiang
#' @references
#'  Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P, 
#'  Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, Kridel R,
#'  Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. 
#'  Double-Hit Trait Expression Signature Defines a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell
#'  Lymphoma. J Clin Oncol. 2018 Dec 3:JCO1801583. doi: 10.1200/JCO.18.01583.
#'  
#' Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A trait expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.
#'  
#' @export


PRPSSLextension = function(PRPSSLObj, newdat, standardization=FALSE,  classProbCut = 0.9,
                           imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  imputeValue = imputeValue[1]
  
  if(is.null(PRPSSLObj)){print("Please input your PRPS self learning object")}
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

