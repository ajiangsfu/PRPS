
#' LPS testing
#' @description This is the function to calculate LPS (Linear Prediction Score) scores for a testing data set 
#' with a given LPS training object. The selected features with their weights, and two classes' LPS score distribution parametes
#' are extracted from the given LPS training object.
#' @details  This is the function to calculate LPS scores and make classification based on Empirical Bayesian probabilities 
#' for a new testing data set, which should be comparable to the training data set as much as possible. 
#' 
#' Within LPStraining and this LPStesting functions, standardization step is included as an option to minimize the difference 
#' between training and testing data sets. Whether or not a user decides to do standardization, this should be consistent between training 
#' and testing data sets, otherwise this current testing function will not work.
#' 
#' Notice that this step is only to make distributions of each selected features comparable within training or testing data sets. 
#' Be aware that this feature-wise standardization cannot make the sample-wise distributions
#' comparable. For example, the training data set must have two classification groups, however, 
#' the proportion of one group might be much smaller than the other group in the testing data set 
#' compared to the training data set, or even worse, the testing data set might only contain one classification 
#' group only. This is the common problem for classification, and feature-wise standardization cannot solve the problem. 
#' 
#' In order to solve the problem, we should make data comparable as much as possbile before classification step. 
#' For example, use the same pre-processing settings and make suitable batch effect correction. 
#' For classification with LPS approach, we also suggest to combine training and testing data together as a full data set 
#' in this LPStesting function, to avoid forcing samples into two groups' classification while there is actual only one group
#' in the testing data set.
#' 
#' LPS calculation is based on Wright 2003. The fomula is straightforward:
#'   \eqn{LPS(X) = \sum a_j x_ij}
#' Here a_j represents the jth selected feature weights, and x_ij is the corresponding feature value for the ith sample.
#'
#' When a Empirical Bayesian probability is calculated, by default, the 1st group in the input mean and sd vectors is treated as the 
#' test group. When we calculate the probabilities, we first calcualte probability that a sample belongs to either group,
#' and then use the following formula to get Empirical Bayesian probability:
#' \eqn{prob(x) = d_test(x)/(d_test(x) + d_ref(x))}
#' Here prob(x) is the Empirical Bayesian probability of a given sample, d_test(x) is the density value assuming that a given sample
#' belongs to the test group, d_ref(x) is the density value assuming that a given sample belongs to the reference group.
#' In the current function, however, we calculate Empirical Bayesian probabilities for both directions.
#' 
#' @param LPStraingObj a LPS training object, which is the output from LPStraining function 
#' @param newdat a new data matrix or data frame, which is comparable to the training data set, 
#'  its columns are for samples and rows are for features
#' @param standardization a logic variable to indicate if standardization is needed before classification 
#'  score calculation
#' @param classProbCut a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
#'  often used values are 0.8 and 0.9, default value is 0.9. The same classProbCut is used for both groups, 
#'  the samples that are not included in either group will be assigned as UNCLASS
#' @param imputeNA a logic variable to indicate if NA imputation is needed, if it is TRUE, NA imputation is 
#'  processed before any other steps, the default is FALSE
#' @param byrow a logic variable to indicate direction for imputation, default is TRUE, 
#'  which will use the row data for imputation
#' @param imputeValue a character variable to indicate which value to be used to replace NA, default is "median", 
#'  the median value of the chosen direction with "byrow" data to be used
#' @return A data frame with LPS score, Empirical Bayesian probabilites for two groups and classification
#' @keywords LPS
#' @author Aixiang Jiang
#' @references Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A gene expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.
#' @export

LPStesting = function(LPStrainObj, newdat, standardization=FALSE, classProbCut = 0.9,  imputeNA = FALSE, byrow = TRUE, imputeValue =c("median","mean")){
  imputeValue = imputeValue[1]
  
  if(is.null(LPStrainObj)){print("Please input your LPS training object")}
  LPS_pars = LPStrainObj$LPS_pars
  weights = LPS_pars$weights
  
  ## imputee NA if imputeNA is true
  if(imputeNA){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for LPS approach, it does not require standardization, however, if standardization = TRUE, do the standardization
  if(standardization){newdat = standardize(newdat)}
  
  # if NA is not imputed, remove it from LPS score calculation
  LPS_score = apply(data.matrix(newdat[rownames(weights),]), 2, getLPSscore, coefs = weights[,1])
  
  # for LPS, 0 is NOT a natural cutoff for two group classification
  # in order to get classification, we need to get two groups' LPS mean and sd
  
  ### get group info
  testGroup = LPStrainObj$classCompare$positive

  refGroup = setdiff(unique(LPStrainObj$LPS_train$true_class),testGroup)
  
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
