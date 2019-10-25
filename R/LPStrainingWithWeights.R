
#' LPS calculation for training data set with given weights
#' @description This is to calculate LPS (Linear Prediction Score) scores based on a given training data set with classification
#' labels and selected features' weights  
#' @details LPS calculation is based on Wright 2003. The fomula is straightforward:
#' \eqn{LPS(X) = \sum a_j x_ij}
#' Here a_j represents the jth selected feature weights, and x_ij is the corresponding feature value 
#'  for the ith sample.
#' In this function, we use "apply" function to get LPS classification scores and Empirical Bayes' probabilites for all samples.
#' When we calculate a Empirical Bayes' probability, the 1st group in the input mean and sd vectors is treated 
#' as the test group. When we calculate the probabilities, we first calcualte probability that a sample belongs
#' to either group, and then use the following formula to get Empirical Bayes' probability:
#' \eqn{prob(x) = p_test(x)/(p_test(x) + p_ref(x))}
#' Here prob(x) is the Empirical Bayes' probability of a given sample, p_test(x) is the probability 
#' that a given sample belongs to the test group, p_ref(x) is the probability that a given sample 
#' belongs to the reference group.
#' Notice that the test and reference group is just the relative grouping, in fact, for this step, 
#' we often need to calculate Empirical Bayes' probabilities for a given sample from two different standing points.
#' c) . This function also give classification for the training group and confusion matrix to compare 
#' LPS classification with original group info for training data set.
#' If NAs are not imputed, they are ignored for feature selection, weight calculation, LPS parameter estimation, 
#' and LPS calculation.
#' @param trainDat training data set, a data matrix or a data frame, samples are in columns, and features/traits are in rows
#' @param weights a numeric vector with selected features (as names of the vector) and their weights
#' @param groupInfo a known group classification, which order should be the same as in colnames of trainDat
#' @param refGroup the code for reference group, default is the 1st item in groupInfo
#' @param classProbCut a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
#'  often used values are 0.8 and 0.9, default value is 0.9. Only one value is used for both groups, 
#'  the samples that are not included in either group will be assigned as UNCLASS
#' @param imputeNA a logic variable to indicate if NA imputation is needed, if it is TRUE, 
#'  NA imputation is processed before any other steps, the default is FALSE
#' @param byrow a logic variable to indicate direction for imputation, default is TRUE, 
#'  which will use the row data for imputation
#' @param imputeValue a character variable to indicate which value to be used to replace NA, default is "median", 
#'  the median value of the chose direction with "byrow" data to be used
#' @param standardization a logic variable to indicate if standardization is needed before classification 
#'  score calculation
#' @keywords LPS training
#' @return A list with four items is returned: LPS parameters for selected features, LPS scores and classifications for training samples, and confusion matrix to compare classification based on LPS scores and original classification.
#' \item{LPS_pars}{a list of 2 items, the 1st item is a data frame with weights and group testing results of each selected features for LPS calculation, and the 2nd item is a numeric vector containing LPS mean and sd for two groups}
#' \item{LPS_train}{a data frame of LPS score, true classification, Empirical Bayesian probabilites for both groups, and its classification for all training samples, notice that the classification is based on probabilities instead of LPS scores, and there is UNCLASS group besdies the given two groups}
#' \item{classCompare}{a confusion matrix list object that compare LPS classification based on selected features and weights compared to input group classification for training data set, notice that the samples with UNCLASS are excluded since confusion matrix can not compare 3 groups to 2 groups}
#' \item{classTable}{a table to display comparison of LPS classification based on selected features and weights compared to input group classification for training data set. Since UNCLASS is excluded from confusion matrix, add this table for full comparison}
#' @author Aixiang Jiang
#' @references
#' Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A trait expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.
#' @export
LPStrainingWithWeights = function(trainDat, weights, groupInfo, refGroup = NULL, 
          classProbCut = 0.9, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median", "mean"), standardization = FALSE){

  groupInfo = as.character(groupInfo)
  if(is.null(refGroup)){
    refGroup = groupInfo[1]
  }
  
  imputeValue = imputeValue[1]
  
  ## impute NA if imputeNA is true
  if(imputeNA){
    trainDat = imputeNAs(dataIn = trainDat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for LPS approach, it does not require standardization, however, if standardization = TRUE, do the standardization
  if(standardization){trainDat = standardize(trainDat)}

  # get LPS scores for all samples
  LPS_score = apply(data.matrix(trainDat[names(weights),]), 2, getLPSscore, coefs = weights)
  
  # and use get prob function to get classification
  
  testGroup = setdiff(unique(groupInfo), refGroup)
  
  refind = which(groupInfo == refGroup)
  refLPS = LPS_score[refind]
  testLPS = LPS_score[-refind]
  
  refLPSmean = mean(refLPS, na.rm = T)
  refLPSsd = sd(refLPS, na.rm = T)
  
  testLPSmean = mean(testLPS, na.rm = T)
  testLPSsd = sd(testLPS, na.rm = T)
  
  LPS_prob_test = getProb(LPS_score, groupMeans = c(testLPSmean, refLPSmean), groupSds = c(testLPSsd, refLPSsd))
  LPS_prob_ref = getProb(LPS_score, groupMeans = c(refLPSmean, testLPSmean), groupSds = c(refLPSsd, testLPSsd))
  
  LPS_class = rep("UNCLASS",length(LPS_score))
  LPS_class[which(LPS_prob_test >= classProbCut)] = testGroup
  LPS_class[which(LPS_prob_ref >= classProbCut)] = refGroup
  
  LPS_score = data.frame(LPS_score)
  true_class = groupInfo
  LPS_train = cbind(LPS_score, true_class, LPS_class, LPS_prob_test, LPS_prob_ref, stringsAsFactors =F)
  
  groupInfo = factor(groupInfo, levels = c(refGroup, testGroup))
  ## in order to get comparison, change UNCLASS to NA, therefore only two groups are considered in LPS_class
  LPS_class2 = ifelse(LPS_class == "UNCLASS", NA, LPS_class)
  LPS_class2 = factor(LPS_class, levels = c(refGroup,  testGroup))
  
  ## notice that confusion matrix does not work if the number of levels are not the same
  
  classCompare = caret::confusionMatrix(LPS_class2, reference = groupInfo, positive = testGroup)

  meansds = c(testLPSmean, refLPSmean, testLPSsd, refLPSsd)
  names(meansds) = c("testLPSmean","refLPSmean","testLPSsd","refLPSsd")
  
  weights = data.frame(weights)
  
  LPS_pars =  list(weights,meansds)
  names(LPS_pars) = c("weights","meansds")
  
  #### since UNCLASS is excluded from confusion matrix, add one more output for full comparison
  classTable = table(groupInfo, LPS_class)
  
  outs = list(LPS_pars, LPS_train, classCompare, classTable)
  names(outs) = c("LPS_pars","LPS_train","classCompare", "classTable")
  return(outs)
  
}


