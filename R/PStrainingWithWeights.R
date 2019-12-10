#' PS training with weights
#' @description  This is the wrap up function to select top features, estimate parameters, 
#'  and calculate PS (Prediction Strength) scores based on a given training data set. 
#' @details PS calculation is based on Golub 1999. In this wrap up function, we use four steps to calculate 
#'  PS scores and classification. The range of PS scores is [-1,1]. Before these four steps, we also give an option
#'  for NA imputation. The four steps are:
#' a) apply "standardize" to standardize input data matrix for each feature;
#' b) apply "getTrainingWeights" to select features and return weights for these features;
#' c) apply "getMeanOfGroupMeans" to get mean of group means for each selected feature;
#' d) use "apply" function to get PS scores for all samples with "getPS1sample", the formula is:
#'   \eqn{PS = (V_win − V_lose)/(V_win + V_lose)}
#' Here, where V_win and V_lose are the vote totals for the winning and losing features/traits for a given sample
#' This function also give classification for the training group and confusion matrix to compare PS classification
#'  with original group info for training data set.
#' If NAs are not imputed, they are ignored for feature selection, weight calculation, PS parameter estimation,
#'  and PS calculation.
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
#' @keywords PS training limma weight
#' @return A list with three items is returned: PS parameters for selected features, PS scores and classifications for training samples, and confusion matrix to compare classification based on PS scores and original classification.
#' \item{PS_pars}{a data frame with all parameters needed for PS calculation for each selected features}
#' \item{PS_train}{a data frame of PS score, true classification and its classification based on scores for all training samples}
#' \item{classCompare}{a confusion matrix list object that compare PS classification based on selected features and weights compared to input group classification for training data set}
#' \item{classTable}{a table to display comparison of PS classification based on selected features and weights compared to input group classification for training data set}
#' @references 
#' Golub TR, Slonim DK, Tamayo P, Huard C, Gaasenbeek M, Mesirov JP, et al. Molecular classification of cancer: 
#' class discovery and class prediction by gene expression monitoring. Science. 1999;286:531–7
#' @export
PStrainingWithWeights = function(trainDat, groupInfo, refGroup = NULL, weights, classProbCut = 0.9,
                      imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  
  groupInfo = as.character(groupInfo)
  if(is.null(refGroup)){
    refGroup = groupInfo[1]
  }
  
  imputeValue = imputeValue[1]
  
  ## impute NA if imputeNA is true
  if(imputeNA){
    trainDat = imputeNAs(dataIn = trainDat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for PS approach, should standardize data first
  trainDat = standardize(trainDat)
  
  # use reference group, we can select top features and calculate their weights
  tmp = intersect(names(weights), rownames(trainDat))
  weights = weights[tmp]
  traindat = trainDat[tmp,]
  
  #  now, get mean of groups' means for each selected top features
  mean_2means = t(apply(trainDat[names(weights),], 1, getMeanOfGroupMeans, groupInfo = groupInfo, refGroup = refGroup))
  colnames(mean_2means) = c("meanOfGroupMeans","refGroupMean","testGroupMean")
  
  # in fact, the above two objects can be combined into one matrix, and put the two most important items into col1 and col2
  PS_pars = cbind( mean_2means[,1],weights, mean_2means[,-1])
  colnames(PS_pars)[1:2] = c(colnames(mean_2means)[1])
  
  # get PS scores for all samples
  PS_score = apply(trainDat[names(weights),], 2, getPS1sample, PSpars = PS_pars[,1:2])
  
  # for PS, 0 is a natural cutoff for two group classification
  testGroup = setdiff(unique(groupInfo), refGroup)
  PS_class0 = ifelse(PS_score >= 0, testGroup, refGroup)
  
  #### 20190905, use the true class to calculate PS group mean and sd, which are used for empirial Bayesian prob calculation
  refind = which(groupInfo == refGroup)
  refPS = PS_score[refind]
  testPS = PS_score[-refind]
  
  refPSmean = mean(refPS, na.rm = T)
  refPSsd = sd(refPS, na.rm = T)
  
  testPSmean = mean(testPS, na.rm = T)
  testPSsd = sd(testPS, na.rm = T)
 
  PS_prob_test = getProb(PS_score, groupMeans = c(testPSmean, refPSmean), groupSds = c(testPSsd, refPSsd))
  PS_prob_ref= getProb(PS_score, groupMeans = c(refPSmean, testPSmean), groupSds = c(refPSsd, testPSsd))
  
  PS_class = rep("UNCLASS",length(PS_score))
  PS_class[which(PS_prob_test >= classProbCut)] = testGroup
  PS_class[which(PS_prob_ref >= classProbCut)] = refGroup
  
  true_class = groupInfo
  PS_score = data.frame(PS_score)
  PS_train = cbind(PS_score, true_class, PS_class, PS_prob_test, PS_prob_ref,PS_class0, stringsAsFactors = F)
  
  groupInfo = factor(groupInfo, levels = c(refGroup, testGroup))
  ## in order to get comparison, change UNCLASS to NA, therefore only two groups are considered in PS_class
  PS_class2 = ifelse(PS_class == "UNCLASS", NA, PS_class)
  PS_class2 = factor(PS_class, levels = c(refGroup,  testGroup))
  ### PS_class2 is for confusion matirx only, we do not export it
  ## notice that confusion matrix does not work if the number of levels are not the same
  
  classCompare = caret::confusionMatrix(PS_class2, reference = groupInfo, positive = testGroup)
  
  meansds = c(testPSmean, refPSmean, testPSsd, refPSsd)
  names(meansds) = c("testPSmean","refPSmean","testPSsd","refPSsd")
  
  weights = data.frame(weights)
  
  PS_pars =  list(weights,meansds, mean_2means)  ### notice that this is different from the PS_pars above
  names(PS_pars) = c("weights","meansds","traits")
  
  #### since UNCLASS is excluded from confusion matrix, add one more output for full comparison
  classTable = table(groupInfo, PS_class)
  
  outs = list(PS_pars, PS_train, classCompare, classTable)
  names(outs) = c("PS_pars","PS_train","classCompare", "classTable")
  return(outs)
 
}
