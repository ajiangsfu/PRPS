

#' Feature selection, parameter estimation, and PS calculation for training data set
#' @description  This is the wrap up function to select top features, estimate parameters, 
#'  and calculate PS (Prediction Strength) scores based on a given training data set. 
#' @details PS calculation is based on Golub 1999. In this warp up function, we use four steps to calculate 
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
#' @param selectedTraits  a selected trait list if available
#' @param groupInfo a known group classification, which order should be the same as in colnames of trainDat
#' @param refGroup the code for reference group, default is 0, but it can be a string or other number, which will be 
#'  changed to 0 within the function
#' @param topN an integer to indicate how many top features to be selected
#' @param FDRcut  a FDR cutoff to select top features, which is only valid when topN is set as defaul NULL, 
#'  all features will be returned if both topN and FDRcut are set as default NULL
#' @param weightMethod  a string to indicate weight calculation method, there are five choices: 
#'  "limma" for for limma linear model based t value,"ttest" for t test based t value, 
#'  "MannWhitneyU" for Mann Whitney U based rank-biserial,"PearsonR" for Pearson correlation coefficient,
#'  "SpearmanR" for Spearman correlation coefficient, and the defualt value is "limma"
#' @param imputeNA a logic variable to indicate if NA imputation is needed, if it is TRUE, 
#'  NA imputation is processed before any other steps, the default is FALSE
#' @param byrow a logic variable to indicate direction for imputation, default is TRUE, 
#'  which will use the row data for imputation
#' @param imputeValue a character variable to indicate which value to be used to replace NA, default is "median", 
#'  the median value of the chose direction with "byrow" data to be used
#' @keywords PS training limma weight
#' @return A list with three items is returned: PS parameters for selected features, PS scores and classifications for training samples, and confusion matrix to compare classification based on PS scores and original classification.
#' \item{PS_pars}{a data frame with all parameters needed for PS calculation for each selected features}
#' \item{PS_train}{a data frame of PS score and its classification for all training samples}
#' \item{classCompare}{a confusion matrix list object that compare PS classification based on selected features and weights compared to input group classification for training data set}
#' \item{classTable}{a table to display comparison of PS classification based on selected features and weights compared to input group classification for training data set}
#' @references 
#' TR Golub, DK Slonim, P Tamayo, C Huard, M Gaasenbeek, JP Mesirov, H Coller, ML Loh, JR Downing, MA Caligiuri, et al.
#' Molecular classification of cancer: class discovery and class prediction by gene expression monitoring
#' Science, 286 (1999), pp. 531-537
#' @export
PStraining = function(trainDat, selectedTraits = NULL, groupInfo, refGroup = 0, topN = NULL, FDRcut = 0.1,
                      weightMethod = "limma", imputeNA = FALSE, byrow = TRUE, imputeValue = "median"){
  ### STEPs
  ### => before anything, do we need imput NAs?
  ### a) standardize
  ### b) getTrainingWeights
  ### c) getMeanOfGroupMeans
  ### d) use apply to get PS for all samples with getPS1sample
  
  ### return a list with several items
  ### a) wts
  ### b) mean_2means
  ### c) PS scores and classification
  ### d) confusion matrix to compare known groupInfo and the PS classification
  
  ## impute NA if imputeNA is true
  if(imputeNA == TRUE | imputeNA == T){
    trainDat = imputeNAs(dataIn = trainDat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for PS approach, should standardize data first
  trainDat = standardize(trainDat)
  
  # use reference group, we can select top features and calculate their weights
  weights = getTrainingWeights(trainDat = trainDat, selectedTraits = selectedTraits, groupInfo = groupInfo,
                               refGroup = refGroup, topN = topN, FDRcut = FDRcut, weightMethod = weightMethod)
  
  #  now, get mean of groups' means for each selected top features
  mean_2means = t(apply(trainDat[rownames(weights),], 1, getMeanOfGroupMeans, groupInfo = groupInfo, refGroup = refGroup))
  colnames(mean_2means) = c("meanOfGroupMeans","refGroupMean","testGroupMean")
  
  # in fact, the above two objects can be combined into one matrix, and put the two most important items into col1 and col2
  PS_pars = cbind( mean_2means[,1],weights[,1], mean_2means[,-1], weights[,-1])
  colnames(PS_pars)[1:2] = c(colnames(mean_2means)[1],colnames(weights)[1])
  
  # get PS scores for all samples
  PS_score = apply(trainDat[rownames(weights),], 2, getPS1sample, PSpars = PS_pars[,1:2])
  
  # for PS, 0 is a natural cutoff for two group classification
  testGroup = setdiff(unique(groupInfo), refGroup)
  
  PS_class = ifelse(PS_score >= 0, testGroup, refGroup)
  
  PS_score = data.frame(PS_score)
  PS_train = cbind(PS_score, PS_class, stringsAsFactors =F)
  
  groupInfo = factor(groupInfo, levels = c(refGroup, testGroup))
  PS_class = factor(PS_class, levels = c(refGroup, testGroup))
  classCompare = confusionMatrix(PS_class, reference = groupInfo, positive = testGroup)

  classTable = table(groupInfo, PS_class)
  
  outs = list(PS_pars, PS_train, classCompare, classTable)
  names(outs) = c("PS_pars","PS_train","classCompare", "classTable")
  
  
  return(outs)
}
