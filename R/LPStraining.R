
#' LPS training
#' @description This is the wrap up function to select top features, estimate parameters, and calculate LPS (Linear Predictor Score) scores 
#'  based on a given training data set.  
#' @details LPS calculation is based on Wright 2003. The formula is straightforward:
#' \eqn{LPS(X) = \sum a_j x_ij}
#' Here a_j represents the jth selected feature weights, and x_ij is the corresponding feature value 
#'  for the ith sample.
#' Before we start LPS training, we provide options for NA imputation and standardization. 
#' There are three steps for LPS training:
#' a) apply "getTrainingWeights" to select features and return weights for these features;
#' b) use "apply" function to get LPS classification scores and Empirical Bayes' probabilites for all samples;
#' When a Empirical Bayesian probability is calculated, by default, the 1st group in the input mean and sd vectors is treated as the 
#' test group. When we calculate the probabilities, we first calcualte probability that a sample belongs to either group,
#' and then use the following formula to get Empirical Bayesian probability:
#' \eqn{prob(x) = d_test(x)/(d_test(x) + d_ref(x))}
#' Here prob(x) is the Empirical Bayesian probability of a given sample, d_test(x) is the density value assuming that a given sample
#' belongs to the test group, d_ref(x) is the density value assuming that a given sample belongs to the reference group.
#' In the current LPStraining function, however, we calculate Empirical Bayesian probabilities for both directions.
#' c) This wrap-up function also gives classification for the training group and confusion matrix to compare 
#' LPS classification with original group info for training data set.
#' If NAs are not imputed, they are ignored for feature selection, weight calculation, LPS parameter estimation, 
#' and LPS calculation.
#' @param trainDat a training data set that is in data matrix or a data frame, samples are in columns, and features/traits are in rows
#' @param standardization a logic variable to indicate if standardization is needed before classification 
#'  score calculation
#' @param selectedTraits  a selected feature/trait list if available
#' @param groupInfo a known group classification, which order should be in the same as in columns of trainDat
#' @param refGroup the code for reference group, default is the 1st item in groupInfo
#' @param topN an integer to indicate how many top features to be selected
#' @param FDRcut  a FDR cutoff to select top features, which is only valid when topN is set as defaul NULL, 
#'  all features will be returned if both topN and FDRcut are set as default NULL
#' @param weightMethod  a string to indicate a weight calculation method or a significant variable selection method, 
#' there are five choices: "limma" for limma linear model based t value,"ttest" for t test based t value, 
#'  "MannWhitneyU" for Mann Whitney U based rank-biserial,"PearsonR" for Pearson correlation coefficient,
#'  "SpearmanR" for Spearman correlation coefficient, and the defualt value is "ttest"
#' @param classProbCut a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
#'  often used values are 0.8 and 0.9, default value is 0.9. The same classProbCut is used for both groups, 
#'  the samples that are not included in either group will be assigned as UNCLASS
#' @param imputeNA a logic variable to indicate if NA imputation is needed, if it is TRUE, 
#'  NA imputation is processed before any other steps, the default is FALSE
#' @param byrow a logic variable to indicate direction for imputation, default is TRUE, 
#'  which will use the row data for imputation
#' @param imputeValue a character variable to indicate which value to be used to replace NA, default is "median", 
#'  the median value of the chose direction with "byrow" data to be used
#' @keywords LPS training limma weight
#' @return A list with four items is returned: LPS parameters for selected features, LPS scores and classifications for training samples, 
#' confusion matrix to compare classification based on LPS scores and original classification, and a simple classification comparison table.
#' \item{LPS_pars}{a list of 2 items, the 1st item is a data frame with weights and group comparison results of each selected features 
#' for LPS calculation, and the 2nd item is a numeric vector containing LPS mean and sd for two groups}
#' \item{LPS_train}{a data frame of LPS score, true classification, Empirical Bayesian probabilites for both groups, 
#' and its classification for all training samples, notice that the classification is based on probabilities instead of LPS scores, 
#' and there is a UNCLASS group besdies the given two groups}
#' \item{classCompare}{a confusion matrix list object that compares LPS classification (based on selected features and 
#' weights) to input group classification for training data set, notice that the samples with UNCLASS 
#' are excluded since confusion matrix can not compare 3 groups to 2 groups}
#' \item{classTable}{a table to display comparison between LPS classification (based on selected features and weights) 
#' and input group classification for the given training data set. Since UNCLASS is excluded from confusion matrix, we add this table for full comparison}
#' @author Aixiang Jiang
#' @references
#' 
#' Bauer, D. F. (1972). Constructing confidence sets using rank statistics. Journal of the American Statistical Association 67, 687–690. 
#' doi: 10.1080/01621459.1972.10481279.
#'
#' Best, D. J. & Roberts, D. E. (1975). Algorithm AS 89: The Upper Tail Probabilities of Spearman's rho. Applied Statistics, 24, 377–379. doi: 10.2307/2347111.
#'
#' Hollander, M. & Wolfe, D.A. (1973), Nonparametric Statistical Methods. New York: John Wiley & Sons. Pages 185–194 (Kendall and Spearman tests).
#' 
#' Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. 
#' Statistical Applications in Genetics and Molecular Biology, 3, No. 1, Article 3. http://www.statsci.org/smyth/pubs/ebayes.pdf
#' 
#' Smyth, G. K., Michaud, J., and Scott, H. (2005). The use of within-array replicate spots for assessing differential expression in 
#' microarray experiments. Bioinformatics 21(9), 2067-2075.
#'
#' Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A gene expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.

#' @export
LPStraining = function(trainDat, standardization = FALSE, selectedTraits = NULL, groupInfo, refGroup = NULL, topN = NULL, FDRcut = NULL,
  weightMethod = c("ttest","limma","PearsonR", "SpearmanR", "MannWhitneyU"), classProbCut = 0.9, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median", "mean")){
 
  groupInfo = as.character(groupInfo)
  if(is.null(refGroup)){
    refGroup = groupInfo[1]
  }
  
  weightMethod = weightMethod[1]
  imputeValue = imputeValue[1]
  
  ## impute NA if imputeNA is true
  if(imputeNA){
    trainDat = imputeNAs(dataIn = trainDat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for LPS approach, it does not require standardization, however, if standardization = TRUE, do the standardization
  if(standardization){trainDat = standardize(trainDat)}
  
  # use reference group, we can select top features and calculate their weights
  weights = getTrainingWeights(trainDat = trainDat, selectedTraits = selectedTraits, groupInfo = groupInfo,
                               refGroup = refGroup, topN = topN, FDRcut = FDRcut, weightMethod = weightMethod)
  
  # get LPS scores for all samples
  LPS_score = apply(data.matrix(trainDat[rownames(weights),]), 2, getLPSscore, coefs = weights[,1])
  
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
  
  LPS_pars =  list(weights,meansds)
  names(LPS_pars) = c("weights","meansds")
  
  #### since UNCLASS is excluded from confusion matrix, add one more output for full comparison
  classTable = table(groupInfo, LPS_class)
  
  outs = list(LPS_pars, LPS_train, classCompare, classTable)
  names(outs) = c("LPS_pars","LPS_train","classCompare", "classTable")
  return(outs)
  
}


