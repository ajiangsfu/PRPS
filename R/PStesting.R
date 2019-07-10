
#' PS score calculation for a testing new data set
#' @description  This is the function to calculate PS (Prediction Strength) scores for a testing data set
#'  with PS training object. The selected feature list, these features' parameters are 
#'  from the given PS training object.
#' @details This is the function to calculate PS scores and make classification for a testing new data set. 
#'  However, this new data set should be comparable to the training data set as much as possible. 
#'  Within PStraining and within this current PStesting functions, standardization step is included
#'  to minimize the difference between training and testing data sets, but this step is only done to
#'  make distributions of each selected features comparable. Be aware that this feature-wise 
#'  standardization cannot make the sample-wise distributions comparable. 
#'  For example, the training data set must have two classification groups, however, the proportion of 
#'  one group sample might be much less than the other group in the testing data set compared to the 
#'  training data set, or even worse, the testing data set might only contain one classification group only. 
#'  This is the common problem for classification and feature-wise standardization cannot solve the problem. 
#'  In order to solve the problem, we should make data comparable as much as possbile before classification step.
#'  For example, use the same pre-processing settings and make suitable batch effect correction. 
#'  For classification with PS approach, we also suggest to combine traing and testing data together as "newdat"
#'   for this PStesting function, to avoid forcing two groups' classification while there is actual only one group
#'   in the testing group.
#' @param PStraingObj a PS training object, which is the output from function PStraining
#' @param newdat a new data matrix or data frame, which is comparable to training data set, 
#'  with columns for samples and rows for features
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
#' @return A data frame with PS score and classification
#' @keywords PS
#' @author Aixiang Jiang
#' @references
#' Golub TR, Slonim DK, Tamayo P, Huard C, Gaasenbeek M, Mesirov JP, et al. Molecular classification of cancer: 
#' class discovery and class prediction by gene expression monitoring. Science. 1999;286:531–7
#' @export

PStesting = function(PStrainObj, newdat, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  
  if(is.null(PStrainObj)){print("Please input your PS training object")}
  PS_pars = PStrainObj$PS_pars
  
  imputeValue = imputeValue[1]
  
  ## impute NA if imputeNA is true
  if(imputeNA){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  ### standardize data first
  newdat = standardize(newdat)
  
  # get PS scores for all samples
  PS_score = apply(newdat[rownames(PS_pars),], 2, getPS1sample, PSpars = PS_pars[,1:2])
  
  testGroup = PStrainObj$classCompare$positive
  refGroup = setdiff(unique(PStrainObj$PS_train$PS_class),testGroup)
  
  # for PS, 0 is a natural cutoff for two group classification
  
  PS_class = ifelse(PS_score >= 0, testGroup, refGroup)
  
  PS_score = data.frame(PS_score)
  PS_test = cbind(PS_score, PS_class, stringsAsFactors =F)
  
  return(PS_test)
  
}

