

##### this need to be changed dramatically, only keep the format of the following


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
#' @param groupInfo a known group classifica
#' 
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








