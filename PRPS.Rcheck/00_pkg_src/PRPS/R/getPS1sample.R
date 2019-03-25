

#' PS score calculation for a given sample
#' @description This is an internal function to calculate PS (Prediction Strength) for a given sample.  
#' @details PS calculation is based on Golub 1999. In this package, we use four steps to calculate PS scores, 
#' and this function is for the very last step. 
#' The four steps are: a) apply "standardize" to standardize input data matrix for each feature; 
#' b) apply "getTrainingWeights" to select features and return weights for these features; 
#' c) apply "getMeanOfGroupMeans" to get mean of group means for each selected feature; 
#' d) use "apply" function to get PS for all samples with this current function "getPS1sample".
#' NA is ignored for PS calculation.
#' @param vdat a numeric vector of one sample with multiple selected features
#' @param PSpars  a PS parameter matrix, which rows are the features with the same order as in vdat, 
#'   and the columns are the mean of two group means and weight for each feature
#' @return PS score for a single sample
#' @keywords PS 
#' @author Aixiang Jiang
#' @references
#' TR Golub, DK Slonim, P Tamayo, C Huard, M Gaasenbeek, JP Mesirov, H Coller, ML Loh, JR Downing, MA Caligiuri, et al.
#' Molecular classification of cancer: class discovery and class prediction by gene expression monitoring
#' Science, 286 (1999), pp. 531-537

#' @export
getPS1sample = function(vdat, PSpars){
  #### vdat is a vector of feature data
  #### PSpars is a  matrix, the rows are genes that are in the same order as in vdat
  ####                       the col 1 is for mean of two group means, the col2 is weight
  
  #### since there might be some NA in the data, should remove them before PS calculation
  tmp = which(is.na(vdat))
  if(length(tmp) > 0){
    vdat = vdat[-tmp]
    PSpars = PSpars[-tmp,]
  }
  xbg = vdat - PSpars[,1]
  vg = xbg * PSpars[,2]
  t1 = sum(vg)
  t2 = sum(abs(vg))
  ps = t1/t2
}

