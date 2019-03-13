

#' Mean of two groups' means calculation
#' @description This is an internal function called by PStraining to calculate mean of two groups' means for a feature,
#'  which is needed for PS appraoch. However, it can be called directly if you do not have PStraining object but want to calcualate
#'  PS scores, in this case, however, please make sure that you have get subset of selected features data only as traitdat, and
#'  provide groupInfo and refGroup, then, the output of this function can be used to get PS scores via getClassScores function.
#' @details When mean of each group feature data is calculated, the NAs are removed. 
#' @param traitdat a vector of data for a single feature
#' @param groupInfo a known group classification, which sample order should be the same as for traitdat
#' @param refGroup the code for reference group, default is 0, but it can be a string or other number, 
#'  which will be changed to 0 within the function
#' @return A numerical vector of mean of two group means, and the two group means
#' @keywords mean sd
#' @author Aixiang Jiang
#' @export
getMeanOfGroupMeans= function(traitdat, groupInfo, refGroup = 0){  
  
  groupInfo = ifelse(groupInfo == refGroup, 0, 1)
  
  ntmp = table(groupInfo)
  g0 = which(groupInfo == 0)
  x0 = traitdat[g0]
  x1 = traitdat[-g0]
  
  m0 = mean(x0, na.rm = TRUE)
  m1 = mean(x1, na.rm = TRUE)
  ### should also output the mean of the two group mean
  mm = mean(c(m0, m1))
  
  return(c(mm, m0, m1))
  
}
