

#' Group mean and standard deviation calculation for all selected features
#' @description This is a wrap-up function to get group mean and sd for two groups and for all selected features. 
#' @details This is specifically designed for PRPS approach, however, it can be applied for other methods as well when 
#' group information is unknown, or mean and sd info from previous data can not be used in current data set 
#' due to data not comparable
#' @param testdat a data matrix with columns for samples and rows for features
#' @param selectedTraits a selected feature list
#' @param selectedtraitWeights a numeric vector with weight for all selected features, 
#'   which are in the same order as selectedtraits, however, only the sign of weight is used in this function
#' @param group1ratioPrior a prior ratio to indicate the ratio of test group over reference group 
#'   regarding mean and sd calculation for each selected feature
#' @return A numeric matrix of two group means and sds for all selected traits
#' @keywords mean sd prior
#' @author Aixiang Jiang
#' @references 
#' Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P,
#'  Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, 
#'  Kridel R, Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. 
#'  Double-Hit Gene Expression Signature Defines a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell Lymphoma.
#'   J Clin Oncol. 2018 Dec 3:JCO1801583. doi: 10.1200/JCO.18.01583.
#' @export
getMeanSdAllTraits = function(testdat, selectedTraits, selectedTraitWeights, group1ratioPrior = 0.5){
  meanSds = t(mapply(FUN=getMeanSdPerTrait, selectedTraits, selectedTraitWeights, MoreArgs = list(fulldat = testdat, ratioPrior =group1ratioPrior)))
  colnames(meanSds) = c("group1mean","group0mean","group1sd","group0sd")
  meanSds = data.frame(meanSds)
  return(meanSds)
}
