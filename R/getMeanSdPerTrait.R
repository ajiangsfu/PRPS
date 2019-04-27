
#' Group mean and standard deviation calculation for a given feature
#' @description This is a function to get group mean and sd for two groups and for a given feature
#' @details This is specifically designed for PRPS approach, however, it can be applied for other methods as well when 
#' group information is unknown, or mean and sd info from previous data can not be used in current data set 
#' due to data not comparable
#' @param atrait a given feature name
#' @param aweight weight for a given feature, only the sign of weight is used in this function
#' @param fulldat a data matrix with columns for samples and rows for features
#' @param ratioPrior a prior ratio to indicate the ratio of test group over reference group 
#'   regarding mean and sd calculation for each selected feature  
#' @return A numeric vector of two group means and sds
#' @keywords mean sd prior
#' @author Aixiang Jiang
#' @references
#'  Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P,
#'  Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, 
#'  Kridel R, Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. 
#'  Double-Hit Gene Expression Signature Defines a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell Lymphoma.
#'   J Clin Oncol. 2018 Dec 3:JCO1801583. doi: 10.1200/JCO.18.01583.
#' @export
getMeanSdPerTrait = function(atrait, aweight, fulldat, ratioPrior = 0.5) {
  traitdat = as.numeric(fulldat[atrait,])
  mm = quantile(traitdat,probs = 1-ratioPrior)
  m1=subset(traitdat,traitdat >= mm)
  m2=subset(traitdat,traitdat < mm)
  if (aweight<0){
    mm = quantile(traitdat,probs = ratioPrior)
    m1=subset(traitdat,traitdat <= mm)
    m2=subset(traitdat,traitdat > mm)
  }
  
  ### notice that in this function, all data points are pushed to one of the two groups for feature level mean, sd calculation
  
  return(c(mean(m1, na.rm = T), mean(m2, na.rm = T), sd(m1, na.rm = T), sd(m2, na.rm = T)))
  
}
