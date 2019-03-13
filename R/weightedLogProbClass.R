
#' PRPS score calculation for a data set
#' @description This is usually as an internal function called by PRPStraining, PRPStesting and getClassScores, which is used to actually
#'  calculate PRPS scores However, it can be also called directly.
#' @details #'  PRPS calculation is based on Ennishi 2018. The fomula is 
#'  \eqn{PRPS(X) = \sum a_j x_ij}
#'  Here, a_j represents the jth selected feature weightss, and x_ij is the corresponding feature value for the ith sample, P1 and P0 are the probabilities that the ith sample belongs to two different group.
#'  When calculate a Empirical Bayes' probability, the 1st group in the input mean and sd vectors is treated as
#'  the test group. 
#'  If there are NAs in the data and not imputed before calling this function, these NAs will be ignored for PRPS score calculation.
#'  If you want to call this function directly, make sure to give all values for its parameters. If you do not have PRPStraining object
#'   but would like to calculate PRPS scores, this is function is much easier to use than PRPStesting, and this function is 
#'   identical if you work on PRPS approach.
#' @param newdat a new data matrix or data frame, with columns for samples and rows for features
#' @param topTraits selected features used for PRPS calculation
#' @param weights feature weights
#' @param classMeans a numeric vector of two group means, the 1st item is for testing group, 
#'  while the 2nd item is for reference group
#' @param classSds a numeric vector of two group standard deviations (sds), the 1st item is for testing group, 
#'  while the 2nd item is for reference group
#' @return A numeric vector of PRPS 
#' @keywords PRPS
#' @author Aixiang Jiang
#' @references Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P, 
#'  Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, Kridel R, 
#'  Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. Double-Hit Gene Expression Signature Defines
#'  a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell Lymphoma. J Clin Oncol. 
#'  2018 Dec 3:JCO1801583. doi: 10.1200/JCO.18.01583.
#' @export
weightedLogProbClass = function(newdat, topTraits, weights, classMeans, classSds) {
  
  genedat = newdat[topTraits,]
  
  names(weights) = rownames(classMeans)
  weights = weights[topTraits]
  
  #### mapply is good for multiple lists, need change format before using the mapply function
  
  classMeans = t(classMeans[topTraits,])
  colnames(classMeans) = rownames(genedat)
  classMeans = data.frame(classMeans)
  
  classSds = t(classSds[topTraits,])
  colnames(classSds) = rownames(genedat)
  classSds = data.frame(classSds)
  
  gendatt=t(genedat)
  gendatt = data.frame(gendatt)
  
  lograt = mapply(FUN = function(xx,yy,zz){
    dd=0.01
    ### or, we can set dd as a parameter before call the whole function,
    ####   which can be related to 5% or 10% quantile or so
    
    t1=(xx-yy[1])/(zz[1]+dd) ### notice that for a single value xx, n=1, so use sd directly as denominator
    t2=(xx-yy[2])/(zz[2]+dd)
    
    #### then convert t1 and t2 to p values
    #2*(1 - pt(tval, df))
    # or: 2 * pt(abs(t_value), df, lower.tail = FALSE)
    p1=2 * pt(abs(t1), df=1, lower.tail = FALSE)
    p2=2 * pt(abs(t2), df=1, lower.tail = FALSE)
    gg=log10(p1) - log10(p2)
  },gendatt,classMeans, classSds)
  
  rownames(lograt) = colnames(genedat)
  
  rm(genedat)
  rm(gendatt)
  gc()
  
  
  ### get final y for decision make
  
  #res=lograt %*% abs(weights) ### this line does not consider NA, should change
  # borrow the function for LPS, since it considers the NA within the function
  res = apply (lograt, 1, getLPSscore, coefs= abs(weights))
  
  return(res)
  
}
