
#' LPS score calculation for a given sample
#' @description  This is an internal function to calculate LPS (Linear Prediction Score),
#'  which is called by other functions but it can also be called directly.
#' @details 
#' LPS calculation is based on Wright 2003, the formula is staightforward:
#'  \eqn{LPS(X) = \sum a_j X_j}
#' If NAs are not imputed, they are ignored for LPS calculation.
#' @param vdat a numeric vector for all selected features for a given sample
#' @param coefs a numeric vector of weights for all selected features, 
#'  which should be the same order as in "vdat"
#' @return A LPS score
#' @keywords LPS
#' @author Aixiang Jiang
#' @references Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A gene expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.
#' @export
#' 
getLPSscore = function(vdat, coefs){
  tmp = which(is.na(vdat))
  if(length(tmp) > 0){
    vdat = vdat[-tmp]
    coefs = coefs[-tmp]
  }
  aScore = mapply(function(x,y){x*y}, x=vdat, y=coefs)
  aScore = sum(unlist(aScore))
  return(aScore)
}
