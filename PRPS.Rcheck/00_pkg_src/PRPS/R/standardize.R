
#' A standardization Function
#'
#' @description This function is to standardize data in order to make features or samples be comparable. 
#' You can to standardize your data either by row or by column.
#' @details 
#' Standardization is the standard way to make data comparable. For this package, 
#' this should be done right after imputeNAs if necessary. In general, 
#' we do not recommend this step before classification, however, if the training and testing data sets
#' are obviously not comparable, we suggest to standardize training and testing data separately for each gene. 
#' It should be aware that this step might inflate sample differences for the low expression genes.

#' Standardization is done as:
#'  \eqn{ X = (X-mean(X))/sd(X) }

#' In this step, NA are removed from mean and sd calculation, which means that this step is safe even imputeNAs is not done before this step when they are NAs in the data.
#' 
#' @param dataIn  a data matrix or a data frame, samples are in columns, and features/genes are in rows
#' @param byrow although we recommend to standardize data along with features, 
#' but this function can be used for either direction, the default value is byrow = TRUE
#' @return Standardized data matrix or data frame
#' @keywords standardization
#' @author Aixiang Jiang
#' @export
standardize = function(dataIn, byrow = TRUE){
  ind = ifelse(byrow == TRUE, 1,2)
  outs = apply(dataIn, ind, FUN = function(xx){
    mm = mean(xx, na.rm = TRUE)
    ss = sd(xx, na.rm = TRUE)
    (xx-mm)/ss
  })

  if(byrow == TRUE){
    outs = t(outs)
  }
  return(outs)
}
