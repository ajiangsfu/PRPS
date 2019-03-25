
#' NA value imputation
#' @description This function is to impute NA values if there are any NA in the data. 
#' Without this step, NAs will cause a lot of troubles in classification score calculation, and therefore, we have to ignore NAs.
#' @details During the classification score calculation process, NA offen causes troubles especially when matrix multiplication is involved. 
#' To avoid this problem, we suggest to impute NA values before any other functions are called if you have NA in your data. 
#' For classification purpose, the median value is believed to be most noninformative value, so we use median to replace NA. Further more, 
#' we prefer to do feature-wise median replacement although this function can be applied in both feature wise and sample wise replacements.
#' @param dataIn a data matrix or a data frame, samples are in columns, and features/genes are in rows
#' @param byrow although we recommend to impute NAs along with features, but this function can be used for either direction,
#'  the default value is byrow = TRUE
#' @param imputeValue for classification, the most non-informatic value is median/mean in a given row/column, 
#'  so there are two choices for the imputed values: median or mean when NAs are excluded, the default value is "median"
#' @return Same data matrix or data frame as input except that the NA values are imputed.
#' @keywords impute NA
#' @author Aixiang Jiang
#' @export

imputeNAs = function(dataIn, byrow = TRUE, imputeValue = c("median","mean")){ ### I usually put feature in the row, and sample in the column
  imputeValue = imputeValue[1]
  if(byrow){
    dataIn = t(apply(dataIn, 1, function(xx){
      if(imputeValue == "mean"){
        mxx = mean(xx, na.rm = T)
      }else{
        mxx = median(xx, na.rm = T)
      }
      nas = which(is.na(xx))
      xx[nas] = mxx
      return(xx)
    }))
    
  }else{
    dataIn = apply(dataIn, 2, function(xx){
      if(imputValue == "mean"){
        mxx = mean(xx, na.rm = T)
      }else{
        mxx = median(xx, na.rm = T)
      }
      nas = which(is.na(xx))
      xx[nas] = mxx
      return(xx)
    })
  }
  return(dataIn)
}
