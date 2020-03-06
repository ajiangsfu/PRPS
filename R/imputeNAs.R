#' A function rather aimed at developers
#' @noRd

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
