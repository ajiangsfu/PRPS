#' @export
standardize = function(dataIn, byrow = TRUE){
  ind = ifelse(byrow == TRUE, 1,2)
  outs = apply(dataIn, ind, FUN = function(xx){
    mm = mean(xx, na.rm = TRUE)
    ss = sd(xx, na.rm = TRUE)
    (xx-mm)/ss ### should I add a small number for stable standardization? not for now
  })

  if(byrow == TRUE){
    outs = t(outs)
  }
  return(outs)
}
