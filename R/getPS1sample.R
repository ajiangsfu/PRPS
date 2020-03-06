#' A function rather aimed at developers
#' @noRd

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
  #### when all data are NA, return NA
  ps = NA
  if(length(vdat) > 0){
    xbg = vdat - PSpars[,1]
    vg = xbg * PSpars[,2]
    t1 = sum(vg)
    t2 = sum(abs(vg))
    ps = t1/t2
  }
  return(ps)
}

