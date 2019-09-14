#' @export
getMeanSdPerTrait = function(atrait, aweight, fulldat, ratioPrior = 0.5) {
  traitdat = unlist(fulldat[atrait,])
  ### remove NA
  tmp = which(is.na(traitdat))
  if(length(tmp)>0){traitdat = traitdat[-tmp]}
  mm = quantile(traitdat,probs = 1-ratioPrior)
  m1=subset(traitdat,traitdat >= mm)
  m2=subset(traitdat,traitdat < mm)
  if (aweight<0){
    mm = quantile(traitdat,probs = ratioPrior)
    m1=subset(traitdat,traitdat <= mm)
    m2=subset(traitdat,traitdat > mm)
  }
  
  ### notice that in this function, all data points are pushed to one of the two groups for feature level mean, sd calculation
  
  return(c(mean(m1), mean(m2), sd(m1), sd(m2)))
  
}
