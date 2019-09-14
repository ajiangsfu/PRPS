#' @export
getMeanOfGroupMeans= function(traitdat, groupInfo, refGroup = 0){  
  
  groupInfo = ifelse(groupInfo == refGroup, 0, 1)
  
  ntmp = table(groupInfo)
  g0 = which(groupInfo == 0)
  x0 = traitdat[g0]
  x1 = traitdat[-g0]
  
  m0 = mean(x0, na.rm = TRUE)
  m1 = mean(x1, na.rm = TRUE)
  ### should also output the mean of the two group mean
  mm = mean(c(m0, m1))
  
  return(c(mm, m0, m1))
  
}
