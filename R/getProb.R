#' @export
getProb = function(inscore, groupMeans, groupSds){
  ### assume groupMeans contain 2 values for 2 group, and the 1st one is for positive group
  ### assume groupSds contain 2 values for 2 group, and the 1st one is for positive group
  p1 = dnorm(inscore,mean = groupMeans[1], sd= groupSds[1])
  p0 = dnorm(inscore,mean = groupMeans[2], sd= groupSds[2])
  return(p1/(p1+p0))
  ### notice that this is assumed normal distribution as traditionally
  ###    assuming that group mean and sd are calculated based on large sample size
  ### however, for small samples sized based mean and sd, t distribution might be a better choice
}

