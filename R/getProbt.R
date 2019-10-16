#' @export
getProbt = function(inscore, groupMeans, groupSds, dfs){
  t1 = (inscore - groupMeans[1])/groupSds[1]
  t0 = (inscore - groupMeans[2])/groupSds[2]
  d1 = dt(t1, df=dfs[1])
  d0 = dt(t0, df=dfs[2])
  return(d1/(d1+d0))
}
