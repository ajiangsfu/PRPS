#' @export
getTraitParsWithEM = function(datin, weights, EMmaxRuns = 50){
  ### first of all, make sure all of the names in weights are also in datin, and make sure that theya are in the same order
  ### I did this step in the current functions that I call this function, however, to avoid any potential problem when is used elsewhere,
  ### do it again
  tmp = intersect(rownames(datin), names(weights))
  datin = datin[tmp,]
  weights = weights[tmp]
  
  ### for each feature, get EM objects 
  res = t(apply(cbind(weights,datin), 1, function(xx){
    x1 = xx[1]
    xx=xx[-1]

    ### update on 20190503
    #require(mclust)
    #require(AdaptGauss)
    clustG = mclust::Mclust(xx, G=2:4)
    bestG = clustG$G
    emcut = AdaptGauss::EMGauss(xx, K = bestG, fast=TRUE, MaxNumberofIterations = EMmaxRuns)
    
    ### no matter how many bestG, only keep the 1st and last one for the following
    xm = c(emcut$Means[1], emcut$Means[bestG])
    sds = c(emcut$SDs[1], emcut$SDs[bestG])

    ## without knowing further info, all I can do is to make the mean of the 1st group has the consistent direction with weights
    if(sign(xm[1]-xm[2])*sign(x1)>0){
      ms = xm
      sds = sds
    }else{
      ms = rev(xm)
      sds = rev(sds)
    }
    mm = mean(ms)
    
    return(c(ms,sds,mm))
  }))
  
  ### make sure that it is a data.frame with appropriate columns
  colnames(res) = c("group1mean","group0mean","group1sd","group0sd","meanOfMeans")
  res = data.frame(res)
  ### here, group1: same direction as in weight
  ### return
  return(res)
}