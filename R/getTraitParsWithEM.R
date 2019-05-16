#' Calculate two group means, sds, and mean of group means for each of features, which are used in PRPS and PS
#' @description This function is to calculate feature level group mean, sd and mean of means, without knowing sample group label.
#' @details In order to get feature level group means, sd sand mean of means, without knowing sample group label, we use EM algorithm to 
#' decompose a vector of values into two normal distributions, and obtain parameters (mean and sd) for each of them, 
#' and then get mean of group means. This is done for each given feature. In order to determine the direction for each 
#' feature, we do need weight for each selected feature.
#' @param datin a new data matrix or data frame, columns for samples and rows for features
#' @param weights a numeric vector with selected features (as names of the vector) and their weights
#' @param EMmaxRuns number of Iterations for EM searching; default=50
#' @return A data frame with group means, sds and mean of group means
#' @keywords EM mean sd
#' @author Aixiang Jiang
#' @references Ultsch, A., Thrun, M.C., Hansen-Goos, O., Loetsch, J.: Identification of Molecular Fingerprints
#' in Human Heat Pain Thresholds by Use of an Interactive Mixture Model R Toolbox(AdaptGauss),
#' International Journal of Molecular Sciences, doi:10.3390/ijms161025897, 2015.
#' 
#' #' Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and 
#' density estimation using Gaussian finite mixture models, The R Journal, 8/1, pp. 205-233.
#' 
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