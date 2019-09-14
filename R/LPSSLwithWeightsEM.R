#' @export
LPSSLWithWeightsEM = function(newdat, weights, standardization=FALSE, classProbCut = 0.9, LPShighGroup = "LPShigh", 
                    LPSlowGroup = "LPSlow", breaks = 50, EMmaxRuns = 50, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  imputeValue = imputeValue[1]
  ## imputee NA if imputeNA is true
  if(imputeNA){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for LPS approach, it does not require standardization, however, if standardization = TRUE, do the standardization
  if(standardization){newdat = standardize(newdat)}
  
  ### in case that some features in weight list are actually not in data, do the following
  tmp = intersect(names(weights), rownames(newdat))
  weights = weights[tmp]

  # if NA is not imputed, remove it from LPS score calculation in the following getLPSscore function
  LPS_score = apply(data.matrix(newdat[names(weights),]), 2, getLPSscore, coefs = weights)
  
  ### now, use EM to define 1st draft of group, and then use it to calculate group mean and sd
  
  emcut = AdaptGauss::EMGauss(LPS_score, K = 2, fast=TRUE, MaxNumberofIterations = EMmaxRuns)
  
  # > emcut
  # $Means
  # [1] -106.36599  -39.04443
  # 
  # $SDs
  # [1] 23.62815 11.70491
  # 
  # $Weights
  # [1] 0.1139536 0.8860464
  #### this is exact what I need!
  
  ### add a plot, hist with two distribution lines, do not need to save, just plot it
  hist(LPS_score, prob = TRUE, breaks = breaks)
  curve(emcut$Weights[1]*dnorm(x, mean=emcut$Means[1], sd=emcut$SDs[1]), 
        col="red", lwd=2, add=TRUE, yaxt="n")
  curve(emcut$Weights[2]*dnorm(x, mean=emcut$Means[2], sd=emcut$SDs[2]), 
        col="green", lwd=2, add=TRUE, yaxt="n")
  
  LPS_prob1 = getProb(LPS_score, groupMeans = emcut$Means, groupSds = emcut$SDs)
  LPS_prob2 = getProb(LPS_score, groupMeans = rev(emcut$Means), groupSds = rev(emcut$SDs))
  
  if(emcut$Means[1]> emcut$Means[2]){
    name1 = LPShighGroup
    name2 = LPSlowGroup
  }else{
    name1 = LPSlowGroup
    name2 = LPShighGroup
  }
  
  LPS_class = rep("UNCLASS",length(LPS_score))
  LPS_class[which(LPS_prob1 >= classProbCut)] = name1
  LPS_class[which(LPS_prob2 >= classProbCut)] = name2
  
  LPS_score = data.frame(LPS_score)
  LPS_test = cbind(LPS_score, LPS_class, LPS_prob1, LPS_prob2, stringsAsFactors =F)
  
  #####################################################
  #### add more item to return
  
  weights = data.frame(weights)
  
  meansds = c(testLPSmean, refLPSmean, testLPSsd, refLPSsd)
  names(meansds) = c("testLPSmean","refLPSmean","testLPSsd","refLPSsd")
  
  LPS_pars =  list(weights,meansds)
  names(LPS_pars) = c("weights","meansds")
  
  outs = list(LPS_pars, LPS_test)
  names(outs) = c("LPS_pars","LPS_test")
  
  return(outs)
  
}
