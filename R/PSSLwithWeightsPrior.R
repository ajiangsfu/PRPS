PSSLwithWeightsPrior = function(newdat, weights, ratioPrior = 1/3, standardization=FALSE,  PShighGroup = "PShigh", 
                    PSlowGroup = "PSlow", EMmaxRuns = 50, breaks = 50, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  imputeValue = imputeValue[1]
  ## imputee NA if imputeNA is true
  if(imputeNA){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for PS approach, it does not require standardization, however, if standardization = TRUE, do the standardization
  if(standardization){newdat = standardize(newdat)}
  
  ### in case that some features in weight list are actually not in data, do the following
  tmp = intersect(names(weights), rownames(newdat))
  weights = weights[tmp]
  newdat = newdat[tmp,]

  
  #### in order to get mean and sd for each feature, use a ratio prior
  meansd = getMeanSdAllTraits(testdat = newdat, selectedTraits=names(weights), selectedTraitWeights=weights,
                              group1ratioPrior = ratioPrior)
  
  #  now, get mean of groups' means for each selected top features
  mean_2means = rowMeans(meansd[,1:2])
  mean_2means = cbind(mean_2means, meansd[,2], meansd[,1])
  colnames(mean_2means) = c("meanOfGroupMeans","refGroupMean","testGroupMean")
  
  # in fact, the above two objects can be combined into one matrix, and put the two most important items into col1 and col2
  PS_pars = cbind( mean_2means[,1],weights, mean_2means[,-1])
  colnames(PS_pars)[1:2] = c(colnames(mean_2means)[1],colnames(weights)[1])
  
  # get PS scores for all samples
  PS_score = apply(trainDat[names(weights),], 2, getPS1sample, PSpars = PS_pars[,1:2])
  
  PS_class0 = ifelse(PS_score > 0,  PShighGroup,  PSlowGroup)  
  
  #### 20190418: add classification calls based on Empirical Bayesian probabilities with EM (Expectation-Maximization)
  emcut = AdaptGauss::EMGauss(PS_score, K = 2,fast=TRUE, MaxNumberofIterations = EMmaxRuns)
  
  ### add a plot, hist with two distribution lines, do not need to save, just plot it
  hist(PS_score, prob = TRUE, breaks = breaks)
  curve(emcut$Weights[1]*dnorm(x, mean=emcut$Means[1], sd=emcut$SDs[1]), 
        col="red", lwd=2, add=TRUE, yaxt="n")
  curve(emcut$Weights[2]*dnorm(x, mean=emcut$Means[2], sd=emcut$SDs[2]), 
        col="green", lwd=2, add=TRUE, yaxt="n")
  abline(v=0, col = "red")
  
  PS_score = data.frame(PS_score)
  PS_test = cbind(PS_score, PS_class0,stringsAsFactors =F)
 
  
  #######################################################
  #### I should check all PS SL related functions, this is the 1st one, stop here on 0717
  ######################################################
  
  
  
  
  weights = data.frame(weights)
  
  PS_pars =  list(weights, meansds = tmp, traitsmeansds = cbind(allmeans, allsds))
  names(PS_pars) = c("weights","meansds","traitsmeansds")
  
  outs = list(PS_pars, PS_test)
  names(outs) = c("PS_pars","PS_test")
  
  return(outs)
 
}
