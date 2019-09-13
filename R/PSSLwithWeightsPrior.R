### note on 20190904
### 1) all data should be standardized, therefore there is no parameter for standardization in this function
### 2) for PS, should I use empirical call in the end? struggled for a while, decide to add this in
### 3) while I try to keep similar output as for PRPSSLwithWeightsPrior, the 3rd item of this function item is different
###     this is because I use different parameter for PS score calculation if extension is needed

PSSLwithWeightsPrior = function(newdat, weights, ratioPrior = 1/3,  PShighGroup = "PShigh", 
                    PSlowGroup = "PSlow", imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  imputeValue = imputeValue[1]
  ## imputee NA if imputeNA is true
  if(imputeNA){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for PS approach, should standardize data first
  newdat = standardize(newdat)
  
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
  
  # #### 20190418: add classification calls based on Empirical Bayesian probabilities with EM (Expectation-Maximization)
  # emcut = AdaptGauss::EMGauss(PS_score, K = 2,fast=TRUE, MaxNumberofIterations = EMmaxRuns)
  # 
  # ### add a plot, hist with two distribution lines, do not need to save, just plot it
  # hist(PS_score, prob = TRUE, breaks = breaks)
  # curve(emcut$Weights[1]*dnorm(x, mean=emcut$Means[1], sd=emcut$SDs[1]), 
  #       col="red", lwd=2, add=TRUE, yaxt="n")
  # curve(emcut$Weights[2]*dnorm(x, mean=emcut$Means[2], sd=emcut$SDs[2]), 
  #       col="green", lwd=2, add=TRUE, yaxt="n")
  # abline(v=0, col = "red")
  # 

  #### 20190904, use the 0 theoretical cutoff to get two groups, which are used for empirial Bayesian prob calculation
  ttmp = PS_score[which(PS_score >= 0)]
  rtmp = PS_score[which(PS_score < 0)]
  testPSmean = mean(ttmp)
  refPSmean = mean(rtmp)
  testPSsd = sd(ttmp)
  refPSsd = sd(rtmp)
  
  PS_prob1 = getProb(PS_score, groupMeans = c(testPSmean, refPSmean), groupSds = c(testPSsd, refPSsd))
  
  PS_prob2= getProb(PS_score, groupMeans = c(refPSmean, testPSmean), groupSds = c(refPSsd, testPSsd))
  
  PS_class = rep("UNCLASS",length(PS_score))
  PS_class[which(PS_prob1 >= classProbCut)] = PShighGroup
  PS_class[which(PS_prob2 >= classProbCut)] = PSlowGroup
  
  # tmp = c(emcut$Means, emcut$SDs)
  # names(tmp) = c("testPSmean_EM","refPSmean_EM","testPSsd_EM","refPSsd_EM")
  
  #####################################################
  ########### after I used all PS_score, then define it as data.frame, otherwise, many of code will get trouble
  PS_score = data.frame(PS_score)
  # PS_test = cbind(PS_score, PS_class, PS_prob1, PS_prob2, PS_class0, PS_class_EM, PS_prob1_EM, PS_prob2_EM,
  #                   stringsAsFactors =F)
  # 
  # weights = data.frame(weights)
  # 
  # PS_pars =  list(weights, meansds = c(testPSmean, refPSmean, testPSsd, refPSsd), traitsmeansds = meansd,
  #                  tmp)
  # names(PS_pars) = c("weights","meansds","traitsmeansds", "meansds_EM")
  # 
  PS_test = cbind(PS_score, PS_class, PS_prob1, PS_prob2, PS_class0, stringsAsFactors =F)
  
  weights = data.frame(weights)
  
  PS_pars =  list(weights, meansds = c(testPSmean, refPSmean, testPSsd, refPSsd), traitsmeans = mean_2means)
  names(PS_pars) = c("weights","meansds","traitsmeans")
  
  outs = list(PS_pars, PS_test)
  names(outs) = c("PS_pars","PS_test")
  
  return(outs)
  
}
