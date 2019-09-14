#' @export
PRPSSLwithWeights = function(newdat, weights, standardization=FALSE, classProbCut = 0.9, PRPShighGroup = "PRPShigh", 
                    PRPSlowGroup = "PRPSlow", breaks = 50, EMmaxRuns = 50, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  require(mclust)
  imputeValue = imputeValue[1]
  ## imputee NA if imputeNA is true
  if(imputeNA){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for PRPS approach, it does not require standardization, however, if standardization = TRUE, do the standardization
  if(standardization){newdat = standardize(newdat)}
  
  ### in case that some features in weight list are actually not in data, do the following
  tmp = intersect(names(weights), rownames(newdat))
  weights = weights[tmp]
  newdat = newdat[tmp,]
  
  rps = seq(0.05, 0.95, by = 0.05)
  
  rpsres = sapply(rps, FUN = function(xx){
    tmp = PRPSSLwithWeightsPrior(newdat=newdat, weights=weights, ratioPrior = xx, PRPShighGroup = PRPShighGroup, PRPSlowGroup = PRPSlowGroup)
    mcls = mclust::Mclust(tmp$PRPS_test$PRPS_score, G=2)
    return(mcls$classification)
  })
  
  rownames(rpsres) = colnames(newdat)
  ### the classification is coded with 1 and 2
  ### in order to find samples that are always 1 and always 2, what should I do?
  ### for always 1, if I extract 1 for each cell, then row sum for the sample should be 0
  ### for always 2, if I extract 2 for each cell, then row sum for the sample should be 0
  
  res1 = rpsres - 1
  res2 = rpsres - 2
  
  grp1 = rownames(res1[which(rowSums(res1) == 0),])
  grp2 = rownames(res1[which(rowSums(res2) == 0),])
  
  ### now, I need to work on group mean and sd for each feature
  datgrp1 = newdat[,grp1]
  datgrp2 = newdat[,grp2]
  
  means1 = rowMeans(datgrp1) 
  means2 = rowMeans(datgrp2)
  
  sds1 = apply(datgrp1,1,sd)
  sds2 = apply(datgrp2,1,sd)
  
  # if(length(grp1) >= length(grp2)){
  #   allmeans = cbind(means2, means1)
  #   allsds = cbind(sds2, sds1)
  # }else{
  #   allmeans = cbind(means1, means2)
  #   allsds = cbind(sds1, sds2)
  # }
  
  #### I need to determine the direction for allmeans, and allsds
  #### maybe I can use the middle prior to make the decision
  xx=median(rps)
  tmp = PRPSSLwithWeightsPrior(newdat=newdat, weights=weights, ratioPrior = xx, PRPShighGroup = PRPShighGroup, PRPSlowGroup = PRPSlowGroup)
  mcls = mclust::Mclust(tmp$PRPS_test$PRPS_score, G=2)
  #### in fact, there are group labels in mcls$mean
  #### get the two groups score means and decide which group is which
  
  mmean = mcls$parameters$mean
  if(mmean[2] > mmean[1]){
    allmeans = cbind(means2, means1)
    allsds = cbind(sds2, sds1)
  }else{
    allmeans = cbind(means1, means2)
    allsds = cbind(sds1, sds2)
  }
  
  PRPS_score = weightedLogProbClass(newdat = newdat, topTraits=names(weights), weights=weights,
                                    classMeans = allmeans, classSds = allsds)
  
  #### 20190503, call plotHistEM 
  emsearch = plotHistEM(PRPS_score, G = 2:4, breaks = breaks, EMmaxRuns = EMmaxRuns, scoreName = "PRPS_score")
  bestG = emsearch$bestG
  emcut = emsearch$emcut
  
  ### no matter how many bestG, only keep the 1st and last one for the following
  gmeans = c(emcut$Means[1], emcut$Means[bestG])
  gsds = c(emcut$SDs[1], emcut$SDs[bestG])
  
  PRPS_prob1 = getProb(PRPS_score, groupMeans = gmeans, groupSds = gsds)
  PRPS_prob2 = getProb(PRPS_score, groupMeans = rev(gmeans), groupSds = rev(gsds))
  
  if(gmeans[1]> gmeans[2]){
    name1 = PRPShighGroup
    name2 = PRPSlowGroup
    scoreMeanSds = c(gmeans, gsds)
  }else{
    name1 = PRPSlowGroup
    name2 = PRPShighGroup
    scoreMeanSds = c(rev(gmeans), rev(gsds))
  }
  
  PRPS_class = rep("UNCLASS",length(PRPS_score))
  PRPS_class[which(PRPS_prob1 >= classProbCut)] = name1
  PRPS_class[which(PRPS_prob2 >= classProbCut)] = name2
  
  names(scoreMeanSds) = c("testPRPSmean","refPRPSmean","testPRPSsd","refPRPSsd")
  
  PRPS_class0 = ifelse(PRPS_score > 0,  PRPShighGroup, PRPSlowGroup)
  
  PRPS_score = data.frame(PRPS_score)
  PRPS_test = cbind(PRPS_score, PRPS_class, PRPS_prob1, PRPS_prob2, PRPS_class0,stringsAsFactors =F)

  weights = data.frame(weights)
  
  PRPS_pars =  list(weights, meansds = scoreMeanSds, traitsmeansds = cbind(allmeans, allsds))
  names(PRPS_pars) = c("weights","meansds","traitsmeansds")
  
  outs = list(PRPS_pars, PRPS_test)
  names(outs) = c("PRPS_pars","PRPS_test")
  return(outs)
  
}
