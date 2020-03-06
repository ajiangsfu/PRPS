#' A function rather aimed at developers
#' @noRd

PRPSSLwithWeightsPrior = function(newdat, weights, standardization=FALSE,  classProbCut = 0.9, ratioPrior = 1/2, PRPShighGroup = "PRPShigh", 
                                  PRPSlowGroup = "PRPSlow", imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
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
  
  #### in order to get mean and sd for each feature, use a ratio prior
  meansd = getMeanSdAllTraits(testdat = newdat, selectedTraits=names(weights), selectedTraitWeights=weights,
                              group1ratioPrior = ratioPrior)
  
  df1 = as.integer(ratioPrior*dim(newdat)[2]) - 1
  ## df0 = dim(newdat)[2] - df1 - 1 ### this seems not correct
  ## change on 20200306
  ## df0 = as.integer((1-ratioPrior)*dim(newdat)[2]) - 1
  ## or:
  df0 = dim(newdat)[2] - df1 - 2
  
  dfs = c(df1, df0)
  
  PRPS_score = weightedLogProbClass(newdat=newdat, topTraits=names(weights), weights=weights,
                                    classMeans=meansd[,1:2], classSds = meansd[,3:4], dfs = dfs)
  
  PRPS_class0 = ifelse(PRPS_score > 0,  PRPShighGroup,  PRPSlowGroup)  
  
  # #### 20190503, call plotHistEM 
  #emsearch = plotHistEM(PRPS_score, G = 2:4, breaks = breaks, EMmaxRuns = EMmaxRuns, scoreName = "PRPS_score")
  # bestG = emsearch$bestG
  # emcut = emsearch$emcut
  # 
  # ### no matter how many bestG, only keep the 1st and last one for the following
  # gmeans = c(emcut$Means[1], emcut$Means[bestG])
  # gsds = c(emcut$SDs[1], emcut$SDs[bestG])
  # 
  # PRPS_prob1_EM = getProb(PRPS_score, groupMeans = gmeans, groupSds = gsds)
  # PRPS_prob2_EM = getProb(PRPS_score, groupMeans = rev(gmeans), groupSds = rev(gsds))
  # 
  # if(gmeans[1]> gmeans[2]){
  #   name1 = PRPShighGroup
  #   name2 = PRPSlowGroup
  # }else{
  #   name1 = PRPSlowGroup
  #   name2 = PRPShighGroup
  # }
  # 
  # PRPS_class_EM = rep("UNCLASS",length(PRPS_score))
  # PRPS_class_EM[which(PRPS_prob1_EM >= classProbCut)] = name1
  # PRPS_class_EM[which(PRPS_prob2_EM >= classProbCut)] = name2
  # 
  #### work on a given prior based classification
  # tmpcut = quantile(PRPS_score, probs = 1-ratioPrior)
  # ttmp = PRPS_score[which(PRPS_score >= tmpcut)]
  # rtmp = PRPS_score[which(PRPS_score < tmpcut)]
  # 
  # testPRPSmean = mean(ttmp)
  # refPRPSmean = mean(rtmp)
  # testPRPSsd = sd(ttmp)
  # refPRPSsd = sd(rtmp)
  
  #### note on 20200306, I am wondering why I do not the ratioPrior to get score parameters?
  ####  -> reason, the ratioPrior is the searching tool for gene level parameter estimation
  ####     once we have score, we should ignore our ratio assumption since it is most likely not corret
  ####  -> the following is exactly what I did for JCO 2018 paper, except here I use prob as final classification instead of PRPS 0 cutoff
  #### 20190904, use the 0 theoretical cutoff to get two groups, which are used for empirial Bayesian prob calculation
  ttmp = PRPS_score[which(PRPS_score >= 0)]
  rtmp = PRPS_score[which(PRPS_score < 0)]
  testPRPSmean = mean(ttmp)
  refPRPSmean = mean(rtmp)
  testPRPSsd = sd(ttmp)
  refPRPSsd = sd(rtmp)
  
  PRPS_prob1 = getProb(PRPS_score, groupMeans = c(testPRPSmean, refPRPSmean), groupSds = c(testPRPSsd, refPRPSsd))
  
  PRPS_prob2= getProb(PRPS_score, groupMeans = c(refPRPSmean, testPRPSmean), groupSds = c(refPRPSsd, testPRPSsd))
  
  PRPS_class = rep("UNCLASS",length(PRPS_score))
  PRPS_class[which(PRPS_prob1 >= classProbCut)] = PRPShighGroup
  PRPS_class[which(PRPS_prob2 >= classProbCut)] = PRPSlowGroup
  
  # tmp = c(emcut$Means, emcut$SDs)
  # names(tmp) = c("testPRPSmean_EM","refPRPSmean_EM","testPRPSsd_EM","refPRPSsd_EM")
  
  #####################################################
  ########### after I used all PRPS_score, then define it as data.frame, otherwise, many of code will get trouble
  PRPS_score = data.frame(PRPS_score)
  # PRPS_test = cbind(PRPS_score, PRPS_class, PRPS_prob1, PRPS_prob2, PRPS_class0, PRPS_class_EM, PRPS_prob1_EM, PRPS_prob2_EM,
  #                   stringsAsFactors =F)
  # 
  # weights = data.frame(weights)
  # 
  # PRPS_pars =  list(weights, meansds = c(testPRPSmean, refPRPSmean, testPRPSsd, refPRPSsd), traitsmeansds = meansd,
  #                  tmp)
  # names(PRPS_pars) = c("weights","meansds","traitsmeansds", "meansds_EM")
  # 
  PRPS_test = cbind(PRPS_score, PRPS_class, PRPS_prob1, PRPS_prob2, PRPS_class0, stringsAsFactors =F)
  
  #weights = data.frame(weights)
  
  PRPS_pars =  list(weights, meansds = c(testPRPSmean, refPRPSmean, testPRPSsd, refPRPSsd), traitsmeansds = meansd, dfs)
  names(PRPS_pars) = c("weights","meansds","traitsmeansds", "dfs")
  
  outs = list(PRPS_pars, PRPS_test)
  names(outs) = c("PRPS_pars","PRPS_test")
  
  return(outs)
  
}
