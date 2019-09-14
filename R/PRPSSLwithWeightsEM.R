#' @export
PRPSSLwithWeightsEM = function(newdat, weights, standardization=FALSE, classProbCut = 0.9, PRPShighGroup = "PRPShigh", 
                    PRPSlowGroup = "PRPSlow", breaks = 50, EMmaxRuns = 50, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
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

  #### in order to get mean and sd for each feature
  #### which I could 
  #### a) re-write the original code to be a new function -> this seems most reasonable solution
  #### b) re-write code and post here directly
  #### c) or change the original function
  
  # if NA is not imputed, remove it from PRPS score calculation in the following getPRPSscore function
  #PRPS_score = apply(data.matrix(newdat[names(weights),]), 2, getPRPSscore, coefs = weights)
   ############## 
  
  ### think again, on 20190329
  ### I can still call: weightedLogProbClass = function(newdat, topTraits, weights, classMeans, classSds)
  ### however, before that, I need to use EM to get classMeans, classSds
  
  ### this is to get mean and sd for each selected feature based on EM, write an independent function and then call here
  ### in order to be used in both PRPS and PS, give mean of group means as well
  ### this is: need: 2 means, 2 sds, mean of 2 means, return should be a data frame, and call by colnames in the following part like here
  ### PRPSpars = apply(data.matrix(newdat[names(weights),]), 1, function(xx){
  ### }, coefs = weights)
  PRPSpars = getTraitParsWithEM(datin = newdat, weights = weights, EMmaxRuns = EMmaxRuns)
  PRPS_score = weightedLogProbClass(newdat=newdat, topTraits=names(weights), weights=weights, 
                                       classMeans=PRPSpars[,1:2], classSds=PRPSpars[,3:4])
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
