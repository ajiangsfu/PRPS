LPSSLwithWeightsPrior = function(newdat, weights, ratioPrior = 1/3, testGroup, refGroup, isTestGroupHighLPS = TRUE,
                      standardization=FALSE, classProbCut = 0.9,  imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  imputeValue = imputeValue[1]
  ## imputee NA if imputeNA is true
  if(imputeNA){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for LPS approach, it does not require standardization, however, if standardization = TRUE, do the standardization
  if(standardization){newdat = standardize(newdat)}
  
  # if NA is not imputed, remove it from LPS score calculation
  LPS_score = apply(data.matrix(newdat[names(weights),]), 2, getLPSscore, coefs = weights)
  
  # for LPS, 0 is a NOT natural cutoff for two group classification
  # in order to get classification, need to get two groups' LPS mean and sd
  
  # if(isTestGroupHighLPS == TRUE | isTestGroupHighLPS == T){
  #   ttmp = quantile(LPS_score, probs = 1-ratioPrior)
  #   rtmp = quantile(LPS_score, probs = ratioPrior)
  #   ttmp = LPS_score[which(LPS_score >= ttmp)]
  #   rtmp = LPS_score[which(LPS_score < rtmp)]
  # } else{
  #   ttmp = quantile(LPS_score, probs = ratioPrior)
  #   rtmp = quantile(LPS_score, probs = 1- ratioPrior)
  #   ttmp = LPS_score[which(LPS_score <= ttmp)]
  #   rtmp = LPS_score[which(LPS_score > rtmp)]
  # }
  # ### notice that both ttmp and rtmp are the two ends, meaning that length(rtmp + ttmp) < length(all scores), it will be equal if ratioRrior is 0.5
  # 
  # ### make changes on 20190425 after comparisons on 20190424
  # ### the following 3 lines are actually consistent to the feature level setting for ratio prior usage
  # ###  more importantly, the following 3 lines are also consistent to JCO paper where the original PRPS was defined and used with ratio prior
  tmpcut = quantile(LPS_score, probs = 1-ratioPrior)
  ttmp = LPS_score[which(LPS_score >= tmpcut)]
  rtmp = LPS_score[which(LPS_score < tmpcut)]
  
  testLPSmean = mean(ttmp)
  refLPSmean = mean(rtmp)
  testLPSsd = sd(ttmp)
  refLPSsd = sd(rtmp)
  
  LPS_prob_test = getProb(LPS_score, groupMeans = c(testLPSmean, refLPSmean), groupSds = c(testLPSsd, refLPSsd))
  
  LPS_prob_ref = getProb(LPS_score, groupMeans = c(refLPSmean, testLPSmean), groupSds = c(refLPSsd, testLPSsd))
  
  LPS_class = rep("UNCLASS",length(LPS_score))
  LPS_class[which(LPS_prob_test >= classProbCut)] = testGroup
  LPS_class[which(LPS_prob_ref >= classProbCut)] = refGroup
  
  LPS_score = data.frame(LPS_score)
  LPS_test = cbind(LPS_score, LPS_class, LPS_prob_test, LPS_prob_ref, stringsAsFactors =F)
  
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
