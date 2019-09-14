#' @export
LPSSLextension = function(LPSSLObj, newdat, standardization=FALSE, classProbCut = 0.9,  imputeNA = FALSE, byrow = TRUE, imputeValue =c("median","mean")){
  imputeValue = imputeValue[1]
  
  if(is.null(LPSSLObj)){print("Please input your LPS self learning object")}
  LPS_pars = LPSSLObj$LPS_pars
  weights = LPS_pars$weights
  
  ## imputee NA if imputeNA is true
  if(imputeNA){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for LPS approach, it does not require standardization, however, if standardization = TRUE, do the standardization
  if(standardization){newdat = standardize(newdat)}
  
  # if NA is not imputed, remove it from LPS score calculation
  LPS_score = apply(data.matrix(newdat[rownames(weights),]), 2, getLPSscore, coefs = weights[,1])
  
  # for LPS, 0 is a NOT natural cutoff for two group classification
  # in order to get classification, need to get two groups' LPS mean and sd
  
  ### get group info
  testres = LPSSLObj$LPS_test
  testres = testres[order(testres[,1], decreasing = T),]
  testGroup = testres[1,2]
  
  refGroup = setdiff(unique(LPSSLObj$LPS_test$PRPS_class),c(testGroup, "UNCLASS"))
  
  # for LPS, 0 is a NOT natural cutoff for two group classification
  # in order to get classification, need to get two groups' LPS mean and sd
  
  testLPSmean = LPS_pars$meansds[1]
  refLPSmean = LPS_pars$meansds[2]
  testLPSsd = LPS_pars$meansds[3]
  refLPSsd = LPS_pars$meansds[4]
 
  LPS_prob_test = getProb(LPS_score, groupMeans = c(testLPSmean, refLPSmean), groupSds = c(testLPSsd, refLPSsd))
  
  LPS_prob_ref = getProb(LPS_score, groupMeans = c(refLPSmean, testLPSmean), groupSds = c(refLPSsd, testLPSsd))
  
  LPS_class = rep("UNCLASS",length(LPS_score))
  LPS_class[which(LPS_prob_test >= classProbCut)] = testGroup
  LPS_class[which(LPS_prob_ref >= classProbCut)] = refGroup
  
  LPS_score = data.frame(LPS_score)
  LPS_test = cbind(LPS_score, LPS_class, LPS_prob_test, LPS_prob_ref, stringsAsFactors =F)
  
  return(LPS_test)
  
}
