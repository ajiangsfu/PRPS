
#' Self learning binary classification with selected features and their weights
#' @description This function is to calculate PS (Prediction Strength) scores and make binary classification calls 
#' for a testing data set without PS training object. This function involves a self learning 
#' process with weights + priors + EM + Bayes but no need to input a group ratio prior. However, we do need selected feature list
#' with the feature weights, and we use self learning method to estimate parameters of the two groups for each selected feature 
#' in order to calculate PS scores. 
#' @details This function is trying to get reasonable PS based classification without training data set, but with 
#' selected features and their weights. The actual steps are as following:
#' 1) assume that we have a pool for group ratio priors such as: seq(0.05, 0.95, by = 0.05), this will give us 19 ratio priors
#' 2) for each prior in 1), call PSSLwithWeightsPrior to achieve PS scores
#' 3) apply EM on PS scores from 2) with Mclust, which includes 2 group classification
#' 4) use the samples that are always in the same groups to get group means for each group and mean of these two means for each feature 
#' 5) calculate PS scores
#' 6) Once we have PS scores, we use the theoretic natual cutoff 0 to make classification calls
#' @param newdat a input data matrix or data frame, columns for samples and rows for features
#' @param weights a numeric vector with selected features (as names of the vector) and their weights
#' @param PShighGroup a string to indicate group name with high PS score
#' @param PSlowGroup a string to indicate group name with low PS score
#' @param breaks a integer to indicate number of bins in histogram, default is 50
#' @param EMmaxRuns number of Iterations for EM searching; default=50
#' @param imputeNA a logic variable to indicate if NA imputation is needed, if it is TRUE, NA imputation is 
#'  processed before any other steps, the default is FALSE
#' @param byrow a logic variable to indicate direction for imputation, default is TRUE, 
#'  which will use the row data for imputation
#' @param imputeValue a character variable to indicate which value to be used to replace NA, default is "median", 
#'  the median value of the chose direction with "byrow" data to be used
#' @return A list with two items is returned: PS parameters for selected features, PS scores and classifications for the given samples.
#' \item{PS_pars}{a list of 3 items, the 1st item is a data frame with weights of each selected features for PS
#'  calculation, the 2nd item is a numeric vector containing PS mean and sd for two groups，the 3rd item is a data frame contains 
#'  group means for each group and mean of these two means for each feature based on stable classes}
#' \item{PS_test}{a data frame of PS score and classification with natural 0 cutoff}
#' @keywords PS EM self-learning
#' @author Aixiang Jiang
#' @references 
#' Golub TR, Slonim DK, Tamayo P, Huard C, Gaasenbeek M, Mesirov JP, et al. Molecular classification of cancer: 
#' class discovery and class prediction by gene expression monitoring. Science. 1999;286:531–7
#' 
#' Ultsch, A., Thrun, M.C., Hansen-Goos, O., Loetsch, J.: Identification of Molecular Fingerprints
#' in Human Heat Pain Thresholds by Use of an Interactive Mixture Model R Toolbox(AdaptGauss),
#' International Journal of Molecular Sciences, doi:10.3390/ijms161025897, 2015.
#' 
#' Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and 
#' density estimation using Gaussian finite mixture models, The R Journal, 8/1, pp. 205-233.
#' 

#' @export

PSstableSLwithWeights = function(newdat, weights, PShighGroup = "PShigh", PSlowGroup = "PSlow", breaks = 50, EMmaxRuns = 50, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  require(mclust)
  imputeValue = imputeValue[1]
  ## imputee NA if imputeNA is true
  if(imputeNA){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for PS approach, standardization is required
  newdat = standardize(newdat)
  
  ### in case that some features in weight list are actually not in data, do the following
  tmp = intersect(names(weights), rownames(newdat))
  weights = weights[tmp]
  newdat = newdat[tmp,]
  
  rps = seq(0.05, 0.95, by = 0.05)
  
  rpsres = sapply(rps, FUN = function(xx){
    tmp = PSSLwithWeightsPrior(newdat=newdat, weights=weights, ratioPrior = xx, PShighGroup = PShighGroup, PSlowGroup = PSlowGroup)
    ### note on 2019-05-08: realize that the EM related output from PSSLwithWeightsPrior is actually not used at all in mcls
    ###                     only the prior based scores are used in the following step
    ###                     question: should I not call PSSLwithWeightsPrior at all? should I just copy the prior PS score part?
    ###                               or should I use the EM results from output of PSSLwithWeightsPrior?
    mcls = mclust::Mclust(tmp$PS_test$PS_score, G=2)
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
  
  ### now, I need to work on group mean and mean of two means for each feature
  datgrp1 = newdat[,grp1]
  datgrp2 = newdat[,grp2]
  
  means1 = rowMeans(datgrp1, na.rm = T) 
  means2 = rowMeans(datgrp2, na.rm = T)
  
  mean_2means = rowMeans(cbind(means1, means2))  ### this is a vector of mean of group means for each feature
  mean_2means = cbind(mean_2means, means1, means2)
  colnames(mean_2means) = c("meanOfGroupMeans","groupMean1","groupMean2")

  PS_pars = cbind( mean_2means[,1],weights, mean_2means[,-1])
  colnames(PS_pars)[1] = colnames(mean_2means)[1]
  
  # get PS scores for all samples
  PS_score = apply(newdat[names(weights),], 2, getPS1sample, PSpars = PS_pars[,1:2])
  
  PS_class0 = ifelse(PS_score > 0,  PShighGroup,  PSlowGroup)  
  
  
  #### 20190503, call plotHistEM 
  emsearch = plotHistEM(PS_score, G = 2:4, breaks = breaks, EMmaxRuns = EMmaxRuns, scoreName = "PS_score")
  # bestG = emsearch$bestG
  # emcut = emsearch$emcut
  # 
  # ### no matter how many bestG, only keep the 1st and last one for the following
  # gmeans = c(emcut$Means[1], emcut$Means[bestG])
  # gsds = c(emcut$SDs[1], emcut$SDs[bestG])
  # 
  # if(gmeans[1]> gmeans[2]){
  #   name1 = PShighGroup
  #   name2 = PSlowGroup
  #   scoreMeanSds_EM = c(gmeans, gsds)
  # }else{
  #   name1 = PSlowGroup
  #   name2 = PShighGroup
  #   scoreMeanSds_EM = c(rev(gmeans), rev(gsds))
  # }
  # 
  # PS_prob1_EM = getProb(PS_score, groupMeans = scoreMeanSds_EM[1:2], groupSds = scoreMeanSds_EM[3:4])
  # PS_prob2_EM = getProb(PS_score, groupMeans = rev(scoreMeanSds_EM[1:2]), groupSds = rev(scoreMeanSds_EM[3:4]))
  # 
  # PS_class_EM = rep("UNCLASS",length(PS_score))
  # PS_class_EM[which(PS_prob1_EM >= classProbCut)] = PShighGroup
  # PS_class_EM[which(PS_prob2_EM >= classProbCut)] = PSlowGroup
  # 
  # names(scoreMeanSds_EM) = c("testPSmean_EM","refPSmean_EM","testPSsd_EM","refPSsd_EM")
  # 
  ##### make changes on 20190502
  ### use my stable two groups across priors to do the last step group mean and sd as well, which might be more reasonable instead of forcing two groups?
  #### also, still keep the histogram plus EM lines
  
  #### 20190913, use the 0 theoretical cutoff to get two groups, which are used for empirial Bayesian prob calculation
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
  
  score1 = PS_score[grp1]
  score2 = PS_score[grp2]
  
  sig1 = PShighGroup
  sig2 = PSlowGroup
  
  if(mean(score1) < mean(score2)){
    sig2 = PShighGroup
    sig1 = PSlowGroup
  }
  
  PS_score = data.frame(PS_score)
  PS_test = cbind(PS_score, PS_class, PS_prob1, PS_prob2, PS_class0,stringsAsFactors =F)
  
  ### 20190503: add the stable classification 
  PS_test$stable_class = "UNCLASS"
  PS_test[grp1,"stable_class"] = sig1
  PS_test[grp2,"stable_class"] = sig2
  
  ###################### remove for now, update description later ############
  #### finally, add EM class into
  #PS_test = cbind(PS_test, PS_class_EM, PS_prob1_EM, PS_prob2_EM, stringsAsFactors =F)

  PS_pars =  list(weights, meansds = c(scoreMeans, scoreSds), traitsmeans = mean_2means)
  names(PS_pars) = c("weights","meansds","traitsmeans")
  
  outs = list(PS_pars, PS_test)
  names(outs) = c("PS_pars","PS_test")
  
  return(outs)
  
}
