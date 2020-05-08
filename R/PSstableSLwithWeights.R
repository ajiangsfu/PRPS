
#' PS stable self-training
#' @description This function is to calculate PS (Prediction Strength) scores and make binary classification calls 
#' for a testing data set without PS training object. It involves a self-training process with given features and their weights.
#' @details This function is trying to get reasonable PS based classification without training data set, but with 
#' selected features and their weights. The actual steps are as following:
#' 1) assume that we have a pool for group ratio priors such as seq(0.05, 0.95, by = 0.05) for default ratioRange = c(0.05, 0.95)
#' 2) With given features and their weights
#'    a) for each prior in 1), call PSSLwithWeightsPrior with given features and weights to achieve PS scores
#'       apply EM on PS scores with Mclust, get 2 group classification
#'    b) define the samples that are always in the same classes across searching range as stable classes
#' 3) repeat step 2) but this time with opposite signs in the given weights, result in another set of stable classes
#' 4) get final stable classes that are common in 2) and 3)
#' 5) use final stable classes to get group means and sds for each feature and for each group
#' 5) calculate PS scores
#' 6) Once we have PS scores, we could use the theoretic natual cutoff 0 to make classification calls, which may or may not appropriate. 
#' Alternatively, with two groups based on stable classes assuming that PS score is a mixture of two normal distributions, 
#' we can get Empirical Bayesian probability and make calls
#' @param newdat a input data matrix or data frame, columns for samples and rows for features
#' @param weights a numeric vector with selected features (as names of the vector) and their weights
#' @param plotName a pdf file name with full path and is ended with ".pdf", which is used to save multiple pages 
#'  of PS histgrams with distribution densities. Default value us NULL, no plot is saved.
#' @param ratioRange a numeric vector with two numbers, which indicates ratio search range. The default is
#'  c(0.1, 0.9)for the current function. If your classification is very
#'  unbalanced such as one group is much smaller than the other, and/or sample variation is quite big,
#'  and/or classification results are far away from what you expect, you might want to change the default values.
#'  c(0.15, 0.85) is recommended as an alternative setting other than default. In an extreme rare situation, c(0.4, 0,6) could a good try.
#' @param classProbCut a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
#'  often used values are 0.8 and 0.9, default value is 0.9. Only one value is used for both groups, 
#'  the samples that are not included in either group will be assigned as UNCLASS
#' @param PShighGroup a string to indicate group name with high PS score
#' @param PSlowGroup a string to indicate group name with low PS score
#' @param breaks a integer to indicate number of bins in histogram, default is 50
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
#' @import mclust
#' @export

PSstableSLwithWeights = function(newdat, weights, plotName = NULL, ratioRange = c(0.1, 0.9), classProbCut = 0.9, PShighGroup = "PShigh", PSlowGroup = "PSlow",
                                 breaks = 50, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
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

  ## change on 20190918
  rps = seq(ratioRange[1], ratioRange[2], by = 0.05)
  
  if(!is.null(plotName)){
    pdf(plotName)
  }
  
  rpsres = sapply(rps, FUN = function(xx){
    tmp = PSSLwithWeightsPrior(newdat=newdat, weights=weights, ratioPrior = xx, PShighGroup = PShighGroup, PSlowGroup = PSlowGroup)
    ### note on 2019-05-08: realize that the EM related output from PSSLwithWeightsPrior is actually not used at all in mcls
    ###                     only the prior based scores are used in the following step
    ###                     question: should I not call PSSLwithWeightsPrior at all? should I just copy the prior PS score part?
    ###                               or should I use the EM results from output of PSSLwithWeightsPrior?
    mcls = mclust::Mclust(tmp$PS_test$PS_score, G=2, warn = TRUE)
    ######  add plot EM step back for this local function on 20191206 ######################
    emsearch = plotHistEM(tmp$PS_test$PS_score, G = 2, breaks = breaks, 
                          scoreName = paste("PS_score with rho = ", xx, sep=""))
    #########################################################################
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
  
  ########### changes on Oct 31, 2019 #######################
  orpsres = sapply(rps, FUN = function(xx){
    tmp = PSSLwithWeightsPrior(newdat=newdat, weights=-weights, ratioPrior = xx, PShighGroup = PSlowGroup, PSlowGroup = PShighGroup)
    mcls = mclust::Mclust(tmp$PS_test$PS_score, G=2, warn = TRUE)
    ######  add plot EM step back for this local function on 20191206 ######################
    emsearch = plotHistEM(tmp$PS_test$PS_score, G = 2, breaks = breaks, 
                          scoreName = paste("Reverse weight PS_score with rho = ", xx, sep=""))
    #########################################################################
    return(mcls$classification)
  })
  
  rownames(orpsres) = colnames(newdat)
  ### the classification is coded with 1 and 2
  ### in order to find samples that are always 1 and always 2, what should I do?
  ### for always 1, if I extract 1 for each cell, then row sum for the sample should be 0
  ### for always 2, if I extract 2 for each cell, then row sum for the sample should be 0
  
  ores1 = orpsres - 1
  ores2 = orpsres - 2
  
  ogrp1 = rownames(ores1[which(rowSums(ores1) == 0),])
  ogrp2 = rownames(ores1[which(rowSums(ores2) == 0),])
  
  grp1 = intersect(grp1, ogrp2)
  grp2 = intersect(grp2, ogrp1)
  
  #########################################################
  
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
  
  ###################################################################
  # #### 20190503, call plotHistEM , make changes of scoreName on 20191206 for this local function
  emsearch = plotHistEM(PS_score, G = 2, breaks = breaks, scoreName = "Final PS_score")
  ####################################################################
  
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
  weights = data.frame(weights)
  scorePars = c(testPSmean, refPSmean, testPSsd, refPSsd)
  names(scorePars) = c("highPSmean","lowPSmean","highPSsd","lowPSsd")
  PS_pars =  list(weights, meansds = scorePars, traitsmeans = mean_2means)
  names(PS_pars) = c("weights","meansds","traitsmeans")
  
  if(!is.null(plotName)){
    dev.off()
  }
  
  outs = list(PS_pars, PS_test)
  names(outs) = c("PS_pars","PS_test")
  
  return(outs)
  
}
