
#' Self learning binary classification with selected features and their weights
#' @description This function is to calculate PS (Prediction Strength) scores and make binary classificationcalls 
#' for a testing data set without PS training object. This function involves a self learning 
#' process with weights + priors + EM + Bayes but no need to input a group ratio prior. However, we do need selected feature list
#' with the feature weights, and we use self learning method to estimate parameters of the two groups for each selected feature 
#' in order to calculate PS scores. 
#' @details This function is trying to get reasonable PS based classification without training data set, but with 
#' selected features and their weights. It combines ideas from PSSLwithWeightsPrior and
#' PSSLwithWeightsEM, which we call it as self learning algorithm, the actual steps are as following:
#' 1) assume that we have a pool for group ratio priors such as: seq(0.05, 0.95, by = 0.05), this will give us 19 ratio priors
#' 2) for each prior in 1), call PSSLwithWeightsPrior to achieve PS scores
#' 3) apply EM on PS scores from 2) with Mclust, which includes 2 group classification
#' 4) use the samples that are always in the same groups to get group means and sds for each feature and for each group
#' 5) calculate PS scores
#' 6) Once we have PS scores, we could use the theoretic natual cutoff 0 to make classification calls, which may or may not appropriate. 
#' Alternatively, we can also apply EM to calcualate mean and sd for the
#' two groups assuming that PS score is a mixture of two normal distributions, followed by Empirical Bayes' probability 
#' calculation and final binary classification calls.
#' PS calculation is based on Ennishi 2018, its formula is:
#' \eqn{PS(X_i) = \sum (|a_j| log(P1(x_ij)/P0(x_ij)))}
#' Here, a_j represents the jth selected feature weights, and x_ij is the corresponding feature value
#'  for the ith sample, 
#' P1 and P0 are the probabilities that the ith sample belongs to two different group.
#' However, in order to calculate P1 and P0, we need to have two group mean and sd for each selected feature. Although there are multiple way
#' to obtain these values, in this function, we design to use EM algorithm to achieve group mean and sd assuming that each selected feature
#' is a mixture of two normal distributions. 
#' After we have PS scores, we provide classification calls based on theoretic 0 cutoff and based on probability, which may or may not appropriate
#' for this function since it is based on self learning. To calculate a Empirical Bayes' 
#' probability to make classification calls, we also need to apply EM to get PS score mean and sd for two groups. 
#' After that, we can calcualte probability that a sample belongs to either group,
#' and then use the following formula to get Empirical Bayes' probability:
#' \eqn{prob(x) = p_test(x)/(p_test(x) + p_ref(x))}
#' Here prob(x) is the Empirical Bayes' probability of a given sample, p_test(x) is the probability that a given sample
#' belongs to the test group, p_ref(x) is the probability that a given sample belongs to the reference group.
#' Notice that the test and reference group is just the relative grouping, in fact, for this step, we often need
#'  to calculate Empirical Bayes' probabilities for a given sample from two different standing points.
#' @param newdat a input data matrix or data frame, columns for samples and rows for features
#' @param weights a numeric vector with selected features (as names of the vector) and their weights
#' @param standardization a logic variable to indicate if standardization is needed before classification 
#'  score calculation
#' @param classProbCut a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
#'  often used values are 0.8 and 0.9, default value is 0.8. Only one value is used for both groups, 
#'  the samples that are not included in either group will be assigned as UNCLASS
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
#' \item{PS_pars}{a list of 4 items, the 1st item is a data frame with weights of each selected features for PS
#'  calculation, the 2nd item is a numeric vector containing PS mean and sd for two groups，the 3rd item is a data frame contains mean and sd
#'   for each group and for each selected feature based on stable classes, and the last item is a numeric vector containing PS mean and sd 
#'   for two groups based on EM}
#' \item{PS_test}{a data frame of PS score, classification and two groups' Empirical Bayesian probabilites based on stable classes, 
#' classification with natural 0 cutoff, classification and two groups' Empirical Bayesian probabilites based on EM}
#' @keywords PS EM 
#' @author Aixiang Jiang
#' @references 
#' #' Golub TR, Slonim DK, Tamayo P, Huard C, Gaasenbeek M, Mesirov JP, et al. Molecular classification of cancer: 
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

PSstableSLwithWeights = function(newdat, weights, standardization=FALSE, classProbCut = 0.9, PShighGroup = "PShigh", 
                    PSlowGroup = "PSlow", breaks = 50, EMmaxRuns = 50, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  require(mclust)
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
  
  ### now, I need to work on group mean and sd for each feature
  datgrp1 = newdat[,grp1]
  datgrp2 = newdat[,grp2]
  
  means1 = rowMeans(datgrp1) 
  means2 = rowMeans(datgrp2)
  
  sds1 = apply(datgrp1,1,sd)
  sds2 = apply(datgrp2,1,sd)
  
  #### I need to determine the direction for allmeans, and allsds
  #### maybe I can use the middle prior to make the decision
  xx=median(rps)
  tmp = PSSLwithWeightsPrior(newdat=newdat, weights=weights, ratioPrior = xx, PShighGroup = PShighGroup, PSlowGroup = PSlowGroup)
  mcls = mclust::Mclust(tmp$PS_test$PS_score, G=2)
  
  mmean = mcls$parameters$mean
  if(mmean[2] > mmean[1]){
    allmeans = cbind(means2, means1)
    allsds = cbind(sds2, sds1)
    tt = grp1
    grp1 = grp2
    grp2 = tt
  }else{
    allmeans = cbind(means1, means2)
    allsds = cbind(sds1, sds2)
  }
  
  PS_score = weightedLogProbClass(newdat = newdat, topTraits=names(weights), weights=weights,
                                    classMeans = allmeans, classSds = allsds)
  
  #### 20190503, call plotHistEM 
  emsearch = plotHistEM(PS_score, G = 2:4, breaks = breaks, EMmaxRuns = EMmaxRuns, scoreName = "PS_score")
  bestG = emsearch$bestG
  emcut = emsearch$emcut

  ### no matter how many bestG, only keep the 1st and last one for the following
  gmeans = c(emcut$Means[1], emcut$Means[bestG])
  gsds = c(emcut$SDs[1], emcut$SDs[bestG])
  
  if(gmeans[1]> gmeans[2]){
    name1 = PShighGroup
    name2 = PSlowGroup
    scoreMeanSds_EM = c(gmeans, gsds)
  }else{
    name1 = PSlowGroup
    name2 = PShighGroup
    scoreMeanSds_EM = c(rev(gmeans), rev(gsds))
  }
  
  PS_prob1_EM = getProb(PS_score, groupMeans = scoreMeanSds_EM[1:2], groupSds = scoreMeanSds_EM[3:4])
  PS_prob2_EM = getProb(PS_score, groupMeans = rev(scoreMeanSds_EM[1:2]), groupSds = rev(scoreMeanSds_EM[3:4]))
  
  PS_class_EM = rep("UNCLASS",length(PS_score))
  PS_class_EM[which(PS_prob1_EM >= classProbCut)] = PShighGroup
  PS_class_EM[which(PS_prob2_EM >= classProbCut)] = PSlowGroup

  names(scoreMeanSds_EM) = c("testPSmean_EM","refPSmean_EM","testPSsd_EM","refPSsd_EM")
  
  ##### make changes on 20190502
  ### use my stable two groups across priors to do the last step group mean and sd as well, which might be more reasonable instead of forcing two groups?
  #### also, still keep the histogram plus EM lines
  
  grp1score = PS_score[grp1]
  grp2score = PS_score[grp2]
  
  scoreMeans = c(mean(grp1score), mean(grp2score))
  scoreSds = c(sd(grp1score), sd(grp2score))
  
  PS_prob1 = getProb(PS_score, groupMeans = scoreMeans, groupSds = scoreSds)
  PS_prob2 = getProb(PS_score, groupMeans = rev(scoreMeans), groupSds = rev(scoreSds))
  
  PS_class = rep("UNCLASS",length(PS_score))
  PS_class[which(PS_prob1 >= classProbCut)] = PShighGroup
  PS_class[which(PS_prob2 >= classProbCut)] = PSlowGroup
  
  PS_class0 = ifelse(PS_score > 0,  PShighGroup, PSlowGroup)
  
  PS_score = data.frame(PS_score)
  PS_test = cbind(PS_score, PS_class, PS_prob1, PS_prob2, PS_class0,stringsAsFactors =F)
  
  ### 20190503: add the stable classification 
  PS_test$stable_class = "UNCLASS"
  PS_test[grp1,"stable_class"] = PShighGroup
  PS_test[grp2,"stable_class"] = PSlowGroup
  
  ###################### remove for now, update description later ############
  #### finally, add EM class into
  #PS_test = cbind(PS_test, PS_class_EM, PS_prob1_EM, PS_prob2_EM, stringsAsFactors =F)
  
  weights = data.frame(weights)
  scoreMeanSds = c(scoreMeans, scoreSds)
  
  #### remove EM part
  #PS_pars =  list(weights, meansds = scoreMeanSds, traitsmeansds = cbind(allmeans, allsds),scoreMeanSds_EM)
  #names(PS_pars) = c("weights","meansds","traitsmeansds", "meanSds_EM")
  
  PS_pars =  list(weights, meansds = scoreMeanSds, traitsmeansds = cbind(allmeans, allsds))
  names(PS_pars) = c("weights","meansds","traitsmeansds")
  
  outs = list(PS_pars, PS_test)
  names(outs) = c("PS_pars","PS_test")
  return(outs)
  
}
