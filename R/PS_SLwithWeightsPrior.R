
#' PS score calculation and binary classification for a testing data set without training data set, but with selected feature weights
#' and a ratio prior
#' @description This is the function to calculate PS (Prediction Strength) scores for a testing data set 
#' without training data set. However, we do need selected feature list with their weights. In addtion, we need mean of two group means for each
#' selected feature that we need to apply a ratio prior with default as 1/3 to achieve. Once we have PS scores, we will use the 
#' theoretic natual cutoff 0 of PS scores to make classification calls, as well as classification calls based on Empirical Bayesian probabilities
#' with EM (Expectation-Maximization) and with a ratio prior.
#' @details  This is the function to calculate PS scores and make classification for a testing data set with known weights and a ratio prior of choice. 
#' The range of PS scores is [-1,1]
#' PS calculation is based on Golub 1999 with the following three steps in this function:
#' a) apply "standardize" function to standardize input data matrix for each feature;
#' c) get mean of group means for each selected feature;
#' d) use "apply" function to get PS scores for all samples with "getPS1sample", the formula is:
#'   \eqn{PS = (V_win − V_lose)/(V_win + V_lose)}
#' Here, where V_win and V_lose are the vote totals for the winning and losing features/traits for a given sample
#' 
#' If NAs are not imputed, they are ignored for all of the steps.
#' 
#' In addtion, we design to use EM algorithm to get PS score mean and sd calculation for the two groups assuming that PS score
#' is a mixture of two normal distributions.
#' 
#' When calculate a Empirical Bayes' probability, we first calcualte probability that a sample belongs to either group,
#' and then use the following formula to get Empirical Bayes' probability:
#' \eqn{prob(x) = p_test(x)/(p_test(x) + p_ref(x))}
#' Here prob(x) is the Empirical Bayes' probability of a given sample, p_test(x) is the probability that a given sample
#' belongs to the test group, p_ref(x) is the probability that a given sample belongs to the reference group.
#' Notice that the test and reference group is just the relative grouping, in fact, for this step, we often need
#' to calculate Empirical Bayes' probabilities for a given sample from two different standing points.
#' 
#' @param newdat a new data matrix or data frame, columns for samples and rows for features
#' @param weights a numeric vector with selected features (as names of the vector) and their weights
#' @param ratioPrior a prior ratio to indicate the ratio of test group over reference group 
#'   regarding mean and sd calculation for each selected feature
#' @param standardization a logic variable to indicate if standardization is needed before classification 
#'  score calculation
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
#'  calculation, the 2nd item is a numeric vector containing PS mean of two group means，and the 3rd item 
#'  is a data frame contains mean and sd for each group and for each selected feature}
#' \item{PS_test}{a data frame of PS score, Empirical Bayesian probabilites for both groups and classification based on EM
#' grouping, classification based on 0 as natural cutoff,  Empirical Bayesian probabilites for both groups and classification based on
#' ratio prior. Notice that the ratio prior is used for two ends with the same prior in this very last step}
#' @keywords PS 
#' @author Aixiang Jiang
#' @references 
#' TR Golub, DK Slonim, P Tamayo, C Huard, M Gaasenbeek, JP Mesirov, H Coller, ML Loh, JR Downing, MA Caligiuri, et al.
#' Molecular classification of cancer: class discovery and class prediction by gene expression monitoring
#' Science, 286 (1999), pp. 531-537
#' @export
PS_SLwithWeightsPrior = function(newdat, weights, ratioPrior = 1/3, standardization=FALSE,  PShighGroup = "PShigh", 
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
  
  PS_score = weightedLogProbClass(newdat=newdat, topTraits=names(weights), weights=weights,
                                 classMeans=meansd[,1:2], classSds = meansd[,3:4])

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
  
  PS_prob1 = getProb(PS_score, groupMeans = emcut$Means, groupSds = emcut$SDs)
  PS_prob2 = getProb(PS_score, groupMeans = rev(emcut$Means), groupSds = rev(emcut$SDs))
  
  if(emcut$Means[1]> emcut$Means[2]){
    name1 = PShighGroup
    name2 = PSlowGroup
  }else{
    name1 = PSlowGroup
    name2 = PShighGroup
  }
  
  PS_class = rep("UNCLASS",length(PS_score))
  PS_class[which(PS_prob1 >= classProbCut)] = name1
  PS_class[which(PS_prob2 >= classProbCut)] = name2
  
  PS_score = data.frame(PS_score)
  PS_test = cbind(PS_score, PS_class, PS_prob1, PS_prob2, PS_class0,stringsAsFactors =F)
  
  #### change on 20190418
  #### return a list withn 3 items: , PS_par, PS_test, stable_samples (sample name with grouping info)
  #### it is better to get consistent PS_par as in PStraining
  #### therefore: think again, maybe I can combine this function with PStrain?, maybe not, since the output is not the same
  ####     and, the stable classification is not the truth, but it is the optimal based on given data
  #### more questions
  #### 1) should I keep the result for 0 cutoff? ->  maybe not...
  #### 2) should I use the group mean and sd PS scores of the stable samples for Bayes prob and classification
  ####        as I did in PStraining, or should I keep the current method with EM + Bayes? -> I will keep current approach
  
  tmp = c(emcut$Means, emcut$SDs)
  names(tmp) = c("testPSmean","refPSmean","testPSsd","refPSsd")
  
  
  ##############################################################
  
  #### add prior based classification as well. Notice that one prior is used for both ends in this function
  #### while in previous step, we use all data for feature level
  ttmp = quantile(PS_score, probs = 1-ratioPrior)
  rtmp = quantile(PS_score, probs = ratioPrior)
  ttmp = PS_score[which(PS_score >= ttmp)]
  rtmp = PS_score[which(PS_score < rtmp)]

  testPSmean = mean(ttmp)
  refPSmean = mean(rtmp)
  testPSsd = sd(ttmp)
  refPSsd = sd(rtmp)
  
  PS_prob_high = getProb(PS_score, groupMeans = c(testPSmean, refPSmean), groupSds = c(testPSsd, refPSsd))
  
  PS_prob_low = getProb(PS_score, groupMeans = c(refPSmean, testPSmean), groupSds = c(refPSsd, testPSsd))
  
  PS_class_prior = rep("UNCLASS",length(PS_score))
  PS_class[which(PS_prob_high >= classProbCut)] = PShighGroup
  PS_class[which(PS_prob_low >= classProbCut)] = PSlowGroup
  
  #####################################################
  PS_test = cbind(PS_test, PS_class_prior, PS_prob_high, PS_prob_low,  stringsAsFactors =F)
  
  weights = data.frame(weights)
  
  PS_pars =  list(weights, meansds = tmp, traitsmeansds = cbind(allmeans, allsds))
  names(PS_pars) = c("weights","meansds","traitsmeansds")
  
  outs = list(PS_pars, PS_test)
  names(outs) = c("PS_pars","PS_test")
  
  return(outs)
 
}
