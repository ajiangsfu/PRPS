
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
#' \item{PS_test}{a data frame of PS score, and classification based on 0 as natural cutoff}
#' @keywords PS 
#' @author Aixiang Jiang
#' @references 
#' TR Golub, DK Slonim, P Tamayo, C Huard, M Gaasenbeek, JP Mesirov, H Coller, ML Loh, JR Downing, MA Caligiuri, et al.
#' Molecular classification of cancer: class discovery and class prediction by gene expression monitoring
#' Science, 286 (1999), pp. 531-537
#' @export
PSSLwithWeightsPrior = function(newdat, weights, ratioPrior = 1/3, standardization=FALSE,  PShighGroup = "PShigh", 
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
  
  #### 20190418: add classification calls based on Empirical Bayesian probabilities with EM (Expectation-Maximization)
  emcut = AdaptGauss::EMGauss(PS_score, K = 2,fast=TRUE, MaxNumberofIterations = EMmaxRuns)
  
  ### add a plot, hist with two distribution lines, do not need to save, just plot it
  hist(PS_score, prob = TRUE, breaks = breaks)
  curve(emcut$Weights[1]*dnorm(x, mean=emcut$Means[1], sd=emcut$SDs[1]), 
        col="red", lwd=2, add=TRUE, yaxt="n")
  curve(emcut$Weights[2]*dnorm(x, mean=emcut$Means[2], sd=emcut$SDs[2]), 
        col="green", lwd=2, add=TRUE, yaxt="n")
  abline(v=0, col = "red")
  
  PS_score = data.frame(PS_score)
  PS_test = cbind(PS_score, PS_class0,stringsAsFactors =F)
 
  
  #######################################################
  #### I should check all PS SL related functions, this is the 1st one, stop here on 0717
  ######################################################
  
  
  
  
  weights = data.frame(weights)
  
  PS_pars =  list(weights, meansds = tmp, traitsmeansds = cbind(allmeans, allsds))
  names(PS_pars) = c("weights","meansds","traitsmeansds")
  
  outs = list(PS_pars, PS_test)
  names(outs) = c("PS_pars","PS_test")
  
  return(outs)
 
}
