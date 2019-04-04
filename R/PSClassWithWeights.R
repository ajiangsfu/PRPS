
#' PS score calculation and binary classification for a testing data set without PStraining output object, but with selected feature weights
#' @description This is the function to calculate PS (Prediction Strength) scores for a testing data set 
#' without PS training object. However, we do need selected feature list with their weights, and mean of two group means for each
#' selected feature that we need to apply EM (Expectation-Maximization) to achieve. Once we have PS scores, we could use the 
#' theoretic natual cutoff 0 to make classification calls. Alternatively, we can also apply EM to calcualate mean and sd for the two groups 
#' assuming that PS score is a mixture of two normal distributions, followed by Empirical Bayes' probability calculation and final binary classification calls.
#' @details  This is the function to calculate PS scores and make classification based on
#' Empirical Bayesian probabilities for a testing new data set. Before we do anything, we need to do standardization,
#' which is requred before PS calculation. 
#' PS calculation is based on:
#' \eqn{PS = (V_win − V_lose)/(V_win + V_lose)}
#' Here, where V_win and V_lose are the vote totals for the winning and losing features/traits for a given sample
#' tp get these vote, however, we need to know mean of two group means for each selected features. In oder to get achieve them,
#' there are multiple ways to do, in this function, we apply EM algorithm on each input feature data to obtain, and further
#' to get PS score for each sample. Once we have PS scores, we can use PS natual therotic cutoff 0 to make classification calls.
#' Alternatively, we can make calls based on Empirical Bayes' probability. To do that, 
#' we also need to apply EM to get PS score mean and sd for two groups. 
#' After that, we can calcualte probability that a sample belongs to either group,
#' and then use the following formula to get Empirical Bayes' probability:
#' \eqn{prob(x) = p_test(x)/(p_test(x) + p_ref(x))}
#' Here prob(x) is the Empirical Bayes' probability of a given sample, p_test(x) is the probability that a given sample
#' belongs to the test group, p_ref(x) is the probability that a given sample belongs to the reference group.
#' Notice that the test and reference group is just the relative grouping, in fact, for this step, we often need
#'  to calculate Empirical Bayes' probabilities for a given sample from two different standing points.
#' @param newdat a new data matrix or data frame, columns for samples and rows for features
#' @param weights a numeric vector with selected features (as names of the vector) and their weights
#' @param classProbCut a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
#'  often used values are 0.8 and 0.9, default value is 0.8. Only one value is used for both groups, 
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
#' @return A data frame with PS score, Empirical Bayesian probabilites for two groups and classification
#' @keywords PS EM 
#' @author Aixiang Jiang
#' @references 
#' TR Golub, DK Slonim, P Tamayo, C Huard, M Gaasenbeek, JP Mesirov, H Coller, ML Loh, JR Downing, MA Caligiuri, et al.
#' Molecular classification of cancer: class discovery and class prediction by gene expression monitoring
#' Science, 286 (1999), pp. 531-537
#' 
#' Ultsch, A., Thrun, M.C., Hansen-Goos, O., Loetsch, J.: Identification of Molecular Fingerprints
#' in Human Heat Pain Thresholds by Use of an Interactive Mixture Model R Toolbox(AdaptGauss),
#' International Journal of Molecular Sciences, doi:10.3390/ijms161025897, 2015.

#' @export

PSClassWithWeight = function(newdat, weights, classProbCut = 0.8, PShighGroup = "PShigh", 
                    PSlowGroup = "PSlow", breaks = 50, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  imputeValue = imputeValue[1]
  ## imputee NA if imputeNA is true
  if(imputeNA){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for PS approach, it does require standardization, so we do not have parameter for standardization, but it should be always done
  newdat = standardize(newdat)
  
  ### in case that some features in weight list are actually not in data, do the following
  tmp = intersect(names(weights), rownames(newdat))
  weights = weights[tmp]
  newdat = newdat[tmp,]

  ## call getTraitParsWithEM to get parameters, although we only need one of them (mean of group means per feature) for PS approach
  PSpars = getTraitParsWithEM(datin = newdat, weights = weights)
  PSpars = cbind(PSpars$meanOfMeans, weights)

  # get PS scores for all samples
  PS_score = apply(newdat, 2, getPS1sample, PSpars = PSpars)
  
  # for PS, first of all, use 0 as a natural cutoff for two group classification
  PS_class0 = ifelse(PS_score >= 0, PShighGroup, PSlowGroup)
  
  ### in addition, after do the following
  ### now, use EM to define 1st draft of group, and then use it to calculate group mean and sd
  
  emcut = AdaptGauss::EMGauss(PS_score, K = 2)
  
  ### add a plot, hist with two distribution lines, do not need to save, just plot it
  hist(PS_score, prob = TRUE, breaks = breaks)
  curve(emcut$Weights[1]*dnorm(x, mean=emcut$Means[1], sd=emcut$SDs[1]), 
        col="red", lwd=2, add=TRUE, yaxt="n")
  curve(emcut$Weights[2]*dnorm(x, mean=emcut$Means[2], sd=emcut$SDs[2]), 
        col="green", lwd=2, add=TRUE, yaxt="n")
  
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
  PS_test = cbind(PS_score, PS_class, PS_prob1, PS_prob2, stringsAsFactors =F)
  
  return(PS_test)
  
}
