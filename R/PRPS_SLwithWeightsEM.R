
#' PRPS score calculation and binary classification for a testing data set without training data set, but with selected feature weights
#' @description This is the function to calculate PRPS (Probability ratio based classification predication score) scores for a testing data set 
#' without training data set. However, we do need selected feature list with their weights and mean and sd for two groups for each
#' selected feature, which we apply EM (Expectation-Maximization) to achieve in this function. Once we have PRPS scores, we could use the 
#' theoretic natual cutoff 0 to make classification calls. Alternatively, we can also apply EM to calcualate mean and sd for the
#' two groups assuming that PRPS score is a mixture of two normal distributions, followed by Empirical Bayes' probability 
#' calculation and final binary classification calls.
#' @details  This is the function to calculate PRPS scores and make classification based on
#' Empirical Bayesian probabilities for a testing new data set. 
#' PRPS calculation is based on Ennishi 2018, its formula is:
#' \eqn{PRPS(X_i) = \sum (|a_j| log(P1(x_ij)/P0(x_ij)))}
#' Here, a_j represents the jth selected feature weights, and x_ij is the corresponding feature value
#' for the ith sample, P1 and P0 are the probabilities that the ith sample belongs to two different group.
#' However, in order to calculate P1 and P0, we need to have two group mean and sd for each selected feature. Although there are multiple way
#' to obtain these values, in this function, we design to use EM algorithm to achieve group mean and sd assuming that each selected feature
#' is a mixture of two normal distributions. 
#' 
#' After we have PRPS scores, we need to calculate a Empirical Bayes' probability to make classification calls. To do that, 
#' we also need to apply EM to get PRPS score mean and sd for two groups. 
#' After that, we can calcualte probability that a sample belongs to either group,
#' and then use the following formula to get Empirical Bayes' probability:
#' \eqn{prob(x) = p_test(x)/(p_test(x) + p_ref(x))}
#' Here prob(x) is the Empirical Bayes' probability of a given sample, p_test(x) is the probability that a given sample
#' belongs to the test group, p_ref(x) is the probability that a given sample belongs to the reference group.
#' Notice that the test and reference group is just the relative grouping, in fact, for this step, we often need
#' to calculate Empirical Bayes' probabilities for a given sample from two different standing points.
#' 
#' #' Within this function, standardization step is included as an option to minimize the difference 
#' between training and testing data sets, be aware that standardization step should be consistent between the data used 
#' for weight calculation and the testing data set. 
#' Also, be aware that this step is only done to make distributions of each selected
#' features comparable. Be aware that this feature-wise standardization cannot make the sample-wise distributions
#' comparable. For example, if the proportion of one group over all samples might be much lessin the testing data set 
#' compared to the data that were used for weight generation, or even worse, the testing data set might only contain one classification 
#' group only, in this case, feature-wise standardization can not solve the comparison problem. 
#' This is the common problem for classification and feature-wise standardization cannot solve the problem. 
#' In order to solve the problem, we should make data comparable as much as possbile before classification step.
#' 
#' @param newdat a new data matrix or data frame, columns for samples and rows for features
#' @param weights a numeric vector with selected features (as names of the vector) and their weights
#' @param standardization a logic variable to indicate if standardization is needed before classification 
#'  score calculation
#' @param classProbCut a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
#'  often used values are 0.8 and 0.9, default value is 0.8. Only one value is used for both groups, 
#'  the samples that are not included in either group will be assigned as UNCLASS
#' @param PRPShighGroup a string to indicate group name with high PRPS score
#' @param PRPSlowGroup a string to indicate group name with low PRPS score
#' @param breaks a integer to indicate number of bins in histogram, default is 50
#' @param EMmaxRuns number of Iterations for EM searching; default=50
#' @param imputeNA a logic variable to indicate if NA imputation is needed, if it is TRUE, NA imputation is 
#'  processed before any other steps, the default is FALSE
#' @param byrow a logic variable to indicate direction for imputation, default is TRUE, 
#'  which will use the row data for imputation
#' @param imputeValue a character variable to indicate which value to be used to replace NA, default is "median", 
#'  the median value of the chose direction with "byrow" data to be used
#' @return A list with two items is returned: PRPS parameters for selected features, PRPS scores and classifications for the given samples.
#' \item{PRPS_pars}{a list of 3 items, the 1st item is a data frame with weights of each selected features for PRPS
#'  calculation, the 2nd item is a numeric vector containing PRPS mean and sd for two groups，and the 3rd item 
#'  is a data frame contains mean and sd for each group and for each selected feature}
#' \item{PRPS_test}{a data frame of PRPS score, Empirical Bayesian probabilites for both groups and classification based on EM
#' grouping, classification based on 0 as natural cutoff}
#' @keywords PRPS EM 
#' @author Aixiang Jiang
#' @references Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P, 
#'  Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, Kridel R, 
#'  Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. Double-Hit Gene Expression Signature Defines
#'  a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell Lymphoma. J Clin Oncol. 
#'  2018 Dec 3:JCO1801583. doi: 10.1200/JCO.18.01583.
#' 
#' Ultsch, A., Thrun, M.C., Hansen-Goos, O., Loetsch, J.: Identification of Molecular Fingerprints
#' in Human Heat Pain Thresholds by Use of an Interactive Mixture Model R Toolbox(AdaptGauss),
#' International Journal of Molecular Sciences, doi:10.3390/ijms161025897, 2015.
#' 
#' #' Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A trait expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.

#' @export

PRPS_SLwithWeightsEM = function(newdat, weights, standardization=FALSE, classProbCut = 0.8, PRPShighGroup = "PRPShigh", 
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
  ### after I have PRPS score, do the following
  ### now, use EM to define 1st draft of group, and then use it to calculate group mean and sd
  
  emcut = AdaptGauss::EMGauss(PRPS_score, K = 2,fast=TRUE, MaxNumberofIterations = EMmaxRuns)
  
  # > emcut
  # $Means
  # [1] -106.36599  -39.04443
  # 
  # $SDs
  # [1] 23.62815 11.70491
  # 
  # $Weights
  # [1] 0.1139536 0.8860464
  #### this is exact what I need!
  
  ### add a plot, hist with two distribution lines, do not need to save, just plot it
  hist(PRPS_score, prob = TRUE, breaks = breaks)
  curve(emcut$Weights[1]*dnorm(x, mean=emcut$Means[1], sd=emcut$SDs[1]), 
        col="red", lwd=2, add=TRUE, yaxt="n")
  curve(emcut$Weights[2]*dnorm(x, mean=emcut$Means[2], sd=emcut$SDs[2]), 
        col="green", lwd=2, add=TRUE, yaxt="n")
  abline(v=0, col="red")
  
  PRPS_prob1 = getProb(PRPS_score, groupMeans = emcut$Means, groupSds = emcut$SDs)
  PRPS_prob2 = getProb(PRPS_score, groupMeans = rev(emcut$Means), groupSds = rev(emcut$SDs))
  
  if(emcut$Means[1]> emcut$Means[2]){
    name1 = PRPShighGroup
    name2 = PRPSlowGroup
  }else{
    name1 = PRPSlowGroup
    name2 = PRPShighGroup
  }
  
  PRPS_class = rep("UNCLASS",length(PRPS_score))
  PRPS_class[which(PRPS_prob1 >= classProbCut)] = name1
  PRPS_class[which(PRPS_prob2 >= classProbCut)] = name2
  
  ### add class0 with natural cutoff as well
  PRPS_class0 = ifelse(PRPS_score > 0, PRPShighGroup, PRPSlowGroup)
  
  PRPS_score = data.frame(PRPS_score)
  PRPS_test = cbind(PRPS_score, PRPS_class, PRPS_prob1, PRPS_prob2, PRPS_class0, stringsAsFactors =F)

  weights = data.frame(weights)
  
  PRPS_pars =  list(weights, meansds = tmp, traitsmeansds = cbind(allmeans, allsds))
  names(PRPS_pars) = c("weights","meansds","traitsmeansds")
  
  outs = list(PRPS_pars, PRPS_test)
  names(outs) = c("PRPS_pars","PRPS_test")
  
  return(outs)
  
  
}
