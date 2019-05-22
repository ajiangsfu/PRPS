
#' PRPS score calculation and binary classification for a testing data set without training data set, but with selected feature weights
#' and a ratio prior
#' @description This is the function to calculate PRPS (Probability ratio based classification predication score) scores for a testing data set 
#' without training data set. However, we do need selected feature list with their weights. In addtion, we need mean and sd for two groups for each
#' selected feature that we need to apply a ratio prior with default as 1/3 to achieve. Once we have PRPS scores, we will use the 
#' theoretic natual cutoff 0 of PRPS scores to make classification calls, as well as classification calls based on Empirical Bayesian probabilities
#' with EM (Expectation-Maximization) and with a ratio prior.
#' @details  This is the function to calculate PRPS scores and make classification for a testing data set with known weights and a ratio prior of choice.
#' PRPS calculation is based on Ennishi 2018, its formula is:
#' \eqn{PRPS(X_i) = \sum (|a_j| log(P1(x_ij)/P0(x_ij)))}
#' Here, a_j represents the jth selected feature weights, and x_ij is the corresponding feature value
#'  for the ith sample, 
#' P1 and P0 are the probabilities that the ith sample belongs to two different group.
#' However, in order to calculate P1 and P0, we need to have two group mean and sd for each selected feature. Although there are multiple way
#' to obtain these values, in this function, we design to use a ratio prior to achieve group mean and sd, which was used in Ennishi 2018 originally.
#' After we have have the scores, the natural cutoff of 0 is used to make classification calls as well as empiracal Bayes probability based classification
#' with the ratio prior.
#' 
#' In addtion, we design to use EM algorithm to get PRPS score mean and sd calculation for the two groups assuming that PRPS score
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
#' Within this function, standardization step is included as an option to minimize the difference 
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
#' @param classProbCut  a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
#'  often used values are 0.8 and 0.9, default value is 0.8. Only one value is used for both groups, 
#'  the samples that are not included in either group will be assigned as UNCLASS
#' @param ratioPrior a prior ratio to indicate the ratio of test group over reference group 
#'   regarding mean and sd calculation for each selected feature
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
#' \item{PRPS_pars}{a list of 4 items, the 1st item is a data frame with weights of each selected features for PRPS
#'  calculation, the 2nd item is a numeric vector containing PRPS mean and sd for two groups based on a given prior，the 3rd item 
#'  is a data frame contains mean and sd for each group and for each selected feature based on a given prior, and the 4th
#'  item is a numeric vector containing PRPS mean and sd for two groups based on EM}
#' \item{PRPS_test}{a data frame of PRPS score, classification and two groups' Empirical Bayesian probabilites based on a given prior, 
#' classification based on 0 as natural cutoff, classification and two groups' Empirical Bayesian probabilites based on EM 
#' Notice that the ratio prior is used for two ends with the same prior for prior based calls}
#' @keywords PRPS
#' @author Aixiang Jiang
#' @references Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P, 
#' Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, Kridel R, 
#' Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. Double-Hit Gene Expression Signature Defines
#' a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell Lymphoma. J Clin Oncol. 
#' 2018 Dec 3:JCO1801583. doi: 10.1200/JCO.18.01583.
#' 
#' #' Ultsch, A., Thrun, M.C., Hansen-Goos, O., Loetsch, J.: Identification of Molecular Fingerprints
#' in Human Heat Pain Thresholds by Use of an Interactive Mixture Model R Toolbox(AdaptGauss),
#' International Journal of Molecular Sciences, doi:10.3390/ijms161025897, 2015.
#' 
#' Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and 
#' density estimation using Gaussian finite mixture models, The R Journal, 8/1, pp. 205-233.
#' 
#' 
#' @export

PRPSSLwithWeightsPrior = function(newdat, weights, standardization=FALSE,  classProbCut = 0.9, ratioPrior = 1/3, PRPShighGroup = "PRPShigh", 
                    PRPSlowGroup = "PRPSlow", EMmaxRuns = 50, breaks = 50, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
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

  #### in order to get mean and sd for each feature, use a ratio prior
  meansd = getMeanSdAllTraits(testdat = newdat, selectedTraits=names(weights), selectedTraitWeights=weights,
                              group1ratioPrior = ratioPrior)
  
  PRPS_score = weightedLogProbClass(newdat=newdat, topTraits=names(weights), weights=weights,
                                 classMeans=meansd[,1:2], classSds = meansd[,3:4])

  PRPS_class0 = ifelse(PRPS_score > 0,  PRPShighGroup,  PRPSlowGroup)  
  
  #### 20190503, call plotHistEM 
  emsearch = plotHistEM(PRPS_score, G = 2:4, breaks = breaks, EMmaxRuns = EMmaxRuns, scoreName = "PRPS_score")
  bestG = emsearch$bestG
  emcut = emsearch$emcut
  
  ### no matter how many bestG, only keep the 1st and last one for the following
  gmeans = c(emcut$Means[1], emcut$Means[bestG])
  gsds = c(emcut$SDs[1], emcut$SDs[bestG])
  
  PRPS_prob1_EM = getProb(PRPS_score, groupMeans = gmeans, groupSds = gsds)
  PRPS_prob2_EM = getProb(PRPS_score, groupMeans = rev(gmeans), groupSds = rev(gsds))
  
  if(gmeans[1]> gmeans[2]){
    name1 = PRPShighGroup
    name2 = PRPSlowGroup
  }else{
    name1 = PRPSlowGroup
    name2 = PRPShighGroup
  }
  
  PRPS_class_EM = rep("UNCLASS",length(PRPS_score))
  PRPS_class_EM[which(PRPS_prob1_EM >= classProbCut)] = name1
  PRPS_class_EM[which(PRPS_prob2_EM >= classProbCut)] = name2
  
  #### work on a given prior based classification
  tmpcut = quantile(PRPS_score, probs = 1-ratioPrior)
  ttmp = PRPS_score[which(PRPS_score >= tmpcut)]
  rtmp = PRPS_score[which(PRPS_score < tmpcut)]
  
  testPRPSmean = mean(ttmp)
  refPRPSmean = mean(rtmp)
  testPRPSsd = sd(ttmp)
  refPRPSsd = sd(rtmp)
  
  PRPS_prob1 = getProb(PRPS_score, groupMeans = c(testPRPSmean, refPRPSmean), groupSds = c(testPRPSsd, refPRPSsd))
  
  PRPS_prob2= getProb(PRPS_score, groupMeans = c(refPRPSmean, testPRPSmean), groupSds = c(refPRPSsd, testPRPSsd))
  
  PRPS_class = rep("UNCLASS",length(PRPS_score))
  PRPS_class[which(PRPS_prob1 >= classProbCut)] = PRPShighGroup
  PRPS_class[which(PRPS_prob2 >= classProbCut)] = PRPSlowGroup
  
  tmp = c(emcut$Means, emcut$SDs)
  names(tmp) = c("testPRPSmean_EM","refPRPSmean_EM","testPRPSsd_EM","refPRPSsd_EM")
  
  #####################################################
  ########### after I used all PRPS_score, then define it as data.frame, otherwise, many of code will get trouble
  PRPS_score = data.frame(PRPS_score)
  PRPS_test = cbind(PRPS_score, PRPS_class, PRPS_prob1, PRPS_prob2, PRPS_class0, PRPS_class_EM, PRPS_prob1_EM, PRPS_prob2_EM,
                    stringsAsFactors =F)

  weights = data.frame(weights)
  
  PRPS_pars =  list(weights, meansds = c(testPRPSmean, refPRPSmean, testPRPSsd, refPRPSsd), traitsmeansds = meansd,
                   tmp)
  names(PRPS_pars) = c("weights","meansds","traitsmeansds", "meansds_EM")
  
  outs = list(PRPS_pars, PRPS_test)
  names(outs) = c("PRPS_pars","PRPS_test")
  
  return(outs)
 
}
