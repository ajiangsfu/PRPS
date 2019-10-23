#' Feature selection, parameter estimation, and PRPS calculation for training data set
#' @description This is the wrap up function to select top features, estimate parameters, and calculate PRPS (Probability
#'  ratio based classification predication score) scores based on a given training data set.
#' @details PRPS calculation is based on Ennishi 2018, its formula is:
#' \eqn{PRPS(X_i) = \sum (|a_j| log(P1(x_ij)/P0(x_ij)))}
#' Here, a_j represents the jth selected feature weights, and x_ij is the corresponding feature value
#'  for the ith sample, 
#' P1 and P0 are the probabilities that the ith sample belongs to two different group.
#' In this wrap up function, we use three steps to calculate PRPS scores and classification. 
#' Before these three steps, we also give an option for NA imputation and for standardization for each feature. 
#' The three steps are:
#' a) Apply "getTrainingWeights" to select features and return weights for these features.
#' b) Use "apply" function to get PRPS classification scores and Empirical Bayes' probabilites for all samples.
#' When we calculate a Empirical Bayes' probability, the 1st group in the input mean and sd vectors is treated
#'  as the test group. 
#' When we calculate the probabilities, we first calcualte probability that a sample belongs to either group, 
#'  and then use the 
#' following formula to get Empirical Bayes' probability:
#' \eqn{prob(x) = d_test(x)/(d_test(x) + d_ref(x))}
#' Here prob(x) is the Empirical Bayes' probability of a given sample, d_test(x) is the density value
#'  that a given sample belongs to the test group, d_ref(x) is the density value that a given sample belongs
#'   to the reference group.
#' Notice that the test and reference group is just the relative grouping, in fact, for this step, 
#' we often need to calculate Empirical Bayes' probabilities for a given sample from two different standing points.
#' c) This function also give classification for the training group and confusion matrix to compare PRPS classification
#'  with original group info for training data set.
#' If NAs are not imputed, they are ignored for feature selection, weight calculation, PRPS parameter estimation, 
#' and PRPS calculation.
#' @param trainDat training data set, a data matrix or a data frame, samples are in columns, and features/traits are in rows
#' @param selectedTraits  a selected trait list if available
#' @param groupInfo a known group classification, which order should be the same as in colnames of trainDat
#' @param refGroup the code for reference group, default is 0, but it can be a string or other number, which will be 
#'  changed to 0 within the function
#' @param topN an integer to indicate how many top features to be selected
#' @param FDRcut  a FDR cutoff to select top features, which is only valid when topN is set as defaul NULL, 
#'  all features will be returned if both topN and FDRcut are set as default NULL
#' @param weightMethod  a string to indicate weight calculation method, there are five choices: 
#'  "limma" for for limma linear model based t value,"ttest" for t test based t value, 
#'  "MannWhitneyU" for Mann Whitney U based rank-biserial,"PearsonR" for Pearson correlation coefficient,
#'  "SpearmanR" for Spearman correlation coefficient, and the defualt value is "limma"
#' @param classProbCut a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
#'  often used values are 0.8 and 0.9, default value is 0.9. Only one value is used for both groups, 
#'  the samples that are not included in either group will be assigned as UNCLASS
#' @param imputeNA a logic variable to indicate if NA imputation is needed, if it is TRUE, 
#'  NA imputation is processed before any other steps, the default is FALSE
#' @param byrow a logic variable to indicate direction for imputation, default is TRUE, 
#'  which will use the row data for imputation
#' @param imputeValue a character variable to indicate which value to be used to replace NA, default is "median", 
#'  the median value of the chose direction with "byrow" data to be used
#' @param standardization a logic variable to indicate if standardization is needed before classification 
#'  score calculation
#' @keywords PRPS training limma weight
#' @return A list with three items is returned: PRPS parameters for selected features, PRPS scores and classifications for training samples, and confusion matrix to compare classification based on PRPS scores and original classification.
#' \item{PRPS_pars}{a list of 3 items, the 1st item is a data frame with weights and group testing results of each selected features for PRPS calculation, the 2nd item is a numeric vector containing PRPS mean and sd for two groups，and the 3rd item is a data frame contains mean and sd for each group and for each selected feature}
#' \item{PRPS_train}{a data frame of PRPS score, true classification, Empirical Bayesian probabilites for both groups, and its classification for all training samples, notice that there are two ways for classifications, one is based on probabilities, and there is UNCLASS group besdies the given two groups, alternatively, the other one is based on PRPS scores directly and 0 treated as a natural cutoff}
#' \item{classCompare}{a confusion matrix list object that compare PRPS classification based on selected features and weights compared to input group classification for training data set, notice that the samples with UNCLASS are excluded since confusion matrix can not compare 3 groups to 2 groups}
#' \item{classTable}{a table to display comparison of PRPS classification based on selected features and weights compared to input group classification for training data set. Since UNCLASS is excluded from confusion matrix, add this table for full comparison}
#' @references Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P, 
#'  Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, Kridel R, 
#'  Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. Double-Hit Gene Expression Signature Defines
#'  a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell Lymphoma. J Clin Oncol. 
#'  2018 Dec 3:JCO1801583. doi: 10.1200/JCO.18.01583.
#'  
#'  Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A trait expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.
#' @export

PRPStraining = function(trainDat, standardization = FALSE, selectedTraits = NULL, groupInfo, refGroup = 0, topN = NULL, FDRcut = 0.1,
                        weightMethod = c("ttest","limma","PearsonR", "SpearmanR", "MannWhitneyU"), classProbCut = 0.9, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  weightMethod = weightMethod[1]
  ## impute NA if imputeNA is true
  imputeValue = imputeValue[1]
  
  if(imputeNA){
    trainDat = imputeNAs(dataIn = trainDat, byrow = byrow, imputeValue = imputeValue)
  }
  
  # for PRPS approach, it does not require standardization, however, if standardization = TRUE, do the standardization
  if(standardization){trainDat = standardize(trainDat)}
  
  # use reference group, we can select top features and calculate their weights
  weights = getTrainingWeights(trainDat = trainDat, selectedTraits = selectedTraits, groupInfo = groupInfo,
                               refGroup = refGroup, topN = topN, FDRcut = FDRcut, weightMethod = weightMethod)
  
  # calculate mean and sd for each group, and then get PRPS scores for all samples
  sigdat = trainDat[rownames(weights),]
  refind = which(groupInfo == refGroup)
  
  ### 20191016, add df0 and df1
  df0 = length(refind) -1 
  g0dat = sigdat[,refind]
  g1dat = sigdat[, -refind]
  df1 = dim(g1dat)[2]-1
  
  g0mean = apply(g0dat,1,mean, na.rm=T)
  g1mean = apply(g1dat,1,mean, na.rm=T)
  
  g0sd = apply(g0dat,1,sd, na.rm=T)
  g1sd = apply(g1dat,1,sd, na.rm=T)
  
  traitsmeansds = cbind(g1mean, g0mean, g1sd, g0sd)
  
  colnames(traitsmeansds) = c("testmean","refmean", "testsd","refsd")
  
  dfs = c(df1, df0)
  
  PRPS_score = weightedLogProbClass(newdat = sigdat, topTraits=rownames(weights), weights=weights[,1],
                                    classMeans = traitsmeansds[,1:2], classSds = traitsmeansds[,3:4], dfs = dfs)
  
  # and use get prob function to get classification
  
  testGroup = setdiff(unique(groupInfo), refGroup)
  
  refPRPS = PRPS_score[refind]
  testPRPS = PRPS_score[-refind]
  
  refPRPSmean = mean(refPRPS, na.rm = T)
  refPRPSsd = sd(refPRPS, na.rm = T)
  
  testPRPSmean = mean(testPRPS, na.rm = T)
  testPRPSsd = sd(testPRPS, na.rm = T)
  
  PRPS_prob_test = getProb(PRPS_score, groupMeans = c(testPRPSmean, refPRPSmean), groupSds = c(testPRPSsd, refPRPSsd))
  PRPS_prob_ref = getProb(PRPS_score, groupMeans = c(refPRPSmean, testPRPSmean), groupSds = c(refPRPSsd, testPRPSsd))
  
  PRPS_class = rep("UNCLASS",length(PRPS_score))
  PRPS_class[which(PRPS_prob_test >= classProbCut)] = testGroup
  PRPS_class[which(PRPS_prob_ref >= classProbCut)] = refGroup
  
  ### add one more column on Dec 27, this is classification based on 0 natural cutoff
  PRPS_class0 = ifelse(PRPS_score>0, testGroup, refGroup)
  
  PRPS_score = data.frame(PRPS_score)
  true_class = groupInfo
  PRPS_train = cbind(PRPS_score, true_class, PRPS_class, PRPS_prob_test, PRPS_prob_ref,PRPS_class0 ,stringsAsFactors =F)
  
  groupInfo = factor(groupInfo, levels = c(refGroup, testGroup))
  ## in order to get comparison, change UNCLASS to NA, therefore only two groups are considered in PRPS_class
  PRPS_class2 = ifelse(PRPS_class == "UNCLASS", NA, PRPS_class)
  PRPS_class2 = factor(PRPS_class, levels = c(refGroup,  testGroup))
  
  ## notice that confusion matrix does not work if the number of levels are not the same
  
  classCompare = caret::confusionMatrix(PRPS_class2, reference = groupInfo, positive = testGroup)
  
  meansds = c(testPRPSmean, refPRPSmean, testPRPSsd, refPRPSsd)
  names(meansds) = c("testPRPSmean","refPRPSmean","testPRPSsd","refPRPSsd")
  
  PRPS_pars =  list(weights,meansds, traitsmeansds, dfs)
  names(PRPS_pars) = c("weights","meansds","traitsmeansds","dfs")
  
  #### since UNCLASS is excluded from confusion matrix, add one more output for full comparison
  classTable = table(groupInfo, PRPS_class)
  
  outs = list(PRPS_pars, PRPS_train, classCompare, classTable)
  names(outs) = c("PRPS_pars","PRPS_train","classCompare", "classTable")  ### dfs is in PRPS_pars
  return(outs)
  
}
