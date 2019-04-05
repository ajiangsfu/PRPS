
#' PRPS score calculation and binary classification for a testing data set without PRPStraining output object, but with selected feature weights
#' and a prior of group ratio for the samples in the given testing data set is reasonable.
#' @description This is the function to calculate PRPS (Probability ratio based classification predication score) scores for a testing data set 
#' without PRPS training object. However, we do need selected feature list with their weights, and mean and sd for two groups for each
#' selected feature that we need to apply a ratio prior with default as 1/3 to achieve. Once we have PRPS scores, we will use the 
#' theoretic natual cutoff 0 of PRPS scores to make classification calls. 
#' @details  This is the function to calculate PRPS scores and make classification for a testing data set with known weights and a ratio prior of choice.
#' PRPS calculation is based on Ennishi 2018, its formula is:
#' \eqn{PRPS(X_i) = \sum (|a_j| log(P1(x_ij)/P0(x_ij)))}
#' Here, a_j represents the jth selected feature weights, and x_ij is the corresponding feature value
#'  for the ith sample, 
#' P1 and P0 are the probabilities that the ith sample belongs to two different group.
#' However, in order to calculate P1 and P0, we need to have two group mean and sd for each selected feature. Although there are multiple way
#' to obtain these values, in this function, we design to use a ratio prior to achieve group mean and sd, which was used in Ennishi 2018 originally.
#' After we have have the scores, the natural cutoff of 0 is used to make classification calls.
#' @param newdat a new data matrix or data frame, columns for samples and rows for features
#' @param weights a numeric vector with selected features (as names of the vector) and their weights
#' @param ratioPrior a prior ratio to indicate the ratio of test group over reference group 
#'   regarding mean and sd calculation for each selected feature
#' @param standardization a logic variable to indicate if standardization is needed before classification 
#'  score calculation
#' @param PRPShighGroup a string to indicate group name with high PRPS score
#' @param PRPSlowGroup a string to indicate group name with low PRPS score
#' @param breaks a integer to indicate number of bins in histogram, default is 50
#' @param imputeNA a logic variable to indicate if NA imputation is needed, if it is TRUE, NA imputation is 
#'  processed before any other steps, the default is FALSE
#' @param byrow a logic variable to indicate direction for imputation, default is TRUE, 
#'  which will use the row data for imputation
#' @param imputeValue a character variable to indicate which value to be used to replace NA, default is "median", 
#'  the median value of the chose direction with "byrow" data to be used
#' @return A data frame with PRPS score and classification with 0 as the cutoff
#' @keywords PRPS 
#' @author Aixiang Jiang
#' @references Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P, 
#' Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, Kridel R, 
#' Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. Double-Hit Gene Expression Signature Defines
#' a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell Lymphoma. J Clin Oncol. 
#' 2018 Dec 3:JCO1801583. doi: 10.1200/JCO.18.01583.
#' @export

PRPSClassWithWeightsPrior = function(newdat, weights, ratioPrior = 1/3, standardization=FALSE,  PRPShighGroup = "PRPShigh", 
                    PRPSlowGroup = "PRPSlow", breaks = 50, imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
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
  
  ### add a plot
  hist(PRPS_score, prob = TRUE, breaks = breaks)
  abline(v=0, col="red")

  PRPS_class = ifelse(PRPS_score > 0,  PRPShighGroup,  PRPSlowGroup)  

  PRPS_score = data.frame(PRPS_score)
  PRPS_test = cbind(PRPS_score, PRPS_class, stringsAsFactors =F)
  
  return(PRPS_test)
  
}
