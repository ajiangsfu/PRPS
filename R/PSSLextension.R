
#' PS self learning extension
#' @description This is the function to calculate PS (Prediction Strength) scores and make binary classification calls 
#' for a testing data set with PS_stableSLwithWeights or PSSLwithWeightsPrior self learning object. 
#'  The selected feature list, these features' parameters are from the given PS_stableSLwithWeights or PSSLwithWeightsPrior object.
#' @details  This is the function to calculate PS scores, Empirical Bayesian probabilities and make classification 
#' for a testing new data set. However, this new data set should be comparable to the self learning data set used for PS_stableSLwithWeights or 
#' PSSLwithWeightsPrior as much as possible. 
#'  PS calculation is based on Golub 1999. The formula is:
#'   \eqn{PS = (V_win − V_lose)/(V_win + V_lose)}
#' Here, where V_win and V_lose are the vote totals for the winning and losing features/traits for a given sample
#' This function also give classification for the training group and confusion matrix to compare PS classification
#'  with original group info for training data set.
#' If NAs are not imputed, they are ignored for feature selection, weight calculation, PS parameter estimation,
#'  and PS calculation.
#'  When calculate the probabilities, we first calcualte probability that a sample belongs to either group, 
#'  and then use the following formula to get Empirical Bayes' probability:
#'  \eqn{prob(x) = p_test(x)/(p_test(x) + p_ref(x))}
#'  Here prob(x) is the Empirical Bayes' probability of a given sample, p_test(x) is the probability that a given sample 
#'  belongs to the test group, p_ref(x) is the probability that a given sample belongs to the reference group.
#'  Notice that the test and reference group is just the relative grouping, in fact, for this step, 
#'  we often need to calculate Empirical Bayes' probabilities for a given sample from two different standing points.
#' @param PSSLObj a PS self learning object that is the output from function
#' PS_stableSLwithWeights or PSSLwithWeightsPrior.
#' @param newdat a new data matrix or data frame, which is comparable to self learning data set, 
#'  with columns for samples and rows for features
#' @param classProbCut a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
#'  often used values are 0.8 and 0.9, default value is 0.9. Only one value is used for both groups, 
#'  the samples that are not included in either group will be assigned as UNCLASS
#' @param imputeNA a logic variable to indicate if NA imputation is needed, if it is TRUE, NA imputation is 
#'  processed before any other steps, the default is FALSE
#' @param byrow a logic variable to indicate direction for imputation, default is TRUE, 
#'  which will use the row data for imputation
#' @param imputeValue a character variable to indicate which value to be used to replace NA, default is "median", 
#'  the median value of the chose direction with "byrow" data to be used
#' @return  A data frame with PS scores, Empirical Bayesian probabilites for two groups and classification, and classification based on 0 natural cutoff on PS scores. 
#' @keywords PS
#' @author Aixiang Jiang
#' @references
#' Golub TR, Slonim DK, Tamayo P, Huard C, Gaasenbeek M, Mesirov JP, et al. Molecular classification of cancer: 
#' class discovery and class prediction by gene expression monitoring. Science. 1999;286:531–7
#' 
#' Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A trait expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.
#'  
#' @export

PSSLextension = function(PSSLObj, newdat,classProbCut = 0.9,
                         imputeNA = FALSE, byrow = TRUE, imputeValue = c("median","mean")){
  imputeValue = imputeValue[1]
  
  if(is.null(PSSLObj)){print("Please input your PS self learning object")}
  PS_pars = PSSLObj$PS_pars
  weights = PS_pars$weights
  
  ## impute NA if imputeNA is true
  if(imputeNA){
    newdat = imputeNAs(dataIn = newdat, byrow = byrow, imputeValue = imputeValue)
  }
  
  newdat = standardize(newdat)
  
  ### calculate PS scores
  PS_pars = PSSLObj$PS_pars
  
  pspars = cbind(PS_pars$traitsmeans[,1], PS_pars$weights[,1])
  
  # get PS scores for all samples
  PS_score = apply(newdat[rownames(PS_pars$weights),], 2, getPS1sample, PSpars = pspars)
  
  ### get group info
  testres = PSSLObj$PS_test
  testres = testres[order(testres[,1], decreasing = T),]
  testGroup = testres[1,2]
  
  refGroup = setdiff(unique(PSSLObj$PS_test$PS_class),c(testGroup, "UNCLASS"))
  
  # for PS, 0 is a natural cutoff for two group classification
  # use 0 to get class0, refer the code in the training part
  PS_class0 = ifelse(PS_score>0, testGroup, refGroup)
  
  # alternatively, we can get classification based on prob, and we need to get two groups' PS mean and sd first from training
  # if testing and trainng data sets (more accurately: if the score are comparable) are comparable
  
  testPSmean = PS_pars$meansds[1]
  refPSmean = PS_pars$meansds[2]
  testPSsd = PS_pars$meansds[3]
  refPSsd = PS_pars$meansds[4]
  
  PS_prob_test = getProb(PS_score, groupMeans = c(testPSmean, refPSmean), groupSds = c(testPSsd, refPSsd))
  
  PS_prob_ref = getProb(PS_score, groupMeans = c(refPSmean, testPSmean), groupSds = c(refPSsd, testPSsd))
  
  PS_class = rep("UNCLASS",length(PS_score))
  PS_class[which(PS_prob_test >= classProbCut)] = testGroup
  PS_class[which(PS_prob_ref >= classProbCut)] = refGroup
  
  PS_score = data.frame(PS_score)
  PS_test = cbind(PS_score, PS_class, PS_prob_test, PS_prob_ref, PS_class0, stringsAsFactors =F)
  
  return(PS_test)
  
}


