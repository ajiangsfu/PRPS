
#' Empirical Bayes' probability calculation
#' @description This is a function to calculate Empirical Bayes' probability, 
#'   which is called by other functions but it can also be called directly.
#' @details 
#' The coding is written based on Wright 2003 for LPS (linear predication score) approach, 
#' but it can be used for other appraoches such as PS and PRPS as well. 
#' For LPS, the probability step is required since LPS itself does not have a natural cutoff for classification, 
#' alternatively, if the testing and training data arec comparable, we can also make classification based on LPS cutoffs 
#' but this is not commonly used. For PS and PRPS, 0 is a natural cutoff for two group classification, however, 
#' this probability step is still useful if we allow UNCLASS group in the final classification besides the two types 
#' of classes from the training.
#' When calculate a Empirical Bayes' probability, the 1st group in the input mean and sd vectors is treated as the test group. 
#' When calculate the probabilities, we first calcualte probability that a sample belongs to either group, and then:
#' \eqn{prob(x) = p_test(x)/(p_test(x) + p_ref(x))}
#' Here prob(x) is the Empirical Bayes' probability of a given sample, p_test(x) is the probability that a given sample 
#' belongs to the test group, p_ref(x) is the probability that a given sample belongs to the reference group.
#' Notice that the test and reference group is just the relative grouping, in fact, for this step, 
#' we often need to calculate Empirical Bayes' probabilities for a given sample from two different standing points.
#' If you need to calculate Empirical Bayes' probabilities of a group of samples, "apply" function is needed to get all done.
#' @param inscore a classification score, which can be any types of scores
#' @param groupMeans a numeric vector of two items: two classification score means for two training groups/classes
#' @param groupSds a numeric vector of two items: two classification score standard deviations for two training groups/classes
#' @return A probability for a sample belong to a group
#' @keywords Empirical Bayes' probability
#' @author Aixiang Jiang
#' @references
#' Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A gene expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.
#' @export
getProb = function(inscore, groupMeans, groupSds){
  ### assume groupMeans contain 2 values for 2 group, and the 1st one is for positive group
  ### assume groupSds contain 2 values for 2 group, and the 1st one is for positive group
  p1 = dnorm(inscore,mean = groupMeans[1], sd= groupSds[1])
  p0 = dnorm(inscore,mean = groupMeans[2], sd= groupSds[2])
  return(p1/(p1+p0))
}
