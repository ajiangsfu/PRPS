% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PSSLextension.R
\name{PSSLextension}
\alias{PSSLextension}
\title{PS self-training extension}
\usage{
PSSLextension(
  PSSLObj,
  newdat,
  classProbCut = 0.9,
  imputeNA = FALSE,
  byrow = TRUE,
  imputeValue = c("median", "mean")
)
}
\arguments{
\item{PSSLObj}{a PS self-training object that is the output from function
PS_stableSLwithWeights or PSSLwithWeightsPrior.}

\item{newdat}{a new data matrix or data frame, which is comparable to self-training data set, 
with columns for samples and rows for features}

\item{classProbCut}{a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
often used values are 0.8 and 0.9, default value is 0.9. Only one value is used for both groups, 
the samples that are not included in either group will be assigned as UNCLASS}

\item{imputeNA}{a logic variable to indicate if NA imputation is needed, if it is TRUE, NA imputation is 
processed before any other steps, the default is FALSE}

\item{byrow}{a logic variable to indicate direction for imputation, default is TRUE, 
which will use the row data for imputation}

\item{imputeValue}{a character variable to indicate which value to be used to replace NA, default is "median", 
the median value of the chose direction with "byrow" data to be used}
}
\value{
A data frame with PS scores, Empirical Bayesian probabilites for two groups and classification, and classification based on 0 natural cutoff on PS scores.
}
\description{
This is the function to calculate PS (Prediction Strength) scores and make binary classification calls 
for a testing data set with PS_stableSLwithWeights or PSSLwithWeightsPrior self-training object. 
 The selected feature list, these features' parameters are from the given PS_stableSLwithWeights or PSSLwithWeightsPrior object.
}
\details{
This is the function to calculate PS scores, Empirical Bayesian probabilities and make classification 
for a testing new data set. However, this new data set should be comparable to the self-training data set used for PS_stableSLwithWeights or 
PSSLwithWeightsPrior as much as possible. 
 PS calculation is based on Golub 1999. The formula is:
  \eqn{PS = (V_win − V_lose)/(V_win + V_lose)}
Here, where V_win and V_lose are the vote totals for the winning and losing features/traits for a given sample
This function also give classification for the training group and confusion matrix to compare PS classification
 with original group info for training data set.
If NAs are not imputed, they are ignored for feature selection, weight calculation, PS parameter estimation,
 and PS calculation.
 When calculate the probabilities, we first calcualte probability that a sample belongs to either group, 
 and then use the following formula to get Empirical Bayes' probability:
 \eqn{prob(x) = p_test(x)/(p_test(x) + p_ref(x))}
 Here prob(x) is the Empirical Bayes' probability of a given sample, p_test(x) is the probability that a given sample 
 belongs to the test group, p_ref(x) is the probability that a given sample belongs to the reference group.
 Notice that the test and reference group is just the relative grouping, in fact, for this step, 
 we often need to calculate Empirical Bayes' probabilities for a given sample from two different standing points.
}
\references{
Golub TR, Slonim DK, Tamayo P, Huard C, Gaasenbeek M, Mesirov JP, et al. Molecular classification of cancer: 
class discovery and class prediction by gene expression monitoring. Science. 1999;286:531–7

Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A trait expression-based method
to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
A. 2003 Aug 19;100(17):9991-6.
}
\author{
Aixiang Jiang
}
\keyword{PS}
