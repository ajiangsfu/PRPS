% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LPStrainingWithWeights.R
\name{LPStrainingWithWeights}
\alias{LPStrainingWithWeights}
\title{LPS training with given weights}
\usage{
LPStrainingWithWeights(
  trainDat,
  weights,
  groupInfo,
  refGroup = NULL,
  classProbCut = 0.9,
  imputeNA = FALSE,
  byrow = TRUE,
  imputeValue = c("median", "mean"),
  standardization = FALSE
)
}
\arguments{
\item{trainDat}{training data set, a data matrix or a data frame, samples are in columns, and features/traits are in rows}

\item{weights}{a numeric vector with selected features (as names of the vector) and their weights}

\item{groupInfo}{a known group classification, which order should be the same as in columns of trainDat}

\item{refGroup}{the code for reference group, default is the 1st item in groupInfo}

\item{classProbCut}{a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
often used values are 0.8 and 0.9, default value is 0.9. The same classProbCut is used for both groups, 
the samples that are not included in either group will be assigned as UNCLASS}

\item{imputeNA}{a logic variable to indicate if NA imputation is needed, if it is TRUE, 
NA imputation is processed before any other steps, the default is FALSE}

\item{byrow}{a logic variable to indicate direction for imputation, default is TRUE, 
which will use the row data for imputation}

\item{imputeValue}{a character variable to indicate which value to be used to replace NA, default is "median", 
the median value of the chose direction with "byrow" data to be used}

\item{standardization}{a logic variable to indicate if standardization is needed before classification 
score calculation}
}
\value{
A list with four items is returned: LPS parameters for selected features, LPS scores and classifications for training samples, 
confusion matrix to compare classification based on LPS scores and original classification, and a simple classification comparison table.
\item{LPS_pars}{a list of 2 items, the 1st item is a data frame with feature weights, and the 2nd item is a numeric vector containing LPS mean and sd for two groups}
\item{LPS_train}{a data frame of LPS score, true classification, Empirical Bayesian probabilites for both groups, 
and its classification for all training samples, notice that the classification is based on probabilities instead of LPS scores, 
and there is a UNCLASS group besdies the given two groups}
\item{classCompare}{a confusion matrix list object that compares LPS classification (based on selected features and 
weights) to input group classification for training data set, notice that the samples with UNCLASS 
are excluded since confusion matrix can not compare 3 groups to 2 groups}
\item{classTable}{a table to display comparison between LPS classification (based on selected features and weights) 
and input group classification for the given training data set. Since UNCLASS is excluded from confusion matrix, we add this table for full comparison}
}
\description{
This is to calculate LPS (Linear Predictor Score) scores based on a given training data set with classification
labels and selected features' weights
}
\details{
This is to calculate LPS (Linear Predictor Score) scores based on a given training data set with classification
labels and selected features' weights. Notice that there is no feature selection and weight calculation step in this function 
since we already have selected features and their weights from elsewhere, which is different from LPStraining. 

First of all, we calculate LPS classification scores for all samples based on formula in Wright 2003:
\eqn{LPS(X) = \sum a_j x_ij}
Here a_j represents the jth selected feature weights, and x_ij is the corresponding feature value 
 for the ith sample.
 
Then, we calculate Empirical Bayesian probabilities. By default, the 1st group in the input mean and sd vectors is treated as the 
test group. When we calculate the probabilities, we first calcualte probability that a sample belongs to either group,
and then use the following formula to get Empirical Bayesian probability:
\eqn{prob(x) = d_test(x)/(d_test(x) + d_ref(x))}
Here prob(x) is the Empirical Bayesian probability of a given sample, d_test(x) is the density value assuming that a given sample
belongs to the test group, d_ref(x) is the density value assuming that a given sample belongs to the reference group.
In the current function, however, we calculate Empirical Bayesian probabilities for both directions.

Finally, this wrap-up function also gives classification for the training group and confusion matrix to compare 
LPS classification with original group info for training data set.
If NAs are not imputed, they are ignored for feature selection, weight calculation, LPS parameter estimation, 
and LPS calculation.
}
\references{
Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A trait expression-based method
to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
A. 2003 Aug 19;100(17):9991-6.
}
\author{
Aixiang Jiang
}
\keyword{LPS}
\keyword{training}
