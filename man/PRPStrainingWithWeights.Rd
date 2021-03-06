% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PRPStrainingWithWeights.R
\name{PRPStrainingWithWeights}
\alias{PRPStrainingWithWeights}
\title{PRPS training with given weights}
\usage{
PRPStrainingWithWeights(
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

\item{groupInfo}{a known group classification, which order should be the same as in colnames of trainDat}

\item{refGroup}{the code for reference group, default is the 1st item in groupInfo}

\item{classProbCut}{a numeric variable within (0,1), which is a cutoff of Empirical Bayesian probability, 
often used values are 0.8 and 0.9, default value is 0.9. Only one value is used for both groups, 
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
A list with three items is returned: PRPS parameters for selected features, PRPS scores and classifications for training samples, and confusion matrix to compare classification based on PRPS scores and original classification.
\item{PRPS_pars}{a list of 3 items, the 1st item is a data frame with weights and group testing results of each selected features for PRPS calculation, the 2nd item is a numeric vector containing PRPS mean and sd for two groups，and the 3rd item is a data frame contains mean and sd for each group and for each selected feature}
\item{PRPS_train}{a data frame of PRPS score, true classification, Empirical Bayesian probabilites for both groups, and its classification for all training samples, notice that there are two ways for classifications, one is based on probabilities, and there is UNCLASS group besdies the given two groups, alternatively, the other one is based on PRPS scores directly and 0 treated as a natural cutoff}
\item{classCompare}{a confusion matrix list object that compare PRPS classification based on selected features and weights compared to input group classification for training data set, notice that the samples with UNCLASS are excluded since confusion matrix can not compare 3 groups to 2 groups}
\item{classTable}{a table to display comparison of PRPS classification based on selected features and weights compared to input group classification for training data set. Since UNCLASS is excluded from confusion matrix, add this table for full comparison}
}
\description{
This is the wrap up function to estimate parameters, and calculate PRPS (Probability
 ratio based classification predication score) scores based on a given training data set with group info, features and their weights.
}
\details{
PRPS calculation is based on Ennishi 2018, its formula is:
\eqn{PRPS(X_i) = \sum (|a_j| log(P1(x_ij)/P0(x_ij)))}
Here, a_j represents the jth selected feature weights, and x_ij is the corresponding feature value
 for the ith sample, 
P1 and P0 are the probabilities that the ith sample belongs to two different groups.
The precedure is:
a) Calculate mean and sd for each feature for both groups
b) Use "apply" function to get PRPS classification scores and Empirical Bayes' probabilites for all samples.
When we calculate a Empirical Bayes' probability, the 1st group in the input mean and sd vectors is treated
 as the test group. 
When we calculate the probabilities, we first calcualte probability that a sample belongs to either group, 
and then use the following formula to get Empirical Bayes' probability:
\eqn{prob(x) = d_test(x)/(d_test(x) + d_ref(x))}
Here prob(x) is the Empirical Bayes' probability of a given sample, d_test(x) is the density value
 that a given sample belongs to the test group, d_ref(x) is the density value that a given sample belongs
  to the reference group.
c) Provides classification for the training group and confusion matrix to compare PRPS classification
 with original group info for training data set.
If NAs are not imputed, they are ignored for PRPS parameter estimation, and PRPS calculation.
}
\references{
Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P, 
 Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, Kridel R, 
 Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. Double-Hit Gene Expression Signature Defines
 a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell Lymphoma. J Clin Oncol. 
 2018 Dec 3:JCO1801583. doi: 10.1200/JCO.18.01583.
 
 Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A trait expression-based method
to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
A. 2003 Aug 19;100(17):9991-6.
}
\keyword{PRPS}
\keyword{training}
