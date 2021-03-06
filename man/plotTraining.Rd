% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotTraining.R
\name{plotTraining}
\alias{plotTraining}
\title{Plot function for a training object}
\usage{
plotTraining(trainObj, plotName, xshift = -0.05, yshift = 0.02, breaks = 30)
}
\arguments{
\item{trainObj}{a training object from LPStraining/LPStrainingWithWeights, or PRPStraining/PRPStrainingWithWeights, or PStraining/PStrainingWithWeights}

\item{plotName}{a string variable to indicate the file name to save, it should include path and plot file name (without .pdf part)}

\item{xshift}{a numeric variable to indicate how much value shift along x-axis to move the label position for the natural 0 cutoff point,
default value is -0.05, which is used for ROC plot but not meaningful for LPS training objects}

\item{yshift}{a numeric variable to indicate how much value shift along y-axis to move the label position for the natural 0 cutoff point,
default value is 0.02, which is used for ROC plot but not meaningful for LPS training objects}

\item{breaks}{a integer to indicate how many breaks in the histogram}
}
\description{
This is to plot a histogram, a scatter plot, and a ROC (Receiver Operating Characteristic curve) curve for a given training object
}
\author{
Aixiang Jiang
}
\keyword{ROC}
\keyword{hist,}
\keyword{plot,}
\keyword{scatter}
