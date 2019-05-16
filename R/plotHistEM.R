#' Histogram with EM curves
#' @description Histogram plot with EM curves, and EM output
#' @details In this function, a histogram of the given score vector is plotted with EM curves. The number of EM
#' curves is selected with Mclust. Based on this best selected number of clusters, EMGauss object based on given score
#' vector and the best number of clusters are also provided.
#' @param scoreIn a numeric vector used to make classification
#' @param G an integer vector to indicate searching range for number of clusters
#' @param breaks an integer to indicate number of bins in histogram, default is 50
#' @param EMmaxRuns number of Iterations for EM searching; default=50
#' @param scoreName a string to indicate in the plot, default is "Input score"
#' @return a list of two items: bestG and emcut
#' \item{bestG}{an integer of the best number of clusters}
#' \item{emcut}{an EMGauss object, a list with Means, SDs, and Weights}
#' @keywords EM 
#' @author Aixiang Jiang
#' @references 
#' Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and 
#' density estimation using Gaussian finite mixture models, The R Journal, 8/1, pp. 205-233.
#' 
#' Ultsch, A., Thrun, M.C., Hansen-Goos, O., Loetsch, J.: Identification of Molecular Fingerprints
#' in Human Heat Pain Thresholds by Use of an Interactive Mixture Model R Toolbox(AdaptGauss),
#' International Journal of Molecular Sciences, doi:10.3390/ijms161025897, 2015.
#' @export

plotHistEM = function(scoreIn, G = 1:9, breaks = 50, EMmaxRuns = 50, scoreName = "Input score"){
  #require(mclust)
  #require(AdaptGauss)
  clustG = mclust::Mclust(scoreIn, G=G)
  bestG = clustG$G
  emcut = AdaptGauss::EMGauss(scoreIn, K = bestG, fast=TRUE, MaxNumberofIterations = EMmaxRuns)
  ### add a plot, hist with distribution lines, do not need to save, just plot it
  hist(scoreIn, prob = TRUE, breaks = breaks, main = paste("Histogram of", scoreName, sep = " "), xlab = scoreName)
  for (i in 1:bestG){
    curve(emcut$Weights[i]*dnorm(x, mean=emcut$Means[i], sd=emcut$SDs[i]), 
          col=i+1, lwd=i, add=TRUE, yaxt="n")
  }
  outs = list(bestG, emcut)
  names(outs) = c("bestG", "emcut")
  return(outs)
  
}