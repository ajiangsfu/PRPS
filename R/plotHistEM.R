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