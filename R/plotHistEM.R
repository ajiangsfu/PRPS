#' A function rather aimed at developers
#' @noRd

plotHistEM = function(scoreIn, G = 1:9, breaks = 50, scoreName = "Input score"){
  # I might always use G=2 for PRPS package, however, it is a good idea to keep G as a range seaching
  clustG = mclust::Mclust(scoreIn, G=G, warn = TRUE)
  bestG = clustG$G
  # emcut = AdaptGauss::EMGauss(scoreIn, K = bestG, fast=TRUE, MaxNumberofIterations = EMmaxRuns)
  ### add a plot, hist with distribution lines, do not need to save, just plot it
  hist(scoreIn, prob = TRUE, breaks = breaks, main = paste("Histogram of", scoreName, sep = " "), xlab = scoreName)
  
  #### I need to dig out the parameters I want
  emcut = clustG$parameters
  emcut$Means = emcut$mean
  emcut$SDs = NA
  if(emcut$variance$modelName == "E"){
    emcut$SDs = rep(sqrt(emcut$variance$sigmasq), bestG)
  }else{
    emcut$SDs = sqrt(emcut$variance$sigmasq)
  }
  
  for (i in 1:bestG){
    curve(emcut$pro[i]*dnorm(x, mean=emcut$Means[i], sd=emcut$SDs[i]), 
          col=i+1, lwd=i, add=TRUE, yaxt="n")
  }
  outs = list(bestG, emcut)
  names(outs) = c("bestG", "emcut")
  return(outs)
  
}