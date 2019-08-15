getMeanSdAllTraits = function(testdat, selectedTraits, selectedTraitWeights, group1ratioPrior = 0.5){
  meanSds = t(mapply(FUN=getMeanSdPerTrait, selectedTraits, selectedTraitWeights, MoreArgs = list(fulldat = testdat, ratioPrior =group1ratioPrior)))
  colnames(meanSds) = c("group1mean","group0mean","group1sd","group0sd")
  meanSds = data.frame(meanSds)
  return(meanSds)
}
