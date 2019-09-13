#' @export
getClassScores = function(testdat, classMethod = c("LPS","PRPS","PS"), weights, classMeans, classSds){
  ### testdat: gene is row, sample is column
  ### weights should have names with it
  classMethod = classMethod[1]
  
  if(classMethod == "PRPS"){
    scoreOuts = weightedLogProbClass(newdat=testdat, topTraits=rownames(classMeans), weights=weights, 
                        classMeans= classMeans, classSds = classSds)
  }else if(classMethod == "PS"){
    testdat = standardize(testdat)
    scoreOuts = apply(testdat[rownames(classMeans),], 2, getPS1sample, PSpars = cbind(rowMeans(classMeans),weights))
    
  }else{  ### default: LPS
    testdat = testdat[names(weights),]
    scoreOuts = apply (testdat, 2, getLPSscore, coefs= weights)
    
  }
  return(scoreOuts)
}




