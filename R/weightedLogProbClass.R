#' @export
weightedLogProbClass = function(newdat, topTraits, weights, classMeans, classSds, dfs) {

  genedat = newdat[topTraits,]
  
  names(weights) = rownames(classMeans)
  weights = weights[topTraits]
  
  #### mapply is good for multiple lists, need change format before using the mapply function
  
  classMeans = t(classMeans[topTraits,])
  colnames(classMeans) = rownames(genedat)
  classMeans = data.frame(classMeans)
  
  classSds = t(classSds[topTraits,])
  colnames(classSds) = rownames(genedat)
  classSds = data.frame(classSds)
  
  gendatt=t(genedat)
  gendatt = data.frame(gendatt)
  
  lograt = mapply(FUN = function(xx,yy,zz){
    ### change on 20191112
    p1 = pnorm(q=xx, mean = yy[1], sd = zz[1]) ### normal dist always works
    p2 = pnorm(q=xx, mean = yy[2], sd = zz[2]) ### normal dist always works
    
    gg=log10(p1) - log10(p2) 
  
    },gendatt,classMeans, classSds)
  
  rownames(lograt) = colnames(genedat)
  
  rm(genedat)
  rm(gendatt)
  gc()
  
  
  ### get final y for decision make
  
  #res=lograt %*% abs(weights) ### this line does not consider NA, should change
  # borrow the function for LPS, since it considers the NA within the function
  res = apply (lograt, 1, getLPSscore, coefs= abs(weights))
  
  return(res)
  
}

