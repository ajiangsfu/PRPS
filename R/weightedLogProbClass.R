#' A function rather aimed at developers
#' @noRd

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
  
  ### change on 20191113
  ee = 0.0001
  tmp = which(classSds < ee, arr.ind = T)
  if(dim(tmp)[1] > 0) {
    classSds = classSds[,-tmp[,2]]
    classMeans = classMeans[,-tmp[,2]]
    gendatt  = gendatt[,-tmp[,2]]
    weights = weights[-tmp[,2]]
  }

  lograt = mapply(FUN = function(xx,yy,zz){
    dd=0.01 ### this value is kind of arbitrary, but seems working well for now, keep it to avoid any possible problem

    ### change on 20191016, struggled on if I need to add dd or not
    ### in the end, add dd in the top and bottom to get benefit of avoiding too small sd and not change the t value too much
    t1=(xx-yy[1]+dd)/(zz[1]+dd) ### notice that for a single value xx, n=1, so use sd directly as denominator
    t2=(xx-yy[2]+dd)/(zz[2]+dd)
  
    p1=2 * pt(abs(t1), df= dfs[1], lower.tail = FALSE) 
    p2=2 * pt(abs(t2), df= dfs[2], lower.tail = FALSE)
  
    #### Feb 8, 2021
    #### to avoid inf value of log, change the code to the following format
    gg=log10(p1+0.00000001) - log10(p2+0.00000001)
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

