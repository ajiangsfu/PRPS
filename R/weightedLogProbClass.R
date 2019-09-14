#' @export
weightedLogProbClass = function(newdat, topTraits, weights, classMeans, classSds) {
  
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
    dd=0.01
    ### or, we can set dd as a parameter before call the whole function,
    ####   which can be related to 5% or 10% quantile or so
    
    t1=(xx-yy[1])/(zz[1]+dd) ### notice that for a single value xx, n=1, so use sd directly as denominator
    t2=(xx-yy[2])/(zz[2]+dd)
    
    #### then convert t1 and t2 to p values
    #2*(1 - pt(tval, df))
    # or: 2 * pt(abs(t_value), df, lower.tail = FALSE)
    p1=2 * pt(abs(t1), df=1, lower.tail = FALSE)
    p2=2 * pt(abs(t2), df=1, lower.tail = FALSE)
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
