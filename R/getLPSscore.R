getLPSscore = function(vdat, coefs){
  tmp = which(is.na(vdat))
  if(length(tmp) > 0){
    vdat = vdat[-tmp]
    coefs = coefs[-tmp]
  }
  aScore = mapply(function(x,y){x*y}, x=vdat, y=coefs)
  aScore = sum(unlist(aScore))
  return(aScore)
}
