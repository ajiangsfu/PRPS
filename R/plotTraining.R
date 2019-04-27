
#' Plot function for training objects
#' @description  This is to plot hist for a given training object, and scatter plot if the training object us from LPStraining or PRPStraining
#' @param trainObj a training object from LPStraining, or PRPStraining, or PStraining
#' @param plotNanme a string variable to indicate the file name to save, it should includes path and plot file name (without .pdf part)
#' @param breaks a integer to indicate how many cells in the histogram
#' @keywords hist, scatter plot
#' @author Aixiang Jiang
#' @export
plotTraining = function(trainObj, plotName, breaks = 30){
  
  ### remove ROC on 20190425
  
  ## call hist, and scatter plot
  
  datin = trainObj[[2]]
  datin[,2] = as.character(datin[,2])
  
  ### in order to get direction, get group name and means for each group 
  
  grps = unique(datin[,2])
  
  mean1 = mean(subset(datin, datin[,2] == grps[1])[,1]) 
  mean2 = mean(subset(datin, datin[,2] == grps[2])[,1]) 
  
  if(mean1 > mean2){
    datin$class01 = ifelse(datin[,2] == grps[1], 1, 0)
  }else{
    datin$class01 = ifelse(datin[,2] == grps[2], 1, 0)
  }
  
  pdf(paste(plotName, ".pdf", sep=""))
  
  ### hist
  hist(datin[,1], breaks = breaks, xlab = colnames(datin)[1], main = paste("Histogram of ",colnames(datin)[1], sep=""))
  
  ### scatter plot if prob is available
  tmp = grep("prob", colnames(datin))
  if(length(tmp)>1){
    
    ### order data first before plot
    datin = datin[order(datin[,1]),]
    
    plot(datin[,1], datin[,tmp[1]],xlab = colnames(datin)[1], ylab = "Empirical Bayesian probabilites",
         main = paste("Scatter Plot of Prob vs. ",colnames(datin)[1], sep=""))
    points(datin[,1], datin[,tmp[2]])
    
    lines(datin[,1], datin[,tmp[1]], col="red")
    lines(datin[,1], datin[,tmp[2]], col="green")
    
    ### now, should figure out labels for the two groups
    prob1 = datin[1,tmp[1]]
    prob2 = datin[1,tmp[2]]
    if(prob1>prob2){
      grp1 = datin[1,2]
      grp2 = setdiff(grps,grp1)
    }else{
      grp2 = datin[1,2]
      grp1 = setdiff(grps,grp2)
    }
    legend("left",c(grp1,grp2), col=c("red", "green"), lty=1)
  }
  
  dev.off()
  
  
  
}








