
#' Plot function for testing objects
#' @description  This is to plot hist for a given testing object, and scatter plot if the testing object is from LPStesting or PRPStesting
#' @param testObj a testing object from LPStesting, or PRPStesting, or PStesting
#' @param plotNanme a string variable to indicate the file name to save, it should includes path and plot file name (without .pdf part)
#' @param breaks a integer to indicate how many cells in the histogram
#' @keywords  hist, scatter plot
#' @author Aixiang Jiang
#' @export
plotTesting = function(testObj, plotName, breaks = 30){

  ### in order to get direction, get group name and means for each group 
  
  grps = unique(testObj[,2])  
  ### notice that, now I might have two or 3 classes
  ### however, only two groups are done for prob
  
  grps = setdiff(grps, "UNCLASS")
  
  pdf(paste(plotName, ".pdf", sep=""))
 
  ### hist
  hist(testObj[,1], breaks = breaks, xlab = colnames(testObj)[1], main = paste("Histogram of ",colnames(testObj)[1], sep=""))
  
  ### scatter plot if prob is available
  tmp = grep("prob", colnames(testObj))
  if(length(tmp)>1){
    
    ### order data first before plot
    testObj = testObj[order(testObj[,1]),]
    
    plot(testObj[,1], testObj[,tmp[1]],xlab = colnames(testObj)[1], ylab = "Empirical Bayesian probabilites",
         main = paste("Scatter Plot of Prob vs. ",colnames(datin)[1], sep=""))
    points(testObj[,1], testObj[,tmp[2]])
    
    lines(testObj[,1], testObj[,tmp[1]], col="red")
    lines(testObj[,1], testObj[,tmp[2]], col="green")
    
    ### now, should figure out labels for the two groups
    prob1 = testObj[1,tmp[1]]
    prob2 = testObj[1,tmp[2]]
    if(prob1>prob2){
      grp1 = testObj[1,2]
      grp2 = setdiff(grps,grp1)
    }else{
      grp2 = testObj[1,2]
      grp1 = setdiff(grps,grp2)
    }
    legend("left",c(grp1,grp2), col=c("red", "green"), lty=1)
  }
  
  dev.off()
  
  
  
}








