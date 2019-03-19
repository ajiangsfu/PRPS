
#' Plot function for training objects
#' @description  This is to plot ROC, hist, and scatter plot for a given training object
#' @param trainObj a training object from LPStrain, or PRPStrain, or PStrain
#' @param plotNanme a string variable to indicate the file name to save, it should includes path and plot file name (without .pdf part)
#' @param xshift a numeric variable to indicate how much value shift along x-axis to move the label position for the natural 0 cutoff point,
#'  default value is 
#' @param yshift a numeric variable to indicate how much value shift along y-axis to move the label position for the natural 0 cutoff point,
#'  default value is 
#' @keywords ROC, hist, scatter plot
#' @author Aixiang Jiang
#' @export
plotTraining = function(trainObj, plotName, xshift = -0.05, yshift = 0.02){
  
  ## call plotROC, hist, and scatter plot
  
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
  plotROC(contdat = datin[,1], contname = colnames(datin)[1], catdat = datin$class01, catname = colnames(datin)[2], xshift = xshift, yshift = yshift)
  dev.off()
  
}








