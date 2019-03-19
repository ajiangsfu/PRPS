
#' ROC curve with AUC, Youden best and natural cutoff of 0
#' @description This is an internal function to plot ROC (Receiver operating characteristic) curve showing AUC, Youden best and natural cutoff of 0 points. 
#' @param contdat a vector of a continuous variable
#' @param contname a string to show name of a continuous vector
#' @param catdat a vector of a dichotomous variable
#' @param catname a string to show name of a  dichotomous variable
#' @param xshift a numeric variable to indicate how much value shift along x-axis to move the label position for the natural 0 cutoff point,
#'  default value is 
#' @param yshift a numeric variable to indicate how much value shift along y-axis to move the label position for the natural 0 cutoff point,
#'  default value is 
#' @keywords ROC
#' @author Aixiang Jiang
#' @export
plotROC = function(contdat, contname, catdat, catname, xshift = -0.05, yshift = 0.02){
  require(pROC)
  roc1 = roc(as.factor(catdat), contdat, direction="<", auc=TRUE, ci=TRUE)
  
  titlename = paste("ROC ",contname, " against ",catname,sep="")
  plot(roc1, print.thres="best",print.thres.best.method="youden",print.auc = TRUE, col="blue", lwd=3, main= titlename)
  #### also put "0" on the plot as well
  #### this is to find the point that is the closet to the 0 on the ROC curve
  cut0 = which.min(abs(contdat))
  cut0 = contdat[cut0]
  p0 = coords(roc1,cut0)
  points(p0[2], p0[3], col="red", pch=19, cex=0.8)
  text(p0[2]+xshift, p0[3]+yshift, labels =paste("0 (",format(p0[2], digits = 3),",", format(p0[3],digits = 3), ")",sep=""), col="red")
  
}
