
#' Plot function for a training object
#' @description  This is to plot a histogram, a scatter plot, and a ROC (Receiver Operating Characteristic curve) curve for a given training object 
#' @param trainObj a training object from LPStraining/LPStrainingWithWeights, or PRPStraining/PRPStrainingWithWeights, or PStraining/PStrainingWithWeights
#' @param plotName a string variable to indicate the file name to save, it should include path and plot file name (without .pdf part)
#' @param xshift a numeric variable to indicate how much value shift along x-axis to move the label position for the natural 0 cutoff point,
#'  default value is -0.05, which is used for ROC plot but not meaningful for LPS training objects
#' @param yshift a numeric variable to indicate how much value shift along y-axis to move the label position for the natural 0 cutoff point,
#'  default value is 0.02, which is used for ROC plot but not meaningful for LPS training objects
#' @param breaks a integer to indicate how many breaks in the histogram
#' @keywords hist, scatter plot, ROC
#' @author Aixiang Jiang
#' @export
plotTraining = function(trainObj, plotName, xshift = -0.05, yshift = 0.02, breaks = 30){
  
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
  
  ################################ ROC plot with discussion, ignore for now ##################################################
  #plotROC(contdat = datin[,1], contname = colnames(datin)[1], catdat = datin[,2], catname = colnames(datin)[2], xshift = xshift, yshift = yshift)
  ### 1st column is the score, 2nd column is the true class
  ### the code itself is correct, for an example:
  # > head(tmp)
  # PRPS_score true_class PRPS_class PRPS_prob_test PRPS_prob_ref PRPS_class0
  # LYM018 -139.91397        GCB        GCB   3.202854e-10  1.000000e+00         GCB
  # LYM120   81.31242        ABC        ABC   9.997699e-01  2.301337e-04         ABC
  # LYM180  154.28566        ABC        ABC   9.999996e-01  4.122399e-07         ABC
  # LYM285  -43.14212        GCB        GCB   1.375906e-03  9.986241e-01         GCB
  # LYM233 -128.44472        GCB        GCB   2.311920e-09  1.000000e+00         GCB
  # LYM407 -123.17014        GCB        GCB   5.652463e-09  1.000000e+00         GCB
  # pRPC::roc(tmp$true_class, tmp$PRPS_score,percent=TRUE, plot=TRUE, ci=TRUE)
  # ### Data: tmp$PRPS_score in 17 controls (tmp$true_class ABC) > 23 cases (tmp$true_class GCB).
  ### Area under the curve: 100%
  ### 95% CI: 100%-100% (DeLong)
  ############ this example gave me 100% AUC, which looks strange, but this is correct
  ##### 100% AUC means: an optimal cutoff exists to separately samples into two groups based on scores without any mistake
  #####                 although this optimal cutoff is not necessary consistent to our final classification
  
  ### hist
  ### hist(datin[,1], breaks = breaks, xlab = colnames(datin)[1], main = paste("Histogram of ",colnames(datin)[1], sep=""))
  
  ## change the hist part on 20200306
  plotHistEM(datin[,1], G = 1:2, breaks = breaks, scoreName = colnames(datin)[1])
  
  ### scatter plot if prob is available
  tmp = grep("prob", colnames(datin))
  if(length(tmp)>1){
    
    ### order data first before plot
    datin = datin[order(datin[,1]),]
    
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
    
    datin$true_class
    plot(datin[,1], datin[,tmp[1]],xlab = colnames(datin)[1], ylab = "Empirical Bayesian probabilites",
         main = paste("Scatter Plot of Prob vs. ",colnames(datin)[1], sep=""))
    points(datin[,1], datin[,tmp[2]])
    
    lines(datin[,1], datin[,tmp[1]], col="red")
    lines(datin[,1], datin[,tmp[2]], col="green")
    
    ### change on 20190731
    pgrp1 = paste("prob", grp1, sep="_")
    pgrp2 = paste("prob", grp2, sep="_")
    legend("left",c(pgrp1,pgrp2), col=c("red", "green"), lty=1, bty = "n")
    abline(h=0.1, lty=3)
    abline(h=0.9, lty=3)
    
    rps = subset(datin, datin$true_class == grp1)
    gps = subset(datin, datin$true_class == grp2)
    
    points(rps[,1], rps[,tmp[1]], col= "red", pch = 19)
    points(rps[,1], rps[,tmp[2]], col= "red", pch = 19)
    
    points(gps[,1], gps[,tmp[1]], col= "green", pch = 19)
    points(gps[,1], gps[,tmp[2]], col= "green", pch = 19)
    
    legend("right", c(grp1, grp2), col=c("red", "green"), pch = 19, bty="n")
    
  }
  
  dev.off()
  
}
