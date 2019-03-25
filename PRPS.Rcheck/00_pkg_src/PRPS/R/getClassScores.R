

#' Classification Score Calculation
#' @description This function is to calculate classification prediction score with three method choices: 
#' LPS (defult, Linear Prediction Score) or PRPS (Probability ratio based classification predication score) or PS (Prediction Strength). 
#' @details For the default method: LPS, the linear prediction scores are calculated. When PS is chosen, standardization is included 
#' before score calculation, and we assume that the input of classMeans is also calculated after standardization. 
#' Notice that this function only issue score values without classification. 
#' If there are NAs in the data and not imputed before calling this function, these NAs will be ignored for score calculation.
#' \eqn{LPS(X_i) = \sum a_j x_ij}
#' Here a_j represents the jth selected feature weights, and x_ij is the corresponding feature value for the ith sample.
#' When PRPS method is selected, probability ratio based classification prediction scores are calculated
#' \eqn{PRPS(X_i) = \sum (|a_j| log(P1(x_ij)/P0(x_ij)))}
#' Again, a_j represents the jth selected feature weights, and x_ij is the corresponding feature value for the ith sample, 
#' P1 and P0 are the probabilities that the ith sample belongs to two different group.
#' When PS method is selected, prediction scores are calculated
#' \eqn{PS = (V_win − V_lose)/(V_win + V_lose)}
#' Here, where V_win and V_lose are the vote totals for the winning and losing features/genes for a given sample
#' @param testdat testing data set, a data matrix or a data frame, samples are in columns, and features/genes are in rows
#' @param classMethod three choices with defaul "LPS", the other two methods are "PRPS" and "PS"
#' @param weights  a numeric vector of selected feature weights, which will be used for classfication score calculation,
#'  this is required for all three methods
#' @param classMeans  a two columns' data frame or data matrix, row names are the selected features that can be matched to testdat, 
#'  columns are for group 1 (the 1st column) and group 0 (the 2nd column), the values are the group means for selected features, 
#'  this is needed for PRPS and PS but not LPS
#' @param classSds a two columns' data frame or data matrix, row names are the selected features that can be matched to testdat, 
#'  columns are for group 1 (the 1st column) and group 0 (the 2nd column), the values are the group sds for selected features, 
#'  this is needed for PRPS but not LPS or PS
#' @return A vector of classification predication score
#' @keywords LPS PS PRPS
#' @author Aixiang Jiang
#' @references 
#' Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A gene expression-based method to diagnose
#'  clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S A. 2003;100(17):9991-6.
#' 
#' Golub, T. R. and Slonim, D. K. and Tamayo, P. and Huard, C. and Gaasenbeek, M. and Mesirov, J. P. and
#'  Coller, H. and Loh, M. L. and Downing, J. R. and Caligiuri, M. A. and Bloomfield, C. D. and Lander, E. S.
#' Molecular Classification of Cancer: Class Discovery and Class Prediction by Gene Expression Monitoring. 
#' Science. 1999; 286(5439): 531-7
#' 
#' Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P,
#'  Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, 
#'  Kridel R, Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. 
#'  Double-Hit Gene Expression Signature Defines a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell Lymphoma.
#'   J Clin Oncol. 2018 Dec 3:JCO1801583. doi: 10.1200/JCO.18.01583.
# 
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




