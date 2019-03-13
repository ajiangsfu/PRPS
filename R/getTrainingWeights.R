
#' Feature selection and their weight calculation for PS, LPS and PRPS approaches
#' @description  this is a wrap up function to select traits, and get their weights for each selected feature
#' @details Within the five weight calculation methods, "limma","ttest", "PearsonR","SpearmanR" are often used and nothing 
#' is changed for them from usual setting in this function. We describe more detail about "MannWhitneyU" here. 
#' Mann-Whitney U test is a well known rank based two group comparison test method, but its statistics value U is not easy 
#' to use for weight. U can be calculated in two different ways depending on which group we stand for. 
#' Assuming that there are two groups with sample sizes n1 and n2, and rank sum for the two groups are R1 and R2, 
#' then we can have the following two U values:
#' \eqn{ U_1 = R_1 - n_1(n_1 + 1)/2}
#' \eqn{ U_2 = R_2 - n_2(n_2 + 1)/2}
#' In this function, we use the test group (non reference group) to calculate U value, and then calculate 
#' its rank-biserial correlation with the following formula:
#' \eqn{ r = 1 - 2U/(n_1 * n_2)}
#' P value for "MannWhitneyU", however, is still based on standard Mann - Whitney U test.
#' NA is ignored before testing.
#' @param trainDat training data set, a data matrix or a data frame, samples are in columns, 
#'   and features/traits are in rows
#' @param selectedTraits a selected trait list if available
#' @param groupInfo  a known group classification, which order should be the same as in colnames of trainDat
#' @param refGroup the code for reference group, default is 0, but it can be a string or other number, 
#'   which will be changed to 0 within the function 
#' @param FDRcut FDR cutoff for feature selection
#' @param weightMethod a string to indicate weight calculation method, there are five choices: 
#' "limma" for for limma linear model based t value,"ttest" for t test based t value, 
#' "MannWhitneyU" for Mann Whitney U based rank-biserial,"PearsonR" for Pearson correlation coefficient,
#' "SpearmanR" for Spearman correlation coefficient, and the defualt value is "ttest"
#' @param topN an integer to indicate how many top features to be selected
#' @return A data frame of selected features' parameters such as weights and FDR
#' @keywords feature selection weight
#' @author Aixiang Jiang
#' @references 
#' Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. 
#' Statistical Applications in traittics and Molecular Biology, Volume 3, Article 3. http://www.statsci.org/smyth/pubs/ebayes.pdf

#' Wendt, H.W. (1972). "Dealing with a common problem in social science: A simplified rank-biserial coefficient of correlation 
#' based on the U statistic". European Journal of Social Psychology. 2 (4): 463–465. doi:10.1002/ejsp.2420020412.
#' 
#' @export

getTrainingWeights = function(trainDat,selectedTraits = NULL, groupInfo, refGroup = 0, topN = NULL, FDRcut = 0.1, weightMethod = c("ttest","limma","PearsonR", "SpearmanR", "MannWhitneyU")) {
  
  ### trainDat: training data set, matrix or data.frame, it should contain all samples in the groupInfo, while it can contain more samples
  ### selectedTraits: if this is known info, otherwise, all genes results will be returned
  ### groupInfo: grouping info of samples with names, the info could be 0-1 format or caterigorical chars
  ### however, the order of the trainDat and groupInfo are not necessary to be the same, I will re-order them any way
  ### alternatively, if the order of groupInfo is aleady matches the order of trainDat, it is fine as well
  
  ### notice that MannWhitneyU or Wilcox W, essentially are the the sum of the signed ranks for one sample
  ### for two groups, the max rank sum for both group U values is u1 + u2=n1*n2, here, n1 and n2 are the sample size of two group
  
  ### this is been said, the default of orderMatches = T, no need to re-order the data,
  ### otherwise need to re-order, in this case, I need the names for groupInfo and colnames for trainDat and they should be comparable
  ###   in this situation, I do not need the same length of samples in the two data, and the order does not matter
  
  weightMethod = weightMethod[1]
  
  #### subset training data if we have the gene list already
  #### make changes on Jan 23, to make sure that I only get common genes
  if(!is.null(selectedTraits)){
    selectedTraits = intersect(selectedTraits, rownames(trainDat))
    if(length(selectedTraits) < 1 | is.na(selectedTraits) | is.null(selectedTraits)){
      print("Please make sure that your selectedTraits has the same format as your gene names in your data")
    }
    trainDat = trainDat[selectedTraits,]
  }
  
  ### change on Dec 3, ignore this part completely
  ### and require the order of samples of groupInfo should be the same as in input data
  ### and no require for the names for items in groupInfo
  #if(!orderMatches){
  #trainDat = trainDat[,names(groupInfo)]
  #}
  
  ### convert groupInfo into 0-1 format, no matter what the refGroup is
  groupInfo = ifelse(groupInfo == refGroup, 0, 1)
  
  #### to avoid trouble, remove the genes if variance is 0
  rem0 = apply(trainDat, 1, sd)
  cut0 = 0.0000001
  rem0 = which(rem0 < cut0)
  
  if(length(rem0)>0){trainDat = trainDat[-rem0,]}
  
  ### there might be many NAs!
  
  ntmp = table(groupInfo)
  g0 = which(groupInfo == 0)
  
  if(weightMethod == "limma"){
    require(limma)
    design = cbind(rep(1,length(groupInfo)),groupInfo)
    
    fit = lmFit(trainDat, design)
    fit = eBayes(fit)
    ### to be consistent to other methods, do not do any selections, but
    ### get t and p values for all features
    res = topTable(fit, coef = "groupInfo", number = dim(trainDat)[1])
    res = res[,3:4]
    colnames(res) = c("limma_t","pValue")
    
  }else if (weightMethod == "MannWhitneyU"){ ## since I am working on two sample test, this is actually Mann-Whitney
    
    ### the following step is much slower than t test, try to use several cpu?
    # res = t(apply(trainDat, 1, function(xx){
    #   tmp = wilcox.test(xx ~ groupInfo)
    #   ### what I need: t stat, p value, add FDR afterwards. I only need t values, however, somebody else might want other info as well
    #   ### return W? and p value, W does not have directions, all of them are pos, then how I determine direcitons? and what to use for weight?
    #   #### seems not a good idea to use Wilcox for weight at all?
    #   #### If I do use it, it is still ok, but I need to get sign from group comparison (same direction as other test), and need to divide it with n1*n2
    #   #### U/(n1*n2) shows the probability that one group is more likely to be higher than another group, which is a good candidate for weight
    #   nn = ntmp[1]*ntmp[2]
    #   xx = rank(xx)
    #   g0mean = sum(xx[g0])/ntmp[1]
    #   g1mean = sum(xx[-g0])/ntmp[2]
    #
    #   ### notice that the smaller rank, the smaller value is
    #   ss = sign(g1mean - g0mean)
    #
    #   ### I actually should report u value from the test group
    #   uout = sum(xx[-g0]) - 0.5*ntmp[2]*(ntmp[2]+1)
    #
    #   return(c(ss*uout/nn, tmp$p.value))
    #   ### notice that the 1st item is not W, but W/(n1*n1) and sign of
    #   ### average testgroup rank  - average refGroup rank
    # }))
    
    #### re-write the above part based on https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test
    # There is a formula to compute the rank-biserial from the Mann–Whitney U and the sample sizes of each group:[11]
    #
    # {\displaystyle r=1-{2U \over n_{1}n_{2}}\,} {\displaystyle r=1-{2U \over n_{1}n_{2}}\,}
    ### to be in the same direction as t test or so, still use my above uout as U
    res = t(apply(trainDat, 1, function(xx){
      
      ## remove NAs
      tmp = which(is.na(xx))
      if(length(tmp) > 0){
        xx = xx[-tmp]
        groupInfo = groupInfo[-tmp]
      }
      
      tmp = wilcox.test(xx ~ groupInfo)
      ### what I need: t stat, p value, add FDR afterwards. I only need t values, however, somebody else might want other info as well
      ### return W? and p value, W does not have directions, all of them are pos, then how I determine direcitons? and what to use for weight?
      #### seems not a good idea to use Wilcox for weight at all?
      #### If I do use it, it is still ok, but I need to get sign from group comparison (same direction as other test), and need to divide it with n1*n2
      #### U/(n1*n2) shows the probability that one group is more likely to be higher than another group, which is a good candidate for weight
      nn = ntmp[1]*ntmp[2]
      xx = rank(xx)
      
      #################### no need any more #######
      #g0mean = sum(xx[g0])/ntmp[1]
      #g1mean = sum(xx[-g0])/ntmp[2]
      ### notice that the smaller rank, the smaller value is
      ##ss = sign(g1mean - g0mean)
      #############################################
      
      ### I actually should report u value from the test group
      uout = sum(xx[-g0]) - 0.5*ntmp[2]*(ntmp[2]+1)
      rbc = 1-2*uout/nn  ### Rank-biserial correlation
      
      return(c(rbc, tmp$p.value))
    }))
    
    colnames(res) = c("Rank_biserial_r","pValue")
    
  }else if (weightMethod == "PearsonR"){
    res = t(apply(trainDat, 1, function(xx){
      ## remove NAs
      tmp = which(is.na(xx))
      if(length(tmp) > 0){
        xx = xx[-tmp]
        groupInfo = groupInfo[-tmp]
      }
      
      ##tmp = t.test(xx ~ groupInfo)  ### the direction seems not correct
      tmp = cor.test(xx, groupInfo)  ### default is pearson
      ### the test data is x, and the reference data is y
      ### what I need: t stat, p value, add FDR afterwards. I only need t values, however, somebody else might want other info as well
      return(c(tmp$estimate, tmp$p.value))  ### return t value and p value
    }))
    
    colnames(res) = c("Pearson_r","pValue")
    
  }else if(weightMethod == "SpearmanR"){
    res = t(apply(trainDat, 1, function(xx){
      ## remove NAs
      tmp = which(is.na(xx))
      if(length(tmp) > 0){
        xx = xx[-tmp]
        groupInfo = groupInfo[-tmp]
      }
      
      ##tmp = t.test(xx ~ groupInfo)  ### the direction seems not correct
      tmp = cor.test(xx, groupInfo, method = "spearman")  ### default is pearson
      ### the test data is x, and the reference data is y
      ### what I need: t stat, p value, add FDR afterwards. I only need t values, however, somebody else might want other info as well
      return(c(tmp$estimate, tmp$p.value))  ### return t value and p value
    }))
    
    # many warnings: In cor.test.default(xx, groupInfo, method = "spearman") :
    #   Cannot compute exact p-value with ties
    
    colnames(res) = c("Spearman_r","pValue")
    
  }else{  ### default is ttest
    
    res = t(apply(trainDat, 1, function(xx){
      ##tmp = t.test(xx ~ groupInfo)  ### the direction seems not correct
      x0 = xx[g0]
      x1 = xx[-g0]
      
      ### remove NAs
      x0 = x0[!is.na(x0)]
      x1 = x1[!is.na(x1)]
      outs = c(NA, NA)
      if(length(x0)>1 & length(x1)>1){
        tmp = t.test(x=x1, y=x0)
        ### the test data is x, and the reference data is y
        outs = c(tmp$statistic, tmp$p.value)
      }
      return(outs)  ### return t value and p value
    }))
    
    colnames(res) = c("tValue","pValue")
    
  }
  
  res = res[!is.na(res[,1]),]
  res = data.frame(res)
  #### add FDR/adj.P.Val in, select topN
  res$FDR = p.adjust(res$pValue, method = "BH")
  res = res[order(res$FDR, decreasing = F),]
  #######################################################################
  ##### when we have the selectedTraits, this chunk is not needed any more
  ### select top N
  ### res = head(res[order(res$FDR, decreasing = F),], topN)
  
  if(!is.null(topN)){
    res = head(res[order(res$FDR, decreasing = F),], topN)
  }
  
  if(!is.null(FDRcut)){
    res = subset(res, res$FDR <= FDRcut)
  }
  
  return(res)
  
}