#' @export
getTrainingParas = function(trainDat, selectedTraits = NULL, groupInfo, refGroup = 0, weightMethod = c("ttest","limma","PearsonR", "SpearmanR", "MannWhitneyU"), topN = NULL){
  ### trainDat: training data set, matrix or data.frame, it should contain all samples in the groupInfo, while it can contain more samples
  ### groupInfo: grouping info of samples with names, the info could be 0-1 format or caterigorical chars
  ### however, the order of the trainDat and groupInfo are not necessary to be the same, I will re-order them any way
  ### alternatively, if the order of groupInfo is aleady matches the order of trainDat, it is fine as well
  
  ### notice that MannWhitneyU or Wilcox W, essentially are the the sum of the signed ranks for one sample
  ### for two groups, the max rank sum for both group U values is u1 + u2=n1*n2, here, n1 and n2 are the sample size of two group
  
  weightMethod = weightMethod[1]
  
  ### this is been said, the default of orderMatches = T, no need to re-order the data,
  ### otherwise need to re-order, in this case, I need the names for groupInfo and colnames for trainDat and they should be comparable
  ###   in this situation, I do not need the same length of samples in the two data, and the order does not matter
  if(!is.null(selectedTraits)){
    trainDat = trainDat[selectedTraits,]
  }
  
  trainDat = trainDat[,names(groupInfo)]
  
  ### convert groupInfo into 0-1 format, no matter what the refGroup is
  groupInfo = ifelse(groupInfo == refGroup, 0, 1)
  
  #### to avoid trouble, remove the traits if variance is 0
  rem0 = apply(trainDat, 1, sd, na.rm = T)
  cut0 = 0.0000001
  rem0 = which(rem0 < cut0)
  
  if(length(rem0)>0){trainDat = trainDat[-rem0,]}
  
  ### there might be many NAs!
  
  ntmp = table(groupInfo)
  g0 = which(groupInfo == 0)
  
  if(weightMethod == "limma"){
    ##require(limma)
    design = cbind(rep(1,length(groupInfo)),groupInfo)
    
    fit = limma::lmFit(trainDat, design)
    fit = limma::eBayes(fit)
    ### to be consistent to other methods, do not do any selections, but
    ### get t and p values for all features
    res = limma::topTable(fit, coef = "groupInfo", number = dim(trainDat)[1])
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
      ##tmp = t.test(xx ~ groupInfo)  ### the direction seems not correct
      tmp = cor.test(xx, groupInfo)  ### default is pearson
      ### the test data is x, and the reference data is y
      ### what I need: t stat, p value, add FDR afterwards. I only need t values, however, somebody else might want other info as well
      return(c(tmp$estimate, tmp$p.value))  ### return t value and p value
    }))
    
    colnames(res) = c("Pearson_r","pValue")
    
  }else if(weightMethod == "SpearmanR"){
    res = t(apply(trainDat, 1, function(xx){
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
  ### select top N
  
  if(!is.null(topN)){
    res = head(res[order(res$FDR, decreasing = F),], topN)
  }
  
  #############################################
  ### then add in two groups mean and sd for the traits in res
  ### actually, I can combine these infor with res and output one data.frame together
  ############################################
  
  sigdat = trainDat[rownames(res),]
  g0dat = sigdat[,g0]
  g1dat = sigdat[, -g0]
  
  g0mean = apply(g0dat,1,mean, na.rm=T)
  g1mean = apply(g1dat,1,mean, na.rm=T)
  
  g0sd = apply(g0dat,1,sd, na.rm=T)
  g1sd = apply(g1dat,1,sd, na.rm=T)
  
  res = cbind(res, g1mean, g0mean, g1sd, g0sd)
  colnames(res)[4:7] = c("group1mean","group0mean", "group1sd","group0sd")
  
  return(res)
  
}
####################################################################
