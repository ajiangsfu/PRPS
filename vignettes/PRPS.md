Type: Package

Title: Binary classification with PRPS and several other choices

Version: 0.1.1

Author: Aixiang Jiang

Maintainer: Aixiang Jiang <aijiang@bccrc.ca> <aixiang.jiang@pathology.ubc.ca>

Depends: R (&gt;= 3.5), lattice, caret, limma, e1071, pROC, mclust

Suggests: knitr

VignetteBuilder: knitr

   

I. Introduction
===============

    PRPS R package is designed to be used for binary classification with three different classification calculation formulas: LPS, PRPS and PS. In addition to training + testing approach, this package also includes self-learning algorithm implemented with PRPS and PS.

1. LPS (Linear Predictor Score);
--------------------------------

Classification scores are calculated with LPS formula (Wright et al., 2003), and probability calculation is based on empirical Bayes method assuming that we know LPS score distribution parameters for the two groups, or if we can estimate these parameters from a training data set. We make binary classification calls based on cutoff on probabilities with default 0.9.

#### References:

Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A trait expression-based method to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S A. 2003 Aug 19;100(17):9991-6.

2. PRPS (Probability ratio based classification predication score);
-------------------------------------------------------------------

Classification scores are calculated with PRPS formula (Ennishi et al., 2019), and probability calculation is based on empirical Bayes method assuming that we know PRPS score distribution parameters for the two groups, or if we can estimate these parameters from a training data set. We make binary classification calls based on cutoff on probabilities with default 0.9. In addition, classification based on theoretical cutoff 0 is also provided.

#### References:

Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P, Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, Kridel R, Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. Double-Hit Trait Expression Signature Defines a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell Lymphoma. Journal of Clinical Oncology 37, no. 3 (January 20, 2019) 190-201. DOI: 10.1200/JCO.18.01583

3. PS (Prediction Strength).
----------------------------

    Classification scores are calculated with PS formula (Golub et al., 1999), and probability calculation is based on empirical Bayes method assuming that we know PS score distribution parameters for the two groups, or if we can estimate these parameters from a training data set. We make binary classification calls based on cutoff on probabilities with default 0.9. In addition, classification based on theoretical cutoff 0 is also provided.

#### References:

TR Golub, DK Slonim, P Tamayo, C Huard, M Gaasenbeek, JP Mesirov, H Coller, ML Loh, JR Downing, MA Caligiuri, et al. Molecular classification of cancer: class discovery and class prediction by gene expression monitoring Science, 286 (1999), pp. 531-537

   

4. Self-learning
----------------

    self-learning algorithm (Jiang et al., 2020) is to identify cases in an unlabeled high dimensional data set that can be confidently classified and to use those cases as pseudo training data thereby allowing classification of all cases in the cohort.    Classification scores are calculated with PRPS (Ennishi et al., 2019) or PS formula (Golub et al., 1999), and probability calculation is based on empirical Bayes method with score distribution parameters from the pseudo training data. We make binary classification calls based on cutoff on probabilities with default 0.9. In addition, classification based on theoretical cutoff 0 is also provided.

#### References:

Jiang A, Hilton LK, Tang J, Rushton CK, Grande BM, Scott DW, Morin RD, Self-learning binary classification and its application to cancer gene expression data (in preparation)

Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P, Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, Kridel R, Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. Double-Hit Trait Expression Signature Defines a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell Lymphoma. Journal of Clinical Oncology 37, no. 3 (January 20, 2019) 190-201. DOI: 10.1200/JCO.18.01583

TR Golub, DK Slonim, P Tamayo, C Huard, M Gaasenbeek, JP Mesirov, H Coller, ML Loh, JR Downing, MA Caligiuri, et al. Molecular classification of cancer: class discovery and class prediction by gene expression monitoring Science, 286 (1999), pp. 531-537

II. Package installation
========================

PRPS R package is currently avaiable at GitHub (<https://github.com/ajiangsfu/PRPS>). It will be available at R (<https://www.r-project.org/>) in the future.

PRPS depends on several existing packages, users should have "limma" installed while other packages can be installed automatically along with PRPS installation. If you do not have "limma", please install it before PRPS installation.

There are at least three ways to install PRPS R package from GitHub.

1. Install PRPS with R package "devtools""
------------------------------------------

2. Install PRPS with R package "remotes"
----------------------------------------

3. Install PRPS with archive files
----------------------------------

There are two archive files for PRPS package available under <https://github.com/ajiangsfu/PRPS>. One file is in source package format for the most updated version (current version: PRPS\_0.1.1.tar.gz), and the other file is in binary package format for the most updated version (current version: PRPS\_0.1.1\_R\_x86\_64-pc-linux-gnu.tar.gz). The source package file can be used in any system, but the binary package file built under Linus cannot be installed in Window system.

After you download the PRPS package file(s), you can install PRPS with the following choices:

### 1) Under R

Packages -&gt; Install Package(s) from local files, then you can select one of your downloaded PRPS package files to install.

### 2) Under Rstudio

Tools -&gt; Install packages, select "Package Archive File(.tar.gz)" from dropdown list of the panel, then you can select one of your downloaded PRPS package files to install.

### 3) Comment line

The comment line to install local package achive file is something like:

install.packages(prpsFile, repos=NULL, type="source")

Here, prpsFile is PRPS package name with your local path, such as "C:/Download/PRPS\_0.1.1.tar.gz".

III. Typical path: training + testing
=====================================

Classification is usually processed when we have comparable training and testing data sets. Training data set is used to select variables to build up a classification model, and classfication is made for a comparable testing data set.

1. Typical workflow when training and testing data sets are comparable
----------------------------------------------------------------------

When training and testing data sets are processed and comparable, the typical classification workflow with PRPS package is:

1). select one of classification algorithms: LPS, PRPS, or PS

2). make decisions on parameter setting for traning functions

3). run LPStraining, or PRPStraining or PStraining on user's training data set

4). run LPStesting, or PRPStesting or PStesting on user's testing data set

2. Example data
---------------

There are five example data files attached with PRPS package.

### 1). rosenwald.cli

This data frame contains subset clinic information about rosenwald dataset, which is downloaded from LPS R package: <https://cran.r-project.org/web/packages/LPS/LPS.pdf>. The original rosenwald data set contains 240 Diffuse Large B-Cell Lymphomas samples (<https://llmpp.nih.gov/DLBCL/>), in this subset data set in LPS R package, however, 40 samples were randomly selected from the training set that are not in Type III sub-types, and 20 samples ere randomly selected from the validation set that are not in Type III sub-types. Here, Type III sub-type refers to samples that are not classified as ABC (activated B cell) or GCB (germinal center), more recently Type III is called as UNC (un-classified) instead.

Together there are 60 Diffuse Large B-Cell Lymphomas samples in rosenwald.cli. The column "set" indicates if a sample was in training or validation data set in the Rosenwald paper, column group is for COO (cell of origin) classification, which could be GCB or ABC since Type III is excluded in LPS R package.

In order to see detail of rosenwald.cli, load PRPS first, and then load rosenwald.cli embedded within PRPS package.

``` r
library(PRPS)
load(system.file("extdata", "rosenwald.cli.rda", package = "PRPS"))
str(rosenwald.cli)
#> 'data.frame':    60 obs. of  4 variables:
#>  $ set      : Factor w/ 2 levels "Training","Validation": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ group    : Factor w/ 2 levels "ABC","GCB": 2 1 1 2 2 2 1 2 1 1 ...
#>  $ follow.up: num  2.4 1 10.5 1 1.6 0.2 10.8 2.3 8.4 1 ...
#>  $ status   : Factor w/ 2 levels "Alive","Dead": 2 2 1 2 2 2 1 2 1 2 ...
```

#### References:

Rosenwald A, Wright G, Chan WC, et al. The use of molecular profiling to predict survival after chemotherapy for diffuse large-B-cell lymphoma. N Engl J Med 2002;346(25):1937-1947

<https://github.com/maressyl/R.LPS/tree/master/LPS>

### 2). rosenwald.expr

This data matrix contains Lymphochip microarrays expression data for the 60 Diffuse Large B-Cell Lymphomas samples as described in rosenwald.cli. To see the detail of rosenwald.expr, load it embedded within PRPS package.

``` r
load(system.file("extdata", "rosenwald.expr.rda", package = "PRPS"))
str(rosenwald.expr)
#>  num [1:7399, 1:60] 0.607 0.239 0.828 0.932 0.86 0 0.476 0.814 0.81 0.87 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:7399] "27481" "17013" "24751" "27498" ...
#>   ..$ : chr [1:60] "LYM018" "LYM120" "LYM180" "LYM285" ...
```

#### References:

Rosenwald A, Wright G, Chan WC, et al. The use of molecular profiling to predict survival after chemotherapy for diffuse large-B-cell lymphoma. N Engl J Med 2002;346(25):1937-1947

<https://github.com/maressyl/R.LPS/tree/master/LPS>

### 3). WrightCOO

This data frame contains 158 COO (cell of origin) related genes with both Ensembl annotation ID (row names) and gene symbols (column "Gene"). "The current standard approach for routinely classifying samples using Affymetrix U133 arrays employs 186 probesets (George Wright, personal communication). The 165 Ensembl genes that correspond to these probesets were used for classification by RNA-seq." (Morin et al., 2011)". In our recent research, we matched 158 genes in Ennishi cohort (Ennishi et al., 2019), which are used in self-learning research (Jiang et al., 2020) and are included in PRPS package.

In order to see detail of WrightCOO, load WrightCOO embedded within PRPS package.

``` r
load(system.file("extdata", "WrightCOO.rda", package = "PRPS"))
str(WrightCOO)
#> 'data.frame':    158 obs. of  1 variable:
#>  $ Gene: chr  "NR3C1" "PMM2" "STS" "BCL2" ...
head(WrightCOO)
#>                  Gene
#> ENSG00000113580 NR3C1
#> ENSG00000140650  PMM2
#> ENSG00000101846   STS
#> ENSG00000171791  BCL2
#> ENSG00000168811 IL12A
#> ENSG00000196549   MME
```

Here, we provide both ENSG IDs and gene symbols for this Wright COO list.

#### References:

Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A trait expression-based method to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S A. 2003 Aug 19;100(17):9991-6.

Morin, R. D., Mendez-Lago, M., Mungall, A. J., Goya, R., Mungall, K. L., Corbett, R. D., ... Marra, M. A. (2011). Frequent mutation of histone-modifying genes in non-Hodgkin lymphoma. Nature, 476(7360), 298-303. <https://doi.org/10.1038/nature10351>

Jiang A, Hilton LK, Tang J, Rushton CK, Grande BM, Scott DW, Morin RD, Self-learning binary classification and its application to cancer gene expression data (in preparation)

### 4). DHITsig and DHITsigENSG

DHITsig is a vector containing 104 DHITsig genes and their importance scores (Ennishi et al., 2019). In order to see detail of DHITsig, load DHITsig embedded within PRPS package.

``` r
load(system.file("extdata", "DHITsig.rda", package = "PRPS"))
str(DHITsig)
#>  Named num [1:104] 0.674 0.666 0.618 0.597 0.582 ...
#>  - attr(*, "names")= chr [1:104] "OR13A1" "FAM216A" "MYC" "SLC25A27" ...
head(DHITsig)
#>    OR13A1   FAM216A       MYC  SLC25A27     ALOX5     UQCRH 
#> 0.6742184 0.6662736 0.6180968 0.5973289 0.5822841 0.5545504
```

We can see that each item's name is in gene symbol format. If your data are in ENSG ID format, you need to load DHITsigENSG instead:

``` r
load(system.file("extdata", "DHITsigENSG.rda", package = "PRPS"))
str(DHITsigENSG)
#>  Named num [1:104] 0.674 0.666 0.618 0.597 0.582 ...
#>  - attr(*, "names")= chr [1:104] "ENSG00000256574" "ENSG00000204856" "ENSG00000136997" "ENSG00000153291" ...
head(DHITsigENSG)
#> ENSG00000256574 ENSG00000204856 ENSG00000136997 ENSG00000153291 ENSG00000012779 
#>       0.6742184       0.6662736       0.6180968       0.5973289       0.5822841 
#> ENSG00000173660 
#>       0.5545504
```

DHITsigENSG contains exactly the same genes as in DHITsig and their importance scores, the only difference is that its item names are in ENSG IDs.

#### References:

Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P, Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, Kridel R, Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. Double-Hit Trait Expression Signature Defines a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell Lymphoma. J Clin Oncol. 2018 Dec 3:JCO1801583. doi: 10.1200/JCO.18.01583.

3. Code Examples
----------------

In this section, we are providing code examples to get classification score when training and testing data sets are comparable. The example data set is rosenwald.expr, which is split into training and testing data sets based on set information in rosenwald.cli. Since each person might have his/her own preference, code examples with all three algorithms: LPS, PRPS, and PS, are provided in the following.

### 1) LPS

If LPS is selected, we first apply LPStraining function to find top genes and their weights to distinguish the two groups given in the training data sets, LPS score distribution parameters are also included in the training object. Then, we call LPStesting function to classify samples in the comparable testing data set with training ovject output.

#### Training with LPStraining

First of all, divide data into training and testing data sets.

``` r
dat = rosenwald.expr
trainset = subset(rosenwald.cli, rosenwald.cli$set == "Training")
testset = subset(rosenwald.cli, rosenwald.cli$set == "Validation")
```

Use GCB as a reference group, then LPStraining object can be produced with the follwoing code assuming that we select top 100 genes for COO classification. There are five different ways to define gene weight, here, the default t test is used.

``` r
nin = 100 
trainLPS = LPStraining (trainDat = dat[,rownames(trainset)], groupInfo = trainset$group, refGroup = "GCB", topN = nin,
                      weightMethod = "ttest")
str(trainLPS)
#> List of 4
#>  $ LPS_pars    :List of 2
#>   ..$ weights:'data.frame':  100 obs. of  3 variables:
#>   .. ..$ tValue: num [1:100] 7.94 7.43 7.14 6.86 6.79 ...
#>   .. ..$ pValue: num [1:100] 1.66e-09 6.59e-09 1.59e-08 4.85e-08 6.07e-08 ...
#>   .. ..$ FDR   : num [1:100] 1.23e-05 2.44e-05 3.93e-05 8.98e-05 8.98e-05 ...
#>   ..$ meansds: Named num [1:4] 373 -321 120 157
#>   .. ..- attr(*, "names")= chr [1:4] "testLPSmean" "refLPSmean" "testLPSsd" "refLPSsd"
#>  $ LPS_train   :'data.frame':    40 obs. of  5 variables:
#>   ..$ LPS_score    : num [1:40] -381.8 238.9 394.9 20.9 -306.2 ...
#>   ..$ true_class   : chr [1:40] "GCB" "ABC" "ABC" "GCB" ...
#>   ..$ LPS_class    : chr [1:40] "GCB" "ABC" "ABC" "UNCLASS" ...
#>   ..$ LPS_prob_test: num [1:40] 3.42e-09 9.98e-01 1.00 1.58e-01 1.39e-07 ...
#>   ..$ LPS_prob_ref : num [1:40] 1.00 2.45e-03 2.36e-05 8.42e-01 1.00 ...
#>  $ classCompare:List of 6
#>   ..$ positive: chr "ABC"
#>   ..$ table   : 'table' int [1:2, 1:2] 22 0 0 16
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ Prediction: chr [1:2] "GCB" "ABC"
#>   .. .. ..$ Reference : chr [1:2] "GCB" "ABC"
#>   ..$ overall : Named num [1:7] 1 1 0.907 1 0.579 ...
#>   .. ..- attr(*, "names")= chr [1:7] "Accuracy" "Kappa" "AccuracyLower" "AccuracyUpper" ...
#>   ..$ byClass : Named num [1:11] 1 1 1 1 1 ...
#>   .. ..- attr(*, "names")= chr [1:11] "Sensitivity" "Specificity" "Pos Pred Value" "Neg Pred Value" ...
#>   ..$ mode    : chr "sens_spec"
#>   ..$ dots    : list()
#>   ..- attr(*, "class")= chr "confusionMatrix"
#>  $ classTable  : 'table' int [1:2, 1:3] 0 16 22 0 1 1
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ groupInfo: chr [1:2] "GCB" "ABC"
#>   .. ..$ LPS_class: chr [1:3] "ABC" "GCB" "UNCLASS"
```

The last item in the LPStraining object provides information on how well our LPS prediction model is doing.

``` r
trainLPS$classTable
#>          LPS_class
#> groupInfo ABC GCB UNCLASS
#>       GCB   0  22       1
#>       ABC  16   0       1
```

This means with given COO info, and with top 100 significant genes based on t test, LPS classification can put all GCB samples back to GCB samples as expected, for ABC samples, however, one of them is put into UNCLASS, all other samples are put back into ABC as expected.

More detail information about LPStraining function can be found with the two options:

``` r
?(LPStraining)
help(LPStraining)
```

#### LPStesting

Now, we can use our training object to make COO classification for a comparable testing data set. In the current example, this is simple to achieve since we already have testing data set shown above. We just use the output (training object) from LPStraining together with testing data set to call LPStesting to get result, and compare it to the original COO classification:

``` r
testLPS = LPStesting(LPStrainObj = trainLPS, newdat = dat[,rownames(testset)])
str(testLPS)
#> 'data.frame':    20 obs. of  4 variables:
#>  $ LPS_score    : num  -81.3 -99.7 -217.2 215.3 358.5 ...
#>  $ LPS_class    : chr  "GCB" "GCB" "GCB" "ABC" ...
#>  $ LPS_prob_test: num  3.19e-03 1.48e-03 8.81e-06 9.95e-01 1.00 ...
#>  $ LPS_prob_ref : num  9.97e-01 9.99e-01 1.00 5.25e-03 6.55e-05 ...
table(testLPS$LPS_class)
#> 
#>     ABC     GCB UNCLASS 
#>       6      13       1
table(testLPS$LPS_class, testset$group)
#>          
#>           ABC GCB
#>   ABC       6   0
#>   GCB       0  13
#>   UNCLASS   0   1
```

We can also combine training and testing data sets together for testing:

``` r
fullLPS = LPStesting(LPStrainObj = trainLPS, newdat = dat)
str(fullLPS)
#> 'data.frame':    60 obs. of  4 variables:
#>  $ LPS_score    : num  -381.8 238.9 394.9 20.9 -306.2 ...
#>  $ LPS_class    : chr  "GCB" "ABC" "ABC" "UNCLASS" ...
#>  $ LPS_prob_test: num  3.42e-09 9.98e-01 1.00 1.58e-01 1.39e-07 ...
#>  $ LPS_prob_ref : num  1.00 2.45e-03 2.36e-05 8.42e-01 1.00 ...
table(fullLPS$LPS_class)
#> 
#>     ABC     GCB UNCLASS 
#>      22      35       3
```

Now, we are wondering if the testLPS and fulltest give us consistent results:

``` r
table(testLPS$LPS_class, fullLPS[rownames(testset), "LPS_class"])
#>          
#>           ABC GCB UNCLASS
#>   ABC       6   0       0
#>   GCB       0  13       0
#>   UNCLASS   0   0       1
```

More detail about LPStesting can be found with the following two options:

``` r
?(LPStesting)
help(LPStesting)
```

### 2) PRPS

If PRPS is selected, we first apply PRPStraining function to find top genes and their weights to distinguish the two groups given in the training data sets, distribution parameters of these genes and PRPS score are also included in the training object. Then, we call PRPStesting function to classify samples in the comparable testing data set with training ovject output.

#### PRPStraining

Use the same training data set as LPStraining, set GCB as the reference group, PRPS training object can be got as following. Notice that now we are calling PRPStraining but all setting is the same as for LPStraining shown above.

``` r
nin = 100 
trainPRPS = PRPStraining (trainDat = dat[,rownames(trainset)], groupInfo = trainset$group, refGroup = "GCB", topN = nin,
                      weightMethod = "ttest")
str(trainPRPS)
#> List of 4
#>  $ PRPS_pars   :List of 4
#>   ..$ weights      :'data.frame':    100 obs. of  3 variables:
#>   .. ..$ tValue: num [1:100] 7.94 7.43 7.14 6.86 6.79 ...
#>   .. ..$ pValue: num [1:100] 1.66e-09 6.59e-09 1.59e-08 4.85e-08 6.07e-08 ...
#>   .. ..$ FDR   : num [1:100] 1.23e-05 2.44e-05 3.93e-05 8.98e-05 8.98e-05 ...
#>   ..$ meansds      : Named num [1:4] 366 -508 126 179
#>   .. ..- attr(*, "names")= chr [1:4] "testPRPSmean" "refPRPSmean" "testPRPSsd" "refPRPSsd"
#>   ..$ traitsmeansds: num [1:100, 1:4] 0.482 0.389 0.434 0.715 0.351 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:100] "17708" "27631" "24796" "26697" ...
#>   .. .. ..$ : chr [1:4] "testmean" "refmean" "testsd" "refsd"
#>   ..$ dfs          : num [1:2] 16 22
#>  $ PRPS_train  :'data.frame':    40 obs. of  6 variables:
#>   ..$ PRPS_score    : num [1:40] -548 237 426 -187 -533 ...
#>   ..$ true_class    : chr [1:40] "GCB" "ABC" "ABC" "GCB" ...
#>   ..$ PRPS_class    : chr [1:40] "GCB" "ABC" "ABC" "GCB" ...
#>   ..$ PRPS_prob_test: num [1:40] 5.48e-12 1.00 1.00 4.72e-04 1.26e-11 ...
#>   ..$ PRPS_prob_ref : num [1:40] 1.00 2.11e-04 9.79e-07 1.00 1.00 ...
#>   ..$ PRPS_class0   : chr [1:40] "GCB" "ABC" "ABC" "GCB" ...
#>  $ classCompare:List of 6
#>   ..$ positive: chr "ABC"
#>   ..$ table   : 'table' int [1:2, 1:2] 23 0 0 16
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ Prediction: chr [1:2] "GCB" "ABC"
#>   .. .. ..$ Reference : chr [1:2] "GCB" "ABC"
#>   ..$ overall : Named num [1:7] 1 1 0.91 1 0.59 ...
#>   .. ..- attr(*, "names")= chr [1:7] "Accuracy" "Kappa" "AccuracyLower" "AccuracyUpper" ...
#>   ..$ byClass : Named num [1:11] 1 1 1 1 1 ...
#>   .. ..- attr(*, "names")= chr [1:11] "Sensitivity" "Specificity" "Pos Pred Value" "Neg Pred Value" ...
#>   ..$ mode    : chr "sens_spec"
#>   ..$ dots    : list()
#>   ..- attr(*, "class")= chr "confusionMatrix"
#>  $ classTable  : 'table' int [1:2, 1:3] 0 16 23 0 0 1
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ groupInfo : chr [1:2] "GCB" "ABC"
#>   .. ..$ PRPS_class: chr [1:3] "ABC" "GCB" "UNCLASS"
```

Still, use default t test to get top 100 gene weights. The last item in the output of PRPStraining provides information on how our PRPS prediction model is doing

``` r
trainPRPS$classTable
#>          PRPS_class
#> groupInfo ABC GCB UNCLASS
#>       GCB   0  23       0
#>       ABC  16   0       1
```

This means with given COO info, and with top 100 significant genes based on t test, PRPS classification is 100% matching to given COO.

More detail information about LPStraining function can be found with the two choices:

``` r
?(PRPStraining)
help(PRPStraining)
```

#### PRPStesting

With the exact same testing data set as above, use the output from PRPStraining, we can apply PRPStesting as following and compare the classification with the truth in the testing data info:

``` r
testPRPS = PRPStesting(PRPStrainObj = trainPRPS, newdat = dat[,rownames(testset)])
str(testPRPS)
#> 'data.frame':    20 obs. of  5 variables:
#>  $ PRPS_score    : num  -187 -251 -411 145 333 ...
#>  $ PRPS_class    : chr  "GCB" "GCB" "GCB" "ABC" ...
#>  $ PRPS_prob_test: num  4.62e-04 2.49e-05 9.28e-09 9.96e-01 1.00 ...
#>  $ PRPS_prob_ref : num  0.999538 0.999975 1 0.004249 0.000012 ...
#>  $ PRPS_class0   : chr  "GCB" "GCB" "GCB" "ABC" ...
table(testPRPS$PRPS_class)
#> 
#>     ABC     GCB UNCLASS 
#>       6      13       1
table(testPRPS$PRPS_class, testset$group)
#>          
#>           ABC GCB
#>   ABC       6   0
#>   GCB       0  13
#>   UNCLASS   0   1
```

We can also combine training and testing data sets together for testing

``` r
fullPRPS = PRPStesting(PRPStrainObj = trainPRPS, newdat = dat)
str(fullPRPS)
#> 'data.frame':    60 obs. of  5 variables:
#>  $ PRPS_score    : num  -548 237 426 -187 -533 ...
#>  $ PRPS_class    : chr  "GCB" "ABC" "ABC" "GCB" ...
#>  $ PRPS_prob_test: num  5.48e-12 1.00 1.00 4.72e-04 1.26e-11 ...
#>  $ PRPS_prob_ref : num  1.00 2.11e-04 9.79e-07 1.00 1.00 ...
#>  $ PRPS_class0   : chr  "GCB" "ABC" "ABC" "GCB" ...
table(fullPRPS$PRPS_class)
#> 
#>     ABC     GCB UNCLASS 
#>      22      36       2
```

Now, we are wondering if the testLPS and fulltest give us consistent results.

``` r
table(testPRPS$PRPS_class, fullPRPS[rownames(testset), "PRPS_class"])
#>          
#>           ABC GCB UNCLASS
#>   ABC       6   0       0
#>   GCB       0  13       0
#>   UNCLASS   0   0       1
```

If we compare PRPS results with LPS results, we can find that the training data itself is 100% matched to the truth for LPS but there is one sample put into UNCLASS based on PRPS. However, for the testing data set, there is one GCB sample was classified to ABC by LPS, while ihe same samppe is classified into UNCLASS by PRPS, which seems more reasonable.

### 3) PS

#### PStraining

Use the same training data set as above, and use GCB as a reference group as in LPS package, then PStraining result can be got as following. Notice that now we are calling PStraining but all setting is the same as for LPStraining

``` r
nin = 100 
trainPS = PStraining (trainDat = dat[,rownames(trainset)], groupInfo = trainset$group, refGroup = "GCB", topN = nin, weightMethod = "ttest")
str(trainPS)
#> List of 4
#>  $ PS_pars     :List of 3
#>   ..$ weights:'data.frame':  100 obs. of  3 variables:
#>   .. ..$ tValue: num [1:100] 7.94 7.43 7.14 6.86 6.79 ...
#>   .. ..$ pValue: num [1:100] 1.66e-09 6.59e-09 1.59e-08 4.85e-08 6.07e-08 ...
#>   .. ..$ FDR   : num [1:100] 1.23e-05 2.44e-05 3.93e-05 8.98e-05 8.98e-05 ...
#>   ..$ meansds: Named num [1:4] 0.813 -0.752 0.169 0.228
#>   .. ..- attr(*, "names")= chr [1:4] "testPSmean" "refPSmean" "testPSsd" "refPSsd"
#>   ..$ traits : num [1:100, 1:3] 0.0999 0.1135 0.1108 0.1107 0.0928 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:100] "17708" "27631" "24796" "26697" ...
#>   .. .. ..$ : chr [1:3] "meanOfGroupMeans" "refGroupMean" "testGroupMean"
#>  $ PS_train    :'data.frame':    40 obs. of  6 variables:
#>   ..$ PS_score    : num [1:40] -0.783 0.53 0.866 -0.27 -0.727 ...
#>   ..$ true_class  : chr [1:40] "GCB" "ABC" "ABC" "GCB" ...
#>   ..$ PS_class    : chr [1:40] "GCB" "ABC" "ABC" "GCB" ...
#>   ..$ PS_prob_test: num [1:40] 7.12e-20 1.00 1.00 1.68e-08 1.54e-18 ...
#>   ..$ PS_prob_ref : num [1:40] 1.00 3.93e-07 8.58e-12 1.00 1.00 ...
#>   ..$ PS_class0   : chr [1:40] "GCB" "ABC" "ABC" "GCB" ...
#>  $ classCompare:List of 6
#>   ..$ positive: chr "ABC"
#>   ..$ table   : 'table' int [1:2, 1:2] 23 0 0 17
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ Prediction: chr [1:2] "GCB" "ABC"
#>   .. .. ..$ Reference : chr [1:2] "GCB" "ABC"
#>   ..$ overall : Named num [1:7] 1 1 0.912 1 0.575 ...
#>   .. ..- attr(*, "names")= chr [1:7] "Accuracy" "Kappa" "AccuracyLower" "AccuracyUpper" ...
#>   ..$ byClass : Named num [1:11] 1 1 1 1 1 1 1 0.425 0.425 0.425 ...
#>   .. ..- attr(*, "names")= chr [1:11] "Sensitivity" "Specificity" "Pos Pred Value" "Neg Pred Value" ...
#>   ..$ mode    : chr "sens_spec"
#>   ..$ dots    : list()
#>   ..- attr(*, "class")= chr "confusionMatrix"
#>  $ classTable  : 'table' int [1:2, 1:2] 0 17 23 0
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ groupInfo: chr [1:2] "GCB" "ABC"
#>   .. ..$ PS_class : chr [1:2] "ABC" "GCB"
```

The last item in the output of PStraining provides information on how our PS prediction model is doing

``` r
trainPS$classTable
#>          PS_class
#> groupInfo ABC GCB
#>       GCB   0  23
#>       ABC  17   0
```

This means with given COO info, and with top 100 significant genes based on t test, PS classification is 100% matching to given COO.

More detail information about LPStraining function can be found with the two choices:

``` r
?(PStraining)
help(PStraining)
```

#### PStesting

With the exact same testing data set as above, use the output from PStraining, we can apply PStesting as following:

``` r
testPS = PStesting(PStrainObj = trainPS, newdat = dat[,rownames(testset)])
str(testPS)
#> 'data.frame':    20 obs. of  5 variables:
#>  $ PS_score    : num  -0.23 -0.367 -0.669 0.609 0.788 ...
#>  $ PS_class    : chr  "GCB" "GCB" "GCB" "ABC" ...
#>  $ PS_prob_test: num  1.31e-03 8.79e-05 1.66e-07 1.00 1.00 ...
#>  $ PS_prob_ref : num  9.99e-01 1.00 1.00 3.71e-04 2.58e-05 ...
#>  $ PS_class0   : chr  "GCB" "GCB" "GCB" "ABC" ...
table(testPS$PS_class)
#> 
#>     ABC     GCB UNCLASS 
#>       6      13       1
table(testPS$PS_class, testset$group)
#>          
#>           ABC GCB
#>   ABC       6   0
#>   GCB       0  13
#>   UNCLASS   0   1
```

We can also combine training and testing data sets together for testing

``` r
fullPS = PStesting(PStrainObj = trainPS, newdat = dat)
str(fullPS)
#> 'data.frame':    60 obs. of  5 variables:
#>  $ PS_score    : num  -0.781 0.51 0.883 -0.224 -0.701 ...
#>  $ PS_class    : chr  "GCB" "ABC" "ABC" "GCB" ...
#>  $ PS_prob_test: num  4.46e-12 1.00 1.00 1.01e-04 5.74e-11 ...
#>  $ PS_prob_ref : num  1.00 6.17e-05 1.83e-08 1.00 1.00 ...
#>  $ PS_class0   : chr  "GCB" "ABC" "ABC" "GCB" ...
table(fullPS$PS_class)
#> 
#>     ABC     GCB UNCLASS 
#>      23      36       1
```

Now, we are wondering if the testLPS and fulltest give us consistent results.

``` r
table(testPS$PS_class, fullPS[rownames(testset), "PS_class"])
#>          
#>           ABC GCB UNCLASS
#>   ABC       6   0       0
#>   GCB       0  13       0
#>   UNCLASS   0   0       1
```

PS performs exactly the same as LPS for the above training and testing data sets.

   

IIII. Special path: classification score calculation without comparable training
================================================================================

In many situations, we might not have comparable training and testing data sets, in this case, we might still calculate classification scores and make classification calls.

1. A similar data set with classificaiton info is available
-----------------------------------------------------------

If there is no comparable training data available, but there is a similar to testing data with classification info available, we might apply calibration or normalization or standardization to make it as our pseudo training data set. In this case, we can do the following:

1.  If we believe that there is only some score shifting between training and testing data sets, calibrate classification scores between testing and traning data sets, and then use the calibrated mean and sd of two groups in the traning data set to calculate Empirical Bayes probabilities for the testing data set, and furthermore get classification calls. Notice that this approach is easy for LPS method but not easy for PROS and PS. This is because PRPS and PS require feature level group info as well, which is more difficult to deal with.

2.  If the assumption in i) does not hold, which means that there are more difference between training and testing data sets, and calibration on classification scores cannot overcome the difference, but if we have house keeping genes, we can use house keeping genes to normalize both training and testing data sets, and use the normalized training and testing data sets and follow instruction in II Typical path: training + testing. This normalization strategy can be applied for all three algorithm in this package: LPS, PRPS and PS. The assumption, however, is that these house keeping genes perform similar across training and testing data sets. When this assumption does not hold, this approach will not work.

3.  If both i) and ii) are not applicable, which means that the calibration on classification scores does not work and we cannot use house keeping genes for normalization, then, we can apply feature-wise standardization for both pseudo training and testing data sets. In this case, when we follow instruction in II Typical path: training + testing, we just need to make sure standardization = TRUE. This approach can be worked for all three algorithms in this package: LPS, PRPS and PS. The assumption is, however, the ture distribution for each feature is similar between training and testing data sets. If not, this approach will not work. For example, if there is a significant higher proportion of one group in the testing data set than in the training data set, normalization method will not work.

``` r
lpstrain = LPStraining (trainDat = dat[,rownames(trainset)], standardization=TRUE, groupInfo = trainset$group, refGroup = "GCB", topN = nin,weightMethod = "ttest")
str(trainLPS)
#> List of 4
#>  $ LPS_pars    :List of 2
#>   ..$ weights:'data.frame':  100 obs. of  3 variables:
#>   .. ..$ tValue: num [1:100] 7.94 7.43 7.14 6.86 6.79 ...
#>   .. ..$ pValue: num [1:100] 1.66e-09 6.59e-09 1.59e-08 4.85e-08 6.07e-08 ...
#>   .. ..$ FDR   : num [1:100] 1.23e-05 2.44e-05 3.93e-05 8.98e-05 8.98e-05 ...
#>   ..$ meansds: Named num [1:4] 373 -321 120 157
#>   .. ..- attr(*, "names")= chr [1:4] "testLPSmean" "refLPSmean" "testLPSsd" "refLPSsd"
#>  $ LPS_train   :'data.frame':    40 obs. of  5 variables:
#>   ..$ LPS_score    : num [1:40] -381.8 238.9 394.9 20.9 -306.2 ...
#>   ..$ true_class   : chr [1:40] "GCB" "ABC" "ABC" "GCB" ...
#>   ..$ LPS_class    : chr [1:40] "GCB" "ABC" "ABC" "UNCLASS" ...
#>   ..$ LPS_prob_test: num [1:40] 3.42e-09 9.98e-01 1.00 1.58e-01 1.39e-07 ...
#>   ..$ LPS_prob_ref : num [1:40] 1.00 2.45e-03 2.36e-05 8.42e-01 1.00 ...
#>  $ classCompare:List of 6
#>   ..$ positive: chr "ABC"
#>   ..$ table   : 'table' int [1:2, 1:2] 22 0 0 16
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ Prediction: chr [1:2] "GCB" "ABC"
#>   .. .. ..$ Reference : chr [1:2] "GCB" "ABC"
#>   ..$ overall : Named num [1:7] 1 1 0.907 1 0.579 ...
#>   .. ..- attr(*, "names")= chr [1:7] "Accuracy" "Kappa" "AccuracyLower" "AccuracyUpper" ...
#>   ..$ byClass : Named num [1:11] 1 1 1 1 1 ...
#>   .. ..- attr(*, "names")= chr [1:11] "Sensitivity" "Specificity" "Pos Pred Value" "Neg Pred Value" ...
#>   ..$ mode    : chr "sens_spec"
#>   ..$ dots    : list()
#>   ..- attr(*, "class")= chr "confusionMatrix"
#>  $ classTable  : 'table' int [1:2, 1:3] 0 16 22 0 1 1
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ groupInfo: chr [1:2] "GCB" "ABC"
#>   .. ..$ LPS_class: chr [1:3] "ABC" "GCB" "UNCLASS"

lpstest = LPStesting(lpstrain, newdat = dat[,rownames(testset)], standardization = TRUE)
str(lpstest)
#> 'data.frame':    20 obs. of  4 variables:
#>  $ LPS_score    : num  -35 -126 -218 218 366 ...
#>  $ LPS_class    : chr  "GCB" "GCB" "GCB" "ABC" ...
#>  $ LPS_prob_test: num  9.22e-03 1.93e-04 3.02e-06 9.94e-01 1.00 ...
#>  $ LPS_prob_ref : num  9.91e-01 1.00 1.00 6.27e-03 4.29e-05 ...
```

In general, this is not a big problem is there is linear shift in classification score level, or even in feature level. But we do have problem to deal with data if distribution is different between training and testing.

2. No similar data set with classificaiton info is available
------------------------------------------------------------

Not only we do not have comparable training and testing data sets, but also we do not have similar training and testing data sets. This means that we do not have any kinds of training data sets at all. In this situation, we should at least have gene features/traits that we would like to distiguish samples. Otherwise no classification score can be calculated.

### 1). Weights for selected features/traits are available

If we know selected features/traits with their weights but no any kind of training data set, and we would like to calculate classification score for a given tesing data set, what can we do?

#### LPS

Without any further information and without any assumption, the only method we can use to calculate classification score is to apply LPS approach. In this case, we can call getClassScores and set method as LPS.

Assuming that we have features and their weights from some where, for example, we can get this information from trainLPS:

``` r
lpswts = trainLPS$LPS_pars$weights
str(lpswts)
#> 'data.frame':    100 obs. of  3 variables:
#>  $ tValue: num  7.94 7.43 7.14 6.86 6.79 ...
#>  $ pValue: num  1.66e-09 6.59e-09 1.59e-08 4.85e-08 6.07e-08 ...
#>  $ FDR   : num  1.23e-05 2.44e-05 3.93e-05 8.98e-05 8.98e-05 ...
```

Here, the tValue is the weight for the top 100 genes selected with LPStraining.

``` r
genes = rownames(lpswts)
lpswts = lpswts[,1]
names(lpswts) = genes
head(lpswts)
#>     17708     27631     24796     26697     27562     17496 
#>  7.938396  7.431466  7.144285  6.859320  6.793850 -6.606511
```
