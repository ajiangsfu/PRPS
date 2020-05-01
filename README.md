PRPS: Binary classification with PRPS, PRPS-ST and more
=======================================================

   

I. Binary classfication methods
-------------------------------

PRPS R package is designed to be used for binary classification with five method choices, they are PRPS (Probability ratio based classification predication score), PRPS-ST (Self-training with Probability ratio based classification predication score), LPS (Linear Predictor Score), PS (Prediction Strength), and PS-ST (Self-training with Prediction Strength).

### 1. PRPS (Probability ratio based classification predication score)

Classification scores are calculated with PRPS formula (Ennishi et al., 2019), and probability calculation is based on empirical Bayes method assuming that we know PRPS score distribution parameters for the two groups, or if we can estimate these parameters from a training data set. We make binary classification calls based on cutoff on probabilities with default 0.9. In addition, classification based on theoretical cutoff 0 is also provided.

#### References:

Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P, Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, Kridel R, Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. Double-Hit Trait Expression Signature Defines a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell Lymphoma. Journal of Clinical Oncology 37, no. 3 (January 20, 2019) 190-201. DOI: 10.1200/JCO.18.01583

### 2. PRPS-ST (Self-training with Probability ratio based classification predication score)

Self-training algorithm (Jiang et al., 2020) is to identify cases in an unlabeled high dimensional data set that can be confidently classified and to use those cases as pseudo training data thereby allowing classification of all cases in the cohort. For PRPS-ST, classification scores are calculated with PRPS (Ennishi et al., 2019), and probability calculation is based on empirical Bayes method with score distribution parameters from the pseudo training data. We make binary classification calls based on cutoff on probabilities with default 0.9. In addition, classification based on theoretical cutoff 0 is also provided.

#### References:

Jiang A, Hilton LK, Tang J, Rushton CK, Grande BM, Scott DW, Morin RD, Robust gene expression-based classification of cancers without normalization (<https://www.biorxiv.org/content/10.1101/2020.04.28.051953v1>)

Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P, Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, Kridel R, Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. Double-Hit Trait Expression Signature Defines a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell Lymphoma. Journal of Clinical Oncology 37, no. 3 (January 20, 2019) 190-201. DOI: 10.1200/JCO.18.01583

TR Golub, DK Slonim, P Tamayo, C Huard, M Gaasenbeek, JP Mesirov, H Coller, ML Loh, JR Downing, MA Caligiuri, et al. Molecular classification of cancer: class discovery and class prediction by gene expression monitoring Science, 286 (1999), pp. 531-537

### 3. LPS (Linear Predictor Score)

Classification scores are calculated with LPS formula (Wright et al., 2003), and probability calculation is based on empirical Bayes method assuming that we know LPS score distribution parameters for the two groups, or if we can estimate these parameters from a training data set. We make binary classification calls based on cutoff on probabilities with default 0.9.

#### References:

Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A trait expression-based method to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S A. 2003 Aug 19;100(17):9991-6.

### 4. PS (Prediction Strength).

Classification scores are calculated with PS formula (Golub et al., 1999), and probability calculation is based on empirical Bayes method assuming that we know PS score distribution parameters for the two groups, or if we can estimate these parameters from a training data set. We make binary classification calls based on cutoff on probabilities with default 0.9. In addition, classification based on theoretical cutoff 0 is also provided.

#### References:

TR Golub, DK Slonim, P Tamayo, C Huard, M Gaasenbeek, JP Mesirov, H Coller, ML Loh, JR Downing, MA Caligiuri, et al. Molecular classification of cancer: class discovery and class prediction by gene expression monitoring Science, 286 (1999), pp. 531-537

### 5. PS-ST (Prediction Strength).

Self-training algorithm (Jiang et al., 2020) is to identify cases in an unlabeled high dimensional data set that can be confidently classified and to use those cases as pseudo training data thereby allowing classification of all cases in the cohort. For PS-ST, classification scores are calculated with PS (Golub et al., 1999), and probability calculation is based on empirical Bayes method with score distribution parameters from the pseudo training data. We make binary classification calls based on cutoff on probabilities with default 0.9. In addition, classification based on theoretical cutoff 0 is also provided.

#### References:

Jiang A, Hilton LK, Tang J, Rushton CK, Grande BM, Scott DW, Morin RD, Robust gene expression-based classification of cancers without normalization (<https://www.biorxiv.org/content/10.1101/2020.04.28.051953v1>)

TR Golub, DK Slonim, P Tamayo, C Huard, M Gaasenbeek, JP Mesirov, H Coller, ML Loh, JR Downing, MA Caligiuri, et al. Molecular classification of cancer: class discovery and class prediction by gene expression monitoring Science, 286 (1999), pp. 531-537

   

II. Package installation
------------------------

PRPS R package is currently avaiable at GitHub (<https://github.com/ajiangsfu/PRPS>). It will be available at R (<https://www.r-project.org/>) in the future.

PRPS depends on several existing packages, users should have "limma" installed while other packages can be installed automatically along with PRPS installation. If you do not have "limma", please install it with the following example code before PRPS installation.

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") BiocManager::install("limma")

There are at least three ways to install PRPS R package from GitHub.

### 1. Install PRPS with R package "devtools""

install.packages("devtools") \#\#\# Run this line only if "devtools"" is not installed, otherwise please ignore this line devtools :: install\_github(repo = "ajiangsfu/PRPS")

### 2. Install PRPS with R package "remotes"

install.packages("remotes") \#\#\# Run this line only if "remotes"" is not installed, otherwise please ignore this line remotes :: install\_github(repo = "ajiangsfu/PRPS")

### 3. Install PRPS with archive files

There are two archive files for PRPS package available under <https://github.com/ajiangsfu/PRPS>. One file is in source package format for the most updated version (current version: PRPS\_1.0.0.tar.gz), and the other file is in binary package format for the most updated version (current version: PRPS\_1.0.0\_R\_x86\_64-pc-linux-gnu.tar.gz). The source package file can be used in any system, but the binary package file built under Linus cannot be installed in Window system.

After you download the PRPS package file(s), you can install PRPS with the following choices:

#### 1) Under R

Packages -&gt; Install Package(s) from local files, then you can select one of your downloaded PRPS package files to install.

#### 2) Under Rstudio

Tools -&gt; Install packages, select "Package Archive File(.tar.gz)" from dropdown list of the panel, then you can select one of your downloaded PRPS package files to install.

#### 3) Comment line

The comment line to install local package achive file is something like:

install.packages(prpsFile, repos=NULL, type="source")

Here, prpsFile is PRPS package name with your local path, such as "C:/Download/PRPS\_1.0.0.tar.gz".

III. Typical path: training + testing
-------------------------------------

Classification is usually processed when we have comparable training and testing data sets. Training data set is used to select variables to build up a classification model, and classfication is made for a comparable testing data set.

### 1. Typical workflow when training and testing data sets are comparable

When training and testing data sets are processed and comparable, the typical classification workflow with PRPS package is:

1). select one of classification algorithms: PRPS, LPS, or PS

2). make decisions on parameter setting for traning functions

3). run PRPStraining, or LPStraining, or PStraining on user's training data set

4). run PRPStesting, LPStesting, or PStesting on user's testing data set

### 2. Example data

There are five example data files attached with PRPS package.

#### 1). rosenwald.cli

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

#### 2). rosenwald.expr

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

#### 3). WrightCOO

This data frame contains 158 COO (cell of origin) related genes with both Ensembl annotation ID (row names) and gene symbols (column "Gene"). "The current standard approach for routinely classifying samples using Affymetrix U133 arrays employs 186 probesets (George Wright, personal communication). The 165 Ensembl genes that correspond to these probesets were used for classification by RNA-seq." (Morin et al., 2011)". In our recent research, we matched 158 genes in Ennishi cohort (Ennishi et al., 2019), which are used in self-training research (Jiang et al., 2020) and are included in PRPS package.

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

Jiang A, Hilton LK, Tang J, Rushton CK, Grande BM, Scott DW, Morin RD, Robust gene expression-based classification of cancers without normalization (<https://www.biorxiv.org/content/10.1101/2020.04.28.051953v1>)

#### 4). DHITsig and DHITsigENSG

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

### 3. Code Examples

In this section, we are providing code examples to get classification score when training and testing data sets are comparable. The example data set is rosenwald.expr, which is split into training and testing data sets based on set information in rosenwald.cli. Since each person might have his/her own preference, code examples with all three algorithms: PRPS, LPS and PS, are provided in the following.

#### 1) PRPS

If PRPS is selected, we first apply PRPStraining function to find top genes and their weights to distinguish the two groups given in the training data sets, distribution parameters of these genes and PRPS score are also included in the training object. Then, we call PRPStesting function to classify samples in the comparable testing data set with training ovject output.

##### PRPStraining

First of all, divide data into training and testing data sets.

``` r
dat = rosenwald.expr
trainset = subset(rosenwald.cli, rosenwald.cli$set == "Training")
testset = subset(rosenwald.cli, rosenwald.cli$set == "Validation")
```

Use GCB as a reference group, then PRPStraining object can be produced with the follwoing code assuming that we select top 100 genes for COO classification. There are five different ways to define gene weight, here, the default t test is used.

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

##### PRPStesting

Now, we can use our training object to make COO classification for a comparable testing data set. In the current example, this is simple to achieve since we already have testing data set shown above. We just use the output (training object) from PRPStraining together with testing data set to call PRPStesting to get result, and compare it to the original COO classification:

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

#### 2) LPS

If LPS is selected, we first apply LPStraining function to find top genes and their weights to distinguish the two groups given in the training data sets, LPS score distribution parameters are also included in the training object. Then, we call LPStesting function to classify samples in the comparable testing data set with training ovject output.

##### Training with LPStraining

Use the same training data set as PRPStraining, set GCB as the reference group, LPStraining object can be produced with the follwoing code assuming that we select top 100 genes for COO classification. There are five different ways to define gene weight, here, the default t test is used.

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

##### LPStesting

With the exact same testing data set as above, use the output from LPStraining, we can apply LPStesting as following and compare the classification with the truth in the testing data info:

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

#### 3) PS

##### PStraining

Use the same training data set as above, and use GCB as a reference group, then PStraining result can be got as following. Notice that now we are calling PStraining but all setting is the same as for LPStraining

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

##### PStesting

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

   

IIII. Special path: classification with self-training
-----------------------------------------------------

In many situations, we might not have comparable training and testing data sets, in this case, we can apply self-training methods for classification.

### 1. Self-training workflow when we do not have comparable training and testing data sets

When we do not have comparable training and testing data sets, the typical self-training classification workflow with PRPS package is:

1). select one of classification algorithms: PRPS-ST or PS-ST

2). make decisions on parameter setting for self-traning functions

3). run PRPSstableSLwithWeights or PRPSstableSLwithWeights on user's self-training data set without classification labels

4). run PRPSSLextension or PSSLextension on user's testing data set that is comparable to self-training data set

### 2. Code examples

In order to show self-training with above example, we will use above training data set as a self-training data set, but we do not use its classification information at all. However, in order to make self-training functions work, we assume that the top 100 genes with weights selected by the above training functions are the genes from previous research and can be applied here. The following example code is to get top 100 genes and their weights from above PRPStraining output.

``` r
wts = trainPRPS$PRPS_pars$weights
genes = rownames(wts)
wts = wts[,1]
names(wts) = genes
```

#### 1). PRPS-ST

If PRPS-ST method is selected, with the above top 100 genes, we can process self-learning, which can be extended to a comparable data set.

##### Self-training with PRPSstableSLwithWeights

``` r
require(mclust)  ### for some reason, this function is needed to be called here although I usually do not need to do so
#> Loading required package: mclust
#> Package 'mclust' version 5.4.6
#> Type 'citation("mclust")' for citing this R package in publications.
prpsst = PRPSstableSLwithWeights(newdat = dat[,rownames(trainset)],  weights = wts, PRPShighGroup = "ABC", PRPSlowGroup = "GCB")
```

![](README_files/figure-markdown_github/unnamed-chunk-28-1.png)![](README_files/figure-markdown_github/unnamed-chunk-28-2.png)![](README_files/figure-markdown_github/unnamed-chunk-28-3.png)![](README_files/figure-markdown_github/unnamed-chunk-28-4.png)![](README_files/figure-markdown_github/unnamed-chunk-28-5.png)![](README_files/figure-markdown_github/unnamed-chunk-28-6.png)![](README_files/figure-markdown_github/unnamed-chunk-28-7.png)![](README_files/figure-markdown_github/unnamed-chunk-28-8.png)![](README_files/figure-markdown_github/unnamed-chunk-28-9.png)![](README_files/figure-markdown_github/unnamed-chunk-28-10.png)![](README_files/figure-markdown_github/unnamed-chunk-28-11.png)![](README_files/figure-markdown_github/unnamed-chunk-28-12.png)![](README_files/figure-markdown_github/unnamed-chunk-28-13.png)![](README_files/figure-markdown_github/unnamed-chunk-28-14.png)![](README_files/figure-markdown_github/unnamed-chunk-28-15.png)![](README_files/figure-markdown_github/unnamed-chunk-28-16.png)![](README_files/figure-markdown_github/unnamed-chunk-28-17.png)![](README_files/figure-markdown_github/unnamed-chunk-28-18.png)![](README_files/figure-markdown_github/unnamed-chunk-28-19.png)![](README_files/figure-markdown_github/unnamed-chunk-28-20.png)![](README_files/figure-markdown_github/unnamed-chunk-28-21.png)![](README_files/figure-markdown_github/unnamed-chunk-28-22.png)![](README_files/figure-markdown_github/unnamed-chunk-28-23.png)![](README_files/figure-markdown_github/unnamed-chunk-28-24.png)![](README_files/figure-markdown_github/unnamed-chunk-28-25.png)![](README_files/figure-markdown_github/unnamed-chunk-28-26.png)![](README_files/figure-markdown_github/unnamed-chunk-28-27.png)![](README_files/figure-markdown_github/unnamed-chunk-28-28.png)![](README_files/figure-markdown_github/unnamed-chunk-28-29.png)![](README_files/figure-markdown_github/unnamed-chunk-28-30.png)![](README_files/figure-markdown_github/unnamed-chunk-28-31.png)![](README_files/figure-markdown_github/unnamed-chunk-28-32.png)![](README_files/figure-markdown_github/unnamed-chunk-28-33.png)![](README_files/figure-markdown_github/unnamed-chunk-28-34.png)![](README_files/figure-markdown_github/unnamed-chunk-28-35.png)![](README_files/figure-markdown_github/unnamed-chunk-28-36.png)![](README_files/figure-markdown_github/unnamed-chunk-28-37.png)![](README_files/figure-markdown_github/unnamed-chunk-28-38.png)![](README_files/figure-markdown_github/unnamed-chunk-28-39.png)

``` r
str(prpsst)
#> List of 2
#>  $ PRPS_pars:List of 4
#>   ..$ weights      :'data.frame':    100 obs. of  1 variable:
#>   .. ..$ weights: num [1:100] 7.94 7.43 7.14 6.86 6.79 ...
#>   ..$ meansds      : num [1:4] 314.7 -438.7 64.3 68.3
#>   ..$ traitsmeansds: num [1:100, 1:4] 0.538 0.479 0.479 0.809 0.389 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:100] "17708" "27631" "24796" "26697" ...
#>   .. .. ..$ : chr [1:4] "means2" "means1" "sds2" "sds1"
#>   ..$ dfs          : num [1:2] 13 9
#>  $ PRPS_test:'data.frame':   40 obs. of  6 variables:
#>   ..$ PRPS_score  : num [1:40] -421 163 336 -148 -343 ...
#>   ..$ PRPS_class  : chr [1:40] "GCB" "ABC" "ABC" "GCB" ...
#>   ..$ PRPS_prob1  : num [1:40] 4.46e-29 1.00 1.00 4.95e-08 5.41e-23 ...
#>   ..$ PRPS_prob2  : num [1:40] 1.00 2.15e-16 9.91e-29 1.00 1.00 ...
#>   ..$ PRPS_class0 : chr [1:40] "GCB" "ABC" "ABC" "GCB" ...
#>   ..$ stable_class: chr [1:40] "GCB" "ABC" "ABC" "UNCLASS" ...
table(trainset$group, prpsst$PRPS_test$PRPS_class)
#>      
#>       ABC GCB
#>   ABC  17   0
#>   GCB   1  22
```

This gives some idea how we are doing with PRPS self-training.

##### Extension with PRPSstableSLwithWeights

The above self-training object can also be used for another data set's classification as long as it is comparable to the self-training data set.

``` r
prpsext = PRPSSLextension(PRPSSLObj = prpsst, newdat = dat[,rownames(testset)])
str(prpsext)
#> 'data.frame':    20 obs. of  5 variables:
#>  $ PRPS_score    : num  -65.4 -224.1 -340.8 130.1 255.8 ...
#>  $ PRPS_class    : chr  "GCB" "GCB" "GCB" "ABC" ...
#>  $ PRPS_prob_test: num  7.93e-02 8.46e-14 8.15e-23 1.00 1.00 ...
#>  $ PRPS_prob_ref : num  9.21e-01 1.00 1.00 4.81e-14 4.73e-23 ...
#>  $ PRPS_class0   : chr  "GCB" "GCB" "GCB" "ABC" ...
table(testset$group, prpsext$PRPS_class)
#>      
#>       ABC GCB UNCLASS
#>   ABC   6   0       0
#>   GCB   1  11       2
```

This gives some idea how we are doing with extension of PRPS self-training.

#### 2). PS-ST

If PS-ST method is selected, with the above top 100 genes, we can process self-learning, which can be extended to a comparable data set.

##### Self-training with PSstableSLwithWeights

``` r
PSst = PSstableSLwithWeights(newdat = dat[,rownames(trainset)],  weights = wts, PShighGroup = "ABC", PSlowGroup = "GCB")
```

![](README_files/figure-markdown_github/unnamed-chunk-30-1.png)![](README_files/figure-markdown_github/unnamed-chunk-30-2.png)![](README_files/figure-markdown_github/unnamed-chunk-30-3.png)![](README_files/figure-markdown_github/unnamed-chunk-30-4.png)![](README_files/figure-markdown_github/unnamed-chunk-30-5.png)![](README_files/figure-markdown_github/unnamed-chunk-30-6.png)![](README_files/figure-markdown_github/unnamed-chunk-30-7.png)![](README_files/figure-markdown_github/unnamed-chunk-30-8.png)![](README_files/figure-markdown_github/unnamed-chunk-30-9.png)![](README_files/figure-markdown_github/unnamed-chunk-30-10.png)![](README_files/figure-markdown_github/unnamed-chunk-30-11.png)![](README_files/figure-markdown_github/unnamed-chunk-30-12.png)![](README_files/figure-markdown_github/unnamed-chunk-30-13.png)![](README_files/figure-markdown_github/unnamed-chunk-30-14.png)![](README_files/figure-markdown_github/unnamed-chunk-30-15.png)![](README_files/figure-markdown_github/unnamed-chunk-30-16.png)![](README_files/figure-markdown_github/unnamed-chunk-30-17.png)![](README_files/figure-markdown_github/unnamed-chunk-30-18.png)![](README_files/figure-markdown_github/unnamed-chunk-30-19.png)![](README_files/figure-markdown_github/unnamed-chunk-30-20.png)![](README_files/figure-markdown_github/unnamed-chunk-30-21.png)![](README_files/figure-markdown_github/unnamed-chunk-30-22.png)![](README_files/figure-markdown_github/unnamed-chunk-30-23.png)![](README_files/figure-markdown_github/unnamed-chunk-30-24.png)![](README_files/figure-markdown_github/unnamed-chunk-30-25.png)![](README_files/figure-markdown_github/unnamed-chunk-30-26.png)![](README_files/figure-markdown_github/unnamed-chunk-30-27.png)![](README_files/figure-markdown_github/unnamed-chunk-30-28.png)![](README_files/figure-markdown_github/unnamed-chunk-30-29.png)![](README_files/figure-markdown_github/unnamed-chunk-30-30.png)![](README_files/figure-markdown_github/unnamed-chunk-30-31.png)![](README_files/figure-markdown_github/unnamed-chunk-30-32.png)![](README_files/figure-markdown_github/unnamed-chunk-30-33.png)![](README_files/figure-markdown_github/unnamed-chunk-30-34.png)![](README_files/figure-markdown_github/unnamed-chunk-30-35.png)

``` r
str(PSst)
#> List of 2
#>  $ PS_pars:List of 3
#>   ..$ weights    :'data.frame':  100 obs. of  1 variable:
#>   .. ..$ weights: num [1:100] 7.94 7.43 7.14 6.86 6.79 ...
#>   ..$ meansds    : Named num [1:4] 0.834 -0.723 0.155 0.249
#>   .. ..- attr(*, "names")= chr [1:4] "highPSmean" "lowPSmean" "highPSsd" "lowPSsd"
#>   ..$ traitsmeans: num [1:100, 1:3] 0.0502 0.0676 0.0553 0.0395 0.0289 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:100] "17708" "27631" "24796" "26697" ...
#>   .. .. ..$ : chr [1:3] "meanOfGroupMeans" "groupMean1" "groupMean2"
#>  $ PS_test:'data.frame': 40 obs. of  6 variables:
#>   ..$ PS_score    : num [1:40] -0.756 0.55 0.882 -0.217 -0.692 ...
#>   ..$ PS_class    : chr [1:40] "GCB" "ABC" "ABC" "GCB" ...
#>   ..$ PS_prob1    : num [1:40] 2.51e-23 1.00 1.00 1.36e-09 1.55e-21 ...
#>   ..$ PS_prob2    : num [1:40] 1.00 7.26e-06 6.70e-10 1.00 1.00 ...
#>   ..$ PS_class0   : chr [1:40] "GCB" "ABC" "ABC" "GCB" ...
#>   ..$ stable_class: chr [1:40] "GCB" "UNCLASS" "ABC" "UNCLASS" ...
table(trainset$group, PSst$PS_test$PS_class)
#>      
#>       ABC GCB
#>   ABC  17   0
#>   GCB   0  23
```

This gives some idea how we are doing with PS self-training.

##### Extension with PSstableSLwithWeights

The above self-training object can also be used for another data set's classification as long as it is comparable to the self-training data set.

``` r
PSext = PSSLextension(PSSLObj = PSst, newdat = dat[,rownames(testset)])
str(PSext)
#> 'data.frame':    20 obs. of  5 variables:
#>  $ PS_score    : num  -0.175 -0.335 -0.628 0.655 0.815 ...
#>  $ PS_class    : chr  "GCB" "GCB" "GCB" "ABC" ...
#>  $ PS_prob_test: num  1.16e-08 2.52e-12 8.86e-20 1.00 1.00 ...
#>  $ PS_prob_ref : num  1.00 1.00 1.00 2.88e-07 3.45e-09 ...
#>  $ PS_class0   : chr  "GCB" "GCB" "GCB" "ABC" ...
table(testset$group, PSext$PS_class)
#>      
#>       ABC GCB UNCLASS
#>   ABC   6   0       0
#>   GCB   0  13       1
```

This gives some idea how we are doing with extension of PS self-training.

#### PRPS R package general information

Version: 1.0.0

Author: Aixiang Jiang

Maintainer: Aixiang Jiang <aijiang@bccrc.ca> <aixiang.jiang@pathology.ubc.ca>

Depends: R (&gt;= 3.5), lattice, caret, limma, e1071, pROC, mclust

Suggests: knitr

VignetteBuilder: knitr
