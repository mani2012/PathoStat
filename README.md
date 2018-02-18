PathoStat: Statistical Microbiome Analysis Toolkit
==================================================

PathoStat is a R shiny package, designed for performing Statistical Microbiome Analysis on 
metagenomics results from sequencing data samples. In particular, it supports 
analyses on the [PathoScope](https://github.com/PathoScope/PathoScope) generated report files. 

The package includes:

    1. Data Summary and Filtering 
    2. Relative Abundance plots (Stacked Bar Plot, Heatmap)
    3. Multiple species boxplot visualization
    4. Diversity analysis (Alpha and Beta diversity)
    5. Differential Expression (DEseq2, edgeR)
    6. Dimension Reduction (PCA, PCoA)
    7. Biomarker identification
    8. Pathway analysis
    

To launch PathoStat in R, just enter the command:
```r
runPathoStat()
```



## Installation

Use [devtools](https://github.com/hadley/devtools) to install it as follows:
```r
require(devtools)
install_github("jasonzhao0307/PathoStat")
```


## Troubleshooting with Installation

If you are having issues with the installation, you may have to setup local 
directory, if you do not have permissions to install in the default location 
for R. You may also want to load a version of R 3.3.1 or higher.
```r
export R_LIBS="/my_own_local_directory/R_libs"
module load R/R-3.2.4
```

And do something like the following
```r
install.packages("devtools", repos="http://cran.r-project.org", 
    lib="/my_own_local_directory/R_libs")
```
