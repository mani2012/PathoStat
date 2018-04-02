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
    
## Run Pathostat
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


