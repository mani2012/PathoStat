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
    

`runPathoStat` is the pipeline function that generates the PathoStat report
and launches shiny app when in interactive mode. It combines all the functions 
into one step.

## Run Pathostat
To launch PathoStat in R, just enter the command:
```r
runPathoStat()
```

## Installation

To begin, install [Bioconductor](http://www.bioconductor.org/) and simply
run the following to automatically install PathoStat and all the dependencies 
as follows.

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("PathoStat")
```

If you want to install the latest development version of PathoStat from Github, 
use [devtools](https://github.com/hadley/devtools) to install it as follows:
```r
require(devtools)
install_github("compbiomed/PathoStat")
```

If all went well you should now be able to load PathoStat:
```r
require(PathoStat)
vignette('PathoStat-vignette', package='PathoStat')
runPathoStat()
```
