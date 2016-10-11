PathoStat: PathoScope Statistical Analysis package
==================================================

The purpose of this package is to perform Statistical Microbiome Analysis on 
metagenomics results from sequencing data samples. In particular, it supports 
analyses on the PathoScope generated report files. PathoStat provides various 
functionalities including Relative Abundance charts, Diversity estimates and 
plots, tests of Differential Abundance, Time Series visualization, and 
Core OTU analysis.

The package includes:

    1. Relative Abundance plots (Stacked Bar Plot, Heatmap)
    2. Diversity plots (Alpha and Beta diversity, Exploratory Tree, BiPlot, 
        Co-Occurrence)
    3. Differential Expression (Expression Plots, Limma)
    4. Confidence Region Plots
    5. PCA plots
    6. PCoA plots
    7. Alluvial Plots for longitudinal data
    8. Core OTU analysis

`runPathoStat` is the pipeline function that generates the PathoStat report
and launches shiny app when in interactive mode. It combines all the functions 
into one step.

## Installation

To begin, install [Bioconductor](http://www.bioconductor.org/) and simply
run the following to automatically install PathoStat and all the dependencies, 
except pandoc, which you have to manually install as follows.

```r
source("http://bioconductor.org/biocLite.R")
biocLite("PathoStat")
```
Install 'pandoc' package by following the instructions at the following URL:
http://pandoc.org/installing.html

Rstudio also provides pandoc binaries at the following location for Windows, 
Linux and Mac:
https://s3.amazonaws.com/rstudio-buildtools/pandoc-1.13.1.zip 

If you want to install the latest development version of PathoStat from Github, 
use [devtools](https://github.com/hadley/devtools) to install it as follows:
```r
require(devtools)
install_github("mani2012/PathoStat", build_vignettes=TRUE)
```

If you want to manually install the PathoStat dependencies, run the following:
```r
source("http://bioconductor.org/biocLite.R")
biocLite(c('MCMCpack', 'limma', 'corpcor', 'rmarkdown', 'knitr', 'pander',
'matrixStats', 'reshape2', 'scales', 'ggplot2', 'rentrez', 'BatchQC', 'DT', 
'gtools', 'plyr', 'tidyr', 'dplyr', 'ape', 'phyloseq', 'shiny', 'grDevices', 
'stats', 'methods', 'XML', 'graphics', 'utils', 'alluvial', 'BiocStyle'))
```

If all went well you should now be able to load PathoStat:
```r
require(PathoStat)
vignette('PathoStatIntro', package='PathoStat')
vignette("PathoStatUserManual", package='PathoStat')
vignette("PathoStatAdvanced", package='PathoStat')
runPathoStat()
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
