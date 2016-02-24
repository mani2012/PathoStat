PathoStat: PathoScope Statistical Analysis package
==================================================

The purpose of this package is to perform Statistical Analysis on the PathoScope generated reports file.

The package includes:

1. Relative Abundance plots

`pathoStat` is the pipeline function that generates the PathoStat report. It combines all the functions into one step.

## Installation

To begin, install [Bioconductor](http://www.bioconductor.org/) along with a
few dependencies that PathoStat uses:

```r
source("http://bioconductor.org/biocLite.R")
biocLite(c('MCMCpack', 'limma', 'preprocessCore', 'devtools', 'stringi', 'corpcor', 'matrixStats', 'shiny', 'ggvis', 'd3heatmap', 'reshape2', 'scales', 'rentrez'))
```

Install 'pandoc' package by following the instructions at the following URL:
http://johnmacfarlane.net/pandoc/installing.html

Rstudio provide pandoc binaries at the following location for (windows, linux and mac):
https://s3.amazonaws.com/rstudio-buildtools/pandoc-1.13.1.zip 

Install 'pander' package by following the instructions at the following URL:
http://rapporter.github.io/pander/

```r
require(devtools)
install_github('Rapporter/pander')
```

Next, use [devtools](https://github.com/hadley/devtools) to install the latest
version of PathoStat from Github:
```r
require(devtools)
install_github("mani2012/PathoStat", build_vignettes=TRUE, auth_token="dadf36cdaef71a2f761f193862a8f6f3f36e3966")
```

If all went well you should now be able to load PathoStat:
```r
require(PathoStat)
vignette('PathoStatIntro', package='PathoStat')
```

## Troubleshooting with Installation

If you are having issues with the installation, you may have to setup local directory, if you do not have permissions to install in the default location for R. You may also want to load a version of R 3.1 or higher.
```r
export R_LIBS="/my_own_local_directory/R_libs"
module load R/R-3.1.1
```

And do something like the following
```r
install.packages("devtools", repos="http://cran.r-project.org", lib="/my_own_local_directory/R_libs")
```