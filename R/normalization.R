#' Normalize the given data based on library size
#'
#' @param zcounts Input counts data matrix
#' @return acounts Normalized counts data matrix
#' @export
#' @examples
#' example_data_dir <- system.file("example/data", package = "PathoStat")
#' pathoreport_file_suffix <- "-sam-report.tsv"
#' datlist <- readPathoscopeData(example_data_dir, pathoreport_file_suffix)
#' countdat <- datlist$countdata
#' acounts <- sizeNormalize(countdat)
sizeNormalize <- function(zcounts) {
    minimum <- min(zcounts)
    if (minimum < 0) {
        zcounts <- zcounts - minimum
    }
    avg <- mean(colSums(zcounts))
    acounts <- apply(zcounts, 2, FUN = function(x)  {
        size <- sum(x)
        x <- x*avg/size
    })
    return(acounts)
}

#' Compute Core OTUs for the given data matrix
#'
#' @param zcounts Standardized counts
#' @param otuthreshold Abundance cutoff threshold for the OTU to be picked
#' @param prevalence Prevalence of the OTU at threshold cutoff among samples
#' @return list containing core OTUs
#' @export
#' @examples
#' example_data_dir <- system.file("example/data", package = "PathoStat")
#' pathoreport_file_suffix <- "-sam-report.tsv"
#' datlist <- readPathoscopeData(example_data_dir, pathoreport_file_suffix)
#' countdat <- datlist$countdata
#' coreotus <- coreOTU(countdat)
coreOTU <- function(zcounts, otuthreshold=0.05, prevalence=0.4) {
    pcounts <- apply(zcounts, 2, FUN = function(x)  {
        size <- sum(x)
        x <- x/size
    })
    coreotus <- c()
    for (i in seq_len(dim(zcounts)[1]))  {
        num <- length(which(pcounts[i,]>otuthreshold))
        prev <- num/dim(zcounts)[2]
        if (prev > prevalence)  coreotus <- c(coreotus, i)
    }
    return(coreotus)
}

#' Compute Empirical Bayes OTU Normalized data
#'
#' @param zcounts counts data to be normalized
#' @param wt Weight parameter indicating how much information to borrow
#' across samples using Empirical Bayes 
#' @param threshold Abundance cutoff threshold for the OTU to be picked
#' @param prevalence Prevalence of the OTU at threshold cutoff among samples
#' @return list containing Empirical Bayes coreOTU Normalized data
#' @export
#' @examples
#' example_data_dir <- system.file("example/data", package = "PathoStat")
#' pathoreport_file_suffix <- "-sam-report.tsv"
#' datlist <- readPathoscopeData(example_data_dir, pathoreport_file_suffix)
#' countdat <- datlist$countdata
#' coreotunormdat <- coreOTUNormalize(countdat)
coreOTUNormalize <- function(zcounts, wt=0.25, otuthreshold=0.05, prevalence=0.4) {
    zcounts <- sizeNormalize(zcounts)
    coreotus <- coreOTU(zcounts, otuthreshold, prevalence)
    if (is.null(coreotus) || length(coreotus)<=1)  {
        acounts <- zcounts
    } else  {
        coreotucounts <- zcounts[coreotus,]
        avg <- mean(colSums(coreotucounts))
        acounts <- apply(zcounts, 2, FUN = function(x, coreotus, avg)  {
            size <- sum(x[coreotus])
            #w <- wt*dim(zcounts)[2]
            if (wt < 0) wt <- 0
            if (wt >= 1) wt <- 1-0.0001
            w <- wt/(1-wt)
            asize <- (size + w*avg)/(w+1)
            x <- x*avg/asize
        }, coreotus, avg)
    }
    return(acounts)
}

#' Compute coreOTU Quantile Normalized data
#'
#' @param zcounts counts data to be normalized
#' @param threshold Abundance cutoff threshold for the OTU to be picked
#' @param prevalence Prevalence of the OTU at threshold cutoff among samples
#' @return list containing coreOTU Quantile Normalized data
#' @importFrom preprocessCore normalize.quantiles
#' @export
#' @examples
#' example_data_dir <- system.file("example/data", package = "PathoStat")
#' pathoreport_file_suffix <- "-sam-report.tsv"
#' datlist <- readPathoscopeData(example_data_dir, pathoreport_file_suffix)
#' countdat <- datlist$countdata
#' coreotunormdat <- coreOTUQuantile(countdat)
coreOTUQuantile <- function(zcounts, otuthreshold=0.05, prevalence=0.4) {
    zcounts <- sizeNormalize(zcounts)
    coreotus <- coreOTU(zcounts, otuthreshold, prevalence)
    acounts <- zcounts
    if (!is.null(coreotus) && length(coreotus)>1)  {
        coreotucounts <- zcounts[coreotus,]
        acounts <- zcounts
        acounts[coreotus,] <- normalize.quantiles(coreotucounts)
    }
    return(acounts)
}

## Example code
# thedata = textConnection(
# "7 8 7 6 9 7 8 7 6 9\
#  8 9 7 8 7 8 9 7 8 7\
#  0 0 0 0 0 10 11 0 0 0\
#  0 0 0 0 0 0 0 9 10 0\
#  0 0 0 0 0 0 0 0 0 11\
#  ")
# data.matrix <- matrix(scan(thedata),ncol=10,byrow=TRUE)
# library(preprocessCore)
# norm.data.matrix <- normalize.quantiles(data.matrix)
# otunorm.data <- otuNormalize(data.matrix, 0)

# library(BatchQC)
# ### simulate data
# nbatch <- 3
# ncond <- 2
# npercond <- 10
# data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
#     npercond, basemean=10000, ggstep=50, bbstep=2000, ccstep=800, 
#     basedisp=100, bdispstep=-10, swvar=1000, seed=1234)
# batch <- rep(1:nbatch, each=ncond*npercond)
# condition <- rep(rep(1:ncond, each=npercond), nbatch)
