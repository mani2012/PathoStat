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
#' @param otuthreshold Abundance cutoff threshold for the OTU to be picked
#' @param prevalence Prevalence of the OTU at threshold cutoff among samples
#' @return list containing Empirical Bayes coreOTU Normalized data
#' @export
#' @examples
#' example_data_dir <- system.file("example/data", package = "PathoStat")
#' pathoreport_file_suffix <- "-sam-report.tsv"
#' datlist <- readPathoscopeData(example_data_dir, pathoreport_file_suffix)
#' countdat <- datlist$countdata
#' coreotunormdat <- coreOTUNormalize(countdat)
coreOTUNormalize <- function(zcounts, wt=0.25, otuthreshold=0.05, 
    prevalence=0.4) {
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
            if (asize <= 0) asize <- 1
            x <- x*avg/asize
        }, coreotus, avg)
    }
    return(acounts)
}

#' Compute coreOTU Quantile Normalized data
#'
#' @param zcounts counts data to be normalized
#' @param otuthreshold Abundance cutoff threshold for the OTU to be picked
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
