################################################################################
#' PathoStat class to store PathoStat input data including phyloseq object
#'
#' Contains all currently-supported BatchQC output data classes: 
#' 
#' slots:
#' \describe{
#'     \item{average_count}{a single object of class otu_tableOrNULL}
#'     \item{besthit_count}{a single object of class otu_tableOrNULL}
#'     \item{highconf_count}{a single object of class otu_tableOrNULL}
#'     \item{lowconf_count}{a single object of class otu_tableOrNULL}
#' }
#' 
#' @name PathoStat-class
#' @rdname PathoStat-class
#' @exportClass PathoStat
pathostat1 <- setClass(Class="PathoStat", 
    representation=representation(
        average_count="otu_tableOrNULL",
        besthit_count="otu_tableOrNULL",
        highconf_count="otu_tableOrNULL",
        lowconf_count="otu_tableOrNULL"
    ), contains = "phyloseq",
    prototype=prototype(average_count=NULL, besthit_count=NULL,
        highconf_count=NULL, lowconf_count=NULL)
)
###############################################################################
