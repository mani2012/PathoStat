require(limma)
require(pander)
require(stats)
require(graphics)
require(reshape2)
require(ggplot2)
require(rentrez)
require(ape)
require(phyloseq)

###############################################################################
#' Generates a PathoStat object from the PathoScope reports for further 
#' analysis using the interactive shiny app
#' 
#' @param input_dir Directory where the tsv files from PathoScope are located
#' @param batch Batch covariate 
#' @param condition Covariates or conditions of interest besides batch
#' @return outputfile The output file with all the statistical plots
#' @import pander stats graphics reshape2 ggplot2 rentrez phyloseq
#' @importFrom scales percent_format
#' @importFrom ape rtree
#' @export
#' @examples
#' example_data_dir <- system.file("example/data", package = "PathoStat")
#' pstat <- createPathoStat(input_dir=example_data_dir, 
#'     sample_data_file="sample_data.tsv")
createPathoStat <- function(input_dir=".", sample_data_file="sample_data.tsv",
    pathoreport_file_suffix="-sam-report.tsv") {
    datlist <- readPathoscopeData(input_dir, pathoreport_file_suffix)
    dat <- datlist$data
    countdat <- datlist$countdata
    ids <- rownames(dat)
    tids <- unlist(lapply(ids, FUN = grepTid))
    taxonLevels <- findTaxonomy(tids)
    taxmat <- findTaxonMat(ids, taxonLevels)
    OTU <- otu_table(countdat, taxa_are_rows = TRUE)
    TAX <- tax_table(taxmat)
    physeq <- phyloseq(OTU, TAX)
    sample_file_path <- file.path(input_dir, sample_data_file)
    tbl <- read.table(sample_file_path, header=TRUE, sep="\t", row.names=1,
        stringsAsFactors=FALSE)
    sampledata = sample_data(data.frame(tbl))
    random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=
        taxa_names(physeq))
    physeq1 <- merge_phyloseq(physeq, sampledata, random_tree)
    pstat <- pathostat(physeq1)
    return(pstat)
}

###############################################################################
#' Statistical Analysis of the PathoScope reports and generates a html report 
#' and produces interactive shiny app plots
#' 
#' @param pstat phyloseq extension pathostat object
#' @param report_file Output report file name 
#' @param report_dir Output report directory path 
#' @param report_option_binary 9 bits Binary String representing the plots to 
#'  display and hide in the report 
#' @param view_report when TRUE, opens the report in a browser 
#' @return outputfile The output file with all the statistical plots
#' @import pander stats graphics reshape2 ggplot2 rentrez phyloseq
#' @importFrom scales percent_format
#' @importFrom utils browseURL
#' @importFrom shiny runApp
#' @export
#' @examples
#' runPathoStat()
runPathoStat <- function(pstat=NULL, report_file = "PathoStat_report.html", 
    report_dir = ".", report_option_binary = "111111111", view_report = FALSE, 
    interactive = TRUE) {
    
    if (is.null(pstat))  {
        data_dir <- system.file("data", package = "PathoStat")
        infileName <- "pstat_data.rda"
        pstat <- loadPstat(data_dir, infileName)
    }
    if (report_dir == ".") {
        report_dir = getwd()
    }
    shinyInput <- list(pstat = pstat, report_dir = report_dir)
    setShinyInput(shinyInput)
    rmdfile <- system.file("reports/PathoStat_report.Rmd", package = 
        "PathoStat")
    report_option_vector <- unlist(strsplit(as.character(report_option_binary), 
        ""))
    static_lib_dir <- system.file("reports/libs", package = "PathoStat")
    file.copy(static_lib_dir, report_dir, recursive = TRUE)
    outputfile <- rmarkdown::render(rmdfile, output_file = report_file, 
        output_dir = report_dir)
    shinyInput <- getShinyInput()
    setShinyInputOrig(shinyInput)
    setShinyInputCombat(NULL)
    if (view_report) {
        browseURL(outputfile)
    }
    if (interactive) {
        appDir <- system.file("shiny", "PathoStat", package = "PathoStat")
        if (appDir == "") {
            stop("Could not find shiny directory. Try re-installing PathoStat.",
                call. = FALSE)
        }
        shiny::runApp(appDir, display.mode = "normal")
    }
    return(outputfile)
}

#' Getter function to get the shinyInput option
#' @return shinyInput option
#' @export
#' @examples
#' getShinyInput()
getShinyInput <- function() {
    shinyInput <- getOption("pathostat.shinyInput")
    return(shinyInput)
}
#' Setter function to set the shinyInput option
#' @param x shinyInput option
#' @return shinyInput option
#' @export
#' @examples
#' setShinyInput(NULL)
setShinyInput <- function(x) {
    options(pathostat.shinyInput = x)
}

#' Getter function to get the shinyInputOrig option
#' @return shinyInputOrig option
#' @export
#' @examples
#' getShinyInputOrig()
getShinyInputOrig <- function() {
    shinyInputOrig <- getOption("pathostat.shinyInputOrig")
    return(shinyInputOrig)
}
#' Setter function to set the shinyInputOrig option
#' @param x shinyInputOrig option
#' @return shinyInputOrig option
#' @export
#' @examples
#' setShinyInputOrig(NULL)
setShinyInputOrig <- function(x) {
    options(pathostat.shinyInputOrig = x)
}

#' Getter function to get the shinyInputCombat option
#' @return shinyInputCombat option
#' @export
#' @examples
#' getShinyInputCombat()
getShinyInputCombat <- function() {
    shinyInputCombat <- getOption("pathostat.shinyInputCombat")
    return(shinyInputCombat)
}
#' Setter function to set the shinyInputCombat option
#' @param x shinyInputCombat option
#' @return shinyInputCombat option
#' @export
#' @examples
#' setShinyInputCombat(NULL)
setShinyInputCombat <- function(x) {
    options(pathostat.shinyInputCombat = x)
}
 
