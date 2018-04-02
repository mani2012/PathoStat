#' pathostat object generated from example pathoscope report files 
#'
#' This example data consists of 33 samples from a diet study with 11 subjects 
#' taking 3 different diets in random order
#'
#' @name pstat_data
#' @format pathostat object extension of phyloseq-class experiment-level object:
#' \describe{
#'     \item{otu_table}{OTU table with  41 taxa and 33 samples}
#'     \item{sample_data}{Sample Data with 33 samples by 18 sample variables}
#'     \item{tax_table}{Taxonomy Table with 41 taxa by 9 taxonomic ranks}
#'     \item{sample_data}{Phylogenetic Tree with 41 tips and 40 internal nodes}
#' }
#' @return pathostat object
"pstat"

###############################################################################
#' Statistical Microbiome Analysis on the pathostat input and generates a
#' html report and produces interactive shiny app plots
#'
#' @param pstat phyloseq extension pathostat object
#' @param report_dir Output report directory path
#' @param report_option_binary 9 bits Binary String representing the plots to
#'  display and hide in the report
#' @param interactive when TRUE, opens the interactive shinyApp
#' @return outputfile The output file with all the statistical plots
#' @import stats graphics reshape2 ggplot2 rentrez phyloseq
#' @import corpcor knitr limma matrixStats methods BiocStyle
#' @import edgeR
#' @importFrom DESeq2 DESeqDataSetFromMatrix sizeFactors estimateDispersions
#'  nbinomWaldTest results
#' @importFrom plyr ddply numcolwise
#' @importFrom dplyr add_rownames mutate_each select_ add_rownames
#' @importFrom scales percent_format
#' @importFrom utils browseURL
#' @importFrom shiny runApp
#' @importFrom shiny shinyOptions getShinyOption
#' @export
#' @examples
#' runPathoStat(interactive = FALSE)
runPathoStat <- function(pstat=NULL,
    report_dir=".", report_option_binary="111111111",interactive=TRUE) {

    if (is.null(pstat))  {
        data_dir <- system.file("data", package = "PathoStat")
        infileName <- "pstat_data.rda"
        pstat <- loadPstat(data_dir, infileName)
    }
    if (report_dir == ".") {
        report_dir = getwd()
    }

    #test and fix the constant/zero row
    row.remove.index <- c()
    if (sum(rowSums(as.matrix(pstat@otu_table@.Data)) == 0) > 0){
        row.remove.index <-
        which(rowSums(as.matrix(pstat@otu_table@.Data)) == 0)
        pstat@otu_table@.Data <- pstat@otu_table@.Data[-row.remove.index,]
        pstat@tax_table@.Data <- pstat@tax_table@.Data[-row.remove.index,]
    }


    # remove variables with identical values
    col.remove.index <- c()
    for (i in 1:ncol(pstat@sam_data)){
        if(dim(unique(pstat@sam_data[,i]))[1] < 2){
            col.remove.index <- c(col.remove.index, i)
        }
    }
    if (!is.null(col.remove.index)){
        pstat@sam_data <- pstat@sam_data[,-col.remove.index]
    }

    shinyInput <- list(pstat = pstat, report_dir = report_dir)
    setShinyInput(shinyInput)
    rmdfile <- system.file("reports/pathostat_report.Rmd", package =
        "PathoStat")
    report_option_vector <- unlist(strsplit(as.character(report_option_binary),
        ""))
    static_lib_dir <- system.file("reports/libs", package = "PathoStat")
    file.copy(static_lib_dir, report_dir, recursive = TRUE)
    shinyInput <- getShinyInput()
    setShinyInputOrig(shinyInput)
    setShinyInputCombat(NULL)

    if (interactive) {
        appDir <- system.file("shiny", "PathoStat", package = "PathoStat")
        if (appDir == "") {
            stop("Could not find shiny directory. Try re-installing PathoStat.",
                call. = FALSE)
        }
        shiny::runApp(appDir, display.mode = "normal")
    }
}

#' Getter function to get the shinyInput option
#' @return shinyInput option
#' @export
#' @examples
#' getShinyInput()
getShinyInput <- function() {
    shinyInput <- getShinyOption("pathostat.shinyInput")
    return(shinyInput)
}
#' Setter function to set the shinyInput option
#' @param x shinyInput option
#' @return shinyInput option
#' @export
#' @examples
#' setShinyInput(NULL)
setShinyInput <- function(x) {
    shinyOptions(pathostat.shinyInput = x)
}

#' Getter function to get the shinyInputOrig option
#' @return shinyInputOrig option
#' @export
#' @examples
#' getShinyInputOrig()
getShinyInputOrig <- function() {
    shinyInputOrig <- getShinyOption("pathostat.shinyInputOrig")
    return(shinyInputOrig)
}
#' Setter function to set the shinyInputOrig option
#' @param x shinyInputOrig option
#' @return shinyInputOrig option
#' @export
#' @examples
#' setShinyInputOrig(NULL)
setShinyInputOrig <- function(x) {
    shinyOptions(pathostat.shinyInputOrig = x)
}

#' Getter function to get the shinyInputCombat option
#' @return shinyInputCombat option
#' @export
#' @examples
#' getShinyInputCombat()
getShinyInputCombat <- function() {
    shinyInputCombat <- getShinyOption("pathostat.shinyInputCombat")
    return(shinyInputCombat)
}
#' Setter function to set the shinyInputCombat option
#' @param x shinyInputCombat option
#' @return shinyInputCombat option
#' @export
#' @examples
#' setShinyInputCombat(NULL)
setShinyInputCombat <- function(x) {
    shinyOptions(pathostat.shinyInputCombat = x)
}
