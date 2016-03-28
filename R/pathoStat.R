require(limma)
require(pander)
require(stats)
require(graphics)
require(reshape2)
require(ggplot2)
require(rentrez)

#' Statistical Analysis of the PathoScope reports and generates a html report 
#' and produces interactive shiny app plots
#' 
#' @param input_dir Directory where the tsv files from PathoScope are located
#' @param batch Batch covariate 
#' @param condition Covariates or conditions of interest besides batch
#' @param report_file Output report file name 
#' @param report_dir Output report directory path 
#' @param report_option_binary 9 bits Binary String representing the plots to 
#'  display and hide in the report 
#' @param view_report when TRUE, opens the report in a browser 
#' @return outputfile The output file with all the statistical plots
#' @import pander stats graphics reshape2 ggplot2 rentrez
#' @importFrom scales percent_format
#' @importFrom utils browseURL
#' @importFrom shiny runApp
#' @export
#' @examples
#' nbatch <- 10
#' nperbatch <- 10
#' batch <- rep(1:nbatch, each=nperbatch)
#' PathoStat(input_dir='.', batch)
pathoStat <- function(input_dir = ".", batch, condition = NULL, report_file = 
    "PathoStat_report.html", report_dir = ".", report_option_binary = 
    "111111111", view_report = TRUE, interactive = TRUE) {
    
    pdata <- data.frame(batch, condition)
    mod = model.matrix(~as.factor(condition), data = pdata)
    if (report_dir == ".") {
        report_dir = getwd()
    }
    datlist <- readPathoscopeData(input_dir)
    dat <- datlist$data
    countdat <- datlist$countdata
    ids <- rownames(dat)
    tids <- unlist(lapply(ids, FUN = grepTid))
    taxonLevels <- findTaxonomy(tids)
    shinyInput <- list(data = dat, batch = batch, condition = condition, 
        report_dir = report_dir, input_dir = input_dir, taxonLevels = 
        taxonLevels, countdata = countdat)
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
 
