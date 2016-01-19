require(limma)
require(pander)
require(stats)
require(graphics)
require(reshape2)
require(scales)
require(ggplot2)

#' Statistical Analysis of the PathoScope reports and generates a html report 
#' and produces interactive shiny app plots
#' 
#' @param input_dir Directory where the tsv files from PathoScope are located
#' @param batch Batch covariate 
#' @param condition Covariates or conditions of interest besides batch
#' @param report_file Output report file name 
#' @param report_dir Output report directory path 
#' @param report_option_binary 9 bits Binary String representing the plots to display and hide in the report 
#' @param view_report when TRUE, opens the report in a browser 
#' @return outputfile The output file with all the statistical plots
#' @export
#' @examples
#' nbatch <- 10
#' nperbatch <- 10
#' batch <- rep(1:nbatch, each=nperbatch)
#' PathoStat(input_dir=".", batch)
pathoStat <- function(input_dir=".", batch, condition=NULL, 
                    report_file="PathoStat_report.html", 
                    report_dir=".", report_option_binary="111111111",
                    view_report=TRUE, interactive=TRUE)  {
  
  pdata <- data.frame(batch, condition)
  mod = model.matrix(~as.factor(condition), data=pdata)
  if (report_dir==".") { report_dir=getwd() }
  dat <- readPathoscopeData(input_dir)
  shinyInput <<- list("data"=dat, "batch"=batch, "condition"=condition, 
                      "report_dir"=report_dir, "input_dir"=input_dir)
  
  rmdfile <- system.file("reports/PathoStat_report.Rmd", package = "PathoStat")
  report_option_vector <- unlist(strsplit(as.character(report_option_binary), ""))
  static_lib_dir <- system.file("reports/libs", package = "PathoStat")
  file.copy(static_lib_dir, report_dir, recursive=TRUE)
  outputfile <- rmarkdown::render(rmdfile, output_file=report_file, output_dir=report_dir)
  shinyInputOrig <<- shinyInput
  shinyInputCombat <<- NULL
  if (view_report)  {
    browseURL(outputfile)
  }
  if (interactive)  {
    appDir <- system.file("shiny", "PathoStat", package = "PathoStat")
    if (appDir == "") {
      stop("Could not find shiny directory. Try re-installing PathoStat.", call. = FALSE)
    }
    shiny::runApp(appDir, display.mode = "normal")
  }
  return(outputfile)
}

