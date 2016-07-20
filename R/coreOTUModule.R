
# if(!require(microbiome, quietly = TRUE)) {
#   source("http://www.bioconductor.org/biocLite.R")
#   biocLite("ade4")
#   biocLite("fastcluster")
#   biocLite("devtools")
#   biocLite("ggplot2")
#   biocLite("Matrix")
#   biocLite("minet")
#   biocLite("mixOmics")
#   biocLite("dplyr")
#   biocLite("qvalue")
#   biocLite("reshape2")
#   biocLite("vegan")
#   biocLite("phyloseq")
#   biocLite("rpart")
#   biocLite("GO.db")
#   devtools::install_github("microbiome/microbiome")
# }

# require(microbiome)

makePhyseqFromShinyInput <- function(shinyInput) {
    ids <- rownames(shinyInput$data)
    taxmat <- findTaxonMat(ids, shinyInput$taxonLevels)
    OTU <- otu_table(shinyInput$countdata, taxa_are_rows = TRUE)
    TAX <- tax_table(taxmat)
    physeq <- phyloseq(OTU, TAX)
    sampledata = sample_data(data.frame(condition=as.factor(
      shinyInput$condition), batch=as.factor(shinyInput$batch), 
      row.names=sample_names(physeq), stringsAsFactors=FALSE))
    random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=
                          taxa_names(physeq))
    physeq1 <- merge_phyloseq(physeq, sampledata, random_tree)
    return(physeq1)
}

#' Creates Core Microbiome plot
#' A (very) thin wrapper around microbiome::plot_core
#'
#' @param psobj Phyloseq object
#' @param prevalence.intervals Prevalence intervals
#' @param detection.thresholds Detection thresholds
#' @param plot.type Plot type
#' 
#' @return ggplot object with core microbiome plot
#' 
plotCoreMicrobiota <- function(psobj,
                               prevalence.intervals = seq(10, 100, 10),
                               detection.thresholds = seq(0,10,0.5),
                               plot.type="lineplot") {
    res <- microbiome::plot_core(psobj, 
                                 prevalence.intervals=prevalence.intervals,
                                 detection.thresholds=detection.thresholds,
                                 plot.type=plot.type)
    res$plot + ggplot2::xlab("Abundance (OTU read count)")
}

#' UI function for Core OTU Module
#'
#' @param id Namespace for module
#' @param label Tab label
#' 
#' @return None
#' 
#' @export
coreOTUModuleUI <-
  function(id, label = "Core OTUs") {
    # This is the namespace for the module
    ns <- NS(id)
    
    fluidRow(
      column(6, plotOutput(ns("coreLine"))),
      column(6, plotOutput(ns("coreHeat")))
    )
  }


#' Server function for Core OTU Module
#'
#' @param input Shiny server input object
#' @param output Shiny server output object
#' @param session This is the session
#' 
#' @return None
#' 
#' @export
coreOTUModule <- function(input, output, session) {
    phyobj <- makePhyseqFromShinyInput(getShinyInput())
    
    output$coreLine <- renderPlot({
      plotCoreMicrobiota(phyobj)
    })
      
    output$coreHeat <- renderPlot({
      plotCoreMicrobiota(phyobj, plot.type="heatmap")
    })

}
