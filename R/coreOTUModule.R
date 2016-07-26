
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
                               prevalence.min,
                               
                               prevalence.intervals = seq(5, 100, 5),
                               detection.thresholds = NULL,
                               plot.type="lineplot") {
    prev <- prevalence.intervals
    if(is.null(detection.thresholds))
        det <- 10^seq(log10(1e-3), log10( max(taxa_abundances(psobj)) ), length = 20)

    res <- microbiome::plot_core(psobj, 
                                 prevalence.intervals=prev,
                                 detection.thresholds=det,
                                 plot.type=plot.type)
    res$plot + ggplot2::xlab("Abundance (OTU read count)")
}

get_coremat <- function(phyobj) {
    det <- unlist(lapply(0:7, function(p){c(1,2,5) * 10^p}))
    det <- det[1:min(which(det > max(otu_table(phyobj))))]
    if(nsamples(phyobj)<101) {
      prev <- seq(1,nsamples(phyobj)) / nsamples(phyobj)
    } else { prev <- seq(0.01, 1, by=0.01) 
    v <- seq(0,1,length.out=45)
    v2 <- seq(1,nsamples(phyobj)) / nsamples(phyobj)
    v <- v[2:45] - v2
    coremat <- data.frame(do.call(rbind, lapply(prev, function(p){
        sapply(det, function(d) sum(rowSums(otu_table(phyobj) > d) >= (p * nsamples(phyobj))))
    })))
    colnames(coremat) <- det
    coremat$prev <- prev
    coremat
}

coremat_lineplot <- function(coremat) {
  coremat %>% 
    dplyr::mutate(prev=factor(prev)) %>%
    tidyr::gather(det, count, -prev) %>%
    dplyr::mutate(det=as.numeric(det)) %>%
    ggplot(aes(x=det, y=count, color=prev, group=prev)) + 
    geom_point() + 
    geom_line() +
    scale_x_log10()
}
coremat_heatmap <- function(phyobj) {
  det <- 10^seq(0,log10(max(otu_table(phyobj), na.rm = T)), length = 20)
  # det <- unlist(lapply(0:7, function(p){c(1,2,5) * 10^p}))
  # det <- det[1:min(which(det > max(otu_table(phyobj))))]
  coremat2 <- data.frame(do.call(cbind, lapply(det, function(d){
    rowSums(otu_table(phyobj) > d) / nsamples(phyobj)
  })))
  colnames(coremat2) <- det
  taxorder <- rownames(coremat2[order(-rowSums(coremat2)),])
  coremat2 <- coremat2[taxorder,]
  coremat2$Taxa <- factor(rownames(coremat2),levels=taxorder)

  coremat2 %>%
    tidyr::gather(det, prev, -Taxa) %>%
    dplyr::mutate(det=as.numeric(det)) %>%
    ggplot(aes(x=det, y=Taxa, fill=prev)) + geom_tile() + scale_x_log10() +
    scale_fill_gradientn("Prevalence", 
                         breaks = seq(from = 0, to = 1, by = 0.2),
                         colours = gray(seq(0,1,length=5)))
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
    tabPanel("Core OTUs",
             sidebarLayout(
               sidebarPanel(
                   sliderInput(ns("detThresh"), "detection", value=20, min=1, max=1000),
                   sliderInput(ns("prevThresh"), "prevalence", value=67, min=0, max=100),
                   width = 3
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("2D lineplot",
                            plotOutput(ns("coreLine"))
                   ),
                   tabPanel("Heatmap", 
                            plotOutput(ns("coreHeat"))
                   ),
                   tabPanel("Table", 
                            DT::dataTableOutput(ns("coreTable"), width='95%')
                   )
                 ),
                 width = 9
               )
             )
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
coreOTUModule <- function(input, output, session, phyobj) {
    output$coreTable <- DT::renderDataTable({
        coretax <- microbiome::core(phyobj,
                                    detection.threshold = input$detThresh,
                                    prevalence.threshold = input$prevThresh)
        tax_table(phyobj)[coretax,]
    })
    
    coremat <- get_coremat(phyobj)
    cols.active <- colorRampPalette(c("black","blue"))(length(coremat$prev))
    
    output$coreLine <- renderPlot({
        # cols <- ifelse(coremat$prev < (input$prevThresh/100), cols.active, cols.inactive)
        cols <- ifelse(coremat$prev == (input$prevThresh/100), cols.active, "#00000011")
        coremat_lineplot(coremat) + 
          scale_colour_manual(guide=F, values=cols) +
          geom_vline(xintercept = input$detThresh, color="#990000", linetype=2)
    })
      
    output$coreHeat <- renderPlot({
      coremat_heatmap(phyobj)
      # plotCoreMicrobiota(phyobj, plot.type="heatmap")
    })

}
