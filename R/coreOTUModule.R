#' Select rows of OTU matrix that meet given detection and prevalence thresholds
#' 
#' @param pstat PathoStat object
#' @param detection An integer specifying the minimum value considered to be 
#' "detected"
#' @param prevalence An integer specifying the minimum number of samples that
#' must be detected
#'
#' @return Subsetted PathoStat object
get_core <- function(pstat, detection, prevalence) {
    filter_taxa(pstat, function(x){ sum(x >= detection, na.rm=TRUE) >= 
        prevalence}, prune = TRUE)
}

#' Create core OTU matrix containing number of OTUs detected at varying 
#' detection and prevalence thresholds.
#' 
#' @param pstat PathoStat object
#'
#' @return Data frame containing number of OTUs at varying detection and 
#' prevalence thresholds, with rows corresponding to number of samples and
#' columns corresponding to detection thresholds. An additional column called
#' "prev"contains the sample threshold for each row.
get_coremat <- function(pstat) {
    # Detection values to calculate
    det <- unlist(lapply(0:7, function(p){c(1,2,5) * 10^p}))
    det <- det[1:min(which(det > max(otu_table(pstat))))]
    # Prevalence values to calculate
    prev <- seq(1,nsamples(pstat))
    
    coremat <- data.frame(do.call(rbind, lapply(prev, function(p){
        sapply(det, function(d) sum(rowSums(otu_table(pstat) >= d) >= p))
    })))
    colnames(coremat) <- det
    coremat$prev <- prev
    coremat  
}

#' Create line plot from core OTU matrix
#' 
#' @param coremat Core OTU matrix (data.frame)
#' 
#' @return Line plot with number of OTUs on the x-axis and detection threshold
#' on the y-axis. Lines connect data points with the same number of samples.
#' 
#' @import ggplot2
#' @importFrom dplyr mutate
#' @importFrom tidyr gather %>%
get_coremat_lineplot <- function(coremat) {
    coremat %>% 
    dplyr::mutate(prev=factor(prev)) %>%
    tidyr::gather(det, count, -prev) %>%
    dplyr::mutate(det=as.numeric(det)) %>%
    ggplot(aes(x=det, y=count, group=prev, color=prev)) + 
    geom_point() + 
    geom_line() +
    scale_x_log10() +
    labs(x="Detection threshold", y="OTU count")
}

#' @import grDevices
#' @importFrom tidyr %>%
get_coremat_heatmap <- function(pstat) {
    det <- 10^seq(0,log10(max(otu_table(pstat), na.rm = TRUE)), length = 20)
    # det <- unlist(lapply(0:7, function(p){c(1,2,5) * 10^p}))
    # det <- det[1:min(which(det > max(otu_table(pstat))))]
    coremat2 <- data.frame(do.call(cbind, lapply(det, function(d){
        rowSums(otu_table(pstat) > d) / nsamples(pstat)
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
#' @importFrom shiny uiOutput h4 textOutput mainPanel fluidRow column plotOutput
#' @importFrom DT dataTableOutput
#' @export
coreOTUModuleUI <- function(id, label = "Core OTUs") {
    # This is the namespace for the module
    ns <- NS(id)
    tabPanel("Core OTUs",
        sidebarLayout(
            sidebarPanel(
                sliderInput(ns("detThresh"), "Detection threshold", value=5, 
                    min=1, max=100, step=1),
                uiOutput(ns("prevThreshControl")),
                uiOutput(ns("taxLevelControl")),
                h4(textOutput(ns("coreSummary")), style="color:goldenrod;"),
                width=3
            ),
            mainPanel(
                fluidRow(column(12,
                    tabsetPanel(
                        tabPanel("2D lineplot",
                            plotOutput(ns("coreLine"), click = "plot_click")
                        ),
                        tabPanel("Heatmap", 
                            plotOutput(ns("coreHeat"))
                        )
                    )
                )),
                #fluidRow(column(12,
                #   h3(textOutput(ns("coreSummary")))
                #), style="padding:9px;"),                 
                fluidRow(column(12,
                    DT::dataTableOutput(ns("coreTable"), width='95%')
                ))
            ) # end mainPanel
        ) # end sidebarLayout
    ) # end tabPanel
}

#' Server function for Core OTU Module
#'
#' @param input Shiny server input object
#' @param output Shiny server output object
#' @param session This is the session
#' @param pstat PathoStat object
#' 
#' @return None
#' 
#' @importFrom shiny reactive renderUI sliderInput selectizeInput renderPlot
#' renderText NS tabPanel sidebarLayout sidebarPanel tabsetPanel
#' @export
coreOTUModule <- function(input, output, session, pstat) {
    glom <- reactive({
        if(is.null(input$taxLevel)) return()
        tax_glom(pstat, taxrank=input$taxLevel)
    })
    
    coremat <- reactive({
      v <- glom()
      if(is.null(v)) return()
      get_coremat(v)
    })
    
    lineplot <- reactive({
        v <- coremat()
        if(is.null(v)) return()
        get_coremat_lineplot(v)
    })
    
    output$prevThreshControl <- renderUI({
      ns <- session$ns
      sliderInput(ns("prevThresh"), "Samples detected",
                  value=5, step=1, min=1, max=nsamples(pstat))
    })
    
    output$taxLevelControl <- renderUI({
      ns <- session$ns
      rnames <- rev(rank_names(pstat))
      selectizeInput(ns("taxLevel") , 'Taxonomy Level', choices = rnames, 
                     selected=rnames[1])
    })

    output$coreTable <- DT::renderDataTable({
        cur <- glom()
        if(is.null(cur)) return(invisible())
        formatTaxTable(tax_table(get_core(cur, input$detThresh, input$prevThresh)))
    })

    output$coreLine <- renderPlot({
        cur <- coremat()
        cur_plot <- lineplot()
        if(is.null(cur) | is.null(cur_plot)) return(invisible())
        prevLabel <- paste0("Samples: ", input$prevThresh, " (",
                            sprintf('%.1f', 100*(input$prevThresh / nsamples(pstat)))
                            ,"% prevalence) ")
        cols <- ifelse(cur$prev == input$prevThresh, "#253494", "#bdbdbd33")
        cur_plot + 
            scale_colour_manual(values=cols, guide=FALSE) +
            annotate("text", x=Inf, y=Inf, vjust=1, hjust=1,
                     size=8, color="#253494",
                     label=prevLabel) +
            geom_vline(xintercept = input$detThresh, color="#990000", linetype=2)
    })
    
    output$coreHeat <- renderPlot({
        get_coremat_heatmap(glom())
    })
    
    output$coreSummary <- renderText({
          cur <- glom()
          if(is.null(cur)) return(invisible())
          nt <- ntaxa(get_core(cur, input$detThresh, input$prevThresh))
          paste0(nt, " core OTUs detected.")
    })
}
