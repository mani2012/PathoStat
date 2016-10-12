#' Select rows of OTU matrix that meet given detection and prevalence 
#' thresholds
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
#' @importFrom tidyr %>%
#' @importFrom phyloseq otu_table nsamples
get_coremat <- function(pstat) {
    # Get the OTU counts table
    otu_counts <- otu_table(pstat)
    
    # Detection values to calculate
    det <- unlist(lapply(0:7, function(p){c(1,2,5) * 10^p}))
    det <- det[1:min(which(det > max(otu_counts)))]
    # Prevalence values to calculate
    prev <- seq(1, nsamples(pstat))
    
    # Difficult to parse:
    # coremat <- data.frame(do.call(rbind, lapply(prev, function(p){
    #     sapply(det, function(d) sum(rowSums(otu_counts >= d) >= p))
    # })))
    
    # Number of samples where OTU is detected
    # Rows are OTUs, columns are detection thresholds
    nsamp_det <- sapply(det, function(d) rowSums(otu_counts >= d))
    
    # Number of OTUs for each prevalence and detection threshold
    # Rows are prevalence (number of samples), columns are detection thresholds
    coremat <- sapply(prev, function(p) colSums(nsamp_det >= p))
    coremat <- t(coremat) %>% data.frame
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
#' @importFrom dplyr mutate_
#' @importFrom tidyr gather_ %>%
get_coremat_lineplot <- function(coremat) {
    cols_to_gather <- colnames(coremat)[colnames(coremat) != "prev"]
    coremat %>% 
        dplyr::mutate_(prev="factor(prev)") %>%
        tidyr::gather_("det", "count", cols_to_gather) %>%
        dplyr::mutate_(det="as.numeric(det)") %>%
        ggplot(aes_string(x="det", y="count", group="prev", color="prev")) + 
            geom_point() + 
            geom_line() +
            scale_x_log10() +
            labs(x="Detection threshold", y="OTU count")
}

#' @import grDevices
#' @importFrom tidyr %>%
#' @importFrom phyloseq nsamples
get_coremat_heatmap <- function(pstat) {
    # Get the OTU counts table
    otu_counts <- otu_table(pstat)

    # Detection values to calculate    
    det <- 10^seq(0,log10(max(otu_counts, na.rm = TRUE)), length = 20)

    # Difficult to parse:
    # coremat2 <- data.frame(do.call(cbind, lapply(det, function(d){
    #    rowSums(otu_counts > d) / nsamples(pstat)
    #})))
    
    # Number of samples where OTU is detected
    # Rows are OTUs, columns are detection thresholds
    nsamp_det <- sapply(det, function(d) rowSums(otu_counts >= d))

    # Proportion of samples where OTU is detected
    coremat2 <- nsamp_det / nsamples(pstat) 
    coremat2 <- coremat2 %>% data.frame
    colnames(coremat2) <- det
    
    # Reorder from most to least prevalent 
    taxorder <- rownames(coremat2[order(-rowSums(coremat2)),])
    coremat2 <- coremat2[taxorder, ]
    coremat2$Taxa <- factor(rownames(coremat2), levels=taxorder)

    coremat2 %>%
        tidyr::gather_("det", "prev", as.character(det)) %>% 
        dplyr::mutate(det=as.numeric(det)) %>%
        ggplot(aes_string(x="det", y="Taxa", fill="prev")) + geom_tile() + 
            scale_x_log10() +
            scale_fill_gradientn("Prevalence", 
                breaks = seq(from = 0, to = 1, by = 0.2),
                colours = gray(seq(0,1,length=5)))
}


#' UI function for Core OTU Module
#' 
#' This function creates the UI for the Core OTU tab. The tab panel can be
#' included within a tabsetPanel, thus providing a simple way to add or remove
#' this module from the Shiny app. The first argument, \code{id}, is the ID to 
#' be used for the namespace \emph{and} must match the \code{id} argument
#' provided to \code{\link{coreOTUModule}}.
#' 
#' @param id Namespace for module
#' @param label Tab label
#' 
#' @return A \code{\link[shiny]{tabPanel}} that can be included within a 
#' \code{\link[shiny]{tabsetPanel}}.
#' 
#' @importFrom shiny uiOutput h4 textOutput mainPanel fluidRow column plotOutput
#' @importFrom DT dataTableOutput
#' @export
#' @examples
#' shiny::mainPanel(
#'     shiny::tabsetPanel(
#'         coreOTUModuleUI("coreOTUModule")
#'     )
#' )
#' 
#' @seealso \code{\link{coreOTUModule}} for the server function, 
#'     \code{\link[shiny]{tabPanel}} for the UI component returned by this 
#'     function, or \url{ http://shiny.rstudio.com/articles/modules.html} for
#'     more information about Shiny modules.
#' 
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
#' This function provides the server logic for the Core OTU tab. This function 
#' is not called directly; instead, it should be invoked within the Shiny app's 
#' server function using the \code{\link[shiny]{callModule}} function.
#' See \url{ http://shiny.rstudio.com/articles/modules.html} for information
#' about this design pattern.
#' 
#' The \code{\link[shiny]{callModule}} function should be invoked with this
#' function as the first argument. \code{callModule} is responsible for creating
#' the namespaced \code{input}, \code{output}, and \code{session} arguments.
#' The second argument to \code{callModule} is the ID to be used for the 
#' namespace and \emph{must} match the \code{id} argument provided to 
#' \link{coreOTUModuleUI}. The third argument to \code{callModule} should
#' be a \linkS4class{PathoStat} object from the app's server function, and is
#' passed to this function as the \code{pstat} argument.
#'
#' @param input Shiny server input object created by \code{callModule}
#' @param output Shiny server output object created by \code{callModule}
#' @param session Session created by \code{callModule}
#' @param pstat PathoStat object (third argument to \code{callModule}).
#' 
#' @return None
#' 
#' @importFrom shiny reactive renderUI sliderInput selectizeInput renderPlot
#' renderText NS tabPanel sidebarLayout sidebarPanel tabsetPanel
#' @export
#' @examples
#' # This function is not called directly; instead, it should be invoked within
#' # the app's server function using the shiny::callModule function.
#' \dontrun{
#' shinyServer(function(input, output, session) {
#'     shinyInput <- getShinyInput()
#'     pstat <- shinyInput$pstat
#'     callModule( coreOTUModule, "coreOTUModule", pstat )
#' }
#' }
#' 
#' @seealso \code{\link{coreOTUModuleUI}} for the UI function, 
#'     \code{\link[shiny]{callModule}} to see how to invoke this function, or
#'     \url{ http://shiny.rstudio.com/articles/modules.html} for more 
#'     information about Shiny modules.
#'
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
        formatTaxTable(tax_table(get_core(cur, input$detThresh, 
            input$prevThresh)))
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
                size=8, color="#253494", label=prevLabel) +
            geom_vline(xintercept = input$detThresh, color="#990000", 
                linetype=2)
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
