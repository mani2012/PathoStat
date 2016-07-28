
tax.name <- c('superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family',
              'genus', 'species', 'no rank')

#' UI function for Diversity Module
#'
#' @param id Namespace for module
#' @param label Tab label
#'
#' @return None
#'
#' @export
diversityModuleUI <-
    function(id, label = "Diversity2") {
        # This is the namespace for the module
        ns <- NS(id)
        tabPanel(label,
                 tabsetPanel(
                     tabPanel("Alpha Diversity", plotOutput(ns("AlphaDiversity"),
                                                            height = "550px")),
                     tabPanel("Beta Diversity",
                              checkboxInput(ns("methodBeta"),
                                            "Weigthed Unifrac (Default: Bray-Curtis)", FALSE),
                              plotOutput(ns("BetaDiversity"), height = "500px")),
                     tabPanel("Exploratory Tree", plotOutput(ns("ExploratoryTree"),
                                                             height = "550px")),
                     tabPanel("BiPlot",
                              sidebarLayout(
                                  sidebarPanel(
                                      selectizeInput(ns('colorBiP'), 'Color', choices =
                                                         c(tax.name[-length(tax.name)], 'condition', 'None'),
                                                     selected='genus'),
                                      selectizeInput(ns('shapeBiP'), 'Shape', choices =
                                                         c(tax.name[-length(tax.name)], 'condition', 'None'),
                                                     selected='condition'),
                                      selectizeInput(ns('labelBiP'), 'Label', choices =
                                                         c(tax.name[-length(tax.name)], 'condition', 'None'),
                                                     selected='None'),
                                      selectizeInput(ns('methodBiP'), 'Method',
                                                     choices=c("DCA", "CCA", "RDA", "DPCoA",
                                                               "NMDS", "MDS", "PCoA"), selected='NMDS'),
                                      width=3
                                  ),
                                  mainPanel(
                                      plotOutput(ns("BiPlot"), height = "550px"), width=9
                                  )
                              )
                     ),
                     tabPanel("Co-Occurrence",
                              sidebarLayout(
                                  sidebarPanel(
                                      selectizeInput(ns('colorCo'), 'Color', choices =
                                                         c(tax.name[-length(tax.name)], 'None'),
                                                     selected='genus'),
                                      selectizeInput(ns('shapeCo'), 'Shape', choices =
                                                         c(tax.name[-length(tax.name)], 'None'),
                                                     selected='None'),
                                      selectizeInput(ns('labelCo'), 'Label', choices =
                                                         c(tax.name[-length(tax.name)], 'None'),
                                                     selected='None'),
                                      sliderInput(ns("max.dist"), "Max Dist:",
                                                  min = 0, max = 1, value = 0.5, step= 0.1),
                                      width=3
                                  ),
                                  mainPanel(
                                      plotOutput(ns("CoOccurrence"), height = "550px"), width=9
                                  )
                              )
                     )
                 )
        )
    }

#' Server function for Diversity Module
#'
#' @param input Shiny server input object
#' @param output Shiny server output object
#' @param session This is the session
#' @param pstat PathoStat object
#'
#' @return None
#'
#' @export
diversityModule <- function(input, output, session, pstat) {

    output$AlphaDiversity <- renderPlot({
        # physeq1 <- shinyInput$pstat
        cn <- colnames(pstat@sam_data)
        cn[1] <- "batch"
        cn[2] <- "condition"
        colnames(pstat@sam_data) <- cn
        alpha_meas <- c("Shannon", "Simpson", "InvSimpson")
        # setGgplotTheme()
        (p <- plot_richness(pstat, "condition", "batch", measures=alpha_meas))
        p + ggplot2::geom_boxplot(data=p$data, ggplot2::aes(x=condition,
                                                            y=value, color=NULL), alpha=0.1)
    })
    output$BetaDiversity <- renderPlot({
        # physeq1 <- shinyInput$pstat
        # setGgplotTheme()
        if (input$methodBeta)  {
            dist = distance(pstat, method = "wUniFrac")
            titleString="Beta Diversity Distance: Weigthed Unifrac"
        } else  {
            dist = distance(pstat, method = "bray")
            titleString="Beta Diversity Distance: Bray-Curtis"
        }
        gplots::heatmap.2(as.matrix(dist), col=gplots::bluered(75), scale="row",
                          key=TRUE, symkey=FALSE, density.info="none", trace="none",
                          margins = c(6, 6))
    })
    output$ExploratoryTree <- renderPlot({
        # physeq1 <- shinyInput$pstat
        cn <- colnames(pstat@sam_data)
        cn[1] <- "batch"
        cn[2] <- "condition"
        colnames(pstat@sam_data) <- cn
        # setGgplotTheme()
        plot_tree(pstat, color="condition", label.tips="genus",
                  size="Abundance")
    })
    output$BiPlot <- renderPlot({
        # physeq1 <- shinyInput$pstat
        cn <- colnames(pstat@sam_data)
        cn[1] <- "batch"
        cn[2] <- "condition"
        colnames(pstat@sam_data) <- cn
        # setGgplotTheme()
        if (input$colorBiP=='None') color <- NULL else color <- input$colorBiP
        if (input$shapeBiP=='None') shape <- NULL else shape <- input$shapeBiP
        if (input$labelBiP=='None') label <- NULL else label <- input$labelBiP
        physeq.ord <- ordinate(pstat, method=input$methodBiP, distance="bray")
        p <- plot_ordination(pstat, physeq.ord, type = "biplot",
                             color=color, shape=shape, label=label)
        if (!is.null(label))  {
            p$layers <- p$layers[-2]
            p <- p+ggplot2::geom_text(mapping=
                                          ggplot2::aes(label=get(label)), size=4, vjust=1.0,
                                      check_overlap=TRUE)
        }
        p
    })
    output$CoOccurrence <- renderPlot({
        # physeq1 <- shinyInput$pstat
        cn <- colnames(pstat@sam_data)
        cn[1] <- "batch"
        cn[2] <- "condition"
        colnames(pstat@sam_data) <- cn
        # setGgplotTheme()
        if (input$colorCo=='None') color <- NULL else color <- input$colorCo
        if (input$shapeCo=='None') shape <- NULL else shape <- input$shapeCo
        if (input$labelCo=='None') label <- NULL else label <- input$labelCo
        ig<-make_network(pstat, dist.fun="jaccard", type = "taxa",
                         max.dist=input$max.dist) #max.dist=0.4 default
        p <- plot_network(ig, pstat, line_weight=0.4, type = "taxa",
                          color = color, shape = shape, label = label)
        if (!is.null(label))  {
            p$layers <- p$layers[-2]
            p <- p+ggplot2::geom_text(mapping=ggplot2::aes(label=get(label)),
                                      size=4, vjust=1.0, check_overlap=TRUE)
        }
        p
    })
}