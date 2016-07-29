
tax.name <- c('superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family',
              'genus', 'species', 'no rank')

availableDistances <- c(phyloseq::distanceMethodList$UniFrac,
                        phyloseq::distanceMethodList$vegdist,
                        phyloseq::distanceMethodList$JSD,
                        phyloseq::distanceMethodList$dist
)

#' Try to identify the main covariate from sample data.
#' Assume that the main covariate is a variable with 2-6 levels, inclusive.
#' If multiple variables satisfy this criteria choose the first column on the
#' left. If none are found return the first variable.
getMainCovariate <- function(pstat) {
    # Get the number of levels for each variable
    nlevels <- sapply(sample_variables(pstat), function(v) {
        length(levels(factor(sample_data(pstat)[[v]])))
    })
    v <- sample_variables(pstat)[nlevels > 1 & nlevels <= 6]
    if(length(v)>0) return(v[1])
    sample_variables(pstat)[1]
}

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
                                                            height = "550px")
                     ), # end tabPanel
                     tabPanel("Ordination",
                              sidebarLayout(
                                  sidebarPanel(
                                      uiOutput(ns("ordMethodControl")),
                                      uiOutput(ns("ordDistanceControl")),
                                      uiOutput(ns("ordFormulaControl"))
                                  ),
                                  mainPanel(
                                      plotOutput(ns("ordinationPlot"), height = "500px")
                                  ) # end mainPanel
                              ) # end sidebarLayout
                     ), #end tabPanel
                     tabPanel("Beta Diversity",
                              checkboxInput(ns("methodBeta"),
                                            "Weigthed Unifrac (Default: Bray-Curtis)", FALSE),
                              plotOutput(ns("BetaDiversity"), height = "500px")
                     ), # end tabPanel
                     tabPanel("Exploratory Tree", plotOutput(ns("ExploratoryTree"),
                                                             height = "550px")
                     ), # end tabPanel
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

    output$ordMethodControl <- renderUI({
        ns <- session$ns
        selectizeInput(ns("ordinationMethod"),
                       'Ordination Method',
                       choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA",
                                   "NMDS", "MDS", "PCoA"),
                       selected="DCA")
    })

    output$ordDistanceControl <- renderUI({
        ns <- session$ns
        if(is.null(input$ordinationMethod)) return(invisible())
        if(! input$ordinationMethod %in% c("CAP", "NMDS", "MDS", "PCoA")) return(invisible())
        selectizeInput(ns("distanceMethod"),
                           'Distance Method',
                           choices = availableDistances,
                       selected="bray"
        )
    })

    output$ordFormulaControl <- renderUI({
        ns <- session$ns
        if(is.null(input$ordinationMethod)) return()
        if(! input$ordinationMethod %in% c("CCA", "RDA", "CAP")) return()
        if(input$ordinationMethod == 'CAP')
            return(selectizeInput(ns("ordFormula"),
                           "Formula",
                           choices = sample_variables(pstat),
                           selected = getMainCovariate(pstat),
                           multiple = TRUE))
        selectizeInput(ns("ordFormula"),
                              "Formula",
                              choices = c("Conditioning variables..." = "", sample_variables(pstat)),
                              multiple = TRUE
                       )
        
        
        # selectizeInput("test","test",choices = sample_variables(pstat),multiple=TRUE)
        #                
        #     'e2', '2. Multi-select', choices = state.name, multiple = TRUE
        # ),
        # 
        # if(input$ordinationMethod %in% c("CCA", "RDA"))
        #     return(textInput(ns("ordFormula"),
        #                'Formula',
        #                 placeholder = "Constraining variables, ex. ~Condition"
        #     ))
        # # ordination method is CAP, formula is required
        # textInput(ns("ordFormula"),
        #               'Formula',
        #               value = paste0('~', getMainCovariate(pstat))
        # )
    })

    output$ordinationPlot <- renderPlot({
        # DCA - none
        # CCA - formula
        # RDA - formula
        # CAP - formula, distance
        # DPCoA - none
        # NMDS - distance
        # MDS - distance
        # PCoA - distance
        
        if(is.null(input$ordinationMethod)) return(invisible())
        ordInput  <- input$ordinationMethod
        distInput <- input$distanceMethod
        frmInput  <-  input$ordFormula
        
        cat("Ordination method ", ordInput, '\n')
        cat("Distance method ", distInput, '\n')
        cat("Formula ", frmInput, '\n')

        if(ordInput=='DCA')
            pstat.ord <- phyloseq::ordinate(pstat,method=ordInput)        
        if(ordInput=='DPCoA')
            # Need to create phyloseq object for some reason
            pstat.ord <- phyloseq::ordinate(phyloseq(otu_table(pstat),phy_tree(pstat)), 
                                                     method=ordInput)
        if(ordInput %in% c('NMDS','MDS','PCoA')) {
            if(is.null(distInput)) return(invisible())
            pstat.ord <- phyloseq::ordinate(pstat, method=ordInput,
                                            distance = distInput)
        }
        if(ordInput %in% c('CCA','RDA')) {
            if(is.null(frmInput)) {
                cat("No formula")
                pstat.ord <- phyloseq::ordinate(pstat, method=ordInput)
            } else {
                frmTxt <- paste0('~ ', paste(frmInput,collapse=' + '))
                cat(frmTxt, '\n')
                frmObj <- formula(frmTxt)
                pstat.ord <- phyloseq::ordinate(pstat, method=ordInput,
                                                formula=frmObj)
            }
        }
        if(ordInput == "CAP") {
                if(is.null(frmInput) | is.null(distInput)) return(invisible())
                frmTxt <- paste0('~ ', paste(frmInput,collapse=' + '))
                cat(frmTxt, '\n')
                frmObj <- formula(frmTxt)
                pstat.ord <- phyloseq::ordinate(pstat, method=ordInput,
                                                formula = frmObj,
                                                distance = distInput)
        }
        
        phyloseq::plot_ordination(pstat, pstat.ord)
    })

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


