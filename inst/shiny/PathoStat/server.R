library(shiny)
library(ggvis)
library(d3heatmap)
library(reshape2)
library(phyloseq)
library(ape)
library(PathoStat)

shinyServer(function(input, output, session) {
    # needed information from PathoStat
    shinyInput <- getShinyInput()
    
    setInputs <- function(combatFlag) {
        if (combatFlag) {
            setShinyInput(getShinyInputCombat())
        } else {
            setShinyInput(getShinyInputOrig())
        }
    }
    # setInputs(FALSE)
    findAllTaxData <- function() {
        taxdata <- findTaxonLevelData(shinyInput$data, shinyInput$taxonLevels, 
            input$taxl)
        if (is.null(shinyInput$taxdata)) {
            shinyInput <- c(shinyInput, list(taxdata = taxdata))
        } else {
            shinyInput$taxdata <- taxdata
        }
        taxcountdata <- findTaxonLevelData(shinyInput$countdata, 
            shinyInput$taxonLevels, input$taxl)
        if (is.null(shinyInput$taxcountdata)) {
            shinyInput <- c(shinyInput, list(taxcountdata = taxcountdata))
        } else {
            shinyInput$taxcountdata <- taxcountdata
        }
        setShinyInput(shinyInput)
    }
    findTaxData <- eventReactive(input$taxl, {
        findAllTaxData()
        shinyInput <- getShinyInput()
        shinyInput$taxdata
    })
    
    findTaxCountData <- eventReactive(input$taxl, {
        findAllTaxData()
        shinyInput <- getShinyInput()
        shinyInput$taxcountdata
    })
    
    tax_ra_bp <- reactive({
        if (is.null(input$taxl)) {
            return()
        }
        taxdata <- findTaxData()
        dat <- melt(cbind(taxdata, ind = as.character(rownames(taxdata))), 
            id.vars = c("ind"))
        dat %>% ggvis(x = ~variable, y = ~value, fill = ~as.factor(ind)) %>% 
        layer_bars(stack = TRUE) %>% 
        add_tooltip(function(dat2) {
            paste0("Sample: ", dat2[2], "<br />", "Genome: ", dat2[1], 
                "<br />", "RA: ", round(dat2[4] - dat2[3], 4))
        }, "hover") %>% 
        # add_axis('x', subdivide = 1, 
        #   values = 1:length(colnames(shinyInput$data)),
        add_axis("x", title = "Samples", properties = axis_props(title = 
            list(fontSize = 15), labels = list(text = "", 
            fontSize = 10))) %>% 
        add_axis("y", title = "Relative Abundance (RA)", properties = 
            axis_props(title = list(fontSize = 15), 
            labels = list(fontSize = 10))) %>% add_legend("fill", title = 
            "Genomes", properties = legend_props(title = list(fontSize = 15), 
            labels = list(fontSize = 10))) %>% set_options(width = "auto", 
            height = "auto")
    })
    tax_ra_bp %>% bind_shiny("TaxRelAbundancePlot")
    
    output$TaxRAsummary <- renderPrint({
        summary(findTaxData())
    })
    
    output$TaxRAtable <- renderTable({
        findTaxData()
    })
    
    output$TaxCountTable <- renderTable({
        findTaxCountData()
    }, digits = 0)
    
    output$downloadData <- downloadHandler(filename = function() {
        paste0("sample_data_", input$taxl, ".csv", sep = "")
    }, content = function(file) {
        write.csv(shinyInput$taxdata, file)
    })
    
    output$downloadCountData <- downloadHandler(filename = function() {
        paste0("sample_data_count_", input$taxl, ".csv", sep = "")
    }, content = function(file) {
        write.csv(shinyInput$taxcountdata, file)
    })
    
    findPhyseqData <- function() {
        ids <- rownames(shinyInput$data)
        taxmat <- findTaxonMat(ids, shinyInput$taxonLevels)
        OTU <- otu_table(shinyInput$countdata, taxa_are_rows = TRUE)
        TAX <- tax_table(taxmat)
        physeq <- phyloseq(OTU, TAX)
        sampledata = sample_data(data.frame(condition=as.factor(condition), 
            batch=as.factor(batch), row.names=sample_names(physeq), 
            stringsAsFactors=FALSE))
        random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=
            taxa_names(physeq))
        physeq1 <- merge_phyloseq(physeq, sampledata, random_tree)
        return(physeq1)
    }
    setGgplotTheme <- function() {
        ggplot2::theme_set(ggplot2::theme_bw())
    }
    pal = "Set1"
    scale_colour_discrete <- function(palname = pal, ...) {
        scale_colour_brewer(palette = palname, ...)
    }
    scale_fill_discrete <- function(palname = pal, ...) {
        scale_fill_brewer(palette = palname, ...)
    }
    
    output$AlphaDiversity <- renderPlot({
        physeq1 <- findPhyseqData()
        #alpha_meas <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", 
        #    "InvSimpson")
        #alpha_meas <- c("Observed", "Chao1", "Shannon", "Simpson", 
        #    "InvSimpson")
        alpha_meas <- c("Shannon", "Simpson", "InvSimpson")
        setGgplotTheme()
        (p <- plot_richness(physeq1, "condition", "batch", measures=alpha_meas))
        p + ggplot2::geom_boxplot(data=p$data, ggplot2::aes(x=condition, 
            y=value, color=NULL), alpha=0.1)
    })
    output$ExploratoryTree <- renderPlot({
        physeq1 <- findPhyseqData()
        setGgplotTheme()
        plot_tree(physeq1, color="condition", label.tips="genus", 
            size="Abundance")
    })
})
