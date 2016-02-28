library(shiny)
library(ggvis)
library(reshape2)

shinyServer(function(input, output, session) {
  #needed information from PathoStat
  
  setInputs <- function(combatFlag)  {
    if (combatFlag)  {
      shinyInput <<- shinyInputCombat
    } else  {
      shinyInput <<- shinyInputOrig
    }
  }
  #setInputs(FALSE)
  findTaxData <- eventReactive(input$taxl, {
    taxdata <- findTaxonLevelData(shinyInput$data, shinyInput$taxonLevels, input$taxl)
    if (is.null(shinyInput$taxdata))  {
      shinyInput <<- c(shinyInput, list("taxdata"=taxdata))
    } else  {
      shinyInput$taxdata <<- taxdata
    }
    shinyInput$taxdata
  })
  
  findTaxCountData <- eventReactive(input$taxl, {
    taxcountdata <- findTaxonLevelData(shinyInput$countdata, shinyInput$taxonLevels, input$taxl)
    if (is.null(shinyInput$taxcountdata))  {
      shinyInput <<- c(shinyInput, list("taxcountdata"=taxcountdata))
    } else  {
      shinyInput$taxcountdata <<- taxcountdata
    }
    shinyInput$taxcountdata
  })
  
  tax_ra_bp <- reactive({
    if (is.null(input$taxl))  {
      return()
    }
    taxdata <- findTaxData()
    dat <- melt(cbind(taxdata, ind = as.character(rownames(taxdata))), id.vars = c('ind'))
    dat %>%  
      ggvis(x=~variable, y=~value, fill=~as.factor(ind)) %>% layer_bars(stack = TRUE) %>%
      add_tooltip(function(dat2){paste0("Sample: ", dat2[2], "<br />", "Genome: ", dat2[1], 
                                        "<br />", "RA: ", round(dat2[4]-dat2[3], 4))}, "hover") %>%
      #add_axis("x", subdivide = 1, values = 1:length(colnames(shinyInput$data)), 
      add_axis("x", 
               title = "Samples", 
               properties = axis_props(
                 title = list(fontSize = 15),
                 labels = list(text="", fontSize = 10)
               )) %>%
      add_axis("y", title = "Relative Abundance (RA)", properties = axis_props(
        title = list(fontSize = 15),
        labels = list(fontSize = 10)
      )) %>%
      add_legend("fill", title = "Genomes", properties = legend_props(
        title = list(fontSize = 15),
        labels = list(fontSize = 10)
      )) %>%
      set_options(width = "auto", height = "auto")
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
  }, digits=0)
  
  output$downloadData <- downloadHandler(
    filename = function() { 
      paste0('sample_data_', input$taxl, '.csv', sep='') 
    },
    content = function(file) {
      write.csv(shinyInput$taxdata, file)
    }
  )

  output$downloadCountData <- downloadHandler(
    filename = function() { 
      paste0('sample_data_count_', input$taxl, '.csv', sep='') 
    },
    content = function(file) {
      write.csv(shinyInput$taxcountdata, file)
    }
  )
})

