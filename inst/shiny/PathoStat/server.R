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
  ra_bp <- reactive({
    dat <- melt(cbind(shinyInput$data, ind = rownames(shinyInput$data)), id.vars = c('ind'))
    dat %>%  
      ggvis(x=~variable, y=~value, fill=~as.factor(ind)) %>% layer_bars(stack = TRUE) %>%
      add_tooltip(function(dat2){paste0("Sample: ", dat2[2], "<br />", "Genome: ", dat2[1], 
                                       "<br />", "RA: ", dat2[4]-dat2[3])}, "hover") %>%
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
  ra_bp %>% bind_shiny("RelAbundancePlot")
  output$RAsummary <- renderPrint({
    summary(shinyInput$data)
  })
  
  output$RAtable <- renderTable({
    shinyInput$data
  })
  
  tax_ra_bp <- reactive({
    taxdata <- findTaxonLevelData(shinyInput$data, shinyInput$taxonLevels, input$taxl)
    if (is.null(shinyInput$taxdata))  {
      shinyInput <<- c(shinyInput, list("taxdata"=taxdata))
    } else  {
      shinyInput$taxdata <<- taxdata
    }
    dat <- melt(cbind(taxdata, ind = as.character(rownames(taxdata))), id.vars = c('ind'))
    dat %>%  
      ggvis(x=~variable, y=~value, fill=~as.factor(ind)) %>% layer_bars(stack = TRUE) %>%
      add_tooltip(function(dat2){paste0("Sample: ", dat2[2], "<br />", "Genome: ", dat2[1], 
                                        "<br />", "RA: ", dat2[4]-dat2[3])}, "hover") %>%
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

  raSummary <- eventReactive(input$taxl, {
    summary(shinyInput$taxdata)
  })
  output$TaxRAsummary <- renderPrint({
    raSummary()
  })
  
  raTable <- eventReactive(input$taxl, {
    shinyInput$taxdata
  })
  output$TaxRAtable <- renderTable({
    raTable()
  })
  
})

