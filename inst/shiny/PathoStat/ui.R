library(shiny)
library(ggvis)
library(d3heatmap)

shinyUI(navbarPage("PathoStat", id="PathoStat", fluid=TRUE, 
                   tabPanel("Relative Abundance",
                                tabsetPanel(
                                  tabPanel("Relative Abundance",ggvisOutput("RelAbundancePlot")),
                                  tabPanel("Summary", verbatimTextOutput("RAsummary")),
                                  tabPanel("Table", tableOutput("RAtable"))
                                )
                            )
                   )
)