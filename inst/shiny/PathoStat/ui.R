library(shiny)
library(ggvis)
library(d3heatmap)

tax.name <- c('superkingdom', 'kingdom', 'phylum', 'class', 'order', 'suborder', 'family', 
              'subfamily', 'genus', 'subgenus', 'species', 'no rank')
shinyUI(navbarPage("PathoStat", id="PathoStat", fluid=TRUE, 
                   tabPanel("Relative Abundance",
                            sidebarLayout(
                              sidebarPanel(
                                #selectizeInput('taxl', 'Taxonomy Level', choices = setNames(tax.abb, tax.name)),
                                selectizeInput('taxl', 'Taxonomy Level', choices = tax.name, selected='no rank'),
                                downloadButton('downloadData', 'Download CSV')
                              ),
                              mainPanel(
                                tabsetPanel(
                                  tabPanel("Taxonomy level RA",ggvisOutput("TaxRelAbundancePlot")),
                                  tabPanel("Summary", verbatimTextOutput("TaxRAsummary")),
                                  tabPanel("Table", tableOutput("TaxRAtable"))
                                )
                              )
                            )
                   )
)
)