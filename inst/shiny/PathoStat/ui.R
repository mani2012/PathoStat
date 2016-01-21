library(shiny)
library(ggvis)
library(d3heatmap)

tax.abb <- c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11", "L12")
tax.name <- c('superkingdom', 'kingdom', 'phylum', 'class', 'order', 'suborder', 'family', 
              'subfamily', 'genus', 'subgenus', 'species', 'no rank')
shinyUI(navbarPage("PathoStat", id="PathoStat", fluid=TRUE, 
                   tabPanel("Relative Abundance",
                                tabsetPanel(
                                  tabPanel("Relative Abundance",ggvisOutput("RelAbundancePlot")),
                                  tabPanel("Summary", verbatimTextOutput("RAsummary")),
                                  tabPanel("Table", tableOutput("RAtable"))
                                )
                            ),
                   tabPanel("Taxonomy Analysis",
                            selectizeInput('taxl', 'Taxonomy Level', choices = setNames(tax.abb, tax.name)),
                            tabsetPanel(
                              tabPanel("Taxonomy level RA",ggvisOutput("TaxRelAbundancePlot")),
                              tabPanel("Summary", verbatimTextOutput("TaxRAsummary")),
                              tabPanel("Table", tableOutput("TaxRAtable"))
                            )
                   )
)
)