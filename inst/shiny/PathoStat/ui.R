library(shiny)
library(ggvis)
library(d3heatmap)

tax.name <- c('superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 
    'genus', 'species', 'no rank')
measure.type <- c('Final Guess', 'Final Best Hit', 'Final High Confidence Hit')
shinyUI(navbarPage("PathoStat", id="PathoStat", fluid=TRUE, 
    tabPanel("Relative Abundance",
        sidebarLayout(
            sidebarPanel(
                #selectizeInput('taxl', 'Taxonomy Level', choices = setNames(
                #   tax.abb, tax.name)),
                selectizeInput('taxl', 'Taxonomy Level', choices = tax.name, 
                    selected='no rank'),
                downloadButton('downloadData', 'Download RA CSV'),
                downloadButton('downloadCountData', 'Download Count CSV'),
                br(),
                p(" "),
                radioButtons('sortBy', 'Sort By',
                    c('None'=0, 'Condition'=1,'Batch'=2), 0),
                width=3
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Taxonomy level RA",
                        ggvisOutput("TaxRelAbundancePlot")),
                    tabPanel("Heatmap", plotOutput("Heatmap", height="550px")),
                    tabPanel("Summary", verbatimTextOutput("TaxRAsummary")),
                    tabPanel("Table", tableOutput("TaxRAtable")),
                    tabPanel("Count Table", tableOutput("TaxCountTable"))
                ), width=9
            )
        )
    ),
    tabPanel("Diversity",
        tabsetPanel(
            tabPanel("Alpha Diversity", plotOutput("AlphaDiversity",
                height = "550px")),
            tabPanel("Beta Diversity", 
                checkboxInput("methodBeta", 
                    "Weigthed Unifrac (Default: Bray-Curtis)", FALSE),
                plotOutput("BetaDiversity", height = "500px")),
            tabPanel("Exploratory Tree", plotOutput("ExploratoryTree", 
                height = "550px"))
        )
    ),
    tabPanel("PCA",
        sidebarLayout(
            sidebarPanel(
                numericInput('xcol', 'Principal Component (x-axis)', 1,
                    min = 1, max = 50),
                numericInput('ycol', 'Principal Component (y-axis)', 2,
                    min = 1, max = 50),
                checkboxInput("colbybatchPCA", 
                    "Color By Batch (Default: Color By Condition)", FALSE),
                width=3
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("PCA", ggvisOutput("PCAplot")),
                    tabPanel("Summary", verbatimTextOutput("PCAsummary")),
                    tabPanel("Table",tableOutput("PCAtable")),
                    tabPanel("Explained Variation",
                        tableOutput("PCAExplainedVariation"))
                ), width=9
            )
        )
    ), 
    tabPanel("PCoA",
        sidebarLayout(
            sidebarPanel(
                numericInput('xcolA', 'Principal Coordinate (x-axis)', 1,
                    min = 1, max = 50),
                numericInput('ycolA', 'Principal Coordinate (y-axis)', 2,
                    min = 1, max = 50),
                checkboxInput("colbybatchPCoA", 
                    "Color By Batch (Default: Color By Condition)", FALSE),
                checkboxInput("methodPCoA", 
                    "Weigthed Unifrac (Default: Bray-Curtis)", FALSE),
                width=3
            ),
            mainPanel(
                plotOutput("PCoAplot", height = "550px"), width=9
            )
        )
    )
)
)
