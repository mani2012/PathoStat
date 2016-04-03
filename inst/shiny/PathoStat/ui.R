library(shiny)
library(ggvis)
library(d3heatmap)

tax.name <- c('superkingdom', 'kingdom', 'phylum', 'class', 'order', 'suborder',
    'family', 'subfamily', 'genus', 'subgenus', 'species', 'subspecies', 
    'no rank')
measure.type <- c('Final Guess', 'Final Best Hit', 'Final High Confidence Hit')
sort.by <- c('None', 'Batch', 'Condition')
shinyUI(navbarPage("PathoStat", id="PathoStat", fluid=TRUE, 
    tabPanel("Relative Abundance",
        #selectizeInput('taxl', 'Taxonomy Level', choices = setNames(tax.abb, 
        #   tax.name)),
        selectizeInput('taxl', 'Taxonomy Level', choices = tax.name, 
            selected='no rank'),
        downloadButton('downloadData', 'Download RA CSV'),
        downloadButton('downloadCountData', 'Download Count CSV'),
        br(),
        p(" "),
        tabsetPanel(
            tabPanel("Taxonomy level RA",ggvisOutput("TaxRelAbundancePlot")),
            tabPanel("Summary", verbatimTextOutput("TaxRAsummary")),
            tabPanel("Table", tableOutput("TaxRAtable")),
            tabPanel("Count Table", tableOutput("TaxCountTable"))
        )
    ),
    tabPanel("Diversity",
        tabsetPanel(
            tabPanel("Alpha Diversity", plotOutput("AlphaDiversity",
                height = "550px")),
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
                    "Color By Batch (Default: Color By Condition)", FALSE)
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("PCA", ggvisOutput("PCAplot")),
                    tabPanel("Summary", verbatimTextOutput("PCAsummary")),
                    tabPanel("Table",tableOutput("PCAtable")),
                    tabPanel("Explained Variation",
                        tableOutput("PCAExplainedVariation"))
                )
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
                    "Weigthed Unifrac (Default: Bray-Curtis)", FALSE)
            ),
            mainPanel(
                plotOutput("PCoAplot", height = "550px")
            )
        )
    )
)
)
