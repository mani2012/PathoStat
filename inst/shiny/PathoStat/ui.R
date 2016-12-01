library(shiny)
library(ggvis)
library(d3heatmap)
library(phyloseq)
library(ape)

tax.name <- c('superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 
    'genus', 'species', 'no rank')
norm.methods <- c('EBayes coreOTU Normalization', 
    'Quantile coreOTU Normalization', 'Library Size Scaling')
measure.type <- c('Final Guess', 'Final Best Hit', 'Final High Confidence Hit')
minbatch <- function(batch1){
    batch2 <- as.factor(batch1)
    batch3 <- split(batch1,batch2)
    return(min(unlist(lapply(1:length(batch3), 
        function(x) length(batch3[[x]])))))
}

shinyInput <- getShinyInput()

pstat <- shinyInput$pstat
covariates <- colnames(sample_data(pstat))
maxbatchElems <- minbatch(c(pstat@sam_data[,1])[[1]])
maxcondElems <- minbatch(c(pstat@sam_data[,2])[[1]])
defaultDisp <- 30
defaultGenesDisp <- 10
maxGenes <- dim(pstat@otu_table)[1]
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
                    tabPanel("RA Table(%)", DT::dataTableOutput("TaxRAtable", 
                        width='95%')),
                    tabPanel("Count Table", DT::dataTableOutput("TaxCountTable",
                        width='95%'))
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
                height = "550px")),
            tabPanel("BiPlot", 
                sidebarLayout(
                    sidebarPanel(
                        selectizeInput('colorBiP', 'Color', choices = 
                            c(tax.name[-length(tax.name)], 'condition', 'None'), 
                            selected='genus'),
                        selectizeInput('shapeBiP', 'Shape', choices = 
                            c(tax.name[-length(tax.name)], 'condition', 'None'), 
                            selected='condition'),
                        selectizeInput('labelBiP', 'Label', choices = 
                            c(tax.name[-length(tax.name)], 'condition', 'None'), 
                            selected='None'),
                        selectizeInput('methodBiP', 'Method', 
                            choices=c("DCA", "CCA", "RDA", "DPCoA", 
                            "NMDS", "MDS", "PCoA"), selected='NMDS'),
                        width=3
                    ),
                    mainPanel(
                        plotOutput("BiPlot", height = "550px"), width=9
                    )
                )
            ),
            tabPanel("Co-Occurrence", 
                sidebarLayout(
                    sidebarPanel(
                        selectizeInput('colorCo', 'Color', choices = 
                            c(tax.name[-length(tax.name)], 'None'), 
                            selected='genus'),
                        selectizeInput('shapeCo', 'Shape', choices = 
                            c(tax.name[-length(tax.name)], 'None'), 
                            selected='None'),
                        selectizeInput('labelCo', 'Label', choices = 
                            c(tax.name[-length(tax.name)], 'None'), 
                            selected='None'),
                        sliderInput("max.dist", "Max Dist:", 
                            min = 0, max = 1, value = 0.5, step= 0.1),
                        width=3
                    ),
                    mainPanel(
                        plotOutput("CoOccurrence", height = "550px"), width=9
                    )
                )
            )
        )
    ),
    tabPanel("Differential Abundance",
        sidebarLayout(
            sidebarPanel(
                selectizeInput('taxlde', 'Taxonomy Level', choices = tax.name, 
                    selected='no rank'),
                selectizeInput('primary', 'Primary Covariate', 
                    choices = covariates, selected=covariates[2]),
                selectizeInput('secondary', 'Secondary Covariate', 
                    choices = covariates, selected=covariates[1]),
                selectizeInput('norm', 'Normalization', choices=norm.methods, 
                    selected='EBayes coreOTU Normalization'),
                actionButton("apply", "Apply"),
                sliderInput("otuthreshold", "OTU cutoff threshold:", 
                    min = 0, max = 1, value = 0.05, step= 0.01),
                sliderInput("prevalence", "OTU cutoff prevalence:", 
                    min = 0, max = 1, value = 0.4, step= 0.01),
                sliderInput("ebweight", "Empirical Bayes Weight:", 
                    min = 0, max = 1, value = 0.25, step= 0.01),
                numericInput('ncSamples', 
                    'No. of Sample(s) Per Primary Covariate', 
                    if (maxcondElems>defaultDisp) defaultDisp 
                    else maxcondElems, min = 1, max = maxcondElems),
                numericInput('noSamples', 
                    'No. of Sample(s) Per Secondary Covariate',  
                    if (maxbatchElems>defaultDisp) defaultDisp 
                    else maxbatchElems, min = 1, max = maxbatchElems),
                checkboxInput("sortbybatch", 
"Sort By Secondary Covariate First (Default: Sort By Primary Covariate First)", 
                    FALSE),
                checkboxInput("colbybatch", 
"Color By Secondary Covariate (Default: Color By Primary Covariate)", FALSE),
                numericInput('noTaxons', 
                    'No. of top Differentially Expressed Taxons to display', 
                    if (maxGenes>defaultGenesDisp) defaultGenesDisp 
                    else maxGenes, min = 1, max = maxGenes),
                width=3
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Expression Plots",ggvisOutput("DiffExPlot")), 
                    tabPanel("Summary", verbatimTextOutput("DEsummary")),
                    tabPanel("Table", tableOutput("DEtable")), 
                    tabPanel("LIMMA",tableOutput("LimmaTable")),
                    tabPanel("EdgeR",tableOutput("EdgeRTable")),
                    tabPanel("DeSeq2",tableOutput("DeSeq2Table"))
                ), width=9
            )
            )
    ),
    tabPanel("Confidence Region",
        sidebarLayout(
            sidebarPanel(
                #selectizeInput('taxlcr', 'Taxonomy Level', choices = tax.name, 
                #    selected='no rank'),
                selectizeInput('taxon1', 'Taxon 1', choices=row.names(
                    shinyInput$pstat@otu_table)),
                selectizeInput('taxon2', 'Taxon 2', choices=row.names(
                    shinyInput$pstat@otu_table)),
                selectizeInput('sample', 'Sample', choices=colnames(
                    shinyInput$pstat@otu_table)),
                checkboxInput("uselogit", 
                    "Use Logit Transformation", FALSE),
                width=5
            ),
            mainPanel(
                plotOutput("confRegion", height = "550px"), width=7
            )
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
    ),
    tabPanel("Time Series",
        tabsetPanel(
            tabPanel("Visualization",
                sidebarLayout(
                    sidebarPanel(
                      selectInput(inputId="Allusset", 
                          label="Visualization column", 
                          choices = colnames(shinyInput$pstat@sam_data)),
                      checkboxInput(inputId="Allurar", 
                          label="Rarefaction? (maximum reads of minimal 
                          sample count)"),
                      selectInput(inputId="Alluglom", label="Agglomerate taxa", 
                          choices = colnames(shinyInput$pstat@tax_table)),
                      uiOutput("Allustax"),
                      downloadButton('downloadAlluvialPlot', 
                                     'Download Plot')
                    ),
                    mainPanel(
                      plotOutput("TimePlotVisu",height = "600px")
                    )
                )
            )
        )
    ),
    coreOTUModuleUI("coreOTUModule")
)
)
