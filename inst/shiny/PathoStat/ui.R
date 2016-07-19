library(shiny)
library(ggvis)
library(d3heatmap)

tax.name <- c('superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 
    'genus', 'species', 'no rank')
measure.type <- c('Final Guess', 'Final Best Hit', 'Final High Confidence Hit')
minbatch <- function(batch1){
    batch2 <- as.factor(batch1)
    batch3 <- split(batch1,batch2)
    return(min(unlist(lapply(1:length(batch3), 
        function(x) length(batch3[[x]])))))
}
findphyloseqData <- function() {
  ids <- rownames(shinyInput$data)
  taxmat <- findTaxonMat(ids, shinyInput$taxonLevels)
  OTU <- otu_table(shinyInput$countdata, taxa_are_rows = TRUE)
  TAX <- tax_table(taxmat)
  physeq <- phyloseq(OTU, TAX)
  sampledata = sample_data(data.frame(condition=as.factor(
    shinyInput$condition), batch=as.factor(shinyInput$batch), 
    row.names=sample_names(physeq), stringsAsFactors=FALSE))
  random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=
                        taxa_names(physeq))
  physeq1 <- merge_phyloseq(physeq, sampledata, random_tree)
  return(physeq1)
}
shinyInput <- getShinyInput()
phyloseq1 <- findphyloseqData()
covariates <- colnames(sample_data(phyloseq1))
maxbatchElems <- minbatch(shinyInput$batch)
maxcondElems <- minbatch(shinyInput$condition)
defaultDisp <- 30
defaultGenesDisp <- 10
maxGenes <- dim(shinyInput$data)[1]
nbatch <- nlevels(as.factor(shinyInput$batch))
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
                    tabPanel("Table", DT::dataTableOutput("TaxRAtable", width='95%')),
                    tabPanel("Count Table", DT::dataTableOutput("TaxCountTable", width='95%'))
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
    tabPanel("Differential Expression",
        sidebarLayout(
            sidebarPanel(
                selectizeInput('taxlde', 'Taxonomy Level', choices = tax.name, 
                    selected='no rank'),
                selectizeInput('primary', 'Primary Covariate', choices = covariates),
                selectizeInput('secondary', 'Secondary Covariate', choices = covariates),
                numericInput('npSamples', 'No. of Sample(s) Per Primary Covariate', 
                    if (maxcondElems>defaultDisp) defaultDisp 
                    else maxcondElems, min = 1, max = maxcondElems),
                numericInput('nsSamples', 'No. of Sample(s) Per Secondary Covariate', 
                    if (maxbatchElems>defaultDisp) defaultDisp 
                    else maxbatchElems, min = 1, max = maxbatchElems),
                checkboxInput("sortbysecondary", 
                    "Sort By Secondary Covariate First (Default: Sort By Primary Covariate First)", 
                    FALSE),
                checkboxInput("colbysecondary", 
                    "Color By Secondary Covariate (Default: Color By Primary Covariate)", FALSE),
                numericInput('noTaxons', 
                    'No. of top Differentially Expressed Taxa to display', 
                    if (maxGenes>defaultGenesDisp) defaultGenesDisp 
                    else maxGenes, min = 1, max = maxGenes),
                width=3
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Expression Plots",ggvisOutput("DiffExPlot")), 
                    tabPanel("Summary", verbatimTextOutput("DEsummary")),
                    tabPanel("Table", tableOutput("DEtable")), 
                    tabPanel("LIMMA",tableOutput("LimmaTable"))
                ), width=9
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
    )
)
)
