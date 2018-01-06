library(shiny)
library(shinyjs)
library(ggvis)
library(d3heatmap)
library(phyloseq)
library(ape)
library(plotly)
source("helpers.R")




alpha.methods <- c("Shannon", "Simpson", "InvSimpson")
# Weigthed Unifrac, Bray-Curtis
beta.methods <- c("bray", "wUniFrac")

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

# choose the covariates that has less than 8 levels
covariates.colorbar <- c()
for (i in 1:length(covariates)){
  num.levels <- length(unique(sample_data(pstat)[[covariates[i]]]))
  if (num.levels < 8){
    covariates.colorbar <- c(covariates.colorbar, covariates[i])
  }
}

# choose the covariates that has 2 levels
covariates.two.levels <- c()
for (i in 1:length(covariates)){
  num.levels <- length(unique(sample_data(pstat)[[covariates[i]]]))
  if (num.levels == 2){
    covariates.two.levels <- c(covariates.two.levels, covariates[i])
  }
}


maxbatchElems <- minbatch(c(pstat@sam_data[,1])[[1]])
maxcondElems <- minbatch(c(pstat@sam_data[,2])[[1]])
defaultDisp <- 30
defaultGenesDisp <- 10
maxGenes <- dim(pstat@otu_table)[1]


# load ui tabs right before calling shinyUI()
source("ui_01_upload.R", local = TRUE) #creates shiny_panel_upload variable
source("ui_02_filter.R", local = TRUE) #creates shiny_panel_upload variable




shinyUI(navbarPage(paste("PathoStat v", packageVersion("PathoStat"), sep = ""), id="PathoStat", fluid=TRUE, 
    
    theme = "bootstrap.min.css",
    tabPanel("Upload", shiny_panel_upload),
    tabPanel("Data Summary/Filtering", shiny_panel_filter),
    tabPanel("Relative Abundance",
        sidebarLayout(
            sidebarPanel(
                #selectizeInput('taxl', 'Taxonomy Level', choices = setNames(
                #   tax.abb, tax.name)),
                selectizeInput('taxl', 'Taxonomy Level', choices = tax.name, 
                    selected='no rank'),
                br(),
                p(" "),
                selectInput("select_condition", "Select Condition:",
                            covariates),
                width=3
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Taxonomy level RA",
                        ggvisOutput("TaxRelAbundancePlot")),
                    # new heatmap
                    tabPanel("Heatmap", 
                             helpText("Note: Only variables with less than 8 levels could be mapped to color bar."),
                             selectInput("select_heatmap_condition_1", "Add colorbar based on:",
                                         covariates.colorbar),
                             selectInput("select_heatmap_condition_2", "Add second colorbar based on:",
                                         covariates.colorbar),
                             checkboxInput("checkbox_heatmap_scale", "Row scaling", value = TRUE),
                             checkboxInput("checkbox_heatmap", "Add colorbar", value = TRUE),
                             downloadButton('download_heatmap_pdf', 'Download heatmap PDF'),
                             plotOutput("Heatmap", height="550px"))
                ), width=9
            )
        )
    ),
    tabPanel("Diversity",
        tabsetPanel(
            tabPanel("Alpha Diversity", 
                     
                sidebarLayout(
                    sidebarPanel(
                     selectizeInput('taxl.alpha', 'Taxonomy Level', choices = tax.name, 
                                     selected='no rank'),
                     selectInput("select_alpha_div_condition", "Compare between:",
                                 covariates.colorbar),
                     selectInput("select_alpha_div_method", "Choose method:",
                                 alpha.methods)
                    ),
                    mainPanel(
                     
                     tabsetPanel(
                       tabPanel("Boxplot", 
                                plotlyOutput("AlphaDiversity"),
                                actionButton("download_alpha", "Download Alpha diversity pdf"),
                                helpText("Note: Wait for 8-10s after clicking DOWNLOAD, and the figure will be opened externally.")
                       ),
                       tabPanel("Statistical Test", 
                                selectInput("select_alpha_stat_method","Non-parametric Test", c("Mann-Whitney","Kruskal-Wallis")),
                                verbatimTextOutput("alpha.stat.test")
                       ),
                       tabPanel("Alpha Diversity Table", 
                                downloadButton('download_table_alpha', 'Download this table'),
                                DT::dataTableOutput("table.alpha")
                       )
                     )
                    )
                )
          ),
          tabPanel("Beta Diversity", 
                   
                   sidebarLayout(
                     sidebarPanel(
                       selectizeInput('taxl.beta', 'Taxonomy Level', choices = tax.name, 
                                      selected='no rank'),
                       selectInput("select_beta_div_method", "Choose method:",
                                   beta.methods),
                       helpText("Only variables with 2 levels are supported for boxplot and stat test here." ),
                       selectInput("select_beta_condition", "Select condition",
                                   covariates.two.levels)
                     ),
                     mainPanel(
                       
                       tabsetPanel(
                         tabPanel("Heatmap", 
                                  selectInput("select_beta_heatmap_condition_1", "Add colorbar on:",
                                              covariates.colorbar),
                                  selectInput("select_beta_heatmap_condition_2", "Add second colorbar on:",
                                              covariates.colorbar),
                                  checkboxInput("checkbox_beta_heatmap", "Add colorbar", value = TRUE),
                                  plotOutput("BetaDiversityHeatmap"),
                                  downloadButton('download_beta_heatmap_pdf', 'Download heatmap PDF')
                         ),
                         tabPanel("Boxplot", 
                                  plotlyOutput("BetaDiversityBoxplot"),
                                  actionButton("download_beta_boxplot", "Download pdf"),
                                  helpText("Note: Wait for 8-10s after clicking DOWNLOAD, and the figure will be opened externally.")
                         ),
                         tabPanel("Statistical Test", 
                                  selectInput("select_beta_stat_method","Select Test", c("PERMANOVA", "Kruskal-Wallis", "Mann-Whitney")),
                                  numericInput("num.permutation.permanova", "Number of permutations", value = 999, max = 2000),
                                  verbatimTextOutput("beta.stat.test")
                         ),
                         tabPanel("Beta Diversity Table", 
                                  downloadButton('download_table_beta', 'Download this table'),
                                  DT::dataTableOutput("table.beta")
                         )
                       )
                     )
                   )
          ),
          
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
                            "NMDS", "PCoA"), selected='NMDS'),
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
    
    tabPanel("Dimension Reduction",
             sidebarLayout(
                 sidebarPanel(
                     numericInput('xcol.new', 'Principal Component (x-axis)', 1,
                                  min = 1, max = 50),
                     numericInput('ycol.new', 'Principal Component (y-axis)', 2,
                                  min = 1, max = 50),
                     selectizeInput('taxl.pca', 'Taxonomy Level', choices = tax.name, 
                                    selected='no rank'),
                     selectInput("select_pca_color", "Color points by:",
                                 covariates),
                     selectInput("select_pca_shape", "Shape points by:",
                                 covariates.colorbar),
                     width=3
                 ),
                 mainPanel(
                     tabsetPanel(
                         tabPanel("PCA plot", 
                                  # This is a bit different pdf downloading method for plotly,
                                  # as we must buy lisence for that
                                  plotlyOutput("pca.plotly"),
                                  actionButton("download_pca", "Download PCA pdf"),
                                  helpText("Note: Wait for 8-10s after clicking DOWNLOAD, and the figure will be opened externally.")),
                         tabPanel("PCA variance", DT::dataTableOutput("PCAtable")),
                         tabPanel("PCoA plot", 
                                  plotlyOutput("pcoa.plotly"),
                                  selectInput("pcoa.method", "PCoA method:",
                                              beta.methods),
                                  actionButton("download_pcoa", "Download PCoA pdf"),
                                  helpText("Note: Wait for 8-10s after clicking DOWNLOAD, and the figure will be opened externally.")),
                         tabPanel("PCoA variance", DT::dataTableOutput("PCoAtable"))
                     )
                 )
             )
    ),
    tabPanel("Differential Analysis",
        tabsetPanel(
             tabPanel("Deseq2", 
                      
                      sidebarLayout(
                        sidebarPanel(
                          selectizeInput('taxl.da', 'Taxonomy Level', choices = tax.name, 
                                         selected='no rank'),
                          selectizeInput('da.condition', 'Select condition', 
                                         choices = covariates.two.levels),
                          selectizeInput('da.condition.covariate', 'Select (multiple) covariates', 
                                         choices = covariates, multiple = TRUE),
                          helpText("Continuous covariates would be automatically cut into factors with 3 levels."),
                          numericInput('da.count.cutoff', 'Minumum count cut-off', 500,
                                       min = 1, max = 5000),
                          numericInput('da.padj.cutoff', 'Choose padj cut-off', 0.05,
                                       min = 1e-100, max = 1),
                          width=3
                        ),
                        mainPanel(
                            tabPanel("DeSeq2", 
                                     tabsetPanel(
                                       tabPanel("DE output",
                                                DT::dataTableOutput("DeSeq2Table.new"),
                                                downloadButton("download_deseq_tb", "Download this table")
                                       )
                                     )
                            ), width=9
                        )
                      )
             ),
             tabPanel("Statistical Test (presence-absence or count based)",
             sidebarLayout(
               sidebarPanel(
                 selectizeInput('taxl.pa', 'Taxonomy Level', choices = tax.name, 
                                selected='no rank'),
                 selectizeInput('pa.condition', 'Select condition', 
                                choices = covariates.two.levels),
                 numericInput('pa.count.cutoff', 'Minumum count cut-off', 500,
                              min = 1, max = 5000),
                 numericInput('pa.padj.cutoff', 'Choose padj cut-off', 0.05,
                              min = 1e-100, max = 1),
                 width=3
               ),
               mainPanel(
                 tabPanel("Test output", 
                          tabsetPanel(
                            tabPanel("output",
                                     selectizeInput('pa.method', 'Select test method', 
                                                    choices = c("Fisher Exact Test", "Chi-squared Test", "Mann-Whitney Test")),
                                     DT::dataTableOutput("pa.test"),
                                     downloadButton("download_pa_test", "Download this table")
                            )
                          )
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
    )
)
)
