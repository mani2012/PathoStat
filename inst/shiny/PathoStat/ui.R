library(shiny)
library(shinyjs)
library(ggvis)
library(d3heatmap)
library(phyloseq)
library(plotly)
source("helpers.R")




alpha.methods <- c("Shannon", "Simpson", "InvSimpson")
# Weigthed Unifrac, Bray-Curtis
beta.methods <- c("wUniFrac", "bray")

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




#sample name
sample.names.all <- colnames(pstat@otu_table@.Data)


# tags$style(type="text/css",
#            ".shiny-output-error { visibility: hidden; }",
#            ".shiny-output-error:before { visibility: hidden; }"
# )

# load ui tabs right before calling shinyUI()
source("ui_01_upload.R", local = TRUE) #creates shiny_panel_upload variable
source("ui_02_filter.R", local = TRUE) #creates shiny_panel_upload variable
source("ui_07_biomarker.R", local = TRUE) #creates shiny_panel_upload variable



shinyUI(navbarPage(paste("PathoStat v", packageVersion("PathoStat"), sep = ""), id="PathoStat", fluid=TRUE,

    theme = "bootstrap.min.css",
    tabPanel("Upload", shiny_panel_upload),
    tabPanel("Data Summary/Filtering", shiny_panel_filter),
    tabPanel("Relative Abundance",
        sidebarLayout(
            sidebarPanel(
                tags$br(),

                # Sort the samples by a condition
                selectizeInput('select_conditions', 'Order Samples by Conditons:', choices=covariates, multiple=TRUE),

                # Select taxon level
                selectInput("taxl", "Tax Level:", choices=tax.name, selected="family"),

                # Dynamically generate based on tax level
                uiOutput("order_organisms"),

                radioButtons("sort_samples_by", "Sort Samples By:",
                                 c("No Sorting" = "nosort",
                                   "Conditions" = "conditions",
                                   "Organisms" = "organisms"), selected="nosort"),    

                # Legend toggle
                checkboxInput("show_legend", "Show Legend", value=TRUE),
                width=3
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Taxonomy level RA",
                        plotlyOutput("TaxRelAbundancePlot", width="800px", height="600px")),
                    # new heatmap
                    tabPanel("Heatmap",
                             helpText("Note: Only variables with less than 8 levels could be mapped to color bar."),
                             fluidRow(
                                 column(3, selectInput("select_heatmap_condition_1", "Add colorbar based on:",
                                                       covariates.colorbar)),
                                 column(3, selectInput("select_heatmap_condition_2", "Add second colorbar based on:",
                                                                                        covariates.colorbar)),
                                 column(3, checkboxInput("checkbox_heatmap", "Add colorbar", value = TRUE))
                             ),
                             selectInput("ssv_format_new", "Select data format", c("relative abundance", "log10 CPM")),
                             uiOutput("single_species_ui_new"),
                             actionButton("boxplotButtonNew", "Plot"),
                             plotOutput("Heatmap", height="550px"),
                             downloadButton('download_heatmap_pdf', 'Download heatmap PDF'))
                ), width=9
            )
        )
    ),
    tabPanel("Diversity",
        tabsetPanel(
            tabPanel("Alpha Diversity",
                br(),
                sidebarLayout(
                    sidebarPanel(
                        br(),
                     selectizeInput('taxl.alpha', 'Taxonomy Level', choices = tax.name,
                                     selected='no rank'),
                     selectInput("select_alpha_div_condition", "Compare between:",
                                 covariates.colorbar),
                     selectInput("select_alpha_div_method", "Choose method:",
                                 alpha.methods)
                    ),
                    mainPanel(
                     br(),
                     tabsetPanel(
                       tabPanel("Boxplot",
                                plotlyOutput("AlphaDiversity"),
                                actionButton("download_alpha", "Download Alpha diversity pdf"),
                                helpText("Note: Wait for 8-10s after clicking DOWNLOAD, and the figure will be opened externally.")
                       ),
                       tabPanel("Taxa number Barplot",
                                plotlyOutput("AlphaDiversityBarplot")
                       ),
                       tabPanel("Statistical Test",
                                selectInput("select_alpha_stat_method","Non-parametric Test", c("Mann-Whitney","T-test", "Kruskal-Wallis")),
                                tableOutput("alpha.stat.test")
                       ),
                       tabPanel("Alpha Diversity Table",
                                br(),
                                downloadButton('download_table_alpha', 'Download this table'),
                                DT::dataTableOutput("table.alpha")
                       )
                     )
                    )
                )
          ),
          tabPanel("Beta Diversity",
                   br(),
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
                       br(),
                       tabsetPanel(
                         tabPanel("Heatmap",
                                  br(),
                                  fluidRow(
                                      column(3, selectInput("select_beta_heatmap_condition_1", "Add colorbar on:",
                                                            covariates.colorbar)),
                                      column(3, selectInput("select_beta_heatmap_condition_2", "Add second colorbar on:",
                                                            covariates.colorbar))
                                  ),
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
                                  tableOutput("beta.stat.test")
                         ),
                         tabPanel("Beta Diversity Table",
                                  br(),
                                  downloadButton('download_table_beta', 'Download this table'),
                                  DT::dataTableOutput("table.beta")
                         )
                       )
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
                     actionButton("DR_plot", "Plot"),
                     width=3
                 ),
                 mainPanel(
                     tabsetPanel(
                         tabPanel("PCA plot",
                                  # This is a bit different pdf downloading method for plotly,
                                  # as we must buy lisence for that
                                  plotlyOutput("pca.plotly"),
                                  selectInput("select_pca_data_format", "Select data type", c("read count", "log10 CPM", "RA")),
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
             tabPanel("Differential Abundance",
                      selectInput("DAmethod", "Select method",
                        c("DESeq2", "edgeR")),
                      conditionalPanel(condition = "input.DAmethod == 'DESeq2'",
                                       sidebarLayout(
                                         sidebarPanel(
                                           selectizeInput('taxl.da', 'Taxonomy Level', choices = tax.name,
                                                          selected='no rank'),
                                           selectizeInput('da.condition', 'Select condition',
                                                          choices = covariates.colorbar),
                                           conditionalPanel(condition = "output.da_condition_type == 'multiple'",
                                                            helpText("Please select 2 levels to compare"),
                                                            uiOutput("da_condition_options")
                                           ),
                                           selectizeInput('da.condition.covariate', 'Select (multiple) covariates',
                                                          choices = covariates, multiple = TRUE),
                                           helpText("Continuous covariates would be automatically cut into factors with 3 levels."),
                                           numericInput('da.count.cutoff', 'Minumum count cut-off', 500,
                                                        min = 1, max = 5000),
                                           numericInput('da.padj.cutoff', 'Choose padj cut-off', 0.5,
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
                      conditionalPanel(condition = "input.DAmethod == 'edgeR'",
                                       sidebarLayout(
                                         sidebarPanel(
                                           selectizeInput('taxl.edger', 'Taxonomy Level', choices = tax.name,
                                                          selected='no rank'),
                                           selectizeInput('edger.condition', 'Select condition',
                                                          choices = covariates.colorbar),
                                           conditionalPanel(condition = "output.edger_condition_type == 'multiple'",
                                                            helpText("Please select 2 levels to compare"),
                                                            uiOutput("edger_condition_options")
                                           ),
                                           helpText("Continuous covariates would be automatically cut into factors with 3 levels."),
                                           numericInput('edger.padj.cutoff', 'Choose padj cut-off', 0.5,
                                                        min = 1e-100, max = 1),
                                           width=3
                                         ),
                                         mainPanel(
                                           tabPanel("edgeR",
                                                    tabsetPanel(
                                                      tabPanel("DE output",
                                                               DT::dataTableOutput("edgerTable.new"),
                                                               downloadButton("download_edger_tb", "Download this table")
                                                      )
                                                    )
                                           ), width=9
                                         )
                                       )
                      )
             ),

             tabPanel("Statistical Test",
             sidebarLayout(
               sidebarPanel(
                 selectizeInput('taxl.pa', 'Taxonomy Level', choices = tax.name,
                                selected='no rank'),
                 selectizeInput('pa.condition', 'Select condition',
                                choices = covariates.colorbar),
                 conditionalPanel(condition = "output.pa_condition_type == 'multiple'",
                                  helpText("Please select 2 levels to compare"),
                                  uiOutput("pa_condition_options")
                 ),
                 numericInput('pa.count.cutoff', 'Minumum count cut-off', 500,
                              min = 1, max = 5000),
                 numericInput('pa.padj.cutoff', 'Choose padj cut-off', 0.5,
                              min = 1e-100, max = 1),
                 width=3
               ),
               mainPanel(
                 tabPanel("Test output",
                          br(),
                          tabsetPanel(
                            tabPanel("output",
                                     selectInput('pa_method', 'Select test method',
                                                    choices = c("Fisher Exact Test", "Chi-squared Test", "Mann-Whitney Test")),
                                     conditionalPanel(condition = "input.pa_method == 'Mann-Whitney Test'",
                                                      selectInput('pa_mann_data_type', 'Select data type', choices = c("read count", "log10 CPM", "RA"))),
                                     DT::dataTableOutput("pa.test"),
                                     downloadButton("download_pa_test", "Download this table")
                            )
                          )
                 ), width=9
               )
             )
        )
        )


    ),
    tabPanel("Biomarker", shiny_panel_biomarker)
    #tabPanel("Pathway (Under construction by Tyler)")
    # tabPanel("Time Series",
    #     tabsetPanel(
    #         tabPanel("Visualization",
    #             sidebarLayout(
    #                 sidebarPanel(
    #                   selectInput(inputId="Allusset",
    #                       label="Visualization column",
    #                       choices = colnames(shinyInput$pstat@sam_data)),
    #                   checkboxInput(inputId="Allurar",
    #                       label="Rarefaction? (maximum reads of minimal
    #                       sample count)"),
    #                   selectInput(inputId="Alluglom", label="Agglomerate taxa",
    #                       choices = colnames(shinyInput$pstat@tax_table)),
    #                   uiOutput("Allustax"),
    #                   downloadButton('downloadAlluvialPlot',
    #                                  'Download Plot')
    #                 ),
    #                 mainPanel(
    #                   plotOutput("TimePlotVisu",height = "600px")
    #                 )
    #             )
    #         )
    #     )
    # )
)
)
