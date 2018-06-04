# tags$style(type="text/css",
#            ".shiny-output-error { visibility: hidden; }",
#            ".shiny-output-error:before { visibility: hidden; }"
# )

shiny_panel_filter <- fluidPage(

    tabsetPanel(

        tabPanel("Sample Filter",
                 br(),
                 sidebarLayout(
                     sidebarPanel(
                         br(),
                         selectInput(
                           "filter_type", "Select filter type",
                           c("By Name", "By Metadata")),

                         conditionalPanel(condition = "input.filter_type == 'By Name'",
                                          selectizeInput('filterSample', 'Choose sample name(s) for removal', choices = sample.names.all,
                                                         selected=NULL, multiple = TRUE),
                                          withBusyIndicatorUI(
                                            actionButton("filterSampleButton", "Filter")
                                          )
                         ),
                         conditionalPanel(condition = "input.filter_type == 'By Metadata'",
                              selectInput("select_condition_sample_filter_sidebar", "Select a variable",
                                                      c("Read Number", covariates)),
                              conditionalPanel(condition = "output.filter_type == 'num.continuous'",
                                               helpText("Please select range of the feature to keep"),
                                               numericInput("num_filter_min", "Min", 0),
                                               numericInput("num_filter_max", "Max", 0),
                                               withBusyIndicatorUI(
                                                 actionButton("filter_num", "Filter")
                                               )
                                               ),
                              conditionalPanel(condition = "output.filter_type == 'cat'",
                                               helpText("Please select levels of the feature to keep"),
                                               uiOutput("filter_cat_options"),
                                               withBusyIndicatorUI(
                                                 actionButton("filter_cat", "Filter")
                                               )
                              )
                         ),
                         br(),
                         withBusyIndicatorUI(
                           actionButton("resetSampleButton", "Reset")
                         ),
                         br(),
                         br(),
                         downloadButton('download_rda', 'Download pathostat data'),
                         width=3
                     ),
                     mainPanel(
                         br(),
                         tabsetPanel(
                             tabPanel("Sample summary",
                                      tableOutput("contents_summary"),
                                      selectInput("select_condition_sample_filter", "Order by:",
                                                  c("Read Number", covariates)),
                                      plotlyOutput("sampleCountSum")),
                             tabPanel("Sample distribution in metadata",
                                      selectInput("select_condition_sample_distribution", "See distribution in:",
                                                  covariates),
                                      plotlyOutput("sample_metadata_distribution"))
                         ), width=9
                     )
                 )


        ),
        tabPanel("Microbes Filter",
                 br(),
                 sidebarLayout(
                     sidebarPanel(
                         br(),
                         selectInput(
                           "filter_type_micro", "Select filter type",
                           c("By mapped read number", 
                             "By relative abundace",
                             "By prevalence")),
                         conditionalPanel(condition = "input.filter_type_micro == 'By mapped read number'",
                                               helpText("Please select the average minimum read mapped"),
                                               numericInput("read_filter_min_micro", "Min", 0, min = 0, max = 10000),
                                               withBusyIndicatorUI(
                                                 actionButton("filter_read_micro", "Filter")
                                               )
                         ),      
                         conditionalPanel(condition = "input.filter_type_micro == 'By relative abundace'",
                                               helpText("Please select average minimum RA"),
                                               numericInput("ra_filter_min_micro", "Min", 0, min = 0, max = 1),
                                               withBusyIndicatorUI(
                                                 actionButton("filter_ra_micro", "Filter")
                                               )
                         ), 
                         conditionalPanel(condition = "input.filter_type_micro == 'By prevalence'",
                                               helpText("Please select minimum prevalence"),
                                               numericInput("prev_filter_min", "Min", 0, min = 0, max = 1),
                                               withBusyIndicatorUI(
                                                 actionButton("filter_prev_micro", "Filter")
                                               )
                         ), 
                         br(),
                         withBusyIndicatorUI(
                           actionButton("resetSampleButtonMicro", "Reset")
                         ),
                         br(),
                         br(),
                         width=3
                     ),
                     mainPanel(
                         br(),
                         tabsetPanel(
                             tabPanel("Sample summary",
                                      tableOutput("contents_summary_micro"),
                                      selectInput("select_condition_sample_filter_micro", "Order by:",
                                                  c("Taxon elements number", covariates)),
                                      plotlyOutput("sampleTaxon"))                
                         
                         ), width=9
                     )
                 )


        ),
        
        tabPanel("Boxplot visualization",
                 selectizeInput('taxl_single_species', 'Taxonomy Level', choices = tax.name,
                                selected='no rank'),
                 selectInput("select_single_species_condition", "Select condition",
                             covariates.colorbar),
                 selectInput("ssv_format", "Select data format", c("read count", "relative abundance", "log10 CPM")),
                 uiOutput("single_species_ui"),
                 actionButton("boxplotButton", "Plot"),
                 plotlyOutput("single_species_boxplot")

                 #plotlyOutput("single_species_barplot")
        ),
        tabPanel("Read Count & RA",
                 br(),
                 sidebarLayout(
                   sidebarPanel(
                     br(),
                     selectizeInput('taxlTable', 'Taxonomy Level', choices = tax.name,
                                    selected='no rank'),
                     width=3
                   ),
                   mainPanel(
                     tabsetPanel(
                       tabPanel("RA Table(%)",
                                br(),
                                downloadButton('downloadData', 'Download RA CSV'),
                                DT::dataTableOutput("TaxRAtable", width='95%')),
                       tabPanel("Count Table",
                                br(),
                                downloadButton('downloadCountData', 'Download Count CSV'),
                                DT::dataTableOutput("TaxCountTable",width='95%'))
                       #coreOTUModuleUI("coreOTUModule")
                     ), width=9
                   )
                 )
        )
    )

)
