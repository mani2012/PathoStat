shiny_panel_filter <- fluidPage(
    
    tabsetPanel(
        tabPanel("Reads Count & RA",
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
                                      DT::dataTableOutput("TaxCountTable",width='95%')),
                             coreOTUModuleUI("coreOTUModule")
                         ), width=9
                     )
                 )
        ),
        tabPanel("Sample Filter",
                 br(),
                 sidebarLayout(
                     sidebarPanel(
                         br(),
                         selectizeInput('filterSample', 'Choose sample name(s) for removal', choices = sample.names.all, 
                                        selected=NULL, multiple = TRUE),
                         withBusyIndicatorUI(
                           actionButton("filterSampleButton", "Remove")
                         ),
                         width=3
                     ),
                     mainPanel(
                         br(),
                         tabsetPanel(
                             tabPanel("Sample Reads Count Sum",
                                      selectInput("select_condition_sample_filter", "Order by:",
                                                  c("Reads Number", covariates)),
                                      plotlyOutput("sampleCountSum"))
                         ), width=9
                     )
                 )

                 
        )
    )

)