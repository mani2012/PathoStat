shiny_panel_filter <- fluidPage(
    
    tabsetPanel(
        tabPanel("Data summary",
                 plotlyOutput("sampleCountSum"),
                 sidebarLayout(
                     sidebarPanel(
                         selectizeInput('taxlTable', 'Taxonomy Level', choices = tax.name, 
                                        selected='no rank'),
                         width=3
                     ),
                     mainPanel(
                         tabsetPanel(
                             tabPanel("RA Table(%)",                 
                                      downloadButton('downloadData', 'Download RA CSV'),
                                      DT::dataTableOutput("TaxRAtable", width='95%')),
                             tabPanel("Count Table", 
                                      downloadButton('downloadCountData', 'Download Count CSV'),
                                      DT::dataTableOutput("TaxCountTable",width='95%')),
                             coreOTUModuleUI("coreOTUModule")
                         ), width=9
                     )
                 )
        ),
        tabPanel("Data filter",
                 helpText("ohho")
        )
    )

)