shiny_panel_filter <- fluidPage(
         sidebarLayout(
             sidebarPanel(
                 selectInput("select_condition.dsf", "Select Condition:",
                             covariates),
                 selectizeInput('taxlTable', 'Taxonomy Level', choices = tax.name, 
                                selected='no rank'),
                 width=3
             ),
             mainPanel(
                 tabsetPanel(
                     tabPanel("Data summary",
                              helpText("We have nothing here now! :'( ")
                     ),
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
)