is.categorical <- function(v) {
  if (class(v) == "integer" || class(v) == "numeric") {
    return(F)
  } else {
    return(T)
  }
}
sam_temp <- as.data.frame(pstat@sam_data)
num_select <- lapply(covariates, function(x) is.categorical(unlist(sam_temp[,x])))
num_covariates <- covariates[!unlist(num_select)]

tabPanel("Summary and Filter",
  tabsetPanel(
    tabPanel("Sample Filter",
      br(),
      sidebarLayout(
        sidebarPanel(
          br(),
          selectInput("filter_type", "Select filter type", c("By Name", "By Metadata")),
          
          conditionalPanel(condition = "input.filter_type == 'By Name'",
            selectizeInput('filterSample', 'Choose sample name(s) for removal', 
                           choices = sample.names.all, selected=NULL, multiple = TRUE),
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
          downloadButton('download_rds', 'Download .rds file'),
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
              selectInput("select_condition_sample_distribution", "See distribution in:", covariates),
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
          selectInput("filter_type_micro", "Select filter type",
               c("By mapped read number",
                 "By relative abundace",
                 "By prevalence")
          ),
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
            selectInput("select_condition_sample_filter_micro", "Order by:", c("Taxon elements number", covariates)),
            plotlyOutput("sampleTaxon"))
          ), width=9
        )
      )
    ),
    tabPanel("Categorize",
      tags$br(),
      sidebarLayout(
        sidebarPanel(
          selectizeInput('bin_cov', 'Covariate', choices=num_covariates, multiple=FALSE),
          uiOutput("nbins"),
          textInput('bin_breaks', 'Custom Breaks (Comma Delimited)'),
          verbatimTextOutput("bin_to1"),
          textInput('bin_labels', 'Custom Labels (Comma Delimited)'),
          verbatimTextOutput("bin_to2"),

          textInput("new_covariate", "Covariate Label", value = "new_cov"),
          actionButton("create_bins", "Create Bins"),
          width=5
        ),
        mainPanel(
          plotlyOutput("unbin_plot", height="150px"),
          plotlyOutput("bin_plot"),
          width=7
        )
      )
    ),    
    tabPanel("Read Counts & Relative Abundance",
      br(),
      sidebarLayout(
        sidebarPanel(
        br(),
        selectizeInput('taxlTable', 'Taxonomy Level', choices = tax.name, selected='no rank'),
         width=3
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Relative Abundance",
              br(),
              DT::dataTableOutput("TaxRAtable", width='95%'),
              downloadButton('downloadData', 'Download (.csv)')
            ),
            tabPanel("Read Counts",
              br(),
              DT::dataTableOutput("TaxCountTable",width='95%'),
              downloadButton('downloadCountData', 'Download (.csv)')
            )
            #coreOTUModuleUI("coreOTUModule")
          ), width=9
        )
      )
    )
  )
)