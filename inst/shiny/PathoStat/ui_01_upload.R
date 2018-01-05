shiny_panel_upload <- fluidPage(

         useShinyjs(),
         tags$style(appCSS),        
         
         tags$div(
             class = "jumbotron",
             tags$div(
                 class = "container",
                 fluidRow(
                     column(7, h1("PathoStat")),
                     column(2, img(src = "bu_logo.png", height = 60, width = 100)),
                     column(2, img(src = "cbm_logo.jpeg", height = 60, width = 100)),
                     column(1, img(src = "bu_bioinfo_logo.png", height = 60, width = 60))
                     
                 ),
                 
                 p("Statistical Microbiome Analysis on metagenomics")
             )
         ),
         sidebarLayout(
             sidebarPanel(
                 radioButtons("uploadChoice", "Upload:",
                              c("Example data" = "example",
                                "Count File" = "files",
                                "PathoScope Files" = "patho.files"
                              )),
                 br(),
                 p(" "),
                 conditionalPanel(condition = sprintf("input['%s'] == 'files'", "uploadChoice"),
                                  fileInput("countsfile", "Counts .csv file (required):",
                                            accept = c(
                                                "text/csv",
                                                "text/comma-separated-values",
                                                "text/tab-separated-values",
                                                "text/plain",
                                                ".csv",
                                                ".tsv"
                                            )
                                  ),
                                  fileInput("annotfile.count", "Annotation .csv file (required):",
                                            accept = c(
                                                "text/csv",
                                                "text/comma-separated-values",
                                                "text/tab-separated-values",
                                                "text/plain",
                                                ".csv",
                                                ".tsv"
                                            )
                                  ),
                                  # Input: Checkbox if file has header ----
                                  checkboxInput("header.count", "Header", TRUE),
                                  
                                  # Input: Select separator ----
                                  radioButtons("sep.count", "Separator",
                                               choices = c(Tab = "\t",
                                                           Comma = ",",
                                                           Semicolon = ";"
                                               ),
                                               selected = ","),
                                  withBusyIndicatorUI(
                                      actionButton("uploadDataCount", 
                                                   "Upload",
                                                   class = "btn-primary")
                                  )
                 ),
                 conditionalPanel(condition = sprintf("input['%s'] == 'patho.files'", "uploadChoice"),
                                  h5("Upload PathoScope generated .tsv files:"),
                                  fileInput("countsfile.pathoscope", "PathoScope outputs (required):",
                                            multiple = TRUE,
                                            accept = c(
                                                "text/csv",
                                                "text/comma-separated-values",
                                                "text/tab-separated-values",
                                                "text/plain",
                                                ".csv",
                                                ".tsv"
                                            )
                                  ),
                                  fileInput("annotfile.ps", "Annotation .tsv file (required):",
                                            accept = c(
                                                "text/csv",
                                                "text/comma-separated-values",
                                                "text/tab-separated-values",
                                                "text/plain",
                                                ".csv",
                                                ".tsv"
                                            )
                                  ),
                                  # Input: Checkbox if file has header ----
                                  checkboxInput("header.ps", "Header", TRUE),
                                  
                                  # Input: Select separator ----
                                  radioButtons("sep.ps", "Separator",
                                               choices = c(Tab = "\t",
                                                           Comma = ",",
                                                           Semicolon = ";"
                                               ),
                                               selected = "\t"),
                                  withBusyIndicatorUI(
                                      actionButton("uploadDataPs", 
                                                   "Upload",
                                                   class = "btn-primary")
                                  )
                                  
                 )
             ),
             mainPanel(
                 conditionalPanel(condition = sprintf("input['%s'] != 'example'", "uploadChoice"),
                                  helpText("Counts Table"),
                                  tableOutput("contents.count"),
                                  helpText("Annotation table"),
                                  tableOutput("contents.meta")
                 ),
                 conditionalPanel(condition = sprintf("input['%s'] == 'patho.files'", "uploadChoice"),
                                  h4("Please open this app in Chrome for multiple files upload."),
                                  h5("Also, example pathoscope report files and annotation files 
                                     could be found at pathToPathoStat/inst/example/data/pathoscope_example/")
                                  ),
                 conditionalPanel(condition = sprintf("input['%s'] == 'files'", "uploadChoice"),
                                  h5("Example count file and annotation file 
                                     could be found at pathToPathoStat/inst/example/data/count_example/")
                                  )
                 
                 
                 
                 
                 
                 )
             )
         )