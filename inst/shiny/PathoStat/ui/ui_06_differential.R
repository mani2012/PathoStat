

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
                                       actionButton("run_deseq2", "Run"),
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
                                       actionButton("run_edgeR", "Run"),
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
             actionButton("run_stat", "Run"),
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


)
