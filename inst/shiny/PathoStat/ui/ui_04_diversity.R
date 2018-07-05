

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
                             alpha.methods),
                 actionButton("alpha_boxplot", "Run")
                ),
                mainPanel(
                 br(),
                 tabsetPanel(
                   tabPanel("Boxplot",
                            plotlyOutput("AlphaDiversity")
                    
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
                               covariates.two.levels),
                   actionButton("beta_boxplot", "Run")
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
                              plotOutput("BetaDiversityHeatmap")
                     ),
                     tabPanel("Boxplot",
                              plotlyOutput("BetaDiversityBoxplot")
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
)
