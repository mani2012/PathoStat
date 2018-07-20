tabPanel("Diversity",
  tabsetPanel(
    tabPanel("Alpha Diversity",
      br(),
      sidebarLayout(
        sidebarPanel(
          br(),
          selectizeInput('taxl.alpha', 'Taxonomy Level', choices = tax.name, selected='no rank'),
          selectInput("select_alpha_div_condition", "Compare between:", covariates.colorbar),
          selectInput("select_alpha_div_method", "Choose method:", alpha.methods),
          selectInput("select_alpha_stat_method","Non-parametric Test", c("Mann-Whitney","T-test", "Kruskal-Wallis")),
          actionButton("alpha_boxplot", "Run")
        ),
        mainPanel(
          br(),
          tabsetPanel(
            tabPanel("Boxplot",
              plotlyOutput("AlphaDiversity"),
              br(),
              DT::dataTableOutput("alpha.stat.test")
            ),
            tabPanel("Alpha Diversity Table",
              br(),
              DT::dataTableOutput("table.alpha"),
              downloadButton('download_table_alpha', 'Download')
            )
          )
        )
      )
    ),
    tabPanel("Beta Diversity",
      br(),
      sidebarLayout(
        sidebarPanel(
          selectizeInput('taxl.beta', 'Taxonomy Level', choices = tax.name, selected='no rank'),
          selectInput("select_beta_div_method", "Choose method:", beta.methods),
          helpText("Only variables with 2 levels are supported for boxplot and stat test here." ),
          selectInput("select_beta_condition", "Select condition", covariates.two.levels),
          selectInput("select_beta_stat_method","Select Test", c("PERMANOVA", "Kruskal-Wallis", "Mann-Whitney")),
          numericInput("num.permutation.permanova", "Number of permutations", value = 999, max = 2000),
          actionButton("beta_boxplot", "Run")
        ),
        mainPanel(
          br(),
          tabsetPanel(
            tabPanel("Heatmap",
              br(),
              fluidRow(
                selectizeInput('bdhm_select_conditions', 'Color Samples by Condition', choices=covariates.colorbar, multiple=TRUE),
                radioButtons("bdhm_sort_by", "Sort By", c("No Sorting" = "nosort", "Conditions" = "conditions"), selected="nosort")),
                plotlyOutput("BetaDiversityHeatmap")
            ),
            tabPanel("Boxplot",
              plotlyOutput("BetaDiversityBoxplot"),
              DT::dataTableOutput("beta.stat.test")
            ),
            tabPanel("Beta Diversity Table",
              br(),
              DT::dataTableOutput("table.beta"),
              downloadButton('download_table_beta', 'Download')
            )
          )
        )
      )
    )
  )
)