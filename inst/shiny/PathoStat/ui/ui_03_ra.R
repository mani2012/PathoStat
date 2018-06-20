tabPanel("Relative Abundance",
  tabsetPanel(
    tabPanel("Sample Relative Abundance",
      tags$br(),
      sidebarLayout(
        sidebarPanel(

          # Sort the samples by a condition
          conditionalPanel(
            condition = "input.group_samples == false",
            selectizeInput('sra_select_conditions', 'Color Samples by Condition', choices=covariates, multiple=TRUE)
          ),
          conditionalPanel(
            condition = "input.group_samples == true",
            selectizeInput('gra_select_conditions', 'Color Samples by Condition', choices=covariates)
          ),

          # Sample aggregation
          checkboxInput("group_samples", "Group Samples by Condition"),

          # Select taxon level
          selectInput("sra_taxlev", "Tax Level", choices=tax.name, selected="family"),

          # Dynamically generate based on tax level
          uiOutput("sra_order_organisms"),

          # Sort the bars
          radioButtons("sra_sort_by", "Sort By",
          c("No Sorting" = "nosort", "Conditions" = "conditions", "Organisms" = "organisms"), selected="nosort"),

          # Legend toggle
          checkboxInput("sra_show_legend", "Show Legend", value=TRUE),

          actionButton("plot_sra", "Plot"),
          width=3
        ),
        mainPanel(
          plotlyOutput("ra_plot", width="800px", height="600px"),
          width=9
        )
      )
    ),
    tabPanel("Heatmap",
      tags$br(),
      sidebarLayout(
        sidebarPanel(
          # Sort the samples by a condition
          selectizeInput('hmra_select_conditions', 'Color Samples by Condition', choices=covariates, multiple=TRUE),

          # Select taxon level
          selectInput("hmra_taxlev", "Tax Level", choices=tax.name, selected="family"),

          # Dynamically generate based on tax level
          uiOutput("hmra_isolate_organisms"),

          # Sort the bars
          radioButtons("hmra_sort_by", "Sort By",
          c("No Sorting" = "nosort", "Conditions" = "conditions", "Organisms" = "organisms"), selected="nosort"),

          # Legend toggle
          checkboxInput("hmra_logcpm", "log(CPM)", value=FALSE),

          actionButton("plot_hmra", "Plot"),
          width=3
        ),
        mainPanel(
          plotlyOutput("hmra_plot", width="800px", height="600px"),
          width=9        
        )
      )
    )
  )
)
