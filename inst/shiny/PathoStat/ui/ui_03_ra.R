
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
            selectizeInput('select_conditions', 'Order Samples by Conditons:', choices=covariates, multiple=TRUE),

            # Select taxon level
            selectInput("taxl", "Tax Level:", choices=tax.name, selected="family"),
          width=3
        ),
        mainPanel(
          tabPanel("Heatmap",
                   helpText("Note: Only variables with less than 8 levels could be mapped to color bar."),
                   fluidRow(
                       column(3, selectInput("select_heatmap_condition_1", "Add colorbar based on:",
                                             covariates.colorbar)),
                       column(3, selectInput("select_heatmap_condition_2", "Add second colorbar based on:",
                                                                              covariates.colorbar)),
                       column(3, checkboxInput("checkbox_heatmap", "Add colorbar", value = TRUE))
                   ),
                   selectInput("ssv_format_new", "Select data format", c("relative abundance", "log10 CPM")),
                   uiOutput("single_species_ui_new"),
                   actionButton("boxplotButtonNew", "Plot"),
                   plotOutput("Heatmap", height="550px"),
                   downloadButton('download_heatmap_pdf', 'Download heatmap PDF')),
          width=9
        )
      )
    )
  )
)
