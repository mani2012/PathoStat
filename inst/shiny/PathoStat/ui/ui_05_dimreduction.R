tabPanel("Dimension Reduction",
  sidebarLayout(
    sidebarPanel(
      numericInput('xcol.new', 'Principal Component (x-axis)', 1, min = 1, max = 50),
      numericInput('ycol.new', 'Principal Component (y-axis)', 2, min = 1, max = 50),
      selectizeInput('taxl.pca', 'Taxonomy Level', choices = tax.name, selected='no rank'),
      selectInput("select_pca_color", "Color points by:", covariates),
      selectInput("select_pca_shape", "Shape points by:", c("None", covariates.colorbar)),
      actionButton("DR_plot", "Plot"),
      width=3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("PCA",
          fluidRow(
            column(9,
              plotlyOutput("pca.plotly", height="500px"),
              selectInput("select_pca_data_format", "Select data type", c("read count", "log10 CPM", "RA"))
            ),
            column(3,
              DT::dataTableOutput("PCAtable")
            )
          )
        ),
        tabPanel("PCoA",
          fluidRow(
            column(9,
              plotlyOutput("pcoa.plotly", height="500px"),
              selectInput("pcoa.method", "PCoA method:", beta.methods)
            ),
            column(3,
              DT::dataTableOutput("PCoAtable")
            )
          )
        )
      )
    )
  )
)