
tabPanel("Dimension Reduction",
         sidebarLayout(
             sidebarPanel(
                 numericInput('xcol.new', 'Principal Component (x-axis)', 1,
                              min = 1, max = 50),
                 numericInput('ycol.new', 'Principal Component (y-axis)', 2,
                              min = 1, max = 50),
                 selectizeInput('taxl.pca', 'Taxonomy Level', choices = tax.name,
                                selected='no rank'),
                 selectInput("select_pca_color", "Color points by:",
                             covariates),
                 selectInput("select_pca_shape", "Shape points by:",
                             covariates.colorbar),
                 actionButton("DR_plot", "Plot"),
                 width=3
             ),
             mainPanel(
                 tabsetPanel(
                     tabPanel("PCA plot",
                              plotlyOutput("pca.plotly"),
                              selectInput("select_pca_data_format", "Select data type", c("read count", "log10 CPM", "RA"))),
                     tabPanel("PCA variance", DT::dataTableOutput("PCAtable")),
                     tabPanel("PCoA plot",
                              plotlyOutput("pcoa.plotly"),
                              selectInput("pcoa.method", "PCoA method:",
                                          beta.methods)),
                     tabPanel("PCoA variance", DT::dataTableOutput("PCoAtable"))
                 )
             )
         )
)
