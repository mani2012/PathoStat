shiny_panel_biomarker <- fluidPage(

    tabsetPanel(
        tabPanel("Classification Model",
                 br(),
                 sidebarLayout(
                     sidebarPanel(
                         br(),
                         selectizeInput('taxl_biomarker', 'Taxonomy Level', choices = tax.name,
                                        selected='genus'),
                         selectInput("select_target_condition_biomarker", "Select Target Condition:",
                                     covariates.colorbar),
                         conditionalPanel(condition = "output.biomarker_condition_type == 'multiple'",
                                          helpText("Please select 2 levels to compare"),
                                          uiOutput("biomarker_condition_options")
                         ),
                         numericInput("num.cv.nfolds", "Number of CV nfolds", value = 3, max = 20, min = 3),
                         numericInput("num.biomarker.run", "Total runs", value = 100, max = 500, min = 50),
                         selectInput("select_covariate_condition_biomarker", "Select Covarites Conditions:",
                                     covariates, multiple = TRUE),
                         selectInput("select_model_biomarker", "Select Model", c("Lasso Logistic Regression")),
                         actionButton("goButtonBiomarker", "Run!"),
                         width=3
                     ),
                     mainPanel(
                         tabsetPanel(
                             tabPanel("Feature selection",
                                      br(),
                            tableOutput("featureSelectionTmp"))
                         ), width=9
                     )
                 )
        )
    )

)
