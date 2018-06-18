library(shiny)

source(file.path("utils", "old_server_stuff.R"),  local = TRUE)
source(file.path("utils", "old_ui_stuff.R"),  local = TRUE)

ui <- navbarPage(
  title = paste("PathoStat v", packageVersion("PathoStat"), sep = ""),
  id="PathoStat",
  fluid=TRUE,
  theme = "bootstrap.min.css",
  source(file.path("ui", "ui_01_upload.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_02_filter.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_03_ra.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_04_diversity.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_05_dimreduction.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_06_differential.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_07_biomarker.R"),  local = TRUE)$value
)

server <- function(input, output, session) {

  source(file.path("server", "server_01_upload.R"),  local = TRUE)$value
  source(file.path("server", "server_02_filter.R"),  local = TRUE)$value
  source(file.path("server", "server_03_ra.R"),  local = TRUE)$value
  source(file.path("server", "server_04_diversity.R"),  local = TRUE)$value
  source(file.path("server", "server_05_dimreduction.R"),  local = TRUE)$value
  source(file.path("server", "server_06_differential.R"),  local = TRUE)$value
  source(file.path("server", "server_07_biomarker.R"),  local = TRUE)$value
}

shinyApp(ui = ui, server = server)
