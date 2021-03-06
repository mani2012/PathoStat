source(file.path("utils", "old_server_stuff.R"),  local = TRUE)

observeEvent(input$filterSampleButton,{
  withBusyIndicatorServer("filterSampleButton", {
    samples.remove <- input$filterSample
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    samples.remove.index <- which(colnames(pstat@otu_table@.Data) %in% samples.remove)
    pstat@otu_table@.Data <- pstat@otu_table@.Data[,-samples.remove.index]
    pstat@sam_data@.Data <-  lapply(pstat@sam_data@.Data, function(x) {x <- x[-samples.remove.index]})
    pstat@sam_data@row.names <- pstat@sam_data@row.names[-samples.remove.index]
    shinyInput <- list(pstat = pstat)
    vals$shiny.input <- shinyInput
    #updateCovariate()
    updateSample()
  })
})

observeEvent(input$reset_button,{
  withBusyIndicatorServer("reset_button", {
    vals$shiny.input <- vals$shiny.input.backup
    updateSample()
    vals$shiny.input <- vals$shiny.input.backup
    updateTaxLevel()
  })
})

output$download_rds <- downloadHandler(filename = function() {
  paste("pathostat", Sys.Date(), ".rds", sep="")
}, content = function(file) {
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat
  saveRDS(pstat, file=file)
})

output$summary_plot_top <- renderPlotly({
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat 
  tables <- pstat.extraction(pstat)
  OTU_TABLE <- tables$OTU
  SAM_DATA <- tables$SAM
  SAM_DATA[,"Reads"] = colSums(OTU_TABLE[,rownames(SAM_DATA)])
  SAM_DATA[,"Taxnum"] = apply(OTU_TABLE , 2, function(x) sum(x >= 1) )

  if (input$filter_type == "By Microbes") {
    cov <- "Taxnum"
  } else {
    cov <- input$select_condition_sample_filter_sidebar
  }
  
  if (length(unique(unlist(SAM_DATA[,cov]))) > 8) {
    fit <- density(unlist(SAM_DATA[,cov]))
    num.density <- plot_ly(x=fit$x, y=fit$y, 
                           type="scatter", 
                           mode="lines", 
                           fill="tozeroy") %>%
                    layout(title=cov)
    num.density$elementId <- NULL
    return(num.density)
  } else {
    cat.df = data.frame(table(SAM_DATA[,cov]))
    cat.bar <- plot_ly(x = cat.df$Var1,
                       y = cat.df$Freq,
                       type = "bar",
                       showlegend = FALSE) %>%
                layout(title=cov,
                       xaxis = list(title = 'Diet',
                                    tickmode = "array",
                                    showticklabels = TRUE,
                                    categoryorder = 'trace'),
                       yaxis = list(title = 'Frequency'))
    cat.bar$elementId <- NULL
    return(cat.bar)
  }
})

output$summary_plot_bottom <- renderPlotly({
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat
  tables <- pstat.extraction(pstat)
  OTU_TABLE <- tables$OTU
  SAM_DATA <- tables$SAM
  SAM_DATA[,"Reads"] = colSums(OTU_TABLE[,rownames(SAM_DATA)])
  SAM_DATA[,"Taxnum"] = apply(OTU_TABLE , 2, function(x) sum(x >= 1) )

  if (input$filter_type == "By Microbes") {
    cov <- "Taxnum"
  } else {
    cov <- input$select_condition_sample_filter_sidebar
  }

  if (length(unique(unlist(SAM_DATA[,cov]))) > 8) {
    vec <- unlist(SAM_DATA[,cov])
    num.scatter <- plot_ly(y = vec, 
                          jitter = 0.3, 
                          pointpos = -1.8, 
                          boxpoints = 'all',
                          marker = list(color = 'rgb(7,40,89)'),
                          line = list(color = 'rgb(7,40,89)'),
                          name = cov, type="box") %>%
                   layout(title=cov)
    num.scatter$elementId <- NULL
    return(num.scatter)
  } else {
    cat.df = data.frame(table(SAM_DATA[,cov]))
    cat.pie <- plot_ly(cat.df, labels = ~Var1, values = ~Freq,
                       type = "pie",
                       showlegend = FALSE) %>%
                layout(title=cov)
    cat.pie$elementId <- NULL
    return(cat.pie)
  }
})

output$filter_type <- reactive({
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat
  if (input$select_condition_sample_filter_sidebar == "Reads") {
    return("num.continuous")
  }
  variable.vec <- sample_data(pstat)[[
        which(pstat@sam_data@names == input$select_condition_sample_filter_sidebar)]]
  if (is.numeric(variable.vec)){
    # catogorical or not
    num.levels <- length(unique(variable.vec))
    if (num.levels < 6) {
      return("cat")
    } else {
      return("num.continuous")
    }
  } else {
    # non-numeric
    return("cat")
  }
})

outputOptions(output, "filter_type", suspendWhenHidden = FALSE)

observeEvent(input$filter_num,{
  withBusyIndicatorServer("filter_num", {
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    # extract the feature vector

    if (input$select_condition_sample_filter_sidebar == "Reads"){
      feature.selected.vec <- colSums(pstat@otu_table@.Data)
    }else{
      feature.selected.vec <- pstat@sam_data@.Data[[
        which(pstat@sam_data@names == input$select_condition_sample_filter_sidebar)]]
    }

    # get the index of sample keeping
    samples.keep.index <- which(feature.selected.vec >= input$num_filter_min &
                                  feature.selected.vec <= input$num_filter_max)

    pstat@otu_table@.Data <- pstat@otu_table@.Data[,samples.keep.index]
    pstat@sam_data@.Data <- lapply(pstat@sam_data@.Data, function(x) {x <- x[samples.keep.index]})
    pstat@sam_data@row.names <- pstat@sam_data@row.names[samples.keep.index]
    shinyInput <- list(pstat = pstat)
    vals$shiny.input <- shinyInput

    #updateCovariate()
    updateSample()
  })
})

output$filter_cat_options <- renderUI({
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat
  if (input$select_condition_sample_filter_sidebar == "Reads"){
    tagList(
      selectInput("cat_filter_options", "Keep these levels:", choices = c(""))
    )
  } else {
    variable.filter.index <- which(pstat@sam_data@names ==
    input$select_condition_sample_filter_sidebar)
    #print(variable.filter.index)
    variable.vec <- pstat@sam_data@.Data[[variable.filter.index]]
    filter.option.vec <- sort(unique(variable.vec))
    tagList(
      selectInput("cat_filter_options", "Keep these levels:", choices = filter.option.vec, multiple = TRUE)
    )
  }
})

observeEvent(input$filter_cat,{
  withBusyIndicatorServer("filter_cat", {
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    # extract the feature vector
    feature.selected.vec <- pstat@sam_data@.Data[[
      which(pstat@sam_data@names == input$select_condition_sample_filter_sidebar)]]

    # get the index of sample keeping
    samples.keep.index <- which(feature.selected.vec %in% input$cat_filter_options)

    pstat@otu_table@.Data <- pstat@otu_table@.Data[,samples.keep.index]
    pstat@sam_data@.Data <- lapply(pstat@sam_data@.Data, function(x) {x <- x[samples.keep.index]})
    pstat@sam_data@row.names <- pstat@sam_data@row.names[samples.keep.index]
    shinyInput <- list(pstat = pstat)
    vals$shiny.input <- shinyInput

    #updateCovariate()
    updateSample()
  })
})

observeEvent(input$filter_read_micro,{
  withBusyIndicatorServer("filter_read_micro", {
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat


    read.count.mat <- pstat@otu_table@.Data
    read.count.mean <- apply(read.count.mat, 1, mean)
    micro.keep.index <- which(read.count.mean >= input$read_filter_min_micro)

    pstat@otu_table@.Data <- pstat@otu_table@.Data[micro.keep.index,]
    pstat@tax_table@.Data <- pstat@tax_table@.Data[micro.keep.index,]
    shinyInput <- list(pstat = pstat)
    vals$shiny.input <- shinyInput

    updateTaxLevel()
  })
})

observeEvent(input$filter_ra_micro,{
  withBusyIndicatorServer("filter_ra_micro", {
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat


    read.count.mat <- pstat@otu_table@.Data
    ra.mat <- getRelativeAbundance(read.count.mat)
    ra.mean <- apply(ra.mat, 1, mean)
    micro.keep.index <- which(ra.mean >= input$ra_filter_min_micro)

    pstat@otu_table@.Data <- pstat@otu_table@.Data[micro.keep.index,]
    pstat@tax_table@.Data <- pstat@tax_table@.Data[micro.keep.index,]
    shinyInput <- list(pstat = pstat)
    vals$shiny.input <- shinyInput

    updateTaxLevel()
  })
})

observeEvent(input$filter_prev_micro,{
  withBusyIndicatorServer("filter_prev_micro", {
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    read.count.mat <- pstat@otu_table@.Data

    micro.prev <- apply(read.count.mat, 1, function(x) (sum(x >= 1)/ncol(read.count.mat)))
    micro.keep.index <- which(micro.prev >= input$prev_filter_min)

    pstat@otu_table@.Data <- pstat@otu_table@.Data[micro.keep.index,]
    pstat@tax_table@.Data <- pstat@tax_table@.Data[micro.keep.index,]
    shinyInput <- list(pstat = pstat)
    vals$shiny.input <- shinyInput

    updateTaxLevel()
  })
})

# Render summary table
output$contents_summary <- renderTable({
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat
  req(pstat)
  summarizeTable(pstat)
})

# Data input summary
output$contents.count <- DT::renderDataTable({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    if (!is.null(input$countsfile.pathoscope)){
        if (input$uploadChoice == "pathofiles"){
        req(input$countsfile.pathoscope)
        df <- read.csv(input$countsfile.pathoscope[[1, 'datapath']],
                       skip = 1,
                       header = TRUE,
                       sep = input$sep.ps)
        return(df)
        }
    }
},
options = list(
    paging = TRUE, scrollX = TRUE, pageLength = 5
))

output$contents.meta <- DT::renderDataTable({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    if (!is.null(input$annotfile.ps)){
        if (input$uploadChoice == "pathofiles"){
        req(input$countsfile.pathoscope)

        df <- read.csv(input$annotfile.ps$datapath,
                       header = input$header.ps,
                       sep = input$sep.ps)
        return(df)
        }
    }
},
options = list(
    paging = TRUE, scrollX = TRUE, pageLength = 5
))

### data input summary
output$contents.count.2 <- DT::renderDataTable({

  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, head of that data file by default,
  # or all rows if selected, will be shown.

  if (!is.null(input$countsfile)){
    if (input$uploadChoice == "count"){
      req(input$countsfile)
      df <- read.csv(input$countsfile$datapath,
                     header = input$header.count,
                     sep = input$sep.count)
      return(df)
    }
  }
},
options = list(
  paging = TRUE, scrollX = TRUE, pageLength = 5
))

output$contents.meta.2 <- DT::renderDataTable({

  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, head of that data file by default,
  # or all rows if selected, will be shown.

  if (!is.null(input$annotfile.count)){
    if (input$uploadChoice == "count"){
      req(input$annotfile.count)
      df <- read.csv(input$annotfile.count$datapath,
                     header = input$header.count,
                     sep = input$sep.count)
      return(df)
    }
  }
},
options = list(
  paging = TRUE, scrollX = TRUE, pageLength = 5
))

output$contents.taxonomy <- DT::renderDataTable({

  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, head of that data file by default,
  # or all rows if selected, will be shown.


  if (!is.null(input$taxon.table)){
    if (input$uploadChoice == "count"){
      req(input$taxon.table)

      df <- read.csv(input$taxon.table$datapath,
                     header = input$header.count,
                     sep = input$sep.count)
      return(df)
    }
  }
},
options = list(
  paging = TRUE, scrollX = TRUE, pageLength = 5
))

# setInputs(FALSE)
findAllTaxData <- function(taxonLevel) {
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat
  if ((taxonLevel == "no rank"))  {
    taxcountdata <- pstat@otu_table
    taxdata <- findRAfromCount(taxcountdata)
    taxcountdata <- data.frame(taxcountdata)
    taxdata <- data.frame(taxdata)
  } else  {
    physeq2 <- tax_glom(pstat, taxonLevel)
    taxcountdata <- physeq2@otu_table
    taxdata <- findRAfromCount(taxcountdata)
    taxcountdata <- data.frame(taxcountdata)
    taxdata <- data.frame(taxdata)
    labvec <- as(tax_table(physeq2)[, taxonLevel], "character")
    labvec <- unlist(lapply(labvec,
                            function(x){paste0( taxonLevel, "|", x)}))
    labvec[is.na(labvec)] <- ""
    names <- rownames(taxcountdata)
    for (i in 1:length(names))  {
      #tid <- grepTid(names[i])
      #labvec[i] <- paste0( "ti|", tid, "|", labvec[i])
      labvec[i] <- paste0(names[i], labvec[i])
    }
    rownames(taxcountdata) <- labvec
    rownames(taxdata) <- labvec
  }

    vals$taxcountdata <- taxcountdata
    vals$taxdata <- taxdata

}
findTaxData <- reactive({
  findAllTaxData(input$taxl)
  vals$taxdata
})

findTaxCountData <- reactive({
  findAllTaxData(input$taxl)
  vals$taxcountdata
})

findTaxDataTable <- reactive({
    findAllTaxData(input$taxlTable)
    vals$taxdata
})

findTaxCountDataTable <- reactive({
    findAllTaxData(input$taxlTable)
    vals$taxcountdata
})

findTaxCountDataDE <- reactive({
  findAllTaxData(input$taxlde)
  vals$taxcountdata
})

findCountTable <- reactive({
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat
  if (input$taxlTable !="no rank"){
    pstat <- tax_glom(pstat, input$taxlTable)
  }
  df.out <- pstat@otu_table@.Data
  rownames(df.out) <- make.names(TranslateIdToTaxLevel(pstat, rownames(df.out),
                                                       input$taxlTable), unique=TRUE)
  df.out
})

findRATable <- reactive({
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat
  if (input$taxlTable !="no rank"){
    pstat <- tax_glom(pstat, input$taxlTable)
  }
  df.ra <- getRelativeAbundance(pstat@otu_table@.Data)
  rownames(df.ra) <- make.names(TranslateIdToTaxLevel(pstat, rownames(df.ra),
                                                     input$taxlTable), unique=TRUE)
  df.ra
})

output$TaxRAsummary <- renderPrint({
  summary(findTaxData())
})

# These are options for rendering datatables
dtopts <- list(scrollX=TRUE, paging=TRUE)

# Render relative abundance table
output$TaxRAtable <- DT::renderDataTable(findRATable(),
                                         options=dtopts, rownames=T)


# Render count table
output$TaxCountTable <- DT::renderDataTable(findCountTable(),
options=dtopts,
rownames=T)

output$downloadData <- downloadHandler(filename = function() {
  paste0("sample_data_", input$taxlTable, ".csv", sep = "")
}, content = function(file) {
  df.out <- findRATable()
  write.csv(df.out, file)
})

output$downloadCountData <- downloadHandler(filename = function() {
  paste0("sample_data_count_", input$taxlTable, ".csv", sep = "")
}, content = function(file) {
  df.out <- findCountTable()
  write.csv(df.out, file)
})

output$nbins <- renderUI({
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat
  tables <- pstat.extraction(pstat)
  SAM_DATA <- tables$SAM
  vals <-unlist(SAM_DATA[,input$bin_cov,drop=TRUE])
  sliderInput("nbins", label="Number of Bins", min=2, max=length(unique(vals)), value=2, step=1)
})
output$bin_to1 <- renderPrint({
  x <- sort(as.numeric(unlist(strsplit(input$bin_breaks,","))))
  print(x)
})
output$bin_to2 <- renderPrint({
  x <- unlist(strsplit(input$bin_labels,","))
  print(x)
})

output$unbin_plot <- renderPlotly({
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat  
  tables <- pstat.extraction(pstat)
  SAM_DATA <- tables$SAM
  unbinned = as.numeric(unlist(SAM_DATA[,input$bin_cov]))
  fit <- density(unbinned)
  p <- plot_ly(x=fit$x, y=fit$y, type="scatter", mode="lines", fill="tozeroy") %>%
         layout(title="Unbinned Density",
                margin = list(l=0,r=0,t=30,b=30))
  p$elementId <- NULL
  return(p)
})

output$bin_plot <- renderPlotly({
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat
  tables <- pstat.extraction(pstat)
  SAM_DATA <- tables$SAM
  unbinned <- as.numeric(unlist(SAM_DATA[,input$bin_cov]))
  
  # Bins
  nbins <- input$nbins
  n <- input$nbins

  if (!is.null(input$nbins)) {
    # Overrides numnber of bins
    bin_breaks = sort(as.numeric(unlist(strsplit(input$bin_breaks,","))))
    if (length(bin_breaks) > 1) {
      nbins = bin_breaks
      n = length(bin_breaks)-1
    }
    # Overrides labels of bins
    bin_labels = unlist(strsplit(input$bin_labels,","))
    labels <- NULL
    if (length(bin_labels) == n) {
      labels <- bin_labels
    }
    binned <- cut.default(unbinned, nbins, labels=labels)
    p <- plot_ly(y = binned, type = "histogram") %>%
         layout(title = input$bin_cov,
                xaxis = list(title = "Frequency"),
                yaxis = list(title = "Bins"),
                margin = list(l=80)
         )
    p$elementId <- NULL
    return(p)
  }
  return()
})

observeEvent(input$create_bins, {
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat
  tables <- pstat.extraction(pstat)
  SAM_DATA <- tables$SAM
  unbinned <- as.numeric(unlist(SAM_DATA[,input$bin_cov]))
  
  # Bins
  nbins <- input$nbins
  n <- input$nbins

  if (!is.null(input$nbins)) {
    # Overrides numnber of bins
    bin_breaks = sort(as.numeric(unlist(strsplit(input$bin_breaks,","))))
    if (length(bin_breaks) > 1) {
      nbins = bin_breaks
      n = length(bin_breaks)-1
    }
    # Overrides labels of bins
    bin_labels = unlist(strsplit(input$bin_labels,","))
    labels <- NULL
    if (length(bin_labels) == n) {
      labels <- bin_labels
    }
    binned <- cut.default(unbinned, nbins, labels=labels)

    SAM_DATA[,input$new_covariate] <- as.character(binned)
    pstat@sam_data <- SAM_DATA
    shinyInput <- list(pstat = pstat)
    vals$shiny.input <- shinyInput
    updateCovariate()
    updateSample()
  }
})