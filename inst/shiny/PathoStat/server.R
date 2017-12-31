library(shiny)
library(ggvis)
library(d3heatmap)
library(reshape2)
library(limma)
library(edgeR)
library(DESeq2)
library(phyloseq)
library(ape)
library(stats)
library(PathoStat)
library(alluvial)
library(plotly)
library(webshot)
library(vegan)
library(dplyr)

# Converts decimal percentage to string with specified digits
pct2str <- function(v, digits=2) {sprintf(paste0('%.',digits,'f'), v*100)}






shinyServer(function(input, output, session) {
  # needed information from PathoStat
  shinyInput <- getShinyInput()
  pstat <- shinyInput$pstat
  sample_data <- pstat@sam_data
  #pdata <- data.frame(shinyInput$batch, shinyInput$condition)
  #modmatrix = model.matrix(~as.factor(shinyInput$condition), data=pdata)
  pdata <- data.frame(c(sample_data[,1])[[1]], c(sample_data[,2])[[1]])
  modmatrix = model.matrix(~as.factor(c(sample_data[,2])[[1]]), data=pdata)
  pca <- BatchQC::batchqc_pca(pstat@otu_table, batch=c(sample_data[,1])[[1]],
                              mod=modmatrix)
  shinyInput <- c(shinyInput, list("pc"=data.frame(pca$x), "vars"=pca$sdev^2))
  setShinyInput(shinyInput)
  
  setInputs <- function(combatFlag) {
    if (combatFlag) {
      setShinyInput(getShinyInputCombat())
    } else {
      setShinyInput(getShinyInputOrig())
    }
  }
  # setInputs(FALSE)
  findAllTaxData <- function(taxonLevel) {
    #taxcountdata <- findTaxonLevelData(pstat@otu_table,
    #    shinyInput$taxonLevels, taxonLevel)
    shinyInput <- getShinyInput()
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
        tid <- grepTid(names[i])
        labvec[i] <- paste0( "ti|", tid, "|", labvec[i])
      }
      rownames(taxcountdata) <- labvec
      rownames(taxdata) <- labvec
    }
    if (is.null(shinyInput$taxcountdata)) {
      shinyInput <- c(shinyInput, list(taxcountdata = taxcountdata))
    } else {
      shinyInput$taxcountdata <- taxcountdata
    }
    if (is.null(shinyInput$taxdata)) {
      shinyInput <- c(shinyInput, list(taxdata = taxdata))
    } else {
      shinyInput$taxdata <- taxdata
    }
    setShinyInput(shinyInput)
  }
  findTaxData <- eventReactive(input$taxl, {
    findAllTaxData(input$taxl)
    shinyInput <- getShinyInput()
    shinyInput$taxdata
  })
  
  
  findTaxCountData <- eventReactive(input$taxl, {
    findAllTaxData(input$taxl)
    shinyInput <- getShinyInput()
    shinyInput$taxcountdata
  })
  
  findTaxCountDataDE <- eventReactive(input$taxlde, {
    findAllTaxData(input$taxlde)
    shinyInput <- getShinyInput()
    shinyInput$taxcountdata
  })
  
  findNormalizedCount <- function() {
    shinyInput <- getShinyInput()
    if(input$norm == 'Quantile coreOTU Normalization')  {
      dat <- coreOTUQuantile(shinyInput$taxcountdata, input$otuthreshold,
                             input$prevalence)
    } else if(input$norm == 'Library Size Scaling')  {
      dat <- sizeNormalize(shinyInput$taxcountdata)
    } else  {
      dat <- coreOTUNormalize(shinyInput$taxcountdata, input$ebweight,
                              input$otuthreshold, input$prevalence)
    }
    dat
  }
  
  tax_ra_bp <- reactive({
    if (is.null(input$taxl)) {
      return()
    }
    taxdata <- findTaxData()
    dat <- melt(cbind(taxdata, ind = as.character(rownames(taxdata))),
                id.vars = c("ind"))
    
    if ((input$taxl == "no rank")){
      covariates.tmp <- colnames(sample_data(shinyInput$pstat))
      dat$condition.select <- rep(shinyInput$pstat@sam_data@.Data[[
        which(covariates.tmp %in% input$select_condition)]], 
        each = dim(taxdata)[1])
      dat <- dat[order(dat$condition.select),]
      dat$condition.select.id <- paste(dat$condition.select,
                                       as.character(dat$variable), sep = "-")
    }else{
      pstat.new <- tax_glom(shinyInput$pstat, input$taxl)
      covariates.tmp <- colnames(sample_data(pstat.new))
      dat$condition.select <- rep(pstat.new@sam_data@.Data[[
        which(covariates.tmp %in% input$select_condition)]], 
        each = dim(taxdata)[1])
      dat <- dat[order(dat$condition.select),]
      dat$condition.select.id <- paste(dat$condition.select,
                                       as.character(dat$variable), sep = "-")
    }
    
    # sort by selecting variables
    
    
    dat %>% ggvis(x = ~condition.select.id, y = ~value, fill = ~as.factor(ind)) %>%
      layer_bars(stack = TRUE) %>%
      add_tooltip(function(dat2) {
        paste0("Sample: ", dat2[2], "<br />", "Genome: ", dat2[1],
               "<br />", "RA: ", round(dat2[4] - dat2[3], 4))
      }, "hover") %>%
      # add_axis('x', subdivide = 1,
      #   values = 1:length(colnames(shinyInput$data)),
      add_axis("x", title = "", properties = axis_props(title =
                                                          list(fontSize = 15), labels = list(angle = 90,
                                                                                             align = "left", baseline = "middle"))) %>%
      add_axis("y", title = "Relative Abundance (RA)", properties =
                 axis_props(title = list(fontSize = 15),
                            labels = list(fontSize = 10))) %>% add_legend("fill", title =
                                                                            "Genomes", properties = legend_props(title = list(fontSize = 15),
                                                                                                                 labels = list(fontSize = 10))) %>% set_options(width = "auto",
                                                                                                                                                                height = "auto")
  })
  tax_ra_bp %>% bind_shiny("TaxRelAbundancePlot")
  
  output$TaxRAsummary <- renderPrint({
    summary(findTaxData())
  })
  
  # These are options for rendering datatables
  dtopts <- list(scrollX=TRUE, paging=TRUE)
  
  # Format relative abundance table
  # Converts percents to strings and expands taxa name
  format_RA_table <- function(tmp) {
    tmp %>% dplyr::add_rownames("fullname") %>%
      dplyr::mutate_each(dplyr::funs(pct2str), -fullname) %>%
      tidyr::separate(fullname, c('fi1', 'taxid', 'fi2', input$taxl),
                      sep='\\|') %>%
      dplyr::select_(.dots=c(as.name(input$taxl), 'taxid', colnames(tmp)))
  }
  
  # Render relative abundance table
  output$TaxRAtable <- DT::renderDataTable(format_RA_table(findTaxData()),
                                           options=dtopts, rownames=F)
  
  # Format count table:
  # Just expands taxa name
  format_Count_table <- function(tmp) {
    tmp %>% dplyr::add_rownames("fullname") %>%
      tidyr::separate(fullname, c('fi1', 'taxid', 'fi2', input$taxl),
                      sep='\\|') %>%
      dplyr::select_(.dots=c(as.name(input$taxl), 'taxid', colnames(tmp)))
  }
  # Render count table
  output$TaxCountTable <- DT::renderDataTable(format_Count_table(
    findTaxCountData()), options=dtopts, rownames=F)
  
  output$downloadData <- downloadHandler(filename = function() {
    paste0("sample_data_", input$taxl, ".csv", sep = "")
  }, content = function(file) {
    shinyInput <- getShinyInput()
    write.csv(shinyInput$taxdata, file)
  })
  
  output$downloadCountData <- downloadHandler(filename = function() {
    paste0("sample_data_count_", input$taxl, ".csv", sep = "")
  }, content = function(file) {
    shinyInput <- getShinyInput()
    write.csv(shinyInput$taxcountdata, file)
  })
  
  findPhyseqData <- function() {
    ids <- rownames(shinyInput$data)
    taxmat <- findTaxonMat(ids, shinyInput$taxonLevels)
    OTU <- otu_table(shinyInput$countdata, taxa_are_rows = TRUE)
    TAX <- tax_table(taxmat)
    physeq <- phyloseq(OTU, TAX)
    sampledata = sample_data(data.frame(condition=as.factor(
      c(sample_data[,2])[[1]]), batch=as.factor(c(sample_data[,1])[[1]]),
      row.names=sample_names(physeq), stringsAsFactors=FALSE))
    random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=
                          taxa_names(physeq))
    physeq1 <- merge_phyloseq(physeq, sampledata, random_tree)
    return(physeq1)
  }
  setGgplotTheme <- function() {
    ggplot2::theme_set(ggplot2::theme_bw())
  }
  pal = "Set1"
  scale_colour_discrete <- function(palname = pal, ...) {
    scale_colour_brewer(palette = palname, ...)
  }
  scale_fill_discrete <- function(palname = pal, ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  
  
  # independent heatmap plotting function in the server using specific data from input
  plotHeatmapColorServer <- function(){
    physeq1 <- shinyInput$pstat
    
    if (input$taxl=="no rank")  {
      if (input$checkbox_heatmap){
        add.colorbar <- "auto"
      } else{
        add.colorbar <- NULL
      }
      return(plotHeatmapColor(physeq1@otu_table@.Data, 
                              do.scale = input$checkbox_heatmap_scale,
                              physeq1@sam_data[[input$select_heatmap_condition]], 
                              annotationColors = add.colorbar,
                              columnTitle = paste("Heatmap with colorbar representing", 
                                                  input$select_heatmap_condition, sep = " ")))
    } else  {
      physeq2 <- tax_glom(physeq1, input$taxl)
      if (input$checkbox_heatmap){
        add.colorbar <- "auto"
      } else{
        add.colorbar <- NULL
      }
      return(plotHeatmapColor(physeq2@otu_table@.Data, 
                              do.scale = input$checkbox_heatmap_scale,
                              physeq2@sam_data[[input$select_heatmap_condition]], 
                              annotationColors = add.colorbar,
                              columnTitle = paste("Heatmap with colorbar representing", 
                                                  input$select_heatmap_condition, sep = " ")))
    }
  }
  
  # show heatmap in the shiny app by calling the plotting function
  output$Heatmap <- renderPlot({
    plotHeatmapColorServer()
    
  })
  # download heatmap by calling the plotting function.
  # note: add "print()" function to the plotting function, and add dev.off()
  output$download_heatmap_pdf <- downloadHandler(
    filename = function() {
      paste('heatmap', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      pdf(file)
      #### add "print()" to plotting function to work!!
      print(plotHeatmapColorServer())
      ####
      dev.off()
    }
    
  )
  
  plotAlphaServer <- function(){
    physeq1 <- shinyInput$pstat
    
    if (input$taxl.alpha !="no rank")  {
      physeq1 <- tax_glom(physeq1, input$taxl.alpha)
    }
    
    meta.data <- physeq1@sam_data
    meta.data$sample.name <- rownames(meta.data)
    meta.data$richness <- estimate_richness(physeq = physeq1, split = T, measures = input$select_alpha_div_method)[,1]
    colnames(meta.data)[which(colnames(meta.data) == input$select_alpha_div_condition)] <- "condition"
    g <- ggplot(meta.data, aes(condition, richness, text=sample.name, color = condition)) + 
      geom_point() + geom_boxplot() +
      labs(title = paste("Alpha diversity between ", 
                         input$select_alpha_div_condition, 
                         " (", input$select_alpha_div_method, ")", sep = ""))
      
    ggplotly(g, tooltip="text")
  }
  
  output$AlphaDiversity <- renderPlotly({
    plotAlphaServer()
  })
  
  observeEvent(input$download_alpha,{
    if (!require("webshot")) install.packages("webshot")
    tmpFile <- tempfile(pattern = "Alpha_diversity_", fileext = ".pdf")
    export(plotAlphaServer(), file = tmpFile)
    browseURL(tmpFile)}
  )
  
  output$table.alpha <- DT::renderDataTable({
    physeq1 <- shinyInput$pstat
    if (input$taxl.alpha !="no rank")  {
      physeq1 <- tax_glom(physeq1, input$taxl.alpha)
    }
    meta.data <- physeq1@sam_data
    meta.data$sample.name <- rownames(meta.data)
    meta.data$richness <- estimate_richness(physeq = physeq1, split = T, measures = input$select_alpha_div_method)[,1]
    colnames(meta.data)[which(colnames(meta.data) == input$select_alpha_div_condition)] <- "condition"
    rownames(meta.data) <- 1:nrow(meta.data)
    meta.data <- as_tibble(meta.data)
    meta.data <- meta.data %>% select(sample.name, condition, richness)
    DT::datatable(meta.data)
    
  })
  
  output$download_table_alpha <- downloadHandler(
    filename = function() { paste('Alpha_diversity_table', '.csv', sep='') },
    content = function(file) {
      physeq1 <- shinyInput$pstat
      if (input$taxl.alpha !="no rank")  {
        physeq1 <- tax_glom(physeq1, input$taxl.alpha)
      }
      meta.data <- physeq1@sam_data
      meta.data$sample.name <- rownames(meta.data)
      meta.data$richness <- estimate_richness(physeq = physeq1, split = T, measures = input$select_alpha_div_method)[,1]
      colnames(meta.data)[which(colnames(meta.data) == input$select_alpha_div_condition)] <- "condition"
      rownames(meta.data) <- 1:nrow(meta.data)
      meta.data <- as_tibble(meta.data)
      meta.data <- meta.data %>% select(sample.name, condition, richness)
      write.csv(data.frame(meta.data), file)
    }
  )
  
  output$alpha.stat.test <- renderPrint({ 
    physeq1 <- pstat
    if (input$taxl.alpha !="no rank")  {
      physeq1 <- tax_glom(physeq1, input$taxl.alpha)
    }
    meta.data <- physeq1@sam_data
    meta.data$sample.name <- rownames(meta.data)
    meta.data$richness <- estimate_richness(physeq = physeq1, split = T, measures = input$select_alpha_div_method)[,1]
    colnames(meta.data)[which(colnames(meta.data) == input$select_alpha_div_condition)] <- "condition"
    meta.data <- data.frame(meta.data)
    meta.data$condition <- as.factor(meta.data$condition)
    
    if (length(unique(meta.data$condition)) == 2){
      if (input$select_alpha_stat_method == "Mann-Whitney"){
        wilcox.test(richness ~ condition, data = meta.data)
      } else{
        print("Condition level number is 2, please use Mann-Whitney test.")
      }
      
    } else if (length(unique(meta.data$condition)) > 2){
      if (input$select_alpha_stat_method == "Mann-Whitney"){
        print("Condition level number is larger than 2, pairwise Mann-Whitney test is applied. You could also check Kruskal-Wallis test result.")
        print("---------------------------------------")
        result.list <- list()
        for (i in 1:length(unique(meta.data$condition))){
          meta.data.tmp <- meta.data[which(meta.data$condition != unique(meta.data$condition)[i]),]
          print(paste(unique(meta.data.tmp$condition), collapse = " and "))
          result.list[[i]] <- wilcox.test(richness ~ condition, data = meta.data.tmp)
          print(result.list[[i]])
          print("------------------------")
        }
        
      } else{
        kruskal.test(richness ~ condition, data = meta.data)
      }
      
    } else{
      "Condition level must be at least 2."
    }
    
    
  })
  
  
  # independent heatmap plotting function in the server using specific data from input
  plotBetaHeatmapColorServer <- function(){
    physeq1 <- shinyInput$pstat
    
    if (input$taxl.beta=="no rank")  {
      if (input$checkbox_beta_heatmap){
        add.colorbar <- "auto"
      } else{
        add.colorbar <- NULL
      }
      
      dist.mat = phyloseq::distance(physeq1, method = input$select_beta_div_method)
      dist.mat <- as.matrix(dist.mat)
      return(plotHeatmapColor(dist.mat, 
                              do.scale = input$checkbox_beta_heatmap_scale,
                              physeq1@sam_data[[input$select_beta_heatmap_condition]], 
                              annotationColors = add.colorbar,
                              columnTitle = paste("Heatmap with colorbar representing", 
                                                  input$select_beta_heatmap_condition, sep = " ")))
    } else  {
      physeq2 <- tax_glom(physeq1, input$taxl.beta)
      if (input$checkbox_beta_heatmap){
        add.colorbar <- "auto"
      } else{
        add.colorbar <- NULL
      }
      dist.mat = phyloseq::distance(physeq2, method = input$select_beta_div_method)
      dist.mat <- as.matrix(dist.mat)
      return(plotHeatmapColor(dist.mat,
                              do.scale = input$checkbox_beta_heatmap_scale,
                              physeq2@sam_data[[input$select_beta_heatmap_condition]], 
                              annotationColors = add.colorbar,
                              columnTitle = paste("Heatmap with colorbar representing", 
                                                  input$select_beta_heatmap_condition, sep = " ")))
    }
  }
  

  output$BetaDiversityHeatmap <- renderPlot({
    plotBetaHeatmapColorServer()
  })
  
  
  output$download_beta_heatmap_pdf <- downloadHandler(
    filename = function() {
      paste('heatmap_beta', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      pdf(file)
      #### add "print()" to plotting function to work!!
      print(plotBetaHeatmapColorServer())
      ####
      dev.off()
    }
    
  )
  
  
  plotBetaBoxplotServer <- function(){
    physeq1 <- shinyInput$pstat
    
    if (input$taxl.beta !="no rank")  {
      physeq1 <- tax_glom(physeq1, input$taxl.beta)
    }
    
    meta.data <- physeq1@sam_data
    meta.data$sample.name <- rownames(meta.data)
    colnames(meta.data)[which(colnames(meta.data) == input$select_beta_condition)] <- "condition"
    
    dist.tmp = phyloseq::distance(physeq1, method = input$select_beta_div_method)
    dist.mat <- as.matrix(dist.tmp)
    dist.within.a <- c()
    dist.within.b <- c()
    dist.between <- c()
    for (i in 1:nrow(dist.mat)){
      for (j in 1:nrow(dist.mat)) {
        if (meta.data$condition[i] == unique(meta.data$condition)[1] & 
            meta.data$condition[j] == unique(meta.data$condition)[1]){
          dist.within.a <- c(dist.within.a, dist.mat[i,j])
        } else if (meta.data$condition[i] == unique(meta.data$condition)[2] & 
                   meta.data$condition[j] == unique(meta.data$condition)[2]){
          dist.within.b <- c(dist.within.b, dist.mat[i,j])
        } else{
          dist.between <- c(dist.between, dist.mat[i,j])
        }
        
      }
    }
    
    y.axis <- list(
      title = paste(input$select_beta_div_method, "Distance", sep = " ")
    )
    
    p <- plot_ly(y = ~dist.within.a, type = "box", name = paste("Within", unique(meta.data$condition)[1])) %>%
      add_trace(y = ~dist.within.b, name = paste("Within", unique(meta.data$condition)[2])) %>%
      add_trace(y = ~dist.between, name = "Between 2 conditions") %>%
      layout(yaxis = y.axis)
      
    p
  }
  
  output$BetaDiversityBoxplot <- renderPlotly({
    plotBetaBoxplotServer()
  })
  
  observeEvent(input$download_beta_boxplot,{
    if (!require("webshot")) install.packages("webshot")
    tmpFile <- tempfile(pattern = "Beta_diversity_boxplot", fileext = ".pdf")
    export(plotBetaBoxplotServer(), file = tmpFile)
    browseURL(tmpFile)}
  )
  
  
  
  
  output$beta.stat.test <- renderPrint({ 
    
    if (input$select_beta_stat_method == "PERMANOVA"){
      physeq1 <- pstat
      if (input$taxl.beta !="no rank")  {
        physeq1 <- tax_glom(physeq1, input$taxl.beta)
      }
      meta.data <- physeq1@sam_data
      meta.data$sample.name <- rownames(meta.data)
      colnames(meta.data)[which(colnames(meta.data) == input$select_beta_condition)] <- "condition"
      meta.data <- data.frame(meta.data)
      meta.data$condition <- as.factor(meta.data$condition)
      
      set.seed(99)
      dist.tmp = phyloseq::distance(physeq1, method = input$select_beta_div_method)
      beta.div <- adonis2(dist.tmp~condition, data=meta.data, 
                          permutations = input$num.permutation.permanova, strata="PLOT")
      beta.div
      
    } else {
      physeq1 <- shinyInput$pstat
      
      if (input$taxl.beta !="no rank")  {
        physeq1 <- tax_glom(physeq1, input$taxl.beta)
      }
      
      meta.data <- physeq1@sam_data
      meta.data$sample.name <- rownames(meta.data)
      colnames(meta.data)[which(colnames(meta.data) == input$select_beta_condition)] <- "condition"
      
      dist.tmp = phyloseq::distance(physeq1, method = input$select_beta_div_method)
      dist.mat <- as.matrix(dist.tmp)
      dist.within.a <- c()
      dist.within.b <- c()
      dist.between <- c()
      for (i in 1:nrow(dist.mat)){
        for (j in 1:nrow(dist.mat)) {
          if (meta.data$condition[i] == unique(meta.data$condition)[1] & 
              meta.data$condition[j] == unique(meta.data$condition)[1]){
            dist.within.a <- c(dist.within.a, dist.mat[i,j])
          } else if (meta.data$condition[i] == unique(meta.data$condition)[2] & 
                     meta.data$condition[j] == unique(meta.data$condition)[2]){
            dist.within.b <- c(dist.within.b, dist.mat[i,j])
          } else{
            dist.between <- c(dist.between, dist.mat[i,j])
          }
          
        }
      }
      dist.list <- list(dist.within.a, dist.within.b, dist.between)
      names(dist.list) <- c(unique(meta.data$condition)[1], unique(meta.data$condition)[2], "between")
      
      if (input$select_beta_stat_method == "Mann-Whitney"){
        print("Condition level number is larger than 2, pairwise Mann-Whitney test is applied. 
              You could also check Kruskal-Wallis test result.")
        print("---------------------------------------")
        result.list <- list()
        for (i in 1:length(dist.list)){
          dist.list.tmp <- dist.list[which(names(tmp) != names(tmp)[i])]
          print(paste(names(dist.list.tmp), collapse = " and "))
          result.list[[i]] <- wilcox.test(dist.list.tmp[[1]], dist.list.tmp[[2]])
          print(result.list[[i]])
          print("------------------------")
        }
        
      } else{
        kruskal.test(list(dist.within.a, dist.within.b, dist.between))
      }
      
    }

    
    

    
  })
  
  output$table.beta <- DT::renderDataTable({
    physeq1 <- shinyInput$pstat
    
    if (input$taxl.beta=="no rank")  {
      dist.mat = phyloseq::distance(physeq1, method = input$select_beta_div_method)
    } else{
      physeq2 <- tax_glom(physeq1, input$taxl.beta)
      dist.mat = phyloseq::distance(physeq2, method = input$select_beta_div_method)
    }
    dist.mat <- as.matrix(dist.mat)
    DT::datatable(dist.mat)
    
  })
  
  output$download_table_beta <- downloadHandler(
    filename = function() { paste('Beta_diversity_table', '.csv', sep='') },
    content = function(file) {
      physeq1 <- shinyInput$pstat
      
      if (input$taxl.beta=="no rank")  {
        dist.mat = phyloseq::distance(physeq1, method = input$select_beta_div_method)
      } else{
        physeq2 <- tax_glom(physeq1, input$taxl.beta)
        dist.mat = phyloseq::distance(physeq2, method = input$select_beta_div_method)
      }
      dist.mat <- as.matrix(dist.mat)
      write.csv(data.frame(dist.mat), file)
    }
  )
  
  
  
  
  
  output$ExploratoryTree <- renderPlot({
    physeq1 <- shinyInput$pstat
    cn <- colnames(physeq1@sam_data)
    cn[1] <- "batch"
    cn[2] <- "condition"
    colnames(physeq1@sam_data) <- cn
    setGgplotTheme()
    plot_tree(physeq1, color="condition", label.tips="genus",
              size="Abundance")
  })
  output$BiPlot <- renderPlot({
    physeq1 <- shinyInput$pstat
    physeq1 <- phyloseq(otu_table(physeq1), phy_tree(physeq1),
                        tax_table(physeq1), sample_data(physeq1))
    cn <- colnames(physeq1@sam_data)
    cn[1] <- "batch"
    cn[2] <- "condition"
    colnames(physeq1@sam_data) <- cn
    setGgplotTheme()
    if (input$colorBiP=='None') color <- NULL else color <- input$colorBiP
    if (input$shapeBiP=='None') shape <- NULL else shape <- input$shapeBiP
    if (input$labelBiP=='None') label <- NULL else label <- input$labelBiP
    physeq.ord <- ordinate(physeq1, method=input$methodBiP, distance="bray")
    p <- plot_ordination(physeq1, physeq.ord, type = "biplot",
                         color=color, shape=shape, label=label)
    if (!is.null(label))  {
      p$layers <- p$layers[-2]
      p <- p+ggplot2::geom_text(mapping=
                                  ggplot2::aes(label=get(label)), size=4, vjust=1.0,
                                check_overlap=TRUE)
    }
    p
  })
  output$CoOccurrence <- renderPlot({
    physeq1 <- shinyInput$pstat
    cn <- colnames(physeq1@sam_data)
    cn[1] <- "batch"
    cn[2] <- "condition"
    colnames(physeq1@sam_data) <- cn
    setGgplotTheme()
    if (input$colorCo=='None') color <- NULL else color <- input$colorCo
    if (input$shapeCo=='None') shape <- NULL else shape <- input$shapeCo
    if (input$labelCo=='None') label <- NULL else label <- input$labelCo
    ig<-make_network(physeq1, dist.fun="jaccard", type = "taxa",
                     max.dist=input$max.dist) #max.dist=0.4 default
    p <- plot_network(ig, physeq1, line_weight=0.4, type = "taxa",
                      color = color, shape = shape, label = label)
    if (!is.null(label))  {
      p$layers <- p$layers[-2]
      p <- p+ggplot2::geom_text(mapping=ggplot2::aes(label=get(label)),
                                size=4, vjust=1.0, check_overlap=TRUE)
    }
    p
  })
  
  # interactive PCA
  PCA <- reactive({
    shinyInput <- getShinyInput()
    pc <- shinyInput$pc
    data.frame(pc[, c(input$xcol.new, input$ycol.new)])
  })

  # independent heatmap plotting function in the server using specific data from input
  plotPCAPlotlyServer <- function(){
    physeq1 <- shinyInput$pstat
    if (input$taxl.pca=="no rank")  {
      plotPCAPlotly(df.input = physeq1@otu_table@.Data, 
                    condition.vec = physeq1@sam_data[[input$select_pca_condition]],
                    condition.name = input$select_pca_condition,
                    pc.a = paste("PC", input$xcol.new, sep = ""),
                    pc.b = paste("PC", input$ycol.new, sep = ""),
                    columnTitle = paste("PCA with colorbar representing", 
                                        input$select_pca_condition, sep = " "))
    } else  {
      physeq2 <- tax_glom(physeq1, input$taxl.pca)
      plotPCAPlotly(df.input = physeq2@otu_table@.Data, 
                    condition.vec = physeq2@sam_data[[input$select_pca_condition]],
                    condition.name = input$select_pca_condition,
                    pc.a = paste("PC", input$xcol.new, sep = ""),
                    pc.b = paste("PC", input$ycol.new, sep = ""),
                    columnTitle = paste("PCA with colorbar representing", 
                                        input$select_pca_condition, sep = " "))}
  }
  
  # show heatmap in the shiny app by calling the plotting function
  output$pca.plotly <- renderPlotly({
    plotPCAPlotlyServer()
    
  })
  
  
  observeEvent(input$download_pca,{
    if (!require("webshot")) install.packages("webshot")
    tmpFile <- tempfile(pattern = "PCA_", fileext = ".pdf")
    export(plotPCAPlotlyServer(), file = tmpFile)
    browseURL(tmpFile)}
  )

  # interactive PCA table
  output$PCAtable <- DT::renderDataTable({
    physeq1 <- shinyInput$pstat
    if (input$taxl.pca=="no rank")  {
      #test and fix the constant/zero row
      if (sum(rowSums(as.matrix(physeq1@otu_table@.Data)) == 0) > 0){
        physeq1@otu_table@.Data <- data.frame(physeq1@otu_table@.Data[-which
                                                                      (rowSums(as.matrix(physeq1@otu_table@.Data)) == 0),])
      }  
      pca.tmp <- prcomp(t(physeq1@otu_table@.Data), scale = TRUE)
      
    } else  {
      physeq2 <- tax_glom(physeq1, input$taxl.pca)
      if (sum(rowSums(as.matrix(physeq2@otu_table@.Data)) == 0) > 0){
        physeq2@otu_table@.Data <- data.frame(physeq2@otu_table@.Data[-which
                                                                      (rowSums(as.matrix(physeq2@otu_table@.Data)) == 0),])
      }  
      pca.tmp <- prcomp(t(physeq2@otu_table@.Data), scale = TRUE)
      
    }
    
    table.output.pca <- t(summary(pca.tmp)$importance)
    table.output.pca[,2] <- scales::percent(as.numeric(table.output.pca[,2]))
    table.output.pca[,3] <- scales::percent(as.numeric(table.output.pca[,3]))
    DT::datatable(table.output.pca)
    
  })

  ### PCoA
  plotPCoAPlotlyServer <- function(){
    physeq1 <- shinyInput$pstat
    physeq1 <- phyloseq(otu_table(physeq1), phy_tree(physeq1),
                        tax_table(physeq1), sample_data(physeq1))
    if (input$taxl.pca=="no rank")  {
      plotPCoAPlotly(physeq.input = physeq1, 
                     condition.vec = physeq1@sam_data[[input$select_pca_condition]],
                     condition.name = input$select_pca_condition,
                     method = input$pcoa.method,
                     pc.a = paste("Axis", input$xcol.new, sep = "."),
                     pc.b = paste("Axis", input$ycol.new, sep = "."),
                     columnTitle = paste("PCoA with colorbar representing", 
                                         input$select_pca_condition, sep = " "))
    } else  {
      physeq2 <- tax_glom(physeq1, input$taxl.pca)
      plotPCoAPlotly(physeq.input = physeq2, 
                     condition.vec = physeq2@sam_data[[input$select_pca_condition]],
                     condition.name = input$select_pca_condition,
                     method = input$pcoa.method,
                     pc.a = paste("Axis", input$xcol.new, sep = "."),
                     pc.b = paste("Axis", input$ycol.new, sep = "."),
                     columnTitle = paste("PCoA with colorbar representing", 
                                         input$select_pca_condition, sep = " "))}
  }
  
  
  output$pcoa.plotly <- renderPlotly({
    plotPCoAPlotlyServer()
  })

  observeEvent(input$download_pcoa,{
    if (!require("webshot")) install.packages("webshot")
    tmpFile <- tempfile(pattern = "PCoA_", fileext = ".pdf")
    export(plotPCoAPlotlyServer(), file = tmpFile)
    browseURL(tmpFile)}
  )
  
  
  getOrdPCoA <- function(){
    physeq1 <- shinyInput$pstat
    if (input$taxl.pca=="no rank")  {
      #test and fix the constant/zero row
      if (sum(rowSums(as.matrix(physeq1@otu_table@.Data)) == 0) > 0){
        physeq1@otu_table@.Data <- data.frame(physeq1@otu_table@.Data[-which
                                                                      (rowSums(as.matrix(physeq1@otu_table@.Data)) == 0),])
      }
      Dist.tmp <- phyloseq::distance(physeq1, method = input$select_beta_div_method)
      ord.tmp <- ordinate(physeq1, method = "PCoA", distance = Dist.tmp)
      #cat(dim(physeq1@otu_table))
      return(ord.tmp$values)

    } else  {
      physeq2 <- tax_glom(physeq1, input$taxl.pca)
      if (sum(rowSums(as.matrix(physeq2@otu_table@.Data)) == 0) > 0){
        physeq2@otu_table@.Data <- data.frame(physeq2@otu_table@.Data[-which
                                                                      (rowSums(as.matrix(physeq2@otu_table@.Data)) == 0),])
      }
      Dist.tmp <- phyloseq::distance(physeq2, method = input$select_beta_div_method)
      ord.tmp <- ordinate(physeq2, method = "PCoA", distance = Dist.tmp)
      return(ord.tmp$values)
        
    }
  }
  
  
  # interactive PCA table
  output$PCoAtable <- DT::renderDataTable({
    physeq1 <- shinyInput$pstat
    ord <- getOrdPCoA()
    df.output <- ord[,c(1,3,5)]
    colnames(df.output) <- c("eigenvalue", "variance explained", "cumulative variance")
    rownames(df.output) <- paste("Axis", 1:nrow(df.output), sep = ".")
    df.output[,2] <- scales::percent(as.numeric(df.output[,2]))
    df.output[,3] <- scales::percent(as.numeric(df.output[,3]))
    DT::datatable(df.output)
    
  })
  
  
  ## New DA analysis section
  output$DeSeq2Table.new <- DT::renderDataTable({
    physeq1 <- shinyInput$pstat
    if (input$taxl.da !="no rank"){
      physeq1 <- tax_glom(physeq1, input$taxl.da)
    }
    physeq1 <- prune_samples(sample_sums(physeq1) > input$da.count.cutoff, physeq1)
    diagdds = phyloseq_to_deseq2(physeq1, as.formula(paste("~",input$da.condition, sep = " ")))
    # calculate geometric means prior to estimate size factors
    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    geoMeans = apply(counts(diagdds), 1, gm_mean)
    diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
    diagdds = DESeq(diagdds, fitType="local")
    res = results(diagdds)
    res = res[order(res$padj, na.last=NA), ]
    sigtab = res[(res$padj < input$padj.cutoff), ]
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq1)[rownames(sigtab), ], "matrix"))
    sigtab$padj <- as.numeric(formatC(sigtab$padj, format = "e", digits = 2))
    sigtab$log2FoldChange <- as.numeric(formatC(sigtab$log2FoldChange, format = "e", digits = 2))
    DT::datatable(sigtab[,-c(1,3,4,5)])
    
  })

  
  
  
  # interactive Differential Expression boxplot
  BP <- reactive({
    findTaxCountDataDE()
    physeq1 <- shinyInput$pstat
    Inbatch <- sample_data(physeq1)[[input$secondary]]
    shinyInput <- getShinyInput()
    dat <- shinyInput$taxcountdata
    dat <- findNormalizedCount()
    dat <- apply(dat, 1:2, FUN = function(x) {
      ifelse(is.null(x) || is.na(x) || is.nan(x), 0, x)
    })
    # lcpm <- log2CPM(shinyInput$taxcountdata)
    lcpm <- log2CPM(dat)
    lcounts <- lcpm$y
    dat <- lcounts
    batch1 <- as.factor(Inbatch)
    batch2 <- split(which(Inbatch == batch1), batch1)
    batch3 <- unlist(lapply(1:length(batch2),
                            function(x) batch2[[x]][1:input$noSamples]))
    dat1 <- dat[, batch3]
    colnames(dat1) <- seq(1:ncol(dat))[batch3]
    dat1
  })
  DE <- reactive({
    physeq1 <- shinyInput$pstat
    Incondition <- sample_data(physeq1)[[input$primary]]
    findTaxCountDataDE()
    shinyInput <- getShinyInput()
    dat <- shinyInput$taxcountdata
    dat <- findNormalizedCount()
    dat <- apply(dat, 1:2, FUN = function(x) {
      ifelse(is.null(x) || is.na(x) || is.nan(x), 0, x)
    })
    # lcpm <- log2CPM(shinyInput$taxcountdata)
    lcpm <- log2CPM(dat)
    lcounts <- lcpm$y
    dat <- lcounts
    cond1 <- as.factor(Incondition)
    cond2 <- split(which(Incondition == cond1), cond1)
    cond3 <- unlist(lapply(1:length(cond2),
                           function(x) cond2[[x]][1:input$ncSamples]))
    dat1 <- dat[, cond3]
    colnames(dat1) <- seq(1:ncol(dat))[cond3]
    dat1
  })
  diffex_bp <- reactive({
    physeq1 <- shinyInput$pstat
    Incondition <- sample_data(physeq1)[[input$primary]]
    Inbatch <- sample_data(physeq1)[[input$secondary]]
    if (input$sortbybatch) {
      batch4 <- split(Inbatch, as.factor(Inbatch))
      batch5 <- unlist(lapply(1:length(batch4),
                              function(x) batch4[[x]][1:input$noSamples]))
      dat1 <- BP()
      dat2 <- melt(as.data.frame(dat1), measure.var = colnames(dat1))
      dat2$batch <- as.factor(unlist(lapply(1:length(batch5),
                                            function(x) rep(batch5[x], nrow(dat1)))))
      dat2$condition <- as.factor(unlist(lapply(as.numeric(colnames(dat1))
                                                , function(x) rep(Incondition[x], nrow(dat1)))))
      dat2$samples <- unlist(lapply(seq(ncol(dat1)),
                                    function(x) rep(seq(ncol(dat1))[x], nrow(dat1))))
    } else {
      cond4 <- split(Incondition,
                     as.factor(Incondition))
      cond5 <- unlist(lapply(1:length(cond4),
                             function(x) cond4[[x]][1:input$ncSamples]))
      dat1 <- DE()
      dat2 <- melt(as.data.frame(dat1), measure.var = colnames(dat1))
      dat2$condition <- as.factor(unlist(lapply(1:length(cond5),
                                                function(x) rep(cond5[x], nrow(dat1)))))
      dat2$batch <- as.factor(unlist(lapply(as.numeric(colnames(dat1)),
                                            function(x) rep(Inbatch[x], nrow(dat1)))))
      dat2$samples <- unlist(lapply(seq(ncol(dat1)),
                                    function(x) rep(seq(ncol(dat1))[x], nrow(dat1))))
    }
    dat2 %>% group_by(batch) %>% ggvis(~samples, ~value, fill =
                                         if (input$colbybatch) ~batch else ~condition) %>%
      layer_boxplots() %>%
      add_tooltip(function(dat2) { paste0("Sample: ",
                                          colnames(shinyInput$pstat@otu_table)[dat2$samples],
                                          "<br>", if (input$colbybatch) input$secondary
                                          else input$primary, ": ",
                                          if (input$colbybatch) dat2$batch else dat2$condition)
      }, "hover") %>%
      add_axis("x", title = if (input$sortbybatch)
        paste(input$noSamples, "Sample(s) Per ", input$secondary,
              sep = " ")
        else
          paste(input$ncSamples, "Sample(s) Per ", input$primary,
                sep=" "),
        properties = axis_props(title = list(fontSize = 15),
                                labels = list(fontSize = 5, angle = 90))) %>%
      add_axis("y", title = "Transformed Abundance", properties = axis_props(title =
                                                                               list(fontSize = 15),labels = list(fontSize = 10))) %>%
      add_legend("fill", title = if (input$colbybatch)
        input$secondary else input$primary,
        properties = legend_props(title =
                                    list(fontSize = 15), labels = list(fontSize = 10))) %>%
      set_options(width="auto", height="auto")
  })
  diffex_bp %>% bind_shiny("DiffExPlot")
  output$DEsummary <- renderPrint({
    if (input$sortbybatch) {
      summary(BP())
    } else {
      summary(DE())
    }
  })
  
  output$DEtable <- renderTable({
    if (input$sortbybatch) {
      BP()
    } else {
      DE()
    }
  })
  
  DELimnorm <- reactiveValues()
  observeEvent(c(input$taxlde, input$apply), {
    findAllTaxData(input$taxlde)
    dat <- findNormalizedCount()
    lcpm <- log2CPM(dat)
    lcounts <- lcpm$y
    DELimnorm$data <- lcounts
  })
  output$LimmaTable <- renderTable({
    shinyInput <- getShinyInput()
    physeq1 <- shinyInput$pstat
    Incondition <- sample_data(physeq1)[[input$primary]]
    Inbatch <- sample_data(physeq1)[[input$secondary]]
    pdata <- data.frame(Inbatch, Incondition)
    mod <- model.matrix(if (input$primary == input$secondary)
      ~as.factor(Incondition)
      else ~as.factor(Incondition) + ~as.factor(Inbatch), data = pdata)
    dat1 <- DELimnorm
    fit <- lmFit(dat1$data, mod)
    fit2 <- eBayes(fit)
    ncond <- nlevels(as.factor(Incondition))
    limmaTable <- topTable(fit2, coef = 2:ncond, number = input$noTaxons)
    for (j in 2:ncond)  {
      colnames(limmaTable)[j-1] <- paste("Primary Covariate: ",
                                         levels(as.factor(Incondition))[j], " (logFC)", sep='')
    }
    limmaTable
  }, rownames = TRUE)
  
  DEnorm <- reactiveValues()
  observeEvent(c(input$taxlde, input$apply),  {
    findAllTaxData(input$taxlde)
    DEnorm$data <- findNormalizedCount()
  })
  output$EdgeRTable <- renderTable({
    shinyInput <- getShinyInput()
    physeq1 <- shinyInput$pstat
    Incondition <- sample_data(physeq1)[[input$primary]]
    Inbatch <- sample_data(physeq1)[[input$secondary]]
    pdata <- data.frame(Inbatch, Incondition)
    group <- factor(Incondition)
    mod <- model.matrix(if (input$primary == input$secondary)
      ~as.factor(Incondition)
      else ~as.factor(Incondition) + ~as.factor(Inbatch), data = pdata)
    dat1 <- DEnorm
    y <- DGEList(counts=dat1$data,group=group)
    #y <- calcNormFactors(y)
    y <- estimateGLMCommonDisp(y,mod)
    y <- estimateGLMTrendedDisp(y,mod)
    y <- estimateGLMTagwiseDisp(y,mod)
    fit <- glmFit(y,mod)
    ncond <- nlevels(group)
    lrt <- glmLRT(fit,coef=2:ncond)
    edgeRTags <- topTags(lrt, n=input$noTaxons)
    edgeRTable <- edgeRTags$table
    for (j in 2:ncond)  {
      colnames(edgeRTable)[j-1] <- paste("Primary Covariate: ",
                                         levels(group)[j], " (logFC)", sep='')
    }
    edgeRTable
  }, rownames = TRUE)
  
  output$DeSeq2Table <- renderTable({
    shinyInput <- getShinyInput()
    physeq1 <- shinyInput$pstat
    Incondition <- sample_data(physeq1)[[input$primary]]
    Inbatch <- sample_data(physeq1)[[input$secondary]]
    pdata <- data.frame(Inbatch, Incondition)
    group <- factor(Incondition)
    design <- if (input$primary == input$secondary)
      ~as.factor(Incondition)
    else ~as.factor(Incondition) + ~as.factor(Inbatch)
    fCondition <- factor(Incondition)
    fBatch <- factor(Inbatch)
    design <- ~fCondition + ~fBatch
    mod <- model.matrix(design, data = pdata)
    dat1 <- DEnorm
    dat1$data <- apply(as.matrix(dat1$data),1:2, as.integer)
    coldata <- data.frame(row.names=colnames(dat1$data), fCondition, fBatch)
    dds <- DESeqDataSetFromMatrix(countData=dat1$data, colData=coldata,
                                  design=design)
    
    #dds <- DESeq(dds)
    sizeFactors(dds) <- rep(1, dim(dat1$data)[2])
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
    
    cvec <- rep(0, length(resultsNames(dds)))
    cvec[2:ncond] <- 1
    res <- results(dds, contrast=cvec)
    ## Order by adjusted p-value
    DeSeq2table <- res[order(res$padj), ]
    DeSeq2table <- DeSeq2table[1:input$noTaxons,]
    #DeSeq2table <- table(res$padj<0.05)
    # for (j in 2:ncond)  {
    #     colnames(DeSeq2table)[j-1] <- paste("Primary Covariate: ",
    #         levels(group)[j], " (logFC)", sep='')
    # }
    DeSeq2table
  }, rownames = TRUE)
  
  output$confRegion <- renderPlot({
    p1 <- shinyInput$pstat@otu_table[input$taxon1, input$sample]
    if (p1 <= 0) p1 <- 1
    p2 <- shinyInput$pstat@otu_table[input$taxon2, input$sample]
    if (p2 <= 0) p2 <- 1
    size <- sum(shinyInput$pstat@otu_table[,input$sample])
    plotConfRegion(p1, p2, size, uselogit=input$uselogit)
  })
  
  #Time Series
  output$Allustax <- renderUI({
    checkboxGroupInput(inputId="Allustax", label="Taxa of interest ",
                       choices = as.character(unique(unlist((
                         shinyInput$pstat@tax_table)[,input$Alluglom]))))
    
  })
  
  Alluvialdata <- reactive({
    if(input$Allurar==T){
      DR<-rarefy_even_depth(pstat,
                            sample.size =min(colSums(otu_table(pstat))),
                            replace=FALSE, rngseed=T)
    }else{
      DR<-pstat
    }
    
    if(is.null(input$Allusset) ||
       is.null(input$Alluglom) ||
       is.null(input$Allustax)){
      return()
    }
    tryCatch({
      tg<-paste("tax_glom(DR, taxrank='",input$Alluglom,"')",sep ='')
      glom = eval(parse(text=tg))
      tg<-paste("subset_taxa(glom, tax_table(glom)[,",input$Alluglom,"]
                %in% ",input$Allustax ,")",sep ='')
      glom = eval(parse(text=tg))
    },error=function(cond){
      return()
    })
    glom<-transform_sample_counts(glom, function(x) x / sum(x) )
    glom_time <- sample_data(glom)[,input$Allusset]
    glom_time$sample = rownames(glom_time)
    rownames(glom_time) = NULL
    glom_otu <- t(otu_table(glom))
    cbind(glom_time, glom_otu) -> glom_otu_time
    row.names(glom_otu_time) <- NULL
    glom_otu_time$sample <- NULL
    glom_otu_time <- plyr::ddply(glom_otu_time,input$Allusset,
                                 plyr::numcolwise(mean))
    tryCatch({
      colnames(glom_otu_time) <- c(input$Allusset,
                                   tax_table(glom)[as.matrix(colnames(
                                     glom_otu_time)[-1]),input$Alluglom])
      glom_otu_time_melted <- melt(glom_otu_time,
                                   id.vars=c(input$Allusset), measure.vars=input$Allustax,
                                   variable.name="taxa", value.name="proportion")
    },error=function(cond){
      return()
    })
    if(!exists("glom_otu_time_melted")){
      return()
    }
    glom_otu_time_melted$proportion <- glom_otu_time_melted$proportion*100
    glom_otu_time_melted <- glom_otu_time_melted[c(2,1,3)]
    glom_otu_time_melted["proportion"]<-rapply(
      glom_otu_time_melted["proportion"],
      f=function(x) ifelse(is.nan(x),0,x), how="replace" )
    return(glom_otu_time_melted)
  })
  
  output$TimePlotVisu<- renderPlot({
    if(is.null(Alluvialdata())){
      return()
    }
    alluvial_ts(Alluvialdata(), wave = .4, ygap = 2,plotdir = 'centred',
                alpha=.9, rankup = FALSE, grid = TRUE, grid.lwd = 5, xmargin = 0.2,
                lab.cex = 1, xlab = input$Allusset, ylab = '', border = NA,
                axis.cex = 1, title = '')
  })
  
  output$downloadAlluvialPlot <- downloadHandler(
    filename = function() { paste('Alluvialplot', '.pdf', sep='') },
    content = function(file) {
      pdf(file, width = 18, height = 10)
      alluvial_ts(Alluvialdata(), wave = .4, ygap = 2,plotdir = 'centred',
                  alpha=.8, rankup = FALSE, grid = TRUE, grid.lwd = 5,
                  xmargin = 0.2, lab.cex = 1, xlab = input$Allusset, ylab = '',
                  border = NA,axis.cex = 1, title = '')
      dev.off()
    }
  )
  
  callModule( coreOTUModule, "coreOTUModule", pstat )
})
