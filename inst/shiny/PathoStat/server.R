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
    
    vals <- reactiveValues(
        shiny.input = getShinyOption("pathostat.shinyInput"),
        taxdata = NULL,
        taxcountdata = NULL
    )

    observeEvent(input$uploadDataCount,{
        withBusyIndicatorServer("uploadDataCount", {
            
        
        df.input <- read.csv(input$countsfile$datapath,
                             header = input$header.count,
                             row.names = 1,
                             stringsAsFactors = FALSE,
                             sep = input$sep.count)
        cat(dim(df.input))
        df.meta.input <- read.csv(input$annotfile.count$datapath,
                                  header = input$header.count,
                                  sep = input$sep.count,
                                  row.names=1,
                                  stringsAsFactors=FALSE)
        #test and fix the constant/zero row
        row.remove.index <- c()
        if (sum(rowSums(as.matrix(df.input)) == 0) > 0){
            row.remove.index <- which(rowSums(as.matrix(df.input)) == 0)
            df.input <- df.input[-row.remove.index,]
        }
        ids <- rownames(df.input)
        tids <- unlist(lapply(ids, FUN = grepTid))
        taxonLevels <- findTaxonomy(tids)
        taxmat <- findTaxonMat(ids, taxonLevels)
        
        OTU <- otu_table(df.input, taxa_are_rows = TRUE)
        TAX <- tax_table(taxmat)
        physeq <- phyloseq(OTU, TAX)
        
        # change NA/NULL to 0
        # remove variables with identical values
        col.remove.index <- c()
        for (i in 1:ncol(df.meta.input)){
            if(length(unique(df.meta.input[,i])) < 2){
                col.remove.index <- c(col.remove.index, i)
            } 
        }
        if (!is.null(col.remove.index)){
            df.meta.input <- df.meta.input[,-col.remove.index]  
        }
        
        sampledata = sample_data(data.frame(df.meta.input))
        random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=
                                taxa_names(physeq))
        physeq1 <- merge_phyloseq(physeq, sampledata, random_tree)
        pstat <- pathostat1(physeq1)
        shinyInput <- list(pstat = pstat)
        vals$shiny.input <- shinyInput
        })
        
    })
    
    observeEvent(input$uploadDataPs, {
        withBusyIndicatorServer("uploadDataPs", {
 
                df.path.vec <- c()
                df.name.vec <- c()
                for(i in 1:length(input$countsfile.pathoscope[,1])){
                    df.path.vec[i] <- input$countsfile.pathoscope[[i, 'datapath']]
                    df.name.vec[i] <- input$countsfile.pathoscope[[i, 'name']]
                }
                
                cat(df.name.vec)
                datlist <- readPathoscopeData(input_dir, pathoreport_file_suffix, 
                                              use.input.files = TRUE,
                                              input.files.path.vec = df.path.vec,
                                              input.files.name.vec = df.name.vec)
                countdat <- datlist$countdata

                df.meta.input <- read.csv(input$annotfile.ps$datapath,
                                          header = input$header.ps,
                                          sep = input$sep.ps,
                                          row.names=1,
                                          stringsAsFactors=FALSE)
                #cat(dim(df.meta.input))
                #test and fix the constant/zero row
                row.remove.index <- c()
                if (sum(rowSums(as.matrix(countdat)) == 0) > 0){
                    row.remove.index <- which(rowSums(as.matrix(countdat)) == 0)
                    countdat <- countdat[-row.remove.index,]
                }
                
                ids <- rownames(countdat)
                tids <- unlist(lapply(ids, FUN = grepTid))
                taxonLevels <- findTaxonomy(tids)
                taxmat <- findTaxonMat(ids, taxonLevels)
                
                #test and fix the constant/zero row
                if (!is.null(row.remove.index)){
                    taxmat <- taxmat[-row.remove.index,]
                }
                
                
                #cat(rownames(df.meta.input))
                
                OTU <- otu_table(countdat, taxa_are_rows = TRUE)
                TAX <- tax_table(taxmat)
                physeq <- phyloseq(OTU, TAX)

                
                # change NA/NULL to 0
                # remove variables with identical values
                col.remove.index <- c()
                for (i in 1:ncol(df.meta.input)){
                    if(length(unique(df.meta.input[,i])) < 2){
                        col.remove.index <- c(col.remove.index, i)
                    } 
                }
                if (!is.null(col.remove.index)){
                    df.meta.input <- df.meta.input[,-col.remove.index]  
                }
                
                
                
                sampledata = sample_data(data.frame(df.meta.input))
                random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=
                                        taxa_names(physeq))
                physeq1 <- merge_phyloseq(physeq, sampledata, random_tree)
                pstat <- pathostat1(physeq1)
                shinyInput <- list(pstat = pstat)
                #setShinyInput(shinyInput)
                vals$shiny.input <- shinyInput
                
        })
})

    
  # setInputs(FALSE)
  findAllTaxData <- function(taxonLevel) {
    #taxcountdata <- findTaxonLevelData(pstat@otu_table,
    #    shinyInput$taxonLevels, taxonLevel)
    
    #shinyInput <- getShinyInput()
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
        tid <- grepTid(names[i])
        labvec[i] <- paste0( "ti|", tid, "|", labvec[i])
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
  
  findTaxCountDataDE <- reactive({
    findAllTaxData(input$taxlde)
    vals$taxcountdata
  })
  
  findNormalizedCount <- function() {
    if(input$norm == 'Quantile coreOTU Normalization')  {
      dat <- coreOTUQuantile(vals$taxcountdata, input$otuthreshold,
                             input$prevalence)
    } else if(input$norm == 'Library Size Scaling')  {
      dat <- sizeNormalize(vals$taxcountdata)
    } else  {
      dat <- coreOTUNormalize(vals$taxcountdata, input$ebweight,
                              input$otuthreshold, input$prevalence)
    }
    dat
  }
  
  ### data input summary
  
  output$contents.count <- renderTable({
      
      # input$file1 will be NULL initially. After the user selects
      # and uploads a file, head of that data file by default,
      # or all rows if selected, will be shown.
      
      if (!is.null(input$countsfile)){
          req(input$countsfile)
          df <- read.csv(input$countsfile$datapath,
                         header = input$header.count,
                         sep = input$sep.count)
          if (ncol(df) < 4){
              return(head(df)) 
          } else{
              return(head(df[,1:3]))
          }
      } else if (!is.null(input$countsfile.pathoscope)){
          req(input$countsfile.pathoscope)
          df <- read.csv(input$countsfile.pathoscope[[1, 'datapath']],
                         skip = 1,
                         header = TRUE,
                         sep = input$sep.ps)
          if (ncol(df) < 4){
              return(head(df)) 
          } else{
              return(head(df[,1:3]))
          }
      }

  })
  
  output$contents.meta <- renderTable({
      
      # input$file1 will be NULL initially. After the user selects
      # and uploads a file, head of that data file by default,
      # or all rows if selected, will be shown.
      
      if (!is.null(input$annotfile.count)){
          req(input$annotfile.count)
          df <- read.csv(input$annotfile.count$datapath,
                         header = input$header.count,
                         sep = input$sep.count)
          if (ncol(df) < 5){
              return(head(df)) 
          } else{
              return(head(df[,1:4]))
          }
      } else if (!is.null(input$annotfile.ps)){
          req(input$countsfile.pathoscope)
         
          df <- read.csv(input$annotfile.ps$datapath,
                         header = input$header.ps,
                         sep = input$sep.ps)
          if (ncol(df) < 5){
              return(head(df)) 
          } else{
              return(head(df[,1:4]))
          }
      }
  })
  
  
  
  
  
  
  
  
  
  
  tax_ra_bp <- reactive({
    if (is.null(input$taxl)) {
      return()
    }
    if (input$uploadDataPs == TRUE | input$uploadDataCount == TRUE){
        cat("barplot update with new data!")
    }
     
    
    shinyInput <- vals$shiny.input 
    pstat <- shinyInput$pstat
    cat(dim(pstat@otu_table@.Data))
    taxdata <- findTaxData()
    #cat(dim(taxdata))
    dat <- melt(cbind(taxdata, ind = as.character(rownames(taxdata))),
                id.vars = c("ind"))
    
    if ((input$taxl == "no rank")){
      covariates.tmp <- colnames(sample_data(pstat))
      dat$condition.select <- rep(pstat@sam_data@.Data[[
        which(covariates.tmp %in% input$select_condition)]], 
        each = dim(taxdata)[1])
      dat <- dat[order(dat$condition.select),]
      dat$condition.select.id <- paste(dat$condition.select,
                                       as.character(dat$variable), sep = "-")
    }else{
      pstat.new <- tax_glom(pstat, input$taxl)
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
    shinyInput <- vals$shiny.input
    write.csv(shinyInput$taxdata, file)
  })
  
  output$downloadCountData <- downloadHandler(filename = function() {
    paste0("sample_data_count_", input$taxl, ".csv", sep = "")
  }, content = function(file) {
     shinyInput <- vals$shiny.input
    write.csv(shinyInput$taxcountdata, file)
  })
  
  output$confRegion <- renderPlot({
    shinyInput <- vals$shiny.input
    p1 <- shinyInput$pstat@otu_table[input$taxon1, input$sample]
    if (p1 <= 0) p1 <- 1
    p2 <- shinyInput$pstat@otu_table[input$taxon2, input$sample]
    if (p2 <= 0) p2 <- 1
    size <- sum(shinyInput$pstat@otu_table[,input$sample])
    plotConfRegion(p1, p2, size, uselogit=input$uselogit)
  })
  
  #Time Series
  output$Allustax <- renderUI({
      shinyInput <- vals$shiny.input
    checkboxGroupInput(inputId="Allustax", label="Taxa of interest ",
                       choices = as.character(unique(unlist((
                         shinyInput$pstat@tax_table)[,input$Alluglom]))))
    
  })
  
  findPhyseqData <- function() {
      shinyInput <- vals$shiny.input
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
    shinyInput <- vals$shiny.input
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
      shinyInput <- vals$shiny.input
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
      shinyInput <- vals$shiny.input
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
        shinyInput <- vals$shiny.input
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
      shinyInput <- vals$shiny.input
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
      shinyInput <- vals$shiny.input
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
      shinyInput <- vals$shiny.input
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
          dist.list.tmp <- dist.list[which(names(dist.list) != names(dist.list)[i])]
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
      shinyInput <- vals$shiny.input
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
        shinyInput <- vals$shiny.input
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
      shinyInput <- vals$shiny.input
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
      shinyInput <- vals$shiny.input
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
      shinyInput <- vals$shiny.input
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
  

  # independent heatmap plotting function in the server using specific data from input
  plotPCAPlotlyServer <- function(){
      shinyInput <- vals$shiny.input
    physeq1 <- shinyInput$pstat
    if (input$taxl.pca=="no rank")  {
      plotPCAPlotly(df.input = physeq1@otu_table@.Data, 
                    condition.color.vec = physeq1@sam_data[[input$select_pca_color]],
                    condition.color.name = input$select_pca_color,
                    condition.shape.vec = physeq1@sam_data[[input$select_pca_shape]],
                    condition.shape.name = input$select_pca_shape,                    
                    
                    pc.a = paste("PC", input$xcol.new, sep = ""),
                    pc.b = paste("PC", input$ycol.new, sep = ""),
                    columnTitle = paste("PCA with colorbar representing", 
                                        input$select_pca_color, sep = " "))
    } else  {
      physeq2 <- tax_glom(physeq1, input$taxl.pca)
      plotPCAPlotly(df.input = physeq2@otu_table@.Data, 
                    condition.color.vec = physeq2@sam_data[[input$select_pca_color]],
                    condition.color.name = input$select_pca_color,
                    condition.shape.vec = physeq1@sam_data[[input$select_pca_shape]],
                    condition.shape.name = input$select_pca_shape, 
                    pc.a = paste("PC", input$xcol.new, sep = ""),
                    pc.b = paste("PC", input$ycol.new, sep = ""),
                    columnTitle = paste("PCA with colorbar representing", 
                                        input$select_pca_color, sep = " "))}
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
      shinyInput <- vals$shiny.input
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
      shinyInput <- vals$shiny.input
    physeq1 <- shinyInput$pstat
    physeq1 <- phyloseq(otu_table(physeq1), phy_tree(physeq1),
                        tax_table(physeq1), sample_data(physeq1))
    if (input$taxl.pca=="no rank")  {
      plotPCoAPlotly(physeq.input = physeq1, 
                     condition.color.vec = physeq1@sam_data[[input$select_pca_color]],
                     condition.color.name = input$select_pca_color,
                     condition.shape.vec = physeq1@sam_data[[input$select_pca_shape]],
                     condition.shape.name = input$select_pca_shape,                      
                     method = input$pcoa.method,
                     pc.a = paste("Axis", input$xcol.new, sep = "."),
                     pc.b = paste("Axis", input$ycol.new, sep = "."),
                     columnTitle = paste("PCoA with colorbar representing", 
                                         input$select_pca_color, sep = " "))
    } else  {
      physeq2 <- tax_glom(physeq1, input$taxl.pca)
      plotPCoAPlotly(physeq.input = physeq2, 
                     condition.color.vec = physeq2@sam_data[[input$select_pca_color]],
                     condition.color.name = input$select_pca_color,
                     condition.shape.vec = physeq1@sam_data[[input$select_pca_shape]],
                     condition.shape.name = input$select_pca_shape,                      
                     method = input$pcoa.method,
                     pc.a = paste("Axis", input$xcol.new, sep = "."),
                     pc.b = paste("Axis", input$ycol.new, sep = "."),
                     columnTitle = paste("PCoA with colorbar representing", 
                                         input$select_pca_color, sep = " "))}
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
      shinyInput <- vals$shiny.input
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
      shinyInput <- vals$shiny.input
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
      shinyInput <- vals$shiny.input
    physeq1 <- shinyInput$pstat
    if (input$taxl.da !="no rank"){
      physeq1 <- tax_glom(physeq1, input$taxl.da)
    }
    physeq1 <- prune_samples(sample_sums(physeq1) > input$da.count.cutoff, physeq1)
    
    # deal with continuous covariates and multiple covariates
    # target condition is the last one in formula. 
    if (!is.null(input$da.condition.covariate)){
      for (i in 1:length(input$da.condition.covariate)){
        num.levels <- length(unique(sample_data(physeq1)[[input$da.condition.covariate[i]]]))
        if (num.levels >= 8){
          sam.index <- which(physeq1@sam_data@names %in% input$da.condition.covariate[i])
          pstat@sam_data@.Data[[sam.index]] <- cut(pstat@sam_data@.Data[[sam.index]], breaks = 3)
        }
      }
      
      diagdds = phyloseq_to_deseq2(physeq1, 
                                   as.formula(paste("~",
                                                    paste(
                                                      paste(input$da.condition.covariate, 
                                                            collapse = " + "), 
                                                      input$da.condition, 
                                                      sep = " + "), 
                                                    sep = " ")))
    } else{
      diagdds = phyloseq_to_deseq2(physeq1, as.formula(paste("~",input$da.condition, sep = " ")))
    }
    
    # calculate geometric means prior to estimate size factors
    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    geoMeans = apply(counts(diagdds), 1, gm_mean)
    diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
    diagdds = DESeq(diagdds, fitType="local")
    res = results(diagdds)
    res = res[order(res$padj, na.last=NA), ]
    
    if (nrow(res) != 0){
      sigtab = res[(res$padj < input$da.padj.cutoff), ]
      if (nrow(sigtab) == 0){
        DT::datatable(as.matrix("No differentially abundant items found!"))
      } else{
        sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq1)[rownames(sigtab), ], "matrix"))
        sigtab$padj <- as.numeric(formatC(sigtab$padj, format = "e", digits = 2))
        sigtab$log2FoldChange <- as.numeric(formatC(sigtab$log2FoldChange, format = "e", digits = 2))
        output$download_deseq_tb <- downloadHandler(
          filename = function() { paste('download_deseq2_table', '.csv', sep='') },
          content = function(file) {
            dist.mat <- as.matrix(sigtab[,-c(1,3,4,5)])
            write.csv(data.frame(dist.mat), file)
          }
        )
        DT::datatable(sigtab[,-c(1,3,4,5)])
  
      }

    }else{
      DT::datatable(as.matrix("No differentially abundant items found!"))
    }

    
  })

  # Presence-Absence Variance analysis
  output$pa.case <- renderPrint({
      shinyInput <- vals$shiny.input
    physeq1 <- shinyInput$pstat
    paste("Case level is: ", physeq1@sam_data[[input$pa.condition]][1], sep = "")
    
  })
  
  output$pa.test <- DT::renderDataTable({
      shinyInput <- vals$shiny.input
    physeq1 <- shinyInput$pstat
    if (input$taxl.pa !="no rank"){
      physeq1 <- tax_glom(physeq1, input$taxl.pa)
    }
    physeq1 <- prune_samples(sample_sums(physeq1) > input$pa.count.cutoff, physeq1)
    df.pam <- GET_PAM(physeq1@otu_table@.Data)
    if (input$pa.method == "Fisher Exact Test"){
      output.mat <- Fisher_Test_Pam(df.pam, 
                                    physeq1@sam_data[[input$pa.condition]], 
                                    input$pa.padj.cutoff)
      output$download_pa_test <- downloadHandler(
        filename = function() { paste(input$pa.method, '.csv', sep='') },
        content = function(file) {
          write.csv(output.mat, file)
        }
      )
      
      DT::datatable(output.mat)
    } else if(input$pa.method == "Chi-squared Test"){
      output.mat <- Chisq_Test_Pam(df.pam, 
                                   physeq1@sam_data[[input$pa.condition]], 
                                   input$pa.padj.cutoff)
      output$download_pa_test <- downloadHandler(
        filename = function() { paste(input$pa.method, '.csv', sep='') },
        content = function(file) {
          write.csv(output.mat, file)
        }
      )
      DT::datatable(output.mat)
    } else if (input$pa.method == "Mann-Whitney Test"){
      output.mat <- Wilcox_Test_df(physeq1@otu_table@.Data, 
                                   physeq1@sam_data[[input$pa.condition]], 
                                   input$pa.padj.cutoff)
      output$download_pa_test <- downloadHandler(
        filename = function() { paste(input$pa.method, '.csv', sep='') },
        content = function(file) {
          write.csv(output.mat, file)
        }
      )
      DT::datatable(output.mat)
    }
    

    
    
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
