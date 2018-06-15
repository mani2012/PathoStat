library(shiny)
library(ROCR)
library(ggvis)
library(d3heatmap)
library(reshape2)
library(DESeq2)
library(edgeR)
library(phyloseq)
library(stats)
library(PathoStat)
library(plotly)
library(webshot)
library(vegan)
library(dplyr)
library(ape)






# Converts decimal percentage to string with specified digits
pct2str <- function(v, digits=2) {sprintf(paste0('%.',digits,'f'), v*100)}

# get RA from counts
getRelativeAbundance <- function(df){
  ra.out <- apply(df, 2, function(x) round(x/sum(x), digits = 4))
  return(ra.out)
}

getLogCPM <- function(df){
  logCPM.out <- apply(df, 2, function(y) log10(y*1e6/sum(y) + 1))
  return(logCPM.out)
}




shinyServer(function(input, output, session) {

    vals <- reactiveValues(
        shiny.input = getShinyOption("pathostat.shinyInput"),
        shiny.input.backup = getShinyOption("pathostat.shinyInput"),
        taxdata = NULL,
        taxcountdata = NULL
    )

    updateTaxLevel <- function(){
      shinyInput <- vals$shiny.input
      pstat <- shinyInput$pstat

      updateSelectInput(session, "taxl",
                        choices = colnames(pstat@tax_table@.Data))
      updateSelectInput(session, "taxl.alpha",
                        choices = colnames(pstat@tax_table@.Data))
      updateSelectInput(session, "taxl.beta",
                        choices = colnames(pstat@tax_table@.Data))
      updateSelectInput(session, "taxl.pca",
                        choices = colnames(pstat@tax_table@.Data))
      updateSelectInput(session, "taxl.da",
                        choices = colnames(pstat@tax_table@.Data))
      updateSelectInput(session, "taxl.edger",
                        choices = colnames(pstat@tax_table@.Data))
      updateSelectInput(session, "taxl.pa",
                        choices = colnames(pstat@tax_table@.Data))
      updateSelectInput(session, "taxl",
                        choices = colnames(pstat@tax_table@.Data))
      updateSelectInput(session, "taxl",
                        choices = colnames(pstat@tax_table@.Data))
      updateSelectInput(session, "taxl",
                        choices = colnames(pstat@tax_table@.Data))
      updateSelectInput(session, "taxl",
                        choices = colnames(pstat@tax_table@.Data))
    }

    # update samples
    updateSample <- function(){
        shinyInput <- vals$shiny.input
        pstat <- shinyInput$pstat
        updateSelectInput(session, "filterSample",
                          choices = colnames(pstat@otu_table@.Data))
    }



    #Update covariate names
    updateCovariate <- function(){
        shinyInput <- vals$shiny.input
        pstat <- shinyInput$pstat
        covariates <- colnames(sample_data(pstat))
        # choose the covariates that has less than 8 levels
        covariates.colorbar <- c()
        for (i in 1:length(covariates)){
            num.levels <- length(unique(sample_data(pstat)[[covariates[i]]]))
            if (num.levels < 8){
                covariates.colorbar <- c(covariates.colorbar, covariates[i])
            }
        }
        # choose the covariates that has 2 levels
        covariates.two.levels <- c()
        for (i in 1:length(covariates)){
            num.levels <- length(unique(sample_data(pstat)[[covariates[i]]]))
            if (num.levels == 2){
                covariates.two.levels <- c(covariates.two.levels, covariates[i])
            }
        }

        updateSelectInput(session, "select_covariate_condition_biomarker",
                          choices = covariates)
        updateSelectInput(session, "select_single_species_condition",
                          choices = covariates.colorbar)
        updateSelectInput(session, "select_target_condition_biomarker",
                          choices = covariates.colorbar)
        updateSelectInput(session, "select_condition_sample_filter",
                          choices = c("Read Number", covariates))
        updateSelectInput(session, "select_condition_sample_filter_micro",
                          choices = c("Taxon elements number", covariates))
        updateSelectInput(session, "select_condition_sample_filter_sidebar",
                          choices = c("Read Number", covariates))
        updateSelectInput(session, "select_condition_sample_distribution",
                          choices = covariates)
        updateSelectInput(session, "select_condition",
                          choices = covariates)
        updateSelectInput(session, "select_heatmap_condition_1",
                          choices = covariates.colorbar)
        updateSelectInput(session, "select_heatmap_condition_2",
                          choices = covariates.colorbar)
        updateSelectInput(session, "select_alpha_div_condition",
                          choices = covariates.colorbar)
        updateSelectInput(session, "select_beta_condition",
                          choices = covariates.two.levels)
        updateSelectInput(session, "select_beta_heatmap_condition_1",
                          choices = covariates.colorbar)
        updateSelectInput(session, "select_beta_heatmap_condition_2",
                          choices = covariates.colorbar)
        updateSelectInput(session, "select_pca_color",
                          choices = covariates)
        updateSelectInput(session, "select_pca_shape",
                          choices = covariates.colorbar)
        updateSelectInput(session, "da.condition",
                          choices = covariates.colorbar)
        updateSelectInput(session, "edger.condition",
                          choices = covariates.colorbar)
        updateSelectInput(session, "da.condition.covariate",
                          choices = covariates)
        updateSelectInput(session, "pa.condition",
                          choices = covariates.colorbar)
    }

    observeEvent(input$uploadPathoStat,{
      withBusyIndicatorServer("uploadPathoStat", {
        if (input$rdtype == 'rda') {
          load(input$rdfile$datapath)
        }
        if (input$rdtype == 'rds') {
          pstat <- readRDS(input$rdfile$datapath)
        }
        shinyInput <- list(pstat = pstat)
        vals$shiny.input <- shinyInput
        vals$shiny.input.backup <- shinyInput
        # update ui
        updateCovariate()
        updateSample()

      })

    })


    observeEvent(input$uploadDataCount,{
        withBusyIndicatorServer("uploadDataCount", {

        df.input <- read.csv(input$countsfile$datapath,
                             header = input$header.count,
                             row.names = 1,
                             stringsAsFactors = FALSE,
                             sep = input$sep.count)
        #cat(dim(df.input))
        df.taxon.input <- read.csv(input$taxon.table$datapath,
                                  header = input$header.count,
                                  sep = input$sep.count,
                                  row.names= 1,
                                  stringsAsFactors=FALSE)

        df.meta.input <- read.csv(input$annotfile.count$datapath,
                                  header = input$header.count,
                                  sep = input$sep.count,
                                  row.names=input$metadata_sample_name_col_count,
                                  stringsAsFactors=FALSE)

        # choose only the samples in metadata that have counts data as well
        df.meta.input <- df.meta.input[match(colnames(df.input), rownames(df.meta.input)), ]


        #test and fix the constant/zero row
        row.remove.index <- c()
        if (sum(rowSums(as.matrix(df.input)) == 0) > 0){
            row.remove.index <- which(rowSums(as.matrix(df.input)) == 0)
            df.input <- df.input[-row.remove.index,]
        }


        OTU <- otu_table(df.input, taxa_are_rows = TRUE)
        TAX <- tax_table(as.matrix(df.taxon.input))
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
        vals$shiny.input.backup <- shinyInput
        # update ui
        updateCovariate()
        updateSample()
        updateTaxLevel()

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

                #cat(df.name.vec)
                datlist <- readPathoscopeData(input_dir, pathoreport_file_suffix = input$report_suffix,
                                              use.input.files = TRUE,
                                              input.files.path.vec = df.path.vec,
                                              input.files.name.vec = df.name.vec)
                countdat <- datlist$countdata

                df.meta.input <- read.csv(input$annotfile.ps$datapath,
                                          header = input$header.ps,
                                          sep = input$sep.ps,
                                          row.names=input$metadata_sample_name_col,
                                          stringsAsFactors=FALSE)

                # choose only the samples in metadata that have counts data as well
                df.meta.input <- df.meta.input[match(colnames(countdat), rownames(df.meta.input)), ]


                #cat(colnames(countdat))
                #cat(rownames(df.meta.input))
                #test and fix the constant/zero row
                row.remove.index <- c()
                if (sum(rowSums(as.matrix(countdat)) == 0) > 0){
                    row.remove.index <- which(rowSums(as.matrix(countdat)) == 0)
                    countdat <- countdat[-row.remove.index,]
                }

                ids <- rownames(countdat)
                tids <- unlist(lapply(ids, FUN = grepTid))
                #print(tids)
                taxonLevels <- findTaxonomy(tids)
                #print(str(taxonLevels))
                taxmat <- findTaxonMat(ids, taxonLevels)
                #print(taxmat)
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
                vals$shiny.input.backup <- shinyInput
                # update ui
                updateCovariate()
                updateSample()
        })
})


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

        updateCovariate()
        updateSample()
    })
    })

      observeEvent(input$resetSampleButton,{
        withBusyIndicatorServer("resetSampleButton", {
        vals$shiny.input <- vals$shiny.input.backup
        updateCovariate()
        updateSample()
    })
    })

      observeEvent(input$resetSampleButtonMicro,{
        withBusyIndicatorServer("resetSampleButtonMicro", {
        vals$shiny.input <- vals$shiny.input.backup
        updateTaxLevel()
    })
    })      

      output$download_rda <- downloadHandler(filename = function() {
        paste("pathostat", Sys.Date(), ".rda", sep="")
      }, content = function(file) {
        shinyInput <- vals$shiny.input
        pstat <- shinyInput$pstat
        save(pstat, file=file)
      })

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
    rownames(df.out) <- TranslateIdToTaxLevel(pstat, rownames(df.out), 
                                              input$taxlTable)
    df.out
  })
 
  findRATable <- reactive({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    if (input$taxlTable !="no rank"){
      pstat <- tax_glom(pstat, input$taxlTable)
    }
    df.ra <- getRelativeAbundance(pstat@otu_table@.Data)
    rownames(df.ra) <- TranslateIdToTaxLevel(pstat, rownames(df.ra), 
                                              input$taxlTable)
    df.ra
  })   


  output$filter_type <- reactive({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    if (input$select_condition_sample_filter_sidebar == "Read Number"){
      return("num.continuous")
    }

    variable.vec <- sample_data(pstat)[[
          which(pstat@sam_data@names == input$select_condition_sample_filter_sidebar)]]
    if (is.numeric(variable.vec)){
      # catogorical or not
        num.levels <- length(unique(variable.vec))
        if (num.levels < 6){
          return("cat")
        } else{
          return("num.continuous")
        }
    }else{
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

      if (input$select_condition_sample_filter_sidebar == "Read Number"){
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

      updateCovariate()
      updateSample()
    })
  })

  output$filter_cat_options <- renderUI({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    variable.vec <- sample_data(pstat)[[
      which(pstat@sam_data@names == input$select_condition_sample_filter_sidebar)]]
    filter.option.vec <- sort(unique(variable.vec))
    tagList(
      selectInput("cat_filter_options", 
                  "Keep these levels:", 
                  choices = filter.option.vec,
                  multiple = TRUE
                  )
    )
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

      updateCovariate()
      updateSample()
    })
  })





  output$sampleCountSum <- renderPlotly({
      shinyInput <- vals$shiny.input
      pstat <- shinyInput$pstat
      Sample_Name <- colnames(pstat@otu_table@.Data)
      Read_Number <- colSums(pstat@otu_table@.Data)
      data <- data.frame(Sample_Name, Read_Number, stringsAsFactors = FALSE)

      if (input$select_condition_sample_filter == "Read Number"){
          data$Sample_Name <- factor(data$Sample_Name,
                                     levels = unique(data$Sample_Name)[order(data$Read_Number,
                                                                             decreasing = FALSE)])
      } else{
          data$Sample_Name <- paste(as.character(pstat@sam_data@.Data
                                                 [[which(pstat@sam_data@names == input$select_condition_sample_filter)]]
                                                 ),data$Sample_Name, sep = "-")
          data$Sample_Name <- factor(data$Sample_Name,
                                     levels = unique(data$Sample_Name)
                                     [order(pstat@sam_data@.Data[[which(pstat@sam_data@names == input$select_condition_sample_filter)]],
                                            decreasing = FALSE)])
      }

      p <- plot_ly(data, x = ~Sample_Name, y = ~Read_Number, type = "bar", name = 'Sample read count') %>%
          layout(margin = list(b = 160))
      p
  })


  output$sampleTaxon <- renderPlotly({
      shinyInput <- vals$shiny.input
      pstat <- shinyInput$pstat
      Sample_Name <- colnames(pstat@otu_table@.Data)
      Taxon_Num <- apply(pstat@otu_table@.Data , 2, function(x) sum(x >= 1) )
      data <- data.frame(Sample_Name, Taxon_Num, stringsAsFactors = FALSE)



      if (input$select_condition_sample_filter_micro == "Taxon elements number"){
        data$Sample_Name <- factor(data$Sample_Name,
                                 levels = unique(data$Sample_Name)[order(data$Taxon_Num,
                                                                         decreasing = FALSE)])
      } else{
          data$Sample_Name <- paste(as.character(pstat@sam_data@.Data
                                                 [[which(pstat@sam_data@names == input$select_condition_sample_filter_micro)]]
                                                 ),data$Sample_Name, sep = "-")
          data$Sample_Name <- factor(data$Sample_Name,
                                     levels = unique(data$Sample_Name)
                                     [order(pstat@sam_data@.Data[[which(pstat@sam_data@names == input$select_condition_sample_filter_micro)]],
                                            decreasing = FALSE)])
      }
      p <- plot_ly(data, 
                   x = ~Sample_Name, 
                   y = ~Taxon_Num, 
                   type = "bar", 
                   name = 'Sample taxon elements number') %>%
          layout(margin = list(b = 160))
      p
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

  #Render summary table
  output$contents_summary <- renderTable({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    req(pstat)
    summarizeTable(pstat)
  })

  output$contents_summary_micro <- renderTable({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    req(pstat)
    summarizeTable(pstat)
  })
  



  output$sample_metadata_distribution <- renderPlotly({
      shinyInput <- vals$shiny.input
      pstat <- shinyInput$pstat

      Sample_Name <- colnames(pstat@otu_table@.Data)
      condition.target <- pstat@sam_data@.Data[[which(pstat@sam_data@names == input$select_condition_sample_distribution)]]
      data <- data.frame(Sample_Name, condition.target, stringsAsFactors = FALSE)
          data$Sample_Name <- paste(as.character(pstat@sam_data@.Data
                                                 [[which(pstat@sam_data@names == input$select_condition_sample_distribution)]]
          ),data$Sample_Name, sep = "-")
          data$Sample_Name <- factor(data$Sample_Name,
                                     levels = unique(data$Sample_Name)
                                     [order(pstat@sam_data@.Data[[which(pstat@sam_data@names == input$select_condition_sample_distribution)]],
                                            decreasing = FALSE)])


      p <- plot_ly(data, x = ~Sample_Name, y = ~condition.target, type = "bar", name = 'Sample distribution') %>%
          layout(margin = list(b = 160))
      p
  })

  ### data input summary

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

#### plot the single species boxplot
  output$single_species_ui <- renderUI({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    species.name.vec <- TranslateIdToTaxLevel(pstat, rownames(pstat@otu_table@.Data), input$taxl_single_species)
    tagList(
      selectInput("select_single_species_name_plot", "Select names (support multiple)", species.name.vec, selected = species.name.vec[1], multiple = TRUE)
    )
  })

  plotSingleSpeciesBoxplotServer <- eventReactive(input$boxplotButton,{
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat

    if (input$taxl_single_species !="no rank")  {
      pstat <- tax_glom(pstat, input$taxl_single_species)
    }




    if (length(input$select_single_species_name_plot) == 1){
      condition.vec <- pstat@sam_data@.Data[[which(pstat@sam_data@names == input$select_single_species_condition)]]
      # change microbe names to selected taxon level
      species.name.vec <- TranslateIdToTaxLevel(pstat, rownames(pstat@otu_table@.Data), input$taxl_single_species)
      #cat(paste("condition:", length(condition.vec)))
      read.num.condition.1 <- rep(0, length(condition.vec))

      if (input$ssv_format == "read count"){
        for (i in 1:nrow(pstat@otu_table@.Data)){
          if (species.name.vec[i] == input$select_single_species_name_plot){
            read.num.condition.1 <- read.num.condition.1 + pstat@otu_table@.Data[i,]
          }
        }
        df.tmp <- data.frame(condition = condition.vec, read_Number = read.num.condition.1)
        p <- plot_ly(df.tmp, y = ~read_Number, color = ~condition, type = "box")
      }else if(input$ssv_format == "log10 CPM"){
        df.cpm <- getLogCPM(pstat@otu_table@.Data)
        for (i in 1:nrow(df.cpm)){
          if (species.name.vec[i] == input$select_single_species_name_plot){
            read.num.condition.1 <- read.num.condition.1 + df.cpm[i,]
          }
        }
        df.tmp <- data.frame(condition = condition.vec, logCPM = read.num.condition.1)
        p <- plot_ly(df.tmp, y = ~logCPM, color = ~condition, type = "box")
      }else{
        df.ra <- getRelativeAbundance(pstat@otu_table@.Data)
        for (i in 1:nrow(df.ra)){
          if (species.name.vec[i] == input$select_single_species_name_plot){
            read.num.condition.1 <- read.num.condition.1 + df.ra[i,]
          }
        }
        df.tmp <- data.frame(condition = condition.vec, RA = read.num.condition.1)
        p <- plot_ly(df.tmp, y = ~RA, color = ~condition, type = "box")
      }

      ## plot
      p


    } else{
      condition.vec <- pstat@sam_data@.Data[[which(pstat@sam_data@names == input$select_single_species_condition)]]
      condition.vec.all <- c()
      species.vec <- c()
      read.num.condition.all <- c()
      species.name.vec <- TranslateIdToTaxLevel(pstat, rownames(pstat@otu_table@.Data), input$taxl_single_species)
      for (i in 1:length(input$select_single_species_name_plot)){
        condition.vec.all <- c(condition.vec.all, condition.vec)
        species.vec <- c(species.vec, rep(input$select_single_species_name_plot[i], length(condition.vec)))
        read.num.condition.tmp <- rep(0, length(condition.vec))
        for (j in 1:nrow(pstat@otu_table@.Data)){
          if (species.name.vec[j] == input$select_single_species_name_plot[i]){
            if (input$ssv_format == "read count"){
              read.num.condition.tmp <- read.num.condition.tmp + pstat@otu_table@.Data[j,]
            }else if(input$ssv_format == "log10 CPM"){
              df.cpm <- getLogCPM(pstat@otu_table@.Data)
              read.num.condition.tmp <- read.num.condition.tmp + df.cpm[j,]
            }else{
              df.ra <- getRelativeAbundance(pstat@otu_table@.Data)
              read.num.condition.tmp <- read.num.condition.tmp + df.ra[j,]
            }


          }
        }
        read.num.condition.all <- c(read.num.condition.all, read.num.condition.tmp)

      }

      if (input$ssv_format == "read count"){
      df.tmp <- data.frame(condition = condition.vec.all, read_Number = read.num.condition.all, name = species.vec)
      p <- plot_ly(df.tmp, x = ~name, y = ~read_Number, color = ~condition, type = "box") %>%
        layout(boxmode = "group")
      suppressMessages(p)
      }else if(input$ssv_format == "log10 CPM"){
        df.tmp <- data.frame(condition = condition.vec.all, logCPM = read.num.condition.all, name = species.vec)
        p <- plot_ly(df.tmp, x = ~name, y = ~logCPM, color = ~condition, type = "box") %>%
          layout(boxmode = "group")
        suppressMessages(p)
      }else{
        df.tmp <- data.frame(condition = condition.vec.all, RA = read.num.condition.all, name = species.vec)
        p <- plot_ly(df.tmp, x = ~name, y = ~RA, color = ~condition, type = "box") %>%
          layout(boxmode = "group")
        suppressMessages(p)
      }

    }



  })


      output$single_species_boxplot <- renderPlotly({
        # suppress warnings
        storeWarn<- getOption("warn")
        options(warn = -1)
        plott <- plotSingleSpeciesBoxplotServer()
        #restore warnings, delayed so plot is completed
        shinyjs::delay(expr =({
          options(warn = storeWarn)
        }) ,ms = 100)

        plott
      })

  # Used to dynamically generate selectable organisms based on taxlev
  output$order_organisms <- renderUI({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat    
    taxcountdata <- pstat@otu_table
    taxdata <- findRAfromCount(taxcountdata)
    taxcountdata <- data.frame(taxcountdata)
    taxdata <- data.frame(taxdata)
    taxlevs <- as.data.frame(pstat@tax_table)
    selectizeInput('order_organisms', label='Order Samples by Organism:',choices = unique(taxlevs[[input$taxl]]), multiple=TRUE)
  })

  output$TaxRelAbundancePlot <- renderPlotly({
    if (is.null(input$taxl)) {
      return()
    }
    if (input$uploadDataPs == TRUE | input$uploadDataCount == TRUE){
      cat("barplot update with new data!")
    }

    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    # In the original plot.relative.abundance function()
    taxcountdata <- pstat@otu_table
    taxdata <- findRAfromCount(taxcountdata)
    taxcountdata <- data.frame(taxcountdata)
    taxdata <- data.frame(taxdata)
    taxlevs <- as.data.frame(pstat@tax_table)

    # Sum by taxon level
    df = as.data.frame(taxdata)
    df$newlev = taxlevs[[input$taxl]]
    df.melt = melt(df, id.vars = c("newlev"))
    df.melt.agg = aggregate(.~variable+newlev, data=df.melt, FUN=sum)
    df.ra = dcast(df.melt.agg, variable~newlev)
    rownames(df.ra) = df.ra$variable
    df.ra$variable = NULL

    # Reorder by most prominent organisms
    df.ra = df.ra[,order(colSums(df.ra))]

    # Put selected organisms first
    if (!is.null(input$order_organisms)) {
      organisms.order = c(setdiff(colnames(df.ra), input$order_organisms), rev(input$order_organisms))
      df.ra = df.ra[,organisms.order]
    }

    # If any conditions are selected make a side bar
    if (!is.null(input$select_conditions)) {
      samdata <- as.data.frame(pstat@sam_data)
      samdata = samdata[,input$select_conditions]

      # Order samples by conditions if not by organisms
      if (input$sort_samples_by == "conditions") {
        for (i in ncol(samdata):1) {
          samdata = samdata[order(samdata[[i]]),]
        }
        # Reorder stacked barplot
        df.ra = df.ra[order(match(rownames(df.ra), rownames(samdata))),]
      }

      # Retain hover-text information before conditions are factorized
      hover.txt <- c()
      for (i in 1:ncol(samdata)) {
        hover.txt <- cbind(hover.txt, samdata[[i]])
      }

      # Plotly | Heatmap ---------------------------------------
      samdata[] <- lapply(samdata, factor)
      m = data.matrix(samdata)
      m.row.normalized = apply(m, 2, function(x)(x-min(x))/(max(x)-min(x)))
      hm <- plot_ly(
        x = colnames(m),
        y = rownames(m),
        z = m.row.normalized, 
        type = "heatmap",
        showscale=FALSE,
        hoverinfo = "x+y+text",
        text=hover.txt
        ) %>%
        layout(xaxis = list(title = "",
                            tickangle = -45),
               yaxis = list(showticklabels = FALSE,
                            ticks = ""))
    }
    # --------------------------------------------------------------

    # Order samples by organisms if not by conditons
    if (input$sort_samples_by == "organisms") {
      for (i in 1:ncol(df.ra)) {
        df.ra = df.ra[order(df.ra[,i]),]
      }
    }

    # Plotly | Stacked Bar Plots ---------------------------------
    df.plot = df.ra
    df.plot$samples <- rownames(df.plot)
    p <- plot_ly(df.plot,
                 y = ~samples, 
                 x = df.plot[[colnames(df.plot)[1]]], 
                 type = 'bar', orientation = 'h', 
                 name = substr(colnames(df.plot)[1], 1, 40)) %>%
      layout(font = list(size = 10),
             yaxis = list(title = '', 
                          tickmode = "array",
                          tickvals = rownames(df.plot),
                          showticklabels = FALSE,
                          categoryorder = 'trace'),
             xaxis = list(title = 'Relative Abundance'),
             barmode = 'stack',
             showlegend = input$show_legend)
    for (i in 2:(ncol(df.plot)-1)) {
      p <- add_trace(p, x = df.plot[[colnames(df.plot)[i]]], 
                     name = substr(colnames(df.plot)[i], 1, 40))
    }
    # --------------------------------------------------------------

    # Create a multiplot if any conditions are seleceted
    if (!is.null(input$select_conditions)) {
      multi.plot <- subplot(hm, p, widths = c(0.1,  0.9))
      multi.plot
    } else {
      p
    }
  })

  output$TaxRAsummary <- renderPrint({
    summary(findTaxData())
  })

  # These are options for rendering datatables
  dtopts <- list(scrollX=TRUE, paging=TRUE)

  # Render relative abundance table
  output$TaxRAtable <- DT::renderDataTable(findRATable(),
                                           options=dtopts, rownames=F)


  # Render count table
  output$TaxCountTable <- DT::renderDataTable(findCountTable(),
  options=dtopts, 
  rownames=F)

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

  # heatmap
  output$single_species_ui_new <- renderUI({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    species.name.vec <- TranslateIdToTaxLevel(pstat, rownames(pstat@otu_table@.Data), input$taxl)
    tagList(
      selectInput("select_single_species_name_plot_new", "Select names (support multiple)", species.name.vec, selected = species.name.vec[1], multiple = TRUE)
    )
  })

  # get smaller df
  getDfPlot <- function(){
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat

    if (input$taxl !="no rank")  {
      pstat <- tax_glom(pstat, input$taxl)
    }
      # change microbe names to selected taxon level
      species.name.vec <- TranslateIdToTaxLevel(pstat, rownames(pstat@otu_table@.Data), input$taxl)
      index.take <- which(species.name.vec %in% input$select_single_species_name_plot_new)
      df.plot <- pstat@otu_table@.Data
      rownames(df.plot) <- species.name.vec
      if(input$ssv_format_new == "log10 CPM"){
        df.cpm <- getLogCPM(df.plot)
        df.cpm <- df.cpm[index.take,]
        return(df.cpm)
      }else{
        df.ra <- getRelativeAbundance(df.plot)
        df.ra <- df.ra[index.take,]
        return(df.ra)
      }
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
      df.plot <- getDfPlot()
      return(plotHeatmapColor(df.plot,
                              do.scale = FALSE,
                              condition.vec.1 = physeq1@sam_data[[input$select_heatmap_condition_1]],
                              condition.vec.2 = physeq1@sam_data[[input$select_heatmap_condition_2]],
                              condition.1.name = input$select_heatmap_condition_1,
                              condition.2.name = input$select_heatmap_condition_2,
                              annotationColors = add.colorbar,
                              columnTitle = paste("Heatmap with colorbar representing",
                                                  input$select_heatmap_condition, sep = " ")))
    } else  {
      physeq2 <- tax_glom(physeq1, input$taxl)
      df.plot <- getDfPlot()
      if (input$checkbox_heatmap){
        add.colorbar <- "auto"
      } else{
        add.colorbar <- NULL
      }
      return(plotHeatmapColor(df.plot,
                              do.scale = FALSE,
                              condition.vec.1 = physeq2@sam_data[[input$select_heatmap_condition_1]],
                              condition.vec.2 = physeq2@sam_data[[input$select_heatmap_condition_2]],
                              condition.1.name = input$select_heatmap_condition_1,
                              condition.2.name = input$select_heatmap_condition_2,
                              annotationColors = add.colorbar,
                              columnTitle = paste("Heatmap with colorbar representing",
                                                  input$select_heatmap_condition, sep = " ")))
    }
  }

  
  
  
    plotHeatmapColorServerButton <- eventReactive(input$boxplotButtonNew,{
      plotHeatmapColorServer()
    })
  
  # show heatmap in the shiny app by calling the plotting function
  output$Heatmap <- renderPlot({

    plotHeatmapColorServerButton()

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
    meta.data$richness <- suppressWarnings(estimate_richness(physeq = physeq1, split = T, measures = input$select_alpha_div_method)[,1])
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



  output$AlphaDiversityBarplot <- renderPlotly({
      shinyInput <- vals$shiny.input
      pstat <- shinyInput$pstat
      if (input$taxl.alpha !="no rank")  {
          pstat <- tax_glom(pstat, input$taxl.alpha)
      }
      df.pam <- GET_PAM(pstat@otu_table@.Data)

      Sample_Name <- colnames(pstat@otu_table@.Data)
      taxa.num <- as.numeric(colSums(df.pam))
      data <- data.frame(Sample_Name, taxa.num, stringsAsFactors = FALSE)
      data[,3] <- as.character(pstat@sam_data@.Data
                               [[which(pstat@sam_data@names == input$select_alpha_div_condition)]])
      colnames(data)[3] <- input$select_alpha_div_condition
      data$Sample_Name <- paste(as.character(pstat@sam_data@.Data
                                             [[which(pstat@sam_data@names == input$select_alpha_div_condition)]]
      ),data$Sample_Name, sep = "-")
      data$Sample_Name <- factor(data$Sample_Name,
                                 levels = unique(data$Sample_Name)
                                 [order(pstat@sam_data@.Data[[which(pstat@sam_data@names == input$select_alpha_div_condition)]],
                                        decreasing = FALSE)])


      p <- plot_ly(data, x = ~Sample_Name, y = ~taxa.num, type = "bar", color = as.formula(paste("~", input$select_alpha_div_condition, sep = "")), name = 'Sample taxa number') %>%
          layout(margin = list(b = 160))
      p
  })


  output$table.alpha <- DT::renderDataTable({
      shinyInput <- vals$shiny.input
    physeq1 <- shinyInput$pstat
    if (input$taxl.alpha !="no rank")  {
      physeq1 <- tax_glom(physeq1, input$taxl.alpha)
    }
    meta.data <- physeq1@sam_data
    meta.data$sample.name <- rownames(meta.data)
    meta.data$richness <- suppressWarnings(estimate_richness(physeq = physeq1, split = T, measures = input$select_alpha_div_method)[,1])
    colnames(meta.data)[which(colnames(meta.data) == input$select_alpha_div_condition)] <- "condition"
    rownames(meta.data) <- 1:nrow(meta.data)

    DT::datatable(meta.data %>% dplyr::select(sample.name, condition, richness))

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
      meta.data$richness <- suppressWarnings(estimate_richness(physeq = physeq1, split = T, measures = input$select_alpha_div_method)[,1])
      colnames(meta.data)[which(colnames(meta.data) == input$select_alpha_div_condition)] <- "condition"
      rownames(meta.data) <- 1:nrow(meta.data)
      meta.data <- as_tibble(meta.data)
      meta.data <- meta.data %>% select(sample.name, condition, richness)
      write.csv(data.frame(meta.data), file)
    }
  )


  alpha.stat.output <- reactive({
    shinyInput <- vals$shiny.input
    physeq1 <- shinyInput$pstat
    if (input$taxl.alpha !="no rank")  {
      physeq1 <- tax_glom(physeq1, input$taxl.alpha)
    }
    meta.data <- physeq1@sam_data
    meta.data$sample.name <- rownames(meta.data)
    meta.data$richness <- suppressWarnings(estimate_richness(physeq = physeq1, split = T, measures = input$select_alpha_div_method)[,1])
    colnames(meta.data)[which(colnames(meta.data) == input$select_alpha_div_condition)] <- "condition"
    meta.data <- data.frame(meta.data)
    meta.data$condition <- as.factor(meta.data$condition)

    if (length(unique(meta.data$condition)) == 2){
      if (input$select_alpha_stat_method == "Mann-Whitney"){
        tmp <- wilcox.test(richness ~ condition, data = meta.data)
        output <- c(tmp$method, tmp$p.value)
        output.table <- data.frame(output)
        rownames(output.table) <- c("Method", "P-value")
        output.table
      } else if (input$select_alpha_stat_method == "T-test"){
        tmp <- t.test(richness ~ condition, data = meta.data)
        output <- c(tmp$method, tmp$p.value)
        output.table <- data.frame(output)
        rownames(output.table) <- c("Method", "P-value")
        output.table
      } else{
        print("Condition level number is 2, please use Mann-Whitney test.")
      }

    } else if (length(unique(meta.data$condition)) > 2){
      if (input$select_alpha_stat_method == "Mann-Whitney"){
        result.list <- list()
        meta.data.list <- list()
        for (i in 1:length(unique(meta.data$condition))){
          meta.data.list[[i]] <- meta.data[which(meta.data$condition != unique(meta.data$condition)[i]),]
          result.list[[i]] <- wilcox.test(richness ~ condition, data = meta.data.list[[i]])
        }
        output.table <- NULL
        group.name <- c()
        for (i in 1:length(result.list)){
          output.tmp <- c(result.list[[i]]$method, result.list[[i]]$p.value)
          output.table <- cbind(output.table, output.tmp)
          group.name[i] <- paste(unique(meta.data.list[[i]]$condition), collapse = " and ")
        }
        rownames(output.table) <- c("Method", "P-value")
        colnames(output.table) <- group.name
        output.table



      } else if (input$select_alpha_stat_method == "T-test"){
                result.list <- list()
        meta.data.list <- list()
        for (i in 1:length(unique(meta.data$condition))){
          meta.data.list[[i]] <- meta.data[which(meta.data$condition != unique(meta.data$condition)[i]),]
          result.list[[i]] <- t.test(richness ~ condition, data = meta.data.list[[i]])
        }
        output.table <- NULL
        group.name <- c()
        for (i in 1:length(result.list)){
          output.tmp <- c(result.list[[i]]$method, result.list[[i]]$p.value)
          output.table <- cbind(output.table, output.tmp)
          group.name[i] <- paste(unique(meta.data.list[[i]]$condition), collapse = " and ")
        }
        rownames(output.table) <- c("Method", "P-value")
        colnames(output.table) <- group.name
        output.table

      } else{
        tmp <- kruskal.test(richness ~ condition, data = meta.data)
        output <- c(tmp$method, tmp$p.value)
        output.table <- data.frame(output)
        rownames(output.table) <- c("Method", "P-value")
        output.table
      }

    } else{
      "Condition level must be at least 2."
    }
  })
  output$alpha.stat.test <- renderTable({
    alpha.stat.output()
  },include.rownames=TRUE)




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

      if (input$select_beta_div_method == "bray"){
        #First get otu_table and transpose it:
        dist.matrix <- t(data.frame(otu_table(physeq1)))
        #Then use vegdist from vegan to generate a bray distance object:
        dist.mat <- vegdist(dist.matrix, method = "bray")
      }else{
        dist.mat = phyloseq::distance(physeq1, method = input$select_beta_div_method)
      }


      dist.mat <- as.matrix(dist.mat)
      return(plotHeatmapColor(dist.mat,
                              do.scale = FALSE,
                              condition.vec.1 = physeq1@sam_data[[input$select_beta_heatmap_condition_1]],
                              condition.vec.2 = physeq1@sam_data[[input$select_beta_heatmap_condition_2]],
                              condition.1.name = input$select_beta_heatmap_condition_1,
                              condition.2.name = input$select_beta_heatmap_condition_2,
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

      
      if (input$select_beta_div_method == "bray"){
        #First get otu_table and transpose it:
        dist.matrix <- t(data.frame(otu_table(physeq2)))
        #Then use vegdist from vegan to generate a bray distance object:
        dist.mat <- vegdist(dist.matrix, method = "bray")
      }else{
        dist.mat = phyloseq::distance(physeq2, method = input$select_beta_div_method)
      }
      dist.mat <- as.matrix(dist.mat)
      return(plotHeatmapColor(dist.mat,
                              do.scale = FALSE,
                              condition.vec.1 = physeq2@sam_data[[input$select_beta_heatmap_condition_1]],
                              condition.vec.2 = physeq2@sam_data[[input$select_beta_heatmap_condition_2]],
                              condition.1.name = input$select_beta_heatmap_condition_1,
                              condition.2.name = input$select_beta_heatmap_condition_2,
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

    
    
    if (input$select_beta_div_method == "bray"){
        #First get otu_table and transpose it:
        dist.matrix <- t(data.frame(otu_table(physeq1)))
        #Then use vegdist from vegan to generate a bray distance object:
        dist.tmp <- vegdist(dist.matrix, method = "bray")
    }else{
        dist.tmp = phyloseq::distance(physeq1, method = input$select_beta_div_method)
    }
    
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




  beta.stat.output <- reactive({
    shinyInput <- vals$shiny.input
    physeq1 <- shinyInput$pstat
    if (input$select_beta_stat_method == "PERMANOVA"){
      if (input$taxl.beta !="no rank")  {
        physeq1 <- tax_glom(physeq1, input$taxl.beta)
      }
      meta.data <- physeq1@sam_data
      meta.data$sample.name <- rownames(meta.data)
      colnames(meta.data)[which(colnames(meta.data) == input$select_beta_condition)] <- "condition"
      meta.data <- data.frame(meta.data)
      meta.data$condition <- as.factor(meta.data$condition)

      set.seed(99)
      
      if (input$select_beta_div_method == "bray"){
        #First get otu_table and transpose it:
        dist.matrix <- t(data.frame(otu_table(physeq1)))
        #Then use vegdist from vegan to generate a bray distance object:
        dist.tmp <- vegdist(dist.matrix, method = "bray")
      }else{
        dist.tmp = phyloseq::distance(physeq1, method = input$select_beta_div_method)
      }      
      
      beta.div <- adonis2(dist.tmp~condition, data=meta.data,
                          permutations = input$num.permutation.permanova, strata="PLOT")
      beta.div

    } else {

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
        result.list <- list()
        group.name <- c()
        for (i in 1:length(dist.list)){
          dist.list.tmp <- dist.list[which(names(dist.list) != names(dist.list)[i])]

          group.name[i] <- paste(names(dist.list.tmp), collapse = " and ")
          result.list[[i]] <- wilcox.test(dist.list.tmp[[1]], dist.list.tmp[[2]])
        }
        output.table <- NULL
        group.name <- c()
        for (i in 1:length(result.list)){
          output.tmp <- c(result.list[[i]]$method, result.list[[i]]$p.value)
          output.table <- cbind(output.table, output.tmp)
        }
        rownames(output.table) <- c("Method", "P-value")
        colnames(output.table) <- group.name
        output.table


      } else{
        tmp <- kruskal.test(list(dist.within.a, dist.within.b, dist.between))
        output <- c(tmp$method, tmp$p.value)
        output.table <- data.frame(output)
        rownames(output.table) <- c("Method", "P-value")
        output.table
      }

    }


  })

  output$beta.stat.test <- renderTable({
    beta.stat.output()

  },include.rownames=TRUE)

  output$table.beta <- DT::renderDataTable({
      shinyInput <- vals$shiny.input
    physeq1 <- shinyInput$pstat

    if (input$taxl.beta=="no rank")  {
        if (input$select_beta_div_method == "bray"){
        #First get otu_table and transpose it:
        dist.matrix <- t(data.frame(otu_table(physeq1)))
        #Then use vegdist from vegan to generate a bray distance object:
        dist.mat <- vegdist(dist.matrix, method = "bray")
    }else{
        dist.mat = phyloseq::distance(physeq1, method = input$select_beta_div_method)
    }

    } else{
      physeq2 <- tax_glom(physeq1, input$taxl.beta)
        if (input$select_beta_div_method == "bray"){
        #First get otu_table and transpose it:
        dist.matrix <- t(data.frame(otu_table(physeq2)))
        #Then use vegdist from vegan to generate a bray distance object:
        dist.mat <- vegdist(dist.matrix, method = "bray")
        }else{
        dist.mat = phyloseq::distance(physeq2, method = input$select_beta_div_method)
        }
    }
    dist.mat <- as.matrix(dist.mat)
    return(dist.mat)

  },
  options = list(
      paging = TRUE, scrollX = TRUE
  ))

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
      df.plot <- physeq1@otu_table@.Data
      if (input$select_pca_data_format == "log10 CPM"){
        df.plot <- getLogCPM(df.plot)
      }else if (input$select_pca_data_format == "RA"){
        df.plot <- getRelativeAbundance(df.plot)
      }
      plotPCAPlotly(df.input = df.plot,
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
      df.plot <- physeq2@otu_table@.Data
      if (input$select_pca_data_format == "log10 CPM"){
        df.plot <- getLogCPM(df.plot)
      }else if (input$select_pca_data_format == "RA"){
        df.plot <- getRelativeAbundance(df.plot)
      }
      plotPCAPlotly(df.input = df.plot,
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
  
  
  plotPCAPlotlyServerButton <- eventReactive(input$DR_plot,{
      plotPCAPlotlyServer()
  })
    
  output$pca.plotly <- renderPlotly({
    plotPCAPlotlyServerButton()

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
    #hide std
    DT::datatable(table.output.pca[,-1])

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

  
    
  plotPCoAPlotlyServerButton <- eventReactive(input$DR_plot,{
      plotPCoAPlotlyServer()
  })
    
  output$pcoa.plotly <- renderPlotly({
    plotPCoAPlotlyServerButton()

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
      if (input$select_beta_div_method == "bray"){
        #First get otu_table and transpose it:
        dist.matrix <- t(data.frame(otu_table(physeq1)))
        #Then use vegdist from vegan to generate a bray distance object:
        DistBC <- vegdist(dist.matrix, method = "bray")
        ord.tmp <- ordinate(physeq1, method = "PCoA", distance = DistBC)
      } else{
        Dist.tmp <- phyloseq::distance(physeq1, method = input$select_beta_div_method)
        ord.tmp <- ordinate(physeq1, method = "PCoA", distance = Dist.tmp)
      }

      #cat(dim(physeq1@otu_table))
      return(ord.tmp$values)

    } else  {
      physeq2 <- tax_glom(physeq1, input$taxl.pca)
      if (sum(rowSums(as.matrix(physeq2@otu_table@.Data)) == 0) > 0){
        physeq2@otu_table@.Data <- data.frame(physeq2@otu_table@.Data[-which
                                                                      (rowSums(as.matrix(physeq2@otu_table@.Data)) == 0),])
      }
      if (input$select_beta_div_method == "bray"){
        #First get otu_table and transpose it:
        dist.matrix <- t(data.frame(otu_table(physeq2)))
        #Then use vegdist from vegan to generate a bray distance object:
        DistBC <- vegdist(dist.matrix, method = "bray")
        ord.tmp <- ordinate(physeq2, method = "PCoA", distance = DistBC)
      } else{
        Dist.tmp <- phyloseq::distance(physeq2, method = input$select_beta_div_method)
        ord.tmp <- ordinate(physeq2, method = "PCoA", distance = Dist.tmp)
      }
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
    # hide eigenvalue
    DT::datatable(df.output[,-1])

  })


  ## New DA analysis section

  output$edgerTable.new <- DT::renderDataTable({
      shinyInput <- vals$shiny.input
      pstat <- shinyInput$pstat
      if (input$taxl.edger !="no rank"){
          pstat <- tax_glom(pstat, input$taxl.edger)
      }



      # number of samples in each level of target variable
      target.var.index <- which(pstat@sam_data@names == input$edger.condition)
      label.vec.num <- pstat@sam_data@.Data[[target.var.index]]



      # if selected condition has multiple levels
      if (length(input$edger_condition_options_use) == 2){
        sample.keep.index <- which(label.vec.num %in% input$edger_condition_options_use)
        label.vec.num <- label.vec.num[sample.keep.index]
        pstat@otu_table@.Data <- pstat@otu_table@.Data[,sample.keep.index]
        pstat@sam_data@row.names <- pstat@sam_data@row.names[sample.keep.index]
        for (i in 1:length(pstat@sam_data@names)){
          pstat@sam_data@.Data[[i]] <- pstat@sam_data@.Data[[i]][sample.keep.index]
        }
      }


      label.vec.save <- unique(label.vec.num)

      # transform label into 1 and 0
      label.vec.num[label.vec.num == unique(label.vec.num)[1]] <- 1
      label.vec.num[label.vec.num != 1] <- 0



      dge = phyloseq_to_edgeR(pstat, group=input$edger.condition)
      # Perform binary test
      et = exactTest(dge)
      # Extract values from test results
      tt = topTags(et, n=nrow(dge$table), adjust.method="BH", sort.by="PValue")
      res = tt@.Data[[1]]
      sigtab = res[(res$FDR < input$edger.padj.cutoff), ]
      sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(pstat)[rownames(sigtab), ], "matrix"))

          if (nrow(sigtab) == 0){
              return(as.matrix("No differentially abundant items found!"))
          } else{
              sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(pstat)[rownames(sigtab), ], "matrix"))
              sigtab$FDR <- as.numeric(formatC(sigtab$FDR, format = "e", digits = 2))
              sigtab$logFC <- as.numeric(formatC(sigtab$logFC, format = "e", digits = 2))


              sigtab <- sigtab[,c(which(colnames(sigtab) == input$taxl.edger)[1],
                                  which(colnames(sigtab) == "FDR"),
                                  which(colnames(sigtab) == "logFC"))]
              rownames(sigtab) <- 1:nrow(sigtab)



              # remove "others"
              index.sigtab.save <- c()
              for (i in 1:nrow(sigtab)){
                if (sigtab[i,1] != "others"){
                  index.sigtab.save <- c(index.sigtab.save, i)
                }
              }
              sigtab <- sigtab[index.sigtab.save,]

              num.1 <- c()
              num.2 <- c()
              species.names.tmp <- TranslateIdToTaxLevel(pstat, rownames(pstat@otu_table@.Data), input$taxl.edger)
              for (i in 1:nrow(sigtab)){
                species.index <- which(species.names.tmp == sigtab[i,1])
                num.1 <- c(num.1, sum((pstat@otu_table@.Data[species.index,which(label.vec.num == 1)] > 0)))
                num.2 <- c(num.2, sum((pstat@otu_table@.Data[species.index,which(label.vec.num == 0)] > 0)))
              }

              sigtab <- cbind(sigtab, num.1)
              sigtab <- cbind(sigtab, num.2)


              df.output.prevalence <- percent(round((num.1 + num.2)/ncol(pstat@otu_table@.Data),4))
              sigtab <- cbind(sigtab, df.output.prevalence)


              colnames(sigtab)[ncol(sigtab)-2] <- label.vec.save[1]
              colnames(sigtab)[ncol(sigtab)-1] <- label.vec.save[2]
              colnames(sigtab)[ncol(sigtab)] <- "prevalence"


              output$download_edger_tb <- downloadHandler(
                  filename = function() { paste('download_edger_table', '.csv', sep='') },
                  content = function(file) {
                      dist.mat <- as.matrix(sigtab)
                      write.csv(data.frame(dist.mat), file)
                  }
              )
              foldChange <- c()
              for (i in 1:nrow(sigtab)){
              foldChange[i] <- round((max(as.numeric(c(sigtab[i,5],
                                                 sigtab[i,4]))) /
                       min(as.numeric(c(sigtab[i,5],
                                        sigtab[i,4])))), digits = 2)
              }
              sigtab <- cbind(sigtab, foldChange)
              return(sigtab)

          }


  },
  options = list(
      paging = TRUE
  ))







  output$DeSeq2Table.new <- DT::renderDataTable({
      shinyInput <- vals$shiny.input
    physeq1 <- shinyInput$pstat
    if (input$taxl.da !="no rank"){
      physeq1 <- tax_glom(physeq1, input$taxl.da)
    }
    physeq1 <- prune_samples(sample_sums(physeq1) > input$da.count.cutoff, physeq1)




    # number of samples in each level of target variable
    target.var.index <- which(physeq1@sam_data@names == input$da.condition)
    label.vec.num <- physeq1@sam_data@.Data[[target.var.index]]


    # if selected condition has multiple levels
    if (length(input$da_condition_options_use) == 2){
      sample.keep.index <- which(label.vec.num %in% input$da_condition_options_use)
      label.vec.num <- label.vec.num[sample.keep.index]
      physeq1@otu_table@.Data <- physeq1@otu_table@.Data[,sample.keep.index]
      physeq1@sam_data@row.names <- physeq1@sam_data@row.names[sample.keep.index]
      for (i in 1:length(physeq1@sam_data@names)){
        physeq1@sam_data@.Data[[i]] <- physeq1@sam_data@.Data[[i]][sample.keep.index]
      }
    }




    label.vec.save <- unique(label.vec.num)

    # transform label into 1 and 0
    label.vec.num[label.vec.num == unique(label.vec.num)[1]] <- 1
    label.vec.num[label.vec.num != 1] <- 0



    # deal with continuous covariates and multiple covariates
    # target condition is the last one in formula.
    if (!is.null(input$da.condition.covariate)){
      for (i in 1:length(input$da.condition.covariate)){
        num.levels <- length(unique(sample_data(physeq1)[[input$da.condition.covariate[i]]]))
        if (num.levels >= 8){
          sam.index <- which(physeq1@sam_data@names %in% input$da.condition.covariate[i])
          physeq1@sam_data@.Data[[sam.index]] <- cut(physeq1@sam_data@.Data[[sam.index]], breaks = 3)
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
        return(as.matrix("No differentially abundant items found!"))
      } else{
        sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq1)[rownames(sigtab), ], "matrix"))
        sigtab$padj <- as.numeric(formatC(sigtab$padj, format = "e", digits = 2))
        sigtab$log2FoldChange <- as.numeric(formatC(sigtab$log2FoldChange, format = "e", digits = 2))

        rownames(sigtab) <- 1:nrow(sigtab)






        sigtab <- sigtab[,c(which(colnames(sigtab) == input$taxl.da[1]),
                            which(colnames(sigtab) == "padj"),
                            which(colnames(sigtab) == "log2FoldChange"))]





        # remove "others"
        index.sigtab.save <- c()
        for (i in 1:nrow(sigtab)){
          if (sigtab[i,1] != "others"){
            index.sigtab.save <- c(index.sigtab.save, i)
          }
        }
        sigtab <- sigtab[index.sigtab.save,]

        num.1 <- c()
        num.2 <- c()
        species.names.tmp <- TranslateIdToTaxLevel(physeq1, rownames(physeq1@otu_table@.Data), input$taxl.da)
        for (i in 1:nrow(sigtab)){
          species.index <- which(species.names.tmp == sigtab[i,1])
          num.1 <- c(num.1, sum((physeq1@otu_table@.Data[species.index,which(label.vec.num == 1)] > 0)))
          num.2 <- c(num.2, sum((physeq1@otu_table@.Data[species.index,which(label.vec.num == 0)] > 0)))
        }

        sigtab <- cbind(sigtab, num.1)
        sigtab <- cbind(sigtab, num.2)


        df.output.prevalence <- percent(round((num.1 + num.2)/ncol(physeq1@otu_table@.Data),4))
        sigtab <- cbind(sigtab, df.output.prevalence)


        colnames(sigtab)[ncol(sigtab)-2] <- label.vec.save[1]
        colnames(sigtab)[ncol(sigtab)-1] <- label.vec.save[2]
        colnames(sigtab)[ncol(sigtab)] <- "prevalence"

        output$download_deseq_tb <- downloadHandler(
          filename = function() { paste('download_deseq2_table', '.csv', sep='') },
          content = function(file) {
            dist.mat <- as.matrix(sigtab)
            write.csv(data.frame(dist.mat), file)
          }
        )
        foldChange <- c()
        for (i in 1:nrow(sigtab)){
        foldChange[i] <- round((max(as.numeric(c(sigtab[i,5],
                                                 sigtab[i,4]))) /
                       min(as.numeric(c(sigtab[i,5],
                                        sigtab[i,4])))), digits = 2)
        }
        sigtab <- cbind(sigtab, foldChange)
        return(sigtab)

      }

    }else{
      return(as.matrix("No differentially abundant items found!"))
    }

  },
  options = list(
      paging = TRUE
  ))



  ### check whether selected condition for DA has two levels or more
  output$da_condition_type <- reactive({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    target.var.index <- which(pstat@sam_data@names == input$da.condition)
    label.vec <- pstat@sam_data@.Data[[target.var.index]]
    label.level.num <- length(unique(label.vec))
    if (label.level.num == 2){
      return("binary")
    } else{
      return("multiple")
    }

  })
  outputOptions(output, "da_condition_type", suspendWhenHidden = FALSE)

  # select 2 levels
  output$da_condition_options <- renderUI({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    variable.vec <- sample_data(pstat)[[
      which(pstat@sam_data@names == input$da.condition)]]
    filter.option.vec <- sort(unique(variable.vec))
    tagList(
      selectInput("da_condition_options_use", "Select 2 levels", choices = filter.option.vec, multiple = TRUE)
    )
  })


  ### check whether selected condition for edger has two levels or more
  output$edger_condition_type <- reactive({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    target.var.index <- which(pstat@sam_data@names == input$edger.condition)
    label.vec <- pstat@sam_data@.Data[[target.var.index]]
    label.level.num <- length(unique(label.vec))
    if (label.level.num == 2){
      return("binary")
    } else{
      return("multiple")
    }

  })
  outputOptions(output, "edger_condition_type", suspendWhenHidden = FALSE)

  # select 2 levels
  output$edger_condition_options <- renderUI({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    variable.vec <- sample_data(pstat)[[
      which(pstat@sam_data@names == input$edger.condition)]]
    filter.option.vec <- sort(unique(variable.vec))
    tagList(
      selectInput("edger_condition_options_use", "Select 2 levels", choices = filter.option.vec, multiple = TRUE)
    )
  })







  # Presence-Absence Variance analysis

  output$pa.test <- DT::renderDataTable({
      shinyInput <- vals$shiny.input
    physeq1 <- shinyInput$pstat
    if (input$taxl.pa !="no rank"){
      physeq1 <- tax_glom(physeq1, input$taxl.pa)
    }
    physeq1 <- prune_samples(sample_sums(physeq1) > input$pa.count.cutoff, physeq1)


    target.var.index <- which(physeq1@sam_data@names == input$pa.condition)
    label.vec.num <- physeq1@sam_data@.Data[[target.var.index]]
    # if selected condition has multiple levels
    if (length(input$pa_condition_options_use) == 2){
      sample.keep.index <- which(label.vec.num %in% input$pa_condition_options_use)
      label.vec.num <- label.vec.num[sample.keep.index]
      physeq1@otu_table@.Data <- physeq1@otu_table@.Data[,sample.keep.index]
      physeq1@sam_data@row.names <- physeq1@sam_data@row.names[sample.keep.index]
      for (i in 1:length(physeq1@sam_data@names)){
        physeq1@sam_data@.Data[[i]] <- physeq1@sam_data@.Data[[i]][sample.keep.index]
      }
    }


    df.pam <- GET_PAM(physeq1@otu_table@.Data)


    # change microbe names to selected taxon level
    rownames(df.pam) <- TranslateIdToTaxLevel(physeq1, rownames(df.pam), input$taxl.pa)

    if (input$pa_method == "Fisher Exact Test"){
      output.mat <- Fisher_Test_Pam(df.pam,
                                    physeq1@sam_data[[input$pa.condition]],
                                    input$pa.padj.cutoff)
      output$download_pa_test <- downloadHandler(
        filename = function() { paste(input$pa_method, '.csv', sep='') },
        content = function(file) {
          write.csv(output.mat, file)
        }
      )

      DT::datatable(output.mat)
    } else if(input$pa_method == "Chi-squared Test"){
      output.mat <- Chisq_Test_Pam(df.pam,
                                   physeq1@sam_data[[input$pa.condition]],
                                   input$pa.padj.cutoff)
      output$download_pa_test <- downloadHandler(
        filename = function() { paste(input$pa_method, '.csv', sep='') },
        content = function(file) {
          write.csv(output.mat, file)
        }
      )
      DT::datatable(output.mat, rownames = TRUE)
    } else if (input$pa_method == "Mann-Whitney Test"){
      df.test <- physeq1@otu_table@.Data
      rownames(df.test) <- TranslateIdToTaxLevel(physeq1, rownames(df.test), input$taxl.pa)

      if (input$pa_mann_data_type == "log10 CPM"){
        df.test <- getLogCPM(df.test)
      }else if (input$pa_mann_data_type == "RA"){
        df.test <- getRelativeAbundance(df.test)
      }


      output.mat <- Wilcox_Test_df(df.test,
                                   physeq1@sam_data[[input$pa.condition]],
                                   input$pa.padj.cutoff)
      output$download_pa_test <- downloadHandler(
        filename = function() { paste(input$pa_method, '.csv', sep='') },
        content = function(file) {
          write.csv(output.mat, file)
        }
      )
      DT::datatable(output.mat)
    }


  })


  ### check whether selected condition for pa has two levels or more
  output$pa_condition_type <- reactive({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    target.var.index <- which(pstat@sam_data@names == input$pa.condition)
    label.vec <- pstat@sam_data@.Data[[target.var.index]]
    label.level.num <- length(unique(label.vec))
    if (label.level.num == 2){
      return("binary")
    } else{
      return("multiple")
    }

  })
  outputOptions(output, "pa_condition_type", suspendWhenHidden = FALSE)

  # select 2 levels
  output$pa_condition_options <- renderUI({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    variable.vec <- sample_data(pstat)[[
      which(pstat@sam_data@names == input$pa.condition)]]
    filter.option.vec <- sort(unique(variable.vec))
    tagList(
      selectInput("pa_condition_options_use", "Select 2 levels", choices = filter.option.vec, multiple = TRUE)
    )
  })


### biomarker


  getBiomarker <- function(){
    if (input$select_model_biomarker == "Lasso Logistic Regression"){
      shinyInput <- vals$shiny.input
      physeq1 <- shinyInput$pstat
      if (input$taxl_biomarker !="no rank"){
        physeq1 <- tax_glom(physeq1, input$taxl_biomarker)

      }

      target.var.index <- which(physeq1@sam_data@names == input$select_target_condition_biomarker)
      label.vec.num <- physeq1@sam_data@.Data[[target.var.index]]
      # if selected condition has multiple levels
      if (length(input$biomarker_condition_options_use) == 2){
        sample.keep.index <- which(label.vec.num %in% input$biomarker_condition_options_use)
        label.vec.num <- label.vec.num[sample.keep.index]
        physeq1@otu_table@.Data <- physeq1@otu_table@.Data[,sample.keep.index]
        physeq1@sam_data@row.names <- physeq1@sam_data@row.names[sample.keep.index]
        for (i in 1:length(physeq1@sam_data@names)){
          physeq1@sam_data@.Data[[i]] <- physeq1@sam_data@.Data[[i]][sample.keep.index]
        }
      }


      # use log CPM as normalization.
      df.input <- physeq1@otu_table@.Data
      for(i in 1:ncol(df.input)){
        if (is.numeric(df.input[,i])){
          df.input[,i] <- log10(df.input[,i]*1e6/sum(df.input[,i]) + 0.1)
        }
      }

      # change microbe names to selected taxon level
      rownames(df.input) <- TranslateIdToTaxLevel(physeq1, rownames(df.input), input$taxl_biomarker)

      if (!is.null(input$select_covariate_condition_biomarker)){

        target.tmp <- physeq1@sam_data@.Data[[2]]

        covariate.vec <- input$select_covariate_condition_biomarker
        df.list.tmp <- list()
        for (i in 1:length(covariate.vec)){
          df.list.tmp[[covariate.vec[i]]] <- physeq1@sam_data@.Data[[which(physeq1@sam_data@names == covariate.vec[i])]]
        }
        df.covariate <- data.frame(df.list.tmp)
        rownames(df.covariate) <- colnames(df.input)

        df.input <- cbind.data.frame(t(df.input), df.covariate)
      } else{
        df.input <- t(df.input)
      }


      target.vec <- physeq1@sam_data[[input$select_target_condition_biomarker]]

      output.fs <- getSignatureFromMultipleGlmnet(df.input, target.vec, nfolds = input$num.cv.nfolds, nRun = input$num.biomarker.run)
    }

    # number of samples in each level of target variable
    target.var.index <- which(physeq1@sam_data@names == input$select_target_condition_biomarker)
    label.vec.num <- physeq1@sam_data@.Data[[target.var.index]]
    label.vec.save <- unique(label.vec.num)

    # transform label into 1 and 0
    label.vec.num[label.vec.num == unique(label.vec.num)[1]] <- 1
    label.vec.num[label.vec.num != 1] <- 0

    label.vec.num <- as.numeric(label.vec.num)
    num.1 <- c()
    num.2 <- c()
    species.names.tmp <- TranslateIdToTaxLevel(physeq1, rownames(physeq1@otu_table@.Data), input$taxl_biomarker)
    for (i in 1:length(output.fs$feature)){
      species.index <- which(species.names.tmp == output.fs$feature[i])
      num.1 <- c(num.1, sum((physeq1@otu_table@.Data[species.index,which(label.vec.num == 1)] > 0)))
      num.2 <- c(num.2, sum((physeq1@otu_table@.Data[species.index,which(label.vec.num == 0)] > 0)))
    }

    output.df <- data.frame(biomarker = output.fs$feature,
                            selection_rate = output.fs$selection_rate,
                            average_weights = output.fs$feature_weights,
                            num.1,
                            num.2,
                            percent(round((num.1+num.2)/ncol(physeq1@otu_table@.Data),4)))
    foldChange <- c()
    for (i in 1:nrow(output.df)){
      foldChange[i] <- round((max(as.numeric(c(output.df[i,5], output.df[i,4]))) /
                       min(as.numeric(c(output.df[i,5], output.df[i,4])))), digits = 2)
    }
    output.df <- cbind(output.df, foldChange)
    colnames(output.df)[ncol(output.df)-3] <- label.vec.save[1]
    colnames(output.df)[ncol(output.df)-2] <- label.vec.save[2]
    colnames(output.df)[ncol(output.df)-1] <- "prevalence"
    colnames(output.df)[ncol(output.df)] <- "Fold-Change"
    df.biomarker <- df.input[,which(colnames(df.input) %in% output.fs$feature)]
    return(list(output.df = output.df, df.input = df.biomarker, target.vec = label.vec.num))


  }



  observeEvent(input$goButtonBiomarker, {

     biomarker.vals <- reactiveValues(
      biomarker.list = suppressWarnings(getBiomarker())
     )
     biomarker.vals.2 <- reactiveValues(
       loocv.output.list = LOOAUC_simple_multiple_one_df(biomarker.vals$biomarker.list$df.input,
                                                         biomarker.vals$biomarker.list$target.vec)
     )

      output$featureSelectionTmp <- renderTable({
        biomarker.vals$biomarker.list$output.df
      })

      output$loocv_output_simple <- renderPlot({

        plot(biomarker.vals.2$loocv.output.list$loo.perf.plot, main="ROC curve",col="red",lwd=3, cex=4)
        auc <- performance(biomarker.vals.2$loocv.output.list$loo.pred.plot,"auc")
        auc <- unlist(slot(auc, "y.values"))
        aucRound <- paste("AUC: ", round(auc,3))
        abline(a=0,b=1,lwd=2,lty=2,col="gray")
        legend(0.7,0.5, aucRound)
      })


      output$loocv.violin <- renderPlotly({
        df.tmp <- data.frame(class.vec = biomarker.vals.2$loocv.output.list$testPredictionClassVec,
                             class = c(rep(0, sum(biomarker.vals$biomarker.list$target.vec == 0)),
                                       rep(1, sum(biomarker.vals$biomarker.list$target.vec == 1))))
        p.tmp <- df.tmp %>%
          plot_ly(
            x = ~class,
            y = ~class.vec,
            split = ~class,
            type = 'violin',
            box = list(
              visible = T
            ),
            meanline = list(
              visible = T
            )
          ) %>%
          layout(
            xaxis = list(
              title = "Class"
            ),
            yaxis = list(
              title = "Predicted Probability",
              zeroline = F
            )
          )
        p.tmp

      })

      bootstrap.out <- eventReactive(input$goButtonBiomarkerBoot, {
        suppressWarnings(Bootstrap_LOOCV_LR_AUC(biomarker.vals$biomarker.list$df.input,
                                                biomarker.vals$biomarker.list$target.vec,
                                                nboot = input$num.bootstrap.loocv))
      })

      output$loocv_output <- renderTable({
        bootstrap.out()
      })

  })




  ### check whether selected condition for biomarker has two levels or more
  output$biomarker_condition_type <- reactive({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    target.var.index <- which(pstat@sam_data@names == input$select_target_condition_biomarker)
    label.vec <- pstat@sam_data@.Data[[target.var.index]]
    label.level.num <- length(unique(label.vec))
    if (label.level.num == 2){
      return("binary")
    } else{
      return("multiple")
    }

  })
  outputOptions(output, "biomarker_condition_type", suspendWhenHidden = FALSE)

  # select 2 levels
  output$biomarker_condition_options <- renderUI({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    variable.vec <- sample_data(pstat)[[
      which(pstat@sam_data@names == input$select_target_condition_biomarker)]]
    filter.option.vec <- sort(unique(variable.vec))
    tagList(
      selectInput("biomarker_condition_options_use", "Select 2 levels", choices = filter.option.vec, multiple = TRUE)
    )
  })

})