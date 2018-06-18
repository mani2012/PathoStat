

  ## New DA analysis section

  run.edger <- eventReactive(input$run_edgeR, {
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

  })
  output$edgerTable.new <- DT::renderDataTable({
    run.edger()
  },
  options = list(
      paging = TRUE
  ))




  run.deseq2 <- eventReactive(input$run_deseq2, {
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

  })


  output$DeSeq2Table.new <- DT::renderDataTable({
    run.deseq2()
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

  run.pav <- eventReactive(input$run_stat, {
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
  output$pa.test <- DT::renderDataTable({
    run.pav()
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
