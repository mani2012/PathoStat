

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
    label.vec <- label.vec.num
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

    df.output <- data.frame(biomarker = output.fs$feature,
                            selection_rate = output.fs$selection_rate,
                            average_weights = output.fs$feature_weights,
                            num.1,
                            num.2,
                            percent(round((num.1+num.2)/ncol(physeq1@otu_table@.Data),4)))
    foldChange <- c()
    for (i in 1:nrow(df.output)){
      foldChange[i] <- round((max(as.numeric(c((df.output[i,5] / sum(label.vec.num == 0)),
                                               (df.output[i,4] / sum(label.vec.num == 1))))) /
                                min(as.numeric(c((df.output[i,5] / sum(label.vec.num == 0)),
                                                 (df.output[i,4] / sum(label.vec.num == 1)))))), 
                             digits = 2)
    }
    df.output <- cbind(df.output, foldChange)
    colnames(df.output)[ncol(df.output)-3] <- label.vec.save[1]
    colnames(df.output)[ncol(df.output)-2] <- label.vec.save[2]
    colnames(df.output)[ncol(df.output)-1] <- "prevalence"
    colnames(df.output)[ncol(df.output)] <- "group size adjusted fold change"
    df.biomarker <- df.input[,which(colnames(df.input) %in% output.fs$feature)]
    return(list(df.output = df.output, 
                df.input = df.biomarker, 
                target.vec = label.vec.num,
                label.vec = label.vec))


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
        biomarker.vals$biomarker.list$df.output
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
        df.tmp <- data.frame(class.vec = biomarker.vals.2$loocv.output.list$testPredictionProb,
                             class = biomarker.vals$biomarker.list$label.vec)
        #print(df.tmp)
        p.tmp <- df.tmp %>%
          plot_ly(
            color = ~class,
            y = ~class.vec,
            type = 'box',
            boxpoints = "all", 
            jitter = 0.3,
            pointpos = -1.8,
            colors = c("#132B43", "#56B1F7")
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
    if (is.integer0(target.var.index)){
      target.var.index <- 1
    }    
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
