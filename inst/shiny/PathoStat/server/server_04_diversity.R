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
