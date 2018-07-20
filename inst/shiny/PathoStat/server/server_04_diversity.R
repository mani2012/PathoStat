# Alpha diversity boxplots
plotAlphaServer <- function() {
  shinyInput <- vals$shiny.input
  physeq1 <- shinyInput$pstat
  if (input$taxl.alpha != "no rank") {
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
  g <- ggplotly(g, tooltip="text")
  g$elementId <- NULL # To suppress a shiny warning
  return(g)
}

plotAlphaBoxplotButton <- eventReactive(input$alpha_boxplot,{
  plotAlphaServer()
})
output$AlphaDiversity <- renderPlotly({
  plotAlphaBoxplotButton()
})

# Alpha diversity statistics
do_alpha_stats <- function() {
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
}
plotAlphaBoxplotButton3 <- eventReactive(input$alpha_boxplot,{
  do_alpha_stats()
})
output$alpha.stat.test <- DT::renderDataTable({
  plotAlphaBoxplotButton3()
}, options = list(sDom  = '<"top">t<"bottom">ip'))

# Alpha diversity table
do_alpha_table <- function() {
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
}
plotAlphaBoxplotButton2 <- eventReactive(input$alpha_boxplot,{
  do_alpha_table()
})
output$table.alpha <- DT::renderDataTable({
  plotAlphaBoxplotButton2()
})

# Download alpha diversity table
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

do_beta_heatmap <- function(){
  shinyInput <- vals$shiny.input
  physeq1 <- shinyInput$pstat
  tables <- pstat.extraction(shinyInput$pstat)
  SAM_DATA <- tables$SAM

  if (input$taxl.beta != "no rank") {
    physeq2 <- tax_glom(physeq1, input$taxl.beta)
  } else {
    physeq2 <- physeq1
  }
  if (input$select_beta_div_method == "bray"){
    #First get otu_table and transpose it:
    dist.matrix <- t(data.frame(otu_table(physeq2)))
    #Then use vegdist from vegan to generate a bray distance object:
    dist.mat <- vegdist(dist.matrix, method = "bray")
  }else{
    dist.mat = phyloseq::distance(physeq2, method = input$select_beta_div_method)
  }
  dist.mat <- as.data.frame(dist.mat)
  dist.mat <- dist.mat[order(match(rownames(dist.mat), rev(rownames(dist.mat)))),,drop=FALSE]

  if (!is.null(input$bdhm_select_conditions)) {
    df.sam <- SAM_DATA[,input$bdhm_select_conditions,drop=FALSE]
    if (input$bdhm_sort_by == "conditions") {
      for (i in ncol(df.sam):1) {
        df.sam <- df.sam[order(df.sam[[i]]),,drop=FALSE]
      }
      dist.mat <- dist.mat[order(match(rownames(dist.mat), rownames(df.sam))),,drop=FALSE]
      dist.mat <- dist.mat[,order(match(colnames(dist.mat), rownames(df.sam))),drop=FALSE]
    } else {
      df.sam <- df.sam[order(match(rownames(df.sam), rownames(dist.mat))),,drop=FALSE]
    }
  }

  m <- data.matrix(dist.mat)
  hover.txt <- c()
  for (i in 1:ncol(dist.mat)) {
    hover.txt <- cbind(hover.txt, dist.mat[[i]])
  }
  hm.beta <- plot_ly(x = colnames(m), y = rownames(m), z = m,
                     type = "heatmap",
                     colors= "RdPu",
                     hoverinfo = "x+y+z") %>%
    layout(xaxis = list(showticklabels = FALSE, title = "", ticks = "", tickangle = -45),
           yaxis = list(showticklabels = FALSE, type = 'category', ticks = ""))

  if (!is.null(input$bdhm_select_conditions)) {
    hover.txt <- c()
    for (i in 1:ncol(df.sam)) {
      hover.txt <- cbind(hover.txt, df.sam[[i]])
    }
    df.sam[] <- lapply(df.sam, factor)
    
    # Y-axis of subplot
    m <- data.matrix(df.sam)
    m.row.normalized <- apply(m, 2, function(x)(x-min(x))/(max(x)-min(x)))
    hm.sam.y <- plot_ly(x = colnames(m.row.normalized), 
                        y = rownames(m.row.normalized), 
                        z = m.row.normalized, 
                        type = "heatmap",
                        showscale=FALSE,
                        hoverinfo = "x+y+text",
                        transpose=FALSE,
                        text=hover.txt) %>%
      layout(xaxis = list(title = "", tickangle = -45),
             yaxis = list(showticklabels = FALSE, type = 'category', ticks = ""),
             orientation=TRUE)
    
    # X-axis of subplot
    m <- data.matrix(df.sam)
    m.row.normalized <- apply(m, 2, function(x)(x-min(x))/(max(x)-min(x)))
    m.row.normalized = t(m.row.normalized)
    m.row.normalized = m.row.normalized[order(match(rownames(m.row.normalized), rev(rownames(m.row.normalized)))),,drop=FALSE]
    hm.sam.x <- plot_ly(x = colnames(m.row.normalized), 
                        y = rownames(m.row.normalized), 
                        z = m.row.normalized, 
                        type = "heatmap",
                        showscale=FALSE,
                        hoverinfo = "x+y+text",
                        transpose=FALSE,
                        text=t(hover.txt)) %>%
      layout(xaxis = list(showticklabels = FALSE, type = 'category', ticks = ""),
             yaxis = list(title = "", tickangle = -45),
             orientation=TRUE)
  }

  empty <- plotly_empty(type = "scatter")
  if (!is.null(input$bdhm_select_conditions)) {
    hm.sam.beta.top <- subplot(empty, hm.sam.x, widths=c(0.1,  0.9))
    hm.sam.beta.bot <- subplot(hm.sam.y, hm.beta, widths=c(0.1,  0.9))
    hm.sam.beta <- subplot(hm.sam.beta.top, hm.sam.beta.bot, nrows=2, heights=c(0.1,  0.9))
    hm.sam.beta$elementId <- NULL # To suppress a shiny warning
    return(hm.sam.beta)
  } else {
    hm.beta$elementId <- NULL # To suppress a shiny warning
    return(hm.beta)
  }
}
plotBetaBoxplotServerButton3 <- eventReactive(input$beta_boxplot,{
  do_beta_heatmap()
})
output$BetaDiversityHeatmap <- renderPlotly({
  plotBetaBoxplotServerButton3()
})

# Beta diversity boxplots
plotBetaBoxplotServer <- function() {
  shinyInput <- vals$shiny.input
  physeq1 <- shinyInput$pstat

  if (input$taxl.beta !="no rank")  {
    physeq1 <- tax_glom(physeq1, input$taxl.beta)
  }

  meta.data <- physeq1@sam_data
  meta.data$sample.name <- rownames(meta.data)
  colnames(meta.data)[which(colnames(meta.data) == input$select_beta_condition)] <- "condition"

  if (input$select_beta_div_method == "bray") {
      #First get otu_table and transpose it:
      dist.matrix <- t(data.frame(otu_table(physeq1)))
      #Then use vegdist from vegan to generate a bray distance object:
      dist.tmp <- vegdist(dist.matrix, method = "bray")
  } else {
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
  p$elementId <- NULL # To suppress a shiny warning
  p
}
plotBetaBoxplotServerButton <- eventReactive(input$beta_boxplot,{
  plotBetaBoxplotServer()
})
output$BetaDiversityBoxplot <- renderPlotly({
  plotBetaBoxplotServerButton()
})

# Beta diversity table
do_beta_table <- function() {
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
}
plotBetaBoxplotServerButton2 <- eventReactive(input$beta_boxplot,{
  do_beta_table()
})
output$table.beta <- DT::renderDataTable({
  plotBetaBoxplotServerButton2()
}, options=list(paging = TRUE, scrollX = TRUE))

# Download beta diversity table
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

do_beta_stats <- function() {
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
      for (i in 1:length(result.list)){
        output.tmp <- c(result.list[[i]]$method, result.list[[i]]$p.value)
        output.table <- cbind(output.table, output.tmp)
      }
      rownames(output.table) <- c("Method", "P-value")
      colnames(output.table) <- group.name
      output.table
    } else {
      tmp <- kruskal.test(list(dist.within.a, dist.within.b, dist.between))
      output <- c(tmp$method, tmp$p.value)
      output.table <- data.frame(output)
      rownames(output.table) <- c("Method", "P-value")
      output.table
    }
  }
}
plotBetaBoxplotServerButton3 <- eventReactive(input$beta_boxplot, {
  do_beta_stats()
})
output$beta.stat.test <- DT::renderDataTable({
  plotBetaBoxplotServerButton3()
}, options = list(sDom  = '<"top">t<"bottom">ip'))
