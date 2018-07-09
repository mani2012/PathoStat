# independent heatmap plotting function in the server using specific data from input
plotPCAPlotlyServer <- function(){
  shinyInput <- vals$shiny.input
  physeq1 <- shinyInput$pstat

  # Use shapes is optional
  if (input$select_pca_shape == "None") {
    condition.shape.vec <- NULL
  } else {
    condition.shape.vec <- physeq1@sam_data[[input$select_pca_shape]]
  }

  if (input$taxl.pca=="no rank")  {
    df.plot <- physeq1@otu_table@.Data
    if (input$select_pca_data_format == "log10 CPM"){
      df.plot <- getLogCPM(df.plot)
    } else if (input$select_pca_data_format == "RA"){
      df.plot <- getRelativeAbundance(df.plot)
    }
    p <- plotPCAPlotly(df.input = df.plot,
                       condition.color.vec = physeq1@sam_data[[input$select_pca_color]],
                       condition.color.name = input$select_pca_color,
                       condition.shape.vec = condition.shape.vec,
                       condition.shape.name = input$select_pca_shape,
                       pc.a = paste("PC", input$xcol.new, sep = ""),
                       pc.b = paste("PC", input$ycol.new, sep = ""),
                       columnTitle = paste("PCA with colorbar representing", input$select_pca_color, sep = " "))


   #p$condition.shape.vec = physeq1@sam_data[[input$select_pca_shape]]
   #p$condition.shape.name = input$select_pca_shape

  } else  {
    physeq2 <- tax_glom(physeq1, input$taxl.pca)
    df.plot <- physeq2@otu_table@.Data
    if (input$select_pca_data_format == "log10 CPM") {
      df.plot <- getLogCPM(df.plot)
    } else if (input$select_pca_data_format == "RA"){
      df.plot <- getRelativeAbundance(df.plot)
    }
    p <- plotPCAPlotly(df.input = df.plot,
                       condition.color.vec = physeq2@sam_data[[input$select_pca_color]],
                       condition.color.name = input$select_pca_color,
                       condition.shape.vec = condition.shape.vec,
                       condition.shape.name = input$select_pca_shape,
                       pc.a = paste("PC", input$xcol.new, sep = ""),
                       pc.b = paste("PC", input$ycol.new, sep = ""),
                       columnTitle = paste("PCA with colorbar representing", input$select_pca_color, sep = " ")
    )
  }
  p$elementId <- NULL
  return(p)
}
# Show plot after hitting the plot button
plotPCAPlotlyServerButton <- eventReactive(input$DR_plot,{
  plotPCAPlotlyServer()
})
output$pca.plotly <- renderPlotly({
  plotPCAPlotlyServerButton()
})

# interactive PCA table
do_PCA_table <- function() {
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
  colnames(table.output.pca) = c("Standard deviation", "Variance Explained", "Cumulative Variance")
  table.output.pca[,2] <- scales::percent(as.numeric(table.output.pca[,2]))
  table.output.pca[,3] <- scales::percent(as.numeric(table.output.pca[,3]))
  #hide std
  DT::datatable(table.output.pca[,-1], options = list(sDom  = '<"top">t<"bottom">ip'))
}
# Show table after hitting the plot button
tablePCAServerButton <- eventReactive(input$DR_plot,{
  do_PCA_table()
})
output$PCAtable <- DT::renderDataTable({
  tablePCAServerButton()
})

# PCoA
plotPCoAPlotlyServer <- function(){
  shinyInput <- vals$shiny.input
  physeq1 <- shinyInput$pstat
  physeq1 <- phyloseq(otu_table(physeq1), phy_tree(physeq1),
                      tax_table(physeq1), sample_data(physeq1))

  # Use shapes is optional
  if (input$select_pca_shape == "None") {
    condition.shape.vec <- NULL
  } else {
    condition.shape.vec <- physeq1@sam_data[[input$select_pca_shape]]
  }

  if (input$taxl.pca=="no rank")  {
    p <- plotPCoAPlotly(physeq.input = physeq1,
                        condition.color.vec = physeq1@sam_data[[input$select_pca_color]],
                        condition.color.name = input$select_pca_color,
                        condition.shape.vec = condition.shape.vec,
                        condition.shape.name = input$select_pca_shape,
                        method = input$pcoa.method,
                        pc.a = paste("Axis", input$xcol.new, sep = "."),
                        pc.b = paste("Axis", input$ycol.new, sep = "."),
                        columnTitle = paste("PCoA with colorbar representing", input$select_pca_color, sep = " ")
    )
  } else  {
    physeq2 <- tax_glom(physeq1, input$taxl.pca)
    p <- plotPCoAPlotly(physeq.input = physeq2,
                        condition.color.vec = physeq2@sam_data[[input$select_pca_color]],
                        condition.color.name = input$select_pca_color,
                        condition.shape.vec = condition.shape.vec,
                        condition.shape.name = input$select_pca_shape,
                        method = input$pcoa.method,
                        pc.a = paste("Axis", input$xcol.new, sep = "."),
                        pc.b = paste("Axis", input$ycol.new, sep = "."),
                        columnTitle = paste("PCoA with colorbar representing", input$select_pca_color, sep = " ")
    )
  }
  p$elementId <- NULL
  return(p)
}
# Show plot after hitting the plot button
plotPCoAPlotlyServerButton <- eventReactive(input$DR_plot,{
  plotPCoAPlotlyServer()
})
output$pcoa.plotly <- renderPlotly({
  plotPCoAPlotlyServerButton()
})

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

# Interactive PCA table
do_PCoA_table <- function() {
  shinyInput <- vals$shiny.input
  physeq1 <- shinyInput$pstat
  ord <- getOrdPCoA()
  df.output <- ord[,c(1,3,5)]
  colnames(df.output) <- c("eigenvalue", "Variance Explained", "Cumulative Variance")
  rownames(df.output) <- paste("Axis", 1:nrow(df.output), sep = ".")
  df.output[,2] <- scales::percent(as.numeric(df.output[,2]))
  df.output[,3] <- scales::percent(as.numeric(df.output[,3]))
  # hide eigenvalue
  DT::datatable(df.output[,-1], options = list(sDom  = '<"top">t<"bottom">ip'))
}
# Show table after hitting the plot button
tablePCoAServerButton <- eventReactive(input$DR_plot,{
  do_PCoA_table()
})
output$PCoAtable <- DT::renderDataTable({
  tablePCoAServerButton()
})
