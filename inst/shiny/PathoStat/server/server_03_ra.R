log.cpm <- function(df){
  df.cpm <- apply(df, 2, function(y) log10(y*1e6/sum(y) + 1))
  return(df.cpm)
}

count.to.ra <- function(otu_count) {
  otu_relabu <- otu_count
  for (i in 1:ncol(otu_count))  {
    otu_relabu[,i] <- otu_relabu[,i]/sum(otu_relabu[,i])
  }
  return(otu_relabu)
}

# Trying to not use the pstat object directly
pstat.extraction <- function(p) {
  tables <- list()

  # Taxon Levels
  # unranked taxon x classifications
  tables$TAX <- as.data.frame(p@tax_table)

  # Counts
  # unranked taxon x samples
  tables$OTU <- as.data.frame(p@otu_table)

  # Relative Abundance
  # unranked taxon x samples
  tables$RLA <- as.data.frame(count.to.ra(p@otu_table))

  # Sample Data
  # samples x traits
  tables$SAM <- as.data.frame(p@sam_data)
  return(tables)
}

# Takes the relative abundances at an unranked level and
# upsamples them to higher taxonomical levels
upsample.ra <- function(ra.unranked, tax.table, newlev) {
  df = ra.unranked
  df$newlev = tax.table[[newlev]]
  df.melt = melt(df, id.vars = c("newlev"))
  df.melt.agg = aggregate(.~variable+newlev, data=df.melt, FUN=sum)
  df.ra = dcast(df.melt.agg, variable~newlev)
  rownames(df.ra) = df.ra$variable
  df.ra$variable = NULL
  return(df.ra)
}

# Used to dynamically adjust plot height
output$dynamic_ra_plot <- renderUI({
  height = paste(input$plot_sra_height, "px", sep="")
  plotlyOutput("ra_plot", width="800px", height=height)
})

# Used to dynamically generate selectable organisms based on taxlev
output$sra_order_organisms <- renderUI({
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat
  tables <- pstat.extraction(pstat)
  TAX_TABLE <- tables$TAX
  choices <- unique(TAX_TABLE[[input$sra_taxlev]])
  selectizeInput('sra_order_organisms', label='Reorder Organisms', choices=choices, multiple=TRUE)
})

plot_ra <- function() {

  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat

  # Future proofing this plot
  tables <- pstat.extraction(pstat)
  TAX_TABLE <- tables$TAX
  RA_TABLE <- tables$RLA
  SAM_DATA <- tables$SAM

  # Isolate samples of interest
  if (!is.null(input$sra_isolate_samples)) {
    SAM_DATA <- SAM_DATA[rownames(SAM_DATA) %in% input$sra_isolate_samples,,drop=FALSE] 
    RA_TABLE <- RA_TABLE[,colnames(RA_TABLE) %in% input$sra_isolate_samples,drop=FALSE] 
  }

  # Sum by taxon level
  df.ra <- upsample.ra(RA_TABLE, TAX_TABLE, input$sra_taxlev)

  # If grouping is selected
  if (input$group_samples & !is.null(input$gra_select_conditions)) {
    if (input$gra_select_conditions == "All") {
      df.ra$covariate <- rep("All", nrow(df.ra))
    } else {
      df.ra$covariate <- SAM_DATA[[input$gra_select_conditions]] 
    }
    df.ra.melted <- melt(df.ra, id.vars = "covariate")
    df.avg.ra <- aggregate( . ~ variable + covariate , data = df.ra.melted, mean)
    df.avg.ra <- dcast(data = df.avg.ra,formula = covariate~variable)
    rownames(df.avg.ra) <- df.avg.ra$covariate
    df.sam <- df.avg.ra[,"covariate",drop=FALSE]
    colnames(df.sam) <- input$gra_select_conditions
    df.avg.ra$covariate <- NULL
    df.avg.ra <- df.avg.ra[,order(colSums(df.avg.ra))]
    df.ra <- df.avg.ra
  }

  # Reorder by most prominent organisms
  df.ra <- df.ra[,order(colSums(df.ra))]

  # Put selected organisms first
  if (!is.null(input$sra_order_organisms)) {
    organisms.order <- c(setdiff(colnames(df.ra), input$sra_order_organisms), rev(input$sra_order_organisms))
    df.ra <- df.ra[,organisms.order]
  }

  # Order samples by organisms if not by conditons
  if (input$sra_sort_by == "organisms") {
    for (i in 1:ncol(df.ra)) {
      df.ra <- df.ra[order(df.ra[,i]),]
    }
  }

  # If any conditions are selected make a side bar
  if (!is.null(input$sra_select_conditions) || (input$group_samples & input$gra_select_conditions != "All")) {
    
    if (!input$group_samples) {
      df.sam <- SAM_DATA[,input$sra_select_conditions,drop=FALSE]
    }

    # Order samples by conitions if not by organisms
    if (input$sra_sort_by == "conditions") {
      for (i in ncol(df.sam):1) {
        df.sam <- df.sam[order(df.sam[[i]]),,drop=FALSE]
      }
      # Reorder stacked barplot
      df.ra <- df.ra[order(match(rownames(df.ra), rownames(df.sam))),,drop=FALSE]
    } else {
      df.sam <- df.sam[order(match(rownames(df.sam), rownames(df.ra))),,drop=FALSE]
    }

    # Retain hover-text information before conditions are factorized
    hover.txt <- c()
    for (i in 1:ncol(df.sam)) {
      hover.txt <- cbind(hover.txt, df.sam[[i]])
    }

    # Plotly | Heatmap
    df.sam[] <- lapply(df.sam, factor)
    m <- data.matrix(df.sam)
    m.row.normalized <- apply(m, 2, function(x)(x-min(x))/(max(x)-min(x)))
    hm <- plot_ly(x = colnames(m), y = rownames(m), z = m.row.normalized, 
                  type = "heatmap",
                  showscale=FALSE,
                  hoverinfo = "x+y+text",
                  text=hover.txt) %>%
           layout(xaxis = list(title = "", tickangle = -45),
                  yaxis = list(showticklabels = FALSE, type = 'category', ticks = ""))
  }

  # Plotly | Stacked Bar Plots
  df.plot <- df.ra
  df.plot$samples <- rownames(df.plot)
  sbp <- plot_ly(df.plot, y = ~samples, x = df.plot[[colnames(df.plot)[1]]], 
                 type = 'bar', 
                 orientation = 'h', 
                 name = substr(colnames(df.plot)[1], 1, 40)) %>%
          layout(font = list(size = 10),
                 yaxis = list(title = '', type = 'category',
                              tickmode = "array",
                              tickvals = rownames(df.plot),
                              showticklabels = FALSE,
                              categoryorder = 'trace'),
                 xaxis = list(title = 'Relative Abundance'),
                 barmode = 'stack',
                 showlegend = input$sra_show_legend)
  for (i in 2:(ncol(df.plot)-1)) {
    sbp <- add_trace(sbp, x=df.plot[[colnames(df.plot)[i]]], name=substr(colnames(df.plot)[i], 1, 40))
  } 

  # Create a multiplot if any conditions are selected
  if (!is.null(input$sra_select_conditions) || (input$group_samples & input$gra_select_conditions != "All")) {
    hm.sbp <- subplot(hm, sbp, widths=c(0.1,  0.9))
    hm.sbp$elementId <- NULL # To suppress a shiny warning
    return(hm.sbp)
  } else {
    sbp$elementId <- NULL # To suppress a shiny warning
    return(sbp)
  }
}

# Only plots if button is pressed
do_plot_ra <- eventReactive(input$plot_sra, {
  plot_ra()
})
output$ra_plot <- renderPlotly({
  do_plot_ra()
})


# Used to dynamically adjust plot height
output$dynamic_hmra_plot <- renderUI({
  height = paste(input$plot_hmra_height, "px", sep="")
  plotlyOutput("hmra_plot", width="800px", height=height)
})

# Used to dynamically generate selectable organisms based on taxlev
output$hmra_isolate_organisms <- renderUI({
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat
  tables <- pstat.extraction(pstat)
  TAX_TABLE <- tables$TAX
  choices <- unique(TAX_TABLE[[input$hmra_taxlev]])
  selectizeInput('hmra_isolate_organisms', label='Isolate Organisms', choices=choices, multiple=TRUE)
})

plot_hmra <- function() {

  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat

  # Future proofing this plot
  tables <- pstat.extraction(pstat)
  TAX_TABLE <- tables$TAX
  OTU_TABLE <- tables$OTU
  RA_TABLE <- tables$RLA
  SAM_DATA <- tables$SAM

  # Isolate samples of interest
  if (!is.null(input$hmra_isolate_samples)) {
    SAM_DATA <- SAM_DATA[rownames(SAM_DATA) %in% input$hmra_isolate_samples,,drop=FALSE] 
    OTU_TABLE <- OTU_TABLE[,colnames(OTU_TABLE) %in% input$hmra_isolate_samples,drop=FALSE] 
    RA_TABLE <- RA_TABLE[,colnames(RA_TABLE) %in% input$hmra_isolate_samples,drop=FALSE] 
  }

  if (input$hmra_logcpm) {
    df.ra = OTU_TABLE
  } else {
    df.ra = RA_TABLE
  }

  # Sum by taxon level
  df.ra <- upsample.ra(df.ra, TAX_TABLE, input$hmra_taxlev)

  if (input$hmra_logcpm) {
    df.ra = log.cpm(df.ra)
  }

  if (!is.null(input$hmra_isolate_organisms)) {
    df.ra <- df.ra[,input$hmra_isolate_organisms,drop=FALSE]
  }

  # Reorder by most prominent organisms
  df.ra <- df.ra[,order(colSums(df.ra)),drop=FALSE]

  # Order samples by organisms if not by conditons
  if (input$hmra_sort_by == "organisms") {
    for (i in 1:ncol(df.ra)) {
      df.ra <- df.ra[order(df.ra[,i]),,drop=FALSE]
    }
  }

  if (!is.null(input$hmra_select_conditions)) {
    df.sam <- SAM_DATA[,input$hmra_select_conditions,drop=FALSE]
    if (input$hmra_sort_by == "conditions") {
      for (i in ncol(df.sam):1) {
        df.sam <- df.sam[order(df.sam[[i]]),,drop=FALSE]
      }
      # Reorder stacked barplot
      df.ra <- df.ra[order(match(rownames(df.ra), rownames(df.sam))),,drop=FALSE]
    } else {
      df.sam <- df.sam[order(match(rownames(df.sam), rownames(df.ra))),,drop=FALSE]
    }
  }

  m <- data.matrix(df.ra)
  hover.txt <- c()
  for (i in 1:ncol(df.ra)) {
    hover.txt <- cbind(hover.txt, df.ra[[i]])
  }
  hm.ra <- plot_ly(x = colnames(m), y = rownames(m), z = m,
                type = "heatmap",
                colors= "RdPu",
                hoverinfo = "x+y+z+text",
                text=hover.txt) %>%
    layout(xaxis = list(showticklabels = FALSE, title = "", ticks = "", tickangle = -45),
           yaxis = list(showticklabels = FALSE, type = 'category', ticks = ""))


  if (!is.null(input$hmra_select_conditions)) {
    hover.txt <- c()
    for (i in 1:ncol(df.sam)) {
      hover.txt <- cbind(hover.txt, df.sam[[i]])
    }
    df.sam[] <- lapply(df.sam, factor)
    m <- data.matrix(df.sam)
    m.row.normalized <- apply(m, 2, function(x)(x-min(x))/(max(x)-min(x)))
    hm.sam <- plot_ly(x = colnames(m), y = rownames(m), z = m.row.normalized, 
                  type = "heatmap",
                  showscale=FALSE,
                  hoverinfo = "x+y+text",
                  text=hover.txt) %>%
      layout(xaxis = list(title = "", tickangle = -45),
             yaxis = list(showticklabels = FALSE, type = 'category', ticks = ""))
  }

  # Create a multiplot if any conditions are selected
  if (!is.null(input$hmra_select_conditions)) {
    hm.sam.ra <- subplot(hm.sam, hm.ra, widths=c(0.1,  0.9))
    hm.sam.ra$elementId <- NULL # To suppress a shiny warning
    return(hm.sam.ra)
  } else {
    hm.ra$elementId <- NULL # To suppress a shiny warning
    return(hm.ra)
  }
}

# Only plots if button is pressed
do_plot_hmra <- eventReactive(input$plot_hmra, {
  plot_hmra()
})
output$hmra_plot <- renderPlotly({
  do_plot_hmra()
})

###### Boxplots Moved Here

output$single_species_ui <- renderUI({
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat
  species.name.vec <- TranslateIdToTaxLevel(pstat, rownames(pstat@otu_table@.Data), input$taxl_single_species)
  tagList(
    selectInput("select_single_species_name_plot", 
      "Select names (support multiple)", 
      species.name.vec, selected = species.name.vec[1], multiple = TRUE)
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
