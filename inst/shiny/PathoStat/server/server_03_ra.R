


  # Trying to not use the pstat object directly
  pstat.extraction <- function(p) {
    tables <- list()

    # Taxon Levels
    # unranked taxon x classifications
    tables$TAX <- as.data.frame(p@tax_table)

    # Relative Abundance
    # unranked taxon x samples
    tables$RLA <- as.data.frame(findRAfromCount(p@otu_table))

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

  # Used to dynamically generate selectable organisms based on taxlev
  output$sra_order_organisms <- renderUI({
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    tables <- pstat.extraction(pstat)
    TAX_TABLE <- tables$TAX
    choices <- unique(TAX_TABLE[[input$sra_taxlev]])
    selectizeInput('sra_order_organisms', label='Reorder Organisms', choices=choices, multiple=TRUE)
  })

  output$ra_plot <- renderPlotly({

    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat

    # Future proofing this plot
    tables <- pstat.extraction(pstat)
    TAX_TABLE <- tables$TAX
    RA_TABLE <- tables$RLA
    SAM_DATA <- tables$SAM

    # Sum by taxon level
    df.ra <- upsample.ra(RA_TABLE, TAX_TABLE, input$sra_taxlev)

    # If grouping is selected
    if (input$group_samples & !is.null(input$gra_select_conditions)) {
      df.ra$covariate <- SAM_DATA[[input$gra_select_conditions]]
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
    if (!is.null(input$sra_select_conditions) || input$group_samples) {

      if (!input$group_samples) {
        df.sam <- SAM_DATA[,input$sra_select_conditions]
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
    if (!is.null(input$sra_select_conditions) || input$group_samples) {
      hm.sbp <- subplot(hm, sbp, widths=c(0.1,  0.9))
      hm.sbp$elementId <- NULL # To suppress a shiny warning
      return(hm.sbp)
    } else {
      sbp$elementId <- NULL # To suppress a shiny warning
      return(sbp)
    }
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
