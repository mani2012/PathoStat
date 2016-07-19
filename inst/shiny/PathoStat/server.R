library(shiny)
library(ggvis)
library(d3heatmap)
library(reshape2)
library(limma)
library(phyloseq)
library(ape)
library(PathoStat)

# Converts decimal percentage to string with specified digits
pct2str <- function(v, digits=2) {sprintf(paste0('%.',digits,'f'), v*100)}

shinyServer(function(input, output, session) {
    # needed information from PathoStat
    shinyInput <- getShinyInput()
    pdata <- data.frame(shinyInput$batch, shinyInput$condition)
    modmatrix = model.matrix(~as.factor(shinyInput$condition), data=pdata)
    pca <- BatchQC::batchqc_pca(shinyInput$countdata, batch=shinyInput$batch, 
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
        taxdata <- findTaxonLevelData(shinyInput$data, shinyInput$taxonLevels, 
            taxonLevel)
        if (is.null(shinyInput$taxdata)) {
            shinyInput <- c(shinyInput, list(taxdata = taxdata))
        } else {
            shinyInput$taxdata <- taxdata
        }
        taxcountdata <- findTaxonLevelData(shinyInput$countdata, 
            shinyInput$taxonLevels, taxonLevel)
        if (is.null(shinyInput$taxcountdata)) {
            shinyInput <- c(shinyInput, list(taxcountdata = taxcountdata))
        } else {
            shinyInput$taxcountdata <- taxcountdata
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
    
    tax_ra_bp <- reactive({
        if (is.null(input$taxl)) {
            return()
        }
        taxdata <- findTaxData()
        dat <- melt(cbind(taxdata, ind = as.character(rownames(taxdata))), 
            id.vars = c("ind"))
        dat %>% ggvis(x = ~variable, y = ~value, fill = ~as.factor(ind)) %>% 
        layer_bars(stack = TRUE) %>% 
        add_tooltip(function(dat2) {
            paste0("Sample: ", dat2[2], "<br />", "Genome: ", dat2[1], 
                "<br />", "RA: ", round(dat2[4] - dat2[3], 4))
        }, "hover") %>% 
        # add_axis('x', subdivide = 1, 
        #   values = 1:length(colnames(shinyInput$data)),
        add_axis("x", title = "Samples", properties = axis_props(title = 
            list(fontSize = 15), labels = list(text = "", 
            fontSize = 10))) %>% 
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
        tidyr::separate(fullname, c('fi1', 'taxid', 'fi2', input$taxl), sep='\\|') %>%
        dplyr::select_(.dots=c(as.name(input$taxl), 'taxid', colnames(tmp)))
    }
    
    # Render relative abundance table    
    output$TaxRAtable <- DT::renderDataTable(format_RA_table(findTaxData()), options=dtopts, rownames=F)    
    
    # Format count table:
    # Just expands taxa name
    format_Count_table <- function(tmp) {
      tmp %>% dplyr::add_rownames("fullname") %>%
        tidyr::separate(fullname, c('fi1', 'taxid', 'fi2', input$taxl), sep='\\|') %>%
        dplyr::select_(.dots=c(as.name(input$taxl), 'taxid', colnames(tmp)))
    }
    # Render count table
    output$TaxCountTable <- DT::renderDataTable(format_Count_table(findTaxCountData()), options=dtopts, rownames=F)
    
    output$downloadData <- downloadHandler(filename = function() {
        paste0("sample_data_", input$taxl, ".csv", sep = "")
    }, content = function(file) {
        write.csv(shinyInput$taxdata, file)
    })
    
    output$downloadCountData <- downloadHandler(filename = function() {
        paste0("sample_data_count_", input$taxl, ".csv", sep = "")
    }, content = function(file) {
        write.csv(shinyInput$taxcountdata, file)
    })
    
    findPhyseqData <- function() {
        ids <- rownames(shinyInput$data)
        taxmat <- findTaxonMat(ids, shinyInput$taxonLevels)
        OTU <- otu_table(shinyInput$countdata, taxa_are_rows = TRUE)
        TAX <- tax_table(taxmat)
        physeq <- phyloseq(OTU, TAX)
        sampledata = sample_data(data.frame(condition=as.factor(
            shinyInput$condition), batch=as.factor(shinyInput$batch), 
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
    
    output$Heatmap <- renderPlot({
        physeq1 <- findPhyseqData()
        setGgplotTheme()
        #plot_heatmap(physeq1, sample.label="condition")
        if (input$taxl=="no rank")  {
            plot_heatmap(physeq1)
        } else  {
            physeq2 <- tax_glom(physeq1, input$taxl)
            plot_heatmap(physeq2, taxa.label=input$taxl)
        }
    })
    output$AlphaDiversity <- renderPlot({
        physeq1 <- findPhyseqData()
        #alpha_meas <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", 
        #    "InvSimpson")
        #alpha_meas <- c("Observed", "Chao1", "Shannon", "Simpson", 
        #    "InvSimpson")
        alpha_meas <- c("Shannon", "Simpson", "InvSimpson")
        setGgplotTheme()
        (p <- plot_richness(physeq1, "condition", "batch", measures=alpha_meas))
        p + ggplot2::geom_boxplot(data=p$data, ggplot2::aes(x=condition, 
            y=value, color=NULL), alpha=0.1)
    })
    output$BetaDiversity <- renderPlot({
        physeq1 <- findPhyseqData()
        setGgplotTheme()
        if (input$methodBeta)  {
            dist = distance(physeq1, method = "wUniFrac")
            titleString="Beta Diversity Distance: Weigthed Unifrac"
        } else  {
            dist = distance(physeq1, method = "bray")
            titleString="Beta Diversity Distance: Bray-Curtis"
        }
        gplots::heatmap.2(as.matrix(dist), col=gplots::bluered(75), scale="row",
            key=TRUE, symkey=FALSE, density.info="none", trace="none", 
            margins = c(6, 6))
    })
    output$ExploratoryTree <- renderPlot({
        physeq1 <- findPhyseqData()
        setGgplotTheme()
        plot_tree(physeq1, color="condition", label.tips="genus", 
            size="Abundance")
    })
    output$BiPlot <- renderPlot({
        physeq1 <- findPhyseqData()
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
        physeq1 <- findPhyseqData()
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
        data.frame(pc[, c(input$xcol, input$ycol)])
    })
    # interactive PCA plot
    vis_pc <- reactive({
        shinyInput <- getShinyInput()
        pc <- shinyInput$pc
        
        pc$id <- 1:nrow(pc)
        
        all_values <- function(x) {
            if (is.null(x)) 
                return(NULL)
            row <- pc[pc$id == x$id, ]
            paste0(paste0("PC", input$xcol), ": ", signif(with(row, 
                get(names(row)[input$xcol])), 3), "<br>", 
                paste0("PC", input$ycol), ": ", signif(with(row, 
                get(names(row)[input$ycol])), 3), "<br>")
        }
        
        pc %>% 
        ggvis(~get(names(pc)[input$xcol]), ~get(names(pc)[input$ycol]), 
            fill = if (input$colbybatchPCA) ~factor(shinyInput$batch) 
            else ~factor(condition), `:=`(key, ~id)) %>% 
        layer_points(`:=`(size, 75), `:=`(size.hover, 200)) %>% 
        add_tooltip(all_values, "hover") %>% 
        add_axis("x", title = paste0("PC", input$xcol), properties = axis_props(
            title = list(fontSize = 15), labels = list(fontSize = 10))) %>% 
        add_axis("y", title = paste0("PC", input$ycol), properties = axis_props(
            title = list(fontSize = 15), labels = list(fontSize = 10))) %>% 
        add_legend("fill", title = if (input$colbybatchPCA) 
            "Batches" else "Conditions", properties = legend_props(title = 
            list(fontSize = 15), labels = list(fontSize = 10)))%>% 
        set_options(width = "auto", height = "auto")
    })
    # interactive PCA summary
    vis_pc %>% bind_shiny("PCAplot")
    output$PCAsummary <- renderPrint({
        summary(PCA())
    })
    # interactive PCA table
    output$PCAtable <- renderTable({
        PCA()
    })
    output$PCAExplainedVariation <- renderTable({
        PCA()
        shinyInput <- getShinyInput()
        pc <- shinyInput$pc
        pcs <- t(pc)
        explained_variation <- BatchQC::batchqc_pc_explained_variation(pcs, 
            shinyInput$vars, shinyInput$condition, shinyInput$batch)
        explained_variation
    })
    
    output$PCoAplot <- renderPlot({
        physeq1 <- findPhyseqData()
        setGgplotTheme()
        DistBC = distance(physeq1, method = "bray")
        DistUF = distance(physeq1, method = "wUniFrac")
        ordBC = ordinate(physeq1, method = "PCoA", distance = DistBC)
        ordUF = ordinate(physeq1, method = "PCoA", distance = DistUF)
        if (input$colbybatchPCoA) colorstring="batch" 
        else colorstring="condition"
        if (input$methodPCoA)  {
            ordMethod=ordUF
            titleString="PCoA: Weigthed Unifrac"
        } else  {
            ordMethod=ordBC
            titleString="PCoA: Bray-Curtis"
        }
        plot_ordination(physeq1, ordMethod, axes=c(input$xcolA, input$ycolA), 
            color=colorstring) + 
            ggplot2::geom_point(size=3) +
            #ggplot2::geom_point(mapping=ggplot2::aes(size=2)) +
            #ggplot2::geom_point(mapping=ggplot2::aes(shape=factor(batch))) +
            ggplot2::ggtitle(titleString)
    })
    
    # interactive Differential Expression boxplot
    BP <- reactive({
        getShinyInput()
        phyloseq1 <- findPhyseqData()
        lcpm <- log2CPM(otu_table(phyloseq1))
        lcounts <- lcpm$y
        dat <- lcounts
        second.Cov3 <- as.factor(sample_data(phyloseq1)[[input$secondary]])
        second.Cov4 <- split(which(sample_data(phyloseq1)[[input$secondary]] == second.Cov3), second.Cov3)
        second.Cov5 <- unlist(lapply(1:length(second.Cov4), 
            function(x) second.Cov4[[x]][1:input$nsSamples]))
        dat1 <- dat[, second.Cov5]
        colnames(dat1) <- seq(1:ncol(dat))[second.Cov5]
        dat1
    })
    DE <- reactive({
        getShinyInput()
        phyloseq1 <- findPhyseqData()
        lcpm <- log2CPM(otu_table(phyloseq1))
        lcounts <- lcpm$y
        dat <- lcounts
        first.Cov3 <- as.factor(sample_data(phyloseq1)[[input$primary]])
        first.Cov4 <- split(which(sample_data(phyloseq1)[[input$primary]] == first.Cov3), first.Cov3)
        first.Cov5 <- unlist(lapply(1:length(first.Cov4), 
                                     function(x) first.Cov4[[x]][1:input$npSamples]))
        dat1 <- dat[, first.Cov5]
        colnames(dat1) <- seq(1:ncol(dat))[first.Cov5]
        dat1
    })
    #Differential Expression Barplot
    diffex_bp <- reactive({
        phyloseq2 <- findPhyseqData()
        if (input$sortbysecondary) {
            #Define the primary/secondary covariates (set my user), then sort by the secondary
            second.Cov1 <- split(sample_data(phyloseq2)[[input$secondary]], as.factor(sample_data(phyloseq2)[[input$secondary]]))
            second.Cov2 <- unlist(lapply(1:length(second.Cov1),
                function(x) second.Cov1[[x]][1:input$nsSamples]))
            first.Cov1 <- split(sample_data(phyloseq2)[[input$primary]], as.factor(sample_data(phyloseq2)[[input$primary]]))
            first.Cov2 <- unlist(lapply(1:length(first.Cov1),
                                         function(x) first.Cov1[[x]][1:input$npSamples]))
            dat1 <- BP()
            dat2 <- melt(as.data.frame(dat1), measure.var = colnames(dat1))
            dat2$cov2 <- as.factor(unlist(lapply(1:length(second.Cov2),
                function(x) rep(second.Cov2[x], nrow(dat1)))))
            dat2$cov1 <- as.factor(unlist(lapply(as.numeric(colnames(dat1))
                , function(x) rep(second.Cov1[x], nrow(dat1)))))
            dat2$samples <- unlist(lapply(seq(ncol(dat1)),
                function(x) rep(seq(ncol(dat1))[x], nrow(dat1))))
        } else {
            first.Cov1 <- split(sample_data(phyloseq2)[[input$primary]], as.factor(sample_data(phyloseq2)[[input$primary]]))
            first.Cov2 <- unlist(lapply(1:length(first.Cov1),
                                       function(x) first.Cov1[[x]][1:input$npSamples]))
            second.Cov1 <- split(sample_data(phyloseq2)[[input$secondary]], as.factor(sample_data(phyloseq2)[[input$secondary]]))
            second.Cov2 <- unlist(lapply(1:length(second.Cov1),
                                         function(x) second.Cov1[[x]][1:input$nsSamples]))
            dat1 <- DE()
            dat2 <- melt(as.data.frame(dat1), measure.var = colnames(dat1))
            dat2$cov1 <- as.factor(unlist(lapply(1:length(first.Cov2),
                function(x) rep(first.Cov2[x], nrow(dat1)))))
            dat2$cov2 <- as.factor(unlist(lapply(as.numeric(colnames(dat1)),
                function(x) rep(second.Cov2[x], nrow(dat1)))))
            dat2$samples <- unlist(lapply(seq(ncol(dat1)),
                function(x) rep(seq(ncol(dat1))[x], nrow(dat1))))
        }
        dat2 %>% group_by(cov2) %>% ggvis(~samples, ~value, fill =
                if (input$colbysecondary) ~cov2 else ~cov1) %>%
            layer_boxplots() %>%
            add_tooltip(function(dat2) { paste0("Sample: ", 
                colnames(shinyInput$countdata)[dat2$samples],
                "<br>", if (input$colbysecondary) "Secondary covariate: " else "Primary Covariate: ",
                if (input$colbysecondary) dat2$cov2 else dat2$cov1)
            }, "hover") %>%
            add_axis("x", title = if (input$sortbysecondary)
                paste(input$nsSamples, "Sample(s) Per Secondary Covariate", sep = " ")
                else
                    paste(input$npSamples, "Sample(s) Per Primary Covariate", sep=" "),
                properties = axis_props(title = list(fontSize = 15),
                labels = list(fontSize = 5, angle = 90))) %>%
            add_axis("y", title = "Expression", properties = axis_props(title =
                list(fontSize = 15),labels = list(fontSize = 10))) %>%
            add_legend("fill", title = if (input$colbysecondary)
                "Secondary" else "Primary", properties = legend_props(title =
                list(fontSize = 15), labels = list(fontSize = 10)))
    })
    diffex_bp %>% bind_shiny("DiffExPlot")
    output$DEsummary <- renderPrint({
        if (input$sortbysecondary) {
            summary(BP())
        } else {
            summary(DE())
        }
    })
    
    output$DEtable <- renderTable({
        if (input$sortbysecondary) {
            BP()
        } else {
            DE()
        }
    })
    
    output$LimmaTable <- renderTable({
        shinyInput <- getShinyInput()
        pdata <- data.frame(shinyInput$batch, shinyInput$condition)
        mod <- model.matrix(~as.factor(shinyInput$condition) + 
            ~as.factor(shinyInput$batch), data = pdata)
        dat1 <- DE()
        fit <- lmFit(dat1, mod)
        fit2 <- eBayes(fit)
        ncond <- nlevels(as.factor(shinyInput$condition))
        limmaTable <- topTable(fit2, coef = 2:ncond, number = input$noTaxons)
        for (j in 2:ncond)  {
            colnames(limmaTable)[j-1] <- paste("Condition: ", 
                levels(as.factor(shinyInput$condition))[j], " (logFC)", sep='')
        }
        limmaTable
    })
    
})
