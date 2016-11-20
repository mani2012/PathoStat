library(shiny)
library(ggvis)
library(d3heatmap)
library(reshape2)
library(limma)
library(edgeR)
library(phyloseq)
library(ape)
library(stats)
library(PathoStat)
library(alluvial)

# Converts decimal percentage to string with specified digits
pct2str <- function(v, digits=2) {sprintf(paste0('%.',digits,'f'), v*100)}

shinyServer(function(input, output, session) {
    # needed information from PathoStat
    shinyInput <- getShinyInput()
    pstat <- shinyInput$pstat
    sample_data <- pstat@sam_data
    #pdata <- data.frame(shinyInput$batch, shinyInput$condition)
    #modmatrix = model.matrix(~as.factor(shinyInput$condition), data=pdata)
    pdata <- data.frame(c(sample_data[,1])[[1]], c(sample_data[,2])[[1]])
    modmatrix = model.matrix(~as.factor(c(sample_data[,2])[[1]]), data=pdata)
    pca <- BatchQC::batchqc_pca(pstat@otu_table, batch=c(sample_data[,1])[[1]], 
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
        #taxcountdata <- findTaxonLevelData(pstat@otu_table, 
        #    shinyInput$taxonLevels, taxonLevel)
        shinyInput <- getShinyInput()
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
        if (is.null(shinyInput$taxcountdata)) {
            shinyInput <- c(shinyInput, list(taxcountdata = taxcountdata))
        } else {
            shinyInput$taxcountdata <- taxcountdata
        }
        if (is.null(shinyInput$taxdata)) {
            shinyInput <- c(shinyInput, list(taxdata = taxdata))
        } else {
            shinyInput$taxdata <- taxdata
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
    
    findNormalizedCount <- function() {
        shinyInput <- getShinyInput()
        if(input$norm == 'Quantile coreOTU Normalization')  {
            dat <- coreOTUQuantile(shinyInput$taxcountdata, input$otuthreshold, 
                input$prevalence)
        } else if(input$norm == 'Library Size Scaling')  {
            dat <- sizeNormalize(shinyInput$taxcountdata)
        } else  {
            dat <- coreOTUNormalize(shinyInput$taxcountdata, input$ebweight,
                input$otuthreshold, input$prevalence)
        }
        dat
    }
    
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
    
    output$Heatmap <- renderPlot({
        physeq1 <- shinyInput$pstat
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
        physeq1 <- shinyInput$pstat
        cn <- colnames(physeq1@sam_data)
        cn[1] <- "batch"
        cn[2] <- "condition"
        colnames(physeq1@sam_data) <- cn
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
        physeq1 <- shinyInput$pstat
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
        physeq1 <- shinyInput$pstat
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
            fill = if (input$colbybatchPCA) ~factor(c(sample_data[,1])[[1]]) 
            else ~factor(c(sample_data[,2])[[1]]), `:=`(key, ~id)) %>% 
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
            shinyInput$vars, c(sample_data[,2])[[1]], c(sample_data[,1])[[1]])
        explained_variation
    })
    
    output$PCoAplot <- renderPlot({
        physeq1 <- shinyInput$pstat
        cn <- colnames(physeq1@sam_data)
        cn[1] <- "batch"
        cn[2] <- "condition"
        colnames(physeq1@sam_data) <- cn
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
        findTaxCountDataDE()
        physeq1 <- shinyInput$pstat
        Inbatch <- sample_data(physeq1)[[input$secondary]]
        shinyInput <- getShinyInput()
        dat <- shinyInput$taxcountdata
        dat <- findNormalizedCount()
        # lcpm <- log2CPM(shinyInput$taxcountdata)
        lcpm <- log2CPM(dat)
        lcounts <- lcpm$y
        dat <- lcounts
        batch1 <- as.factor(Inbatch)
        batch2 <- split(which(Inbatch == batch1), batch1)
        batch3 <- unlist(lapply(1:length(batch2), 
            function(x) batch2[[x]][1:input$noSamples]))
        dat1 <- dat[, batch3]
        colnames(dat1) <- seq(1:ncol(dat))[batch3]
        dat1
    })
    DE <- reactive({
        physeq1 <- shinyInput$pstat
        Incondition <- sample_data(physeq1)[[input$primary]]
        findTaxCountDataDE()
        shinyInput <- getShinyInput()
        dat <- shinyInput$taxcountdata
        dat <- findNormalizedCount()
        # lcpm <- log2CPM(shinyInput$taxcountdata)
        lcpm <- log2CPM(dat)
        lcounts <- lcpm$y
        dat <- lcounts
        cond1 <- as.factor(Incondition)
        cond2 <- split(which(Incondition == cond1), cond1)
        cond3 <- unlist(lapply(1:length(cond2), 
            function(x) cond2[[x]][1:input$ncSamples]))
        dat1 <- dat[, cond3]
        colnames(dat1) <- seq(1:ncol(dat))[cond3]
        dat1
    })
    diffex_bp <- reactive({
        physeq1 <- shinyInput$pstat
        Incondition <- sample_data(physeq1)[[input$primary]]
        Inbatch <- sample_data(physeq1)[[input$secondary]]
        if (input$sortbybatch) {
            batch4 <- split(Inbatch, as.factor(Inbatch))
            batch5 <- unlist(lapply(1:length(batch4),
                function(x) batch4[[x]][1:input$noSamples]))
            dat1 <- BP()
            dat2 <- melt(as.data.frame(dat1), measure.var = colnames(dat1))
            dat2$batch <- as.factor(unlist(lapply(1:length(batch5),
                function(x) rep(batch5[x], nrow(dat1)))))
            dat2$condition <- as.factor(unlist(lapply(as.numeric(colnames(dat1))
                , function(x) rep(Incondition[x], nrow(dat1)))))
            dat2$samples <- unlist(lapply(seq(ncol(dat1)),
                function(x) rep(seq(ncol(dat1))[x], nrow(dat1))))
        } else {
            cond4 <- split(Incondition,
                as.factor(Incondition))
            cond5 <- unlist(lapply(1:length(cond4),
                function(x) cond4[[x]][1:input$ncSamples]))
            dat1 <- DE()
            dat2 <- melt(as.data.frame(dat1), measure.var = colnames(dat1))
            dat2$condition <- as.factor(unlist(lapply(1:length(cond5),
                function(x) rep(cond5[x], nrow(dat1)))))
            dat2$batch <- as.factor(unlist(lapply(as.numeric(colnames(dat1)),
                function(x) rep(Inbatch[x], nrow(dat1)))))
            dat2$samples <- unlist(lapply(seq(ncol(dat1)),
                function(x) rep(seq(ncol(dat1))[x], nrow(dat1))))
        }
        dat2 %>% group_by(batch) %>% ggvis(~samples, ~value, fill =
                if (input$colbybatch) ~batch else ~condition) %>%
            layer_boxplots() %>%
            add_tooltip(function(dat2) { paste0("Sample: ", 
                colnames(shinyInput$pstat@otu_table)[dat2$samples],
                "<br>", if (input$colbybatch) input$secondary 
                else input$primary, ": ",
                if (input$colbybatch) dat2$batch else dat2$condition)
            }, "hover") %>%
            add_axis("x", title = if (input$sortbybatch)
                paste(input$noSamples, "Sample(s) Per ", input$secondary,
                    sep = " ")
                else
                    paste(input$ncSamples, "Sample(s) Per ", input$primary, 
                        sep=" "),
                properties = axis_props(title = list(fontSize = 15),
                labels = list(fontSize = 5, angle = 90))) %>%
            add_axis("y", title = "Expression", properties = axis_props(title =
                list(fontSize = 15),labels = list(fontSize = 10))) %>%
            add_legend("fill", title = if (input$colbybatch)
                input$secondary else input$primary, 
                properties = legend_props(title =
                list(fontSize = 15), labels = list(fontSize = 10))) %>%
            set_options(width="auto", height="auto")
    })
    diffex_bp %>% bind_shiny("DiffExPlot")
    output$DEsummary <- renderPrint({
        if (input$sortbybatch) {
            summary(BP())
        } else {
            summary(DE())
        }
    })
    
    output$DEtable <- renderTable({
        if (input$sortbybatch) {
            BP()
        } else {
            DE()
        }
    })
    
    DELimnorm <- reactiveValues()
    observeEvent(c(input$taxlde, input$apply), {
        findAllTaxData(input$taxlde)
        dat <- findNormalizedCount()
        lcpm <- log2CPM(dat)
        lcounts <- lcpm$y
        DELimnorm$data <- lcounts
    })
    output$LimmaTable <- renderTable({
        shinyInput <- getShinyInput()
        physeq1 <- shinyInput$pstat
        Incondition <- sample_data(physeq1)[[input$primary]]
        Inbatch <- sample_data(physeq1)[[input$secondary]]
        pdata <- data.frame(Inbatch, Incondition)
        mod <- model.matrix(if (input$primary == input$secondary) 
            ~as.factor(Incondition)
            else ~as.factor(Incondition) + ~as.factor(Inbatch), data = pdata)
        dat1 <- DELimnorm
        fit <- lmFit(dat1$data, mod)
        fit2 <- eBayes(fit)
        ncond <- nlevels(as.factor(Incondition))
        limmaTable <- topTable(fit2, coef = 2:ncond, number = input$noTaxons)
        for (j in 2:ncond)  {
            colnames(limmaTable)[j-1] <- paste("Primary Covariate: ", 
                levels(as.factor(Incondition))[j], " (logFC)", sep='')
        }
        limmaTable
    }, rownames = TRUE)

    DEEdgeRnorm <- reactiveValues()
    observeEvent(c(input$taxlde, input$apply),  {
        findAllTaxData(input$taxlde)
        DEEdgeRnorm$data <- findNormalizedCount()
    })
    output$EdgeRTable <- renderTable({
        shinyInput <- getShinyInput()
        physeq1 <- shinyInput$pstat
        Incondition <- sample_data(physeq1)[[input$primary]]
        Inbatch <- sample_data(physeq1)[[input$secondary]]
        pdata <- data.frame(Inbatch, Incondition)
        group <- factor(Incondition)
        mod <- model.matrix(if (input$primary == input$secondary) 
            ~as.factor(Incondition)
            else ~as.factor(Incondition) + ~as.factor(Inbatch), data = pdata)
        dat1 <- DEEdgeRnorm
        y <- DGEList(counts=dat1$data,group=group)
        y <- calcNormFactors(y)
        y <- estimateGLMCommonDisp(y,mod)
        y <- estimateGLMTrendedDisp(y,mod)
        y <- estimateGLMTagwiseDisp(y,mod)
        fit <- glmFit(y,mod)
        ncond <- nlevels(group)
        lrt <- glmLRT(fit,coef=2:ncond)
        edgeRTags <- topTags(lrt, n=input$noTaxons)
        edgeRTable <- edgeRTags$table
        for (j in 2:ncond)  {
            colnames(edgeRTable)[j-1] <- paste("Primary Covariate: ",
                levels(group)[j], " (logFC)", sep='')
        }
        edgeRTable
    }, rownames = TRUE)

    output$confRegion <- renderPlot({
        p1 <- shinyInput$pstat@otu_table[input$taxon1, input$sample]
        if (p1 <= 0) p1 <- 1
        p2 <- shinyInput$pstat@otu_table[input$taxon2, input$sample]
        if (p2 <= 0) p2 <- 1
        size <- sum(shinyInput$pstat@otu_table[,input$sample])
        plotConfRegion(p1, p2, size, uselogit=input$uselogit)
    })

    #Time Series
    output$Allustax <- renderUI({
      checkboxGroupInput(inputId="Allustax", label="Taxa of interest ", 
                         choices = as.character(unique(unlist((
                           shinyInput$pstat@tax_table)[,input$Alluglom]))))
      
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

