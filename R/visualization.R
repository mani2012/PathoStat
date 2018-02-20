#' Plot heatmap with color bar
#'
#' @param df.input Input data object that contains the data to be plotted. Required
#' @param condition.vec.1 color vector. Required
#' @param condition.vec.2 color vector 2. Required
#' @param condition.1.name color vector 1 name. Required
#' @param condition.2.name color vector 2 name. Required
#' @param clusterRow Cluster the rows. The default is TRUE
#' @param clusterCol Cluster the columns. The default is TRUE
#' @param displayRowLabels Display the row labels on the heatmap. The default
#' is TRUE.
#' @param displayColumnLabels Display the column labels on the heatmap. The
#' default is TRUE
#' @param displayRowDendrograms Display the row dendrograms on the heatmap. The
#' default is TRUE
#' @param displayColumnDendrograms Display the column dendrograms on the
#' heatmap. The default is TRUE.
#' @param annotationColors Set of annotation colors for color bar. If null,
#' no color bar is shown. If "auto", then colors will be
#' added automatically. The default is "auto".
#' @param columnTitle Title to be displayed at top of heatmap.
#'
#' @return ComplexHeatmap object
#' @export
#' @examples
#' plotHeatmapColor(df.input, condition.vec,
#' clusterRow=TRUE, clusterCol=TRUE, displayRowLabels=FALSE,
#' displayColumnLabels=TRUE, displayRowDendrograms=TRUE,
#' displayColumnDendrograms=TRUE, annotationColors = "auto",
#' columnTitle="Title")

plotHeatmapColor <- function(df.input,
                             condition.vec.1,
                             condition.vec.2,
                             condition.1.name,
                             condition.2.name,
                             do.scale = TRUE,
                             clusterRow=TRUE, clusterCol=TRUE, displayRowLabels=TRUE,
                             displayColumnLabels=TRUE, displayRowDendrograms=TRUE,
                             displayColumnDendrograms=TRUE, annotationColors = "auto",
                             columnTitle="Title"){
    #test and fix the constant/zero row
    if (sum(rowSums(as.matrix(df.input)) == 0) > 0){
        df.input <- df.input[-which(rowSums(as.matrix(df.input)) == 0),]
    }


    if (is.null(annotationColors)){
        topha <- NULL
    } else if (annotationColors == "auto") {
        colors <- RColorBrewer::brewer.pal(8, "Set1")
        cond_levels.1 <- unique(condition.vec.1)
        cond_levels.2 <- unique(condition.vec.2)
        if (length(cond_levels.1) > 8 | length(cond_levels.2) > 8){
            stop("Too many levels in condition for auto coloring")
        }
        col <- list()
        col[[paste(condition.1.name)]] <- setNames(colors[1:length(cond_levels.1)], cond_levels.1)
        col[[paste(condition.2.name)]] <- setNames(colors[1:length(cond_levels.2)], cond_levels.2)

        df.plot <- list()
        df.plot[[paste(condition.1.name)]] <- condition.vec.1
        df.plot[[paste(condition.2.name)]] <- condition.vec.2
        topha <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(df.plot),
                                                   col = col)
    }

    # row scale
    if (do.scale == TRUE){
        df.input <- t(scale(t(df.input)))
    }
    heatmap <- ComplexHeatmap::Heatmap(df.input,
                                       name = "Expression",
                                       column_title = columnTitle,
                                       cluster_rows = clusterRow,
                                       cluster_columns = clusterCol,
                                       top_annotation = topha,
                                       show_row_names = displayRowLabels,
                                       show_column_names = displayColumnLabels,
                                       show_row_dend = displayRowDendrograms,
                                       show_column_dend = displayColumnDendrograms)
    return(heatmap)
}


#' Plot PCA
#'
#' @param df.input Input data object that contains the data to be plotted. Required
#' @param condition.color.vec color vector. Required
#' @param condition.color.name color variable name. Required
#' @param condition.shape.vec shape vector. Required
#' @param condition.shape.name shape variable name. Required
#' @param columnTitle Title to be displayed at top of heatmap.
#' @export
#' @examples
#' plotPCAPlotly(df.input, condition.vec, condition.name = "condition",
#' columnTitle = "Title", pc.a = "PC1", pc.b = "PC2")
plotPCAPlotly <- function(df.input,
                          condition.color.vec, condition.color.name = "condition",
                          condition.shape.vec, condition.shape.name = "condition",
                          columnTitle = "Title", pc.a = "PC1", pc.b = "PC2"){

    #test and fix the constant/zero row
    if (sum(rowSums(as.matrix(df.input)) == 0) > 0){
        df.input <- df.input[-which(rowSums(as.matrix(df.input)) == 0),]
    }

    # conduct PCA
    pca.tmp<- prcomp(t(df.input), scale = TRUE)
    tmp.df <- data.frame(pca.tmp$x)
    # add color variable
    tmp.df[[paste(condition.color.name)]] <- condition.color.vec
    # add shape variable
    tmp.df[[paste(condition.shape.name)]] <- condition.shape.vec
    plot_ly(tmp.df,
            x = as.formula(paste("~", pc.a, sep = "")),
            y = as.formula(paste("~", pc.b, sep = "")),
            mode = "markers",
            color = as.formula(paste("~", condition.color.name, sep = "")),
            symbol = as.formula(paste("~", condition.shape.name, sep = "")),
            type = "scatter",
            text = rownames(tmp.df),
            marker = list(size = 10))
}




#' Plot PCoA
#'
#' @param df.input Input data object that contains the data to be plotted. Required
#' @param condition.color.vec color vector. Required
#' @param condition.color.name color variable name. Required
#' @param condition.shape.vec shape vector. Required
#' @param condition.shape.name shape variable name. Required
#' @param columnTitle Title to be displayed at top of heatmap.
#' @export
#' @examples
#' plotPCoAPlotly(physeq.input, condition.vec, condition.name = "condition", method = "bray",
#' columnTitle = "Title", pc.a = "Axis.1", pc.b = "Axis.2")
plotPCoAPlotly <- function(physeq.input,
                           condition.color.vec, condition.color.name = "condition",
                           condition.shape.vec, condition.shape.name = "condition",
                           method = "bray", columnTitle = "Title", pc.a = "Axis.1", pc.b = "Axis.2"){
  # conduct PCoA
  # wUniFrac or bray

  #test and fix the constant/zero row
  if (sum(rowSums(as.matrix(physeq.input@otu_table@.Data)) == 0) > 0){
    physeq.input@otu_table@.Data <- df.input[-which(rowSums(as.matrix(physeq.input@otu_table@.Data)) == 0),]
  }

  if (method == "bray"){
    DistBC = phyloseq::distance(physeq.input, method = method)
    ordBC = ordinate(physeq.input, method = "PCoA", distance = DistBC)
    tmp.df <- data.frame(ordBC$vectors)
  } else {
    DistUF = phyloseq::distance(physeq.input, method = method)
    ordUF = ordinate(physeq.input, method = "PCoA", distance = DistUF)
    tmp.df <- data.frame(ordUF$vectors)
  }

    # add color variable
    tmp.df[[paste(condition.color.name)]] <- condition.color.vec
    # add shape variable
    tmp.df[[paste(condition.shape.name)]] <- condition.shape.vec

  plot_ly(tmp.df,
          x = as.formula(paste("~", pc.a, sep = "")),
          y = as.formula(paste("~", pc.b, sep = "")),
          mode = "markers",
          color = as.formula(paste("~", condition.color.name, sep = "")),
          symbol = as.formula(paste("~", condition.shape.name, sep = "")),
          type = "scatter",
          text = rownames(tmp.df),
          marker = list(size = 10))
}







