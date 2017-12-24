#' Plot heatmap with color bar
#'
#' @param df.input Input data object that contains the data to be plotted. Required
#' @param condition.vec The condition used for plotting the heatmap. Required
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
plotHeatmapColor <- function(df.input, condition.vec,
                             clusterRow=TRUE, clusterCol=TRUE, displayRowLabels=FALSE,
                             displayColumnLabels=TRUE, displayRowDendrograms=TRUE,
                             displayColumnDendrograms=TRUE, annotationColors = "auto",
                             columnTitle="Title"){
  if (is.null(annotationColors)){
    topha <- NULL
  } else if (annotationColors == "auto") {
    colors <- RColorBrewer::brewer.pal(8, "Set1")
    cond_levels <- unique(condition.vec)
    if (length(cond_levels) > 8){
      stop("Too many levels in condition for auto coloring")
    }
    col <- list(type = setNames(colors[1:length(cond_levels)], cond_levels))
    topha <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(type = condition.vec),
                                               col = col)
  } else {
    topha <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(type = condition.vec),
                                               col = annotationColors)
  }
  
  heatmap <- ComplexHeatmap::Heatmap(t(scale(t(df.input))),
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



