#' @title Co-Clustering Heatmap from a Numeric Matrix
#'
#' @description
#' Generates a heatmap with hierarchical clustering applied to rows and columns.
#' It reorders the input matrix based on clustering results and visualizes the
#' data.
#'
#' @param data A numeric matrix for the heatmap.
#' @param labels.row An optional character vector of new labels for the rows.
#' @param labels.col An optional character vector of new labels for the columns.
#' @param colors A character string specifying one of the predefined palettes or a custom vector of colors. Predefined options include:
#' \itemize{
#'  \item "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds", "RdPu", "Purples", "PuBu", "PuBuGn", "OrRd", "Oranges", "Greys", "Greens", "GnBu", "BuPu", "BuGn", "Blues"
#' }
#'    Defaults to "Blues".
#' @param breaks The number of color breaks to use for the palette. Defaults to `NULL`.
#' @param title Plot title. Default: "Clustered Heatmap".
#' @param border.color Color of cell borders on the heatmap. Use "NA" for no border. Default: "white".
#' @param cell.ratio A numeric value defining the cell's height-to-width ratio. A value of 1 creates square cells. A value > 1 creates vertical rectangles, while a value < 1 creates horizontal rectangles. Default is 1.
#' @param cluster.rows A logical value indicating whether to cluster rows. Default: `TRUE`.
#' @param cluster.cols A logical value indicating whether to cluster columns. Default: `TRUE`.
#' @param clustering.distance.rows Distance measure for row clustering. Can be a string ("euclidean", "correlation", etc.) or a pre-calculated distance object. Default: "euclidean".
#' @param clustering.distance.cols Distance measure for column clustering. Can be a string ("euclidean", "correlation", etc.) or a pre-calculated distance object. Default: "euclidean".
#' @param row.clustering.method Agglomeration method for row clustering (see `?hclust`). Default: "complete".
#' @param col.clustering.method Agglomeration method for column clustering (see `?hclust`). Default: "complete".
#' @param show.row.names Logical, display row names. Default: TRUE.
#' @param show.col.names Logical, display column names. Default: TRUE.
#' @param fontsize Base font size for the plot. Default: 8.
#' @param fontsize.row Font size for row names. If `NULL`, defaults to `fontsize`.
#' @param fontsize.col Font size for column names. If `NULL`, defaults to `fontsize`.
#' @param display.numbers Logical determining if the numeric values are also printed to the cells.
#' @param round.digits Numeric value. For `display.numbers = TRUE`, this sets the decimal places for all values. Default: 2.
#' @param numbers.color Color of the numbers. Default: "black".
#' @param numbers.fontsize Font size of the numbers. Default: `fontsize`.
#' @param label.color Color of the row and column labels. Default: "black".
#' @param col.name.angle Angle for column names. Default: 45.
#' @param legend A logical value indicating whether to show the color legend. Default: `TRUE`.
#' @param legend.labels A character vector of custom labels for the legend.
#' @param show.legend.labels Logical, display legend labels. Default: TRUE.
#'
#' @param ... \dots
#'
#'
#' @return A `ggplot2` object.
#' @import ggplot2
#' @importFrom stats hclust dist as.dist cor
#' @importFrom grDevices colorRampPalette
#' @export
cocluster_plot = function(data,
                          labels.row = NULL,
                          labels.col = NULL,
                          colors = "Blues",
                          breaks = NULL,
                          title = "Clustered Heatmap",
                          border.color = "black",
                          cell.ratio = 1,
                          cluster.rows = TRUE,
                          cluster.cols = TRUE,
                          clustering.distance.rows = "euclidean",
                          clustering.distance.cols = "euclidean",
                          row.clustering.method = "complete",
                          col.clustering.method = "complete",
                          show.row.names = TRUE,
                          show.col.names = TRUE,
                          fontsize = 8,
                          fontsize.row = NULL,
                          fontsize.col = NULL,
                          display.numbers = FALSE,
                          round.digits = 2,
                          numbers.color = "black",
                          numbers.fontsize = 0.5,
                          label.color = "black",
                          col.name.angle = 45,
                          legend = TRUE,
                          legend.labels = NULL,
                          show.legend.labels = TRUE,
                          ...) {

  Col = Row = Value = NULL

  ## Input Validation and Parameter Setup
  if (is.null(rownames(data))) {
    rownames(data) = 1:nrow(data)
  }
  if (is.null(colnames(data))) {
    colnames(data) = 1:ncol(data)
  }
  if (!is.null(labels.row)) {
    rownames(data) = labels.row
  }
  if (!is.null(labels.col)) {
    colnames(data) = labels.col
  }

  row.name.size = if (is.null(fontsize.row)) fontsize else fontsize.row
  col.name.size = if (is.null(fontsize.col)) fontsize else fontsize.col
  final.numbers.fontsize = if (is.null(numbers.fontsize)) fontsize else numbers.fontsize

  ## Hierarchical Clustering and Matrix Reordering
  # Row Clustering
  dist.rows = NULL
  if (cluster.rows) {
    if (is.character(clustering.distance.rows)) {
      if (clustering.distance.rows == "correlation") {
        dist.rows = as.dist(1 - cor(t(data), use = "pairwise.complete.obs"))
      } else {
        dist.rows = dist(data, method = clustering.distance.rows)
      }
    } else if (inherits(clustering.distance.rows, "dist")) {
      dist.rows = clustering.distance.rows
    } else {
      stop("Invalid value for 'clustering.distance.rows'. Must be a string or a dist object.")
    }
  }


  # Column Clustering
  dist.cols = NULL
  if (cluster.cols) {
    if (is.character(clustering.distance.cols)) {
      if (clustering.distance.cols == "correlation") {
        dist.cols = as.dist(1 - cor(data, use = "pairwise.complete.obs"))
      } else {
        dist.cols = dist(t(data), method = clustering.distance.cols)
      }
    } else if (inherits(clustering.distance.cols, "dist")) {
      dist.cols = clustering.distance.cols
    } else {
      stop("Invalid value for 'clustering.distance.cols'. Must be a string or a dist object.")
    }
  }

  row.order = 1:nrow(data)
  col.order = 1:ncol(data)
  if (cluster.rows && nrow(data) > 1) {
    row.order = hclust(dist.rows, method = row.clustering.method)$order
  }
  if (cluster.cols && ncol(data) > 1) {
    col.order = hclust(dist.cols, method = col.clustering.method)$order
  }

  reordered.matrix = data[row.order, col.order]

  ## Prepare Data for Plotting with unique names
  plot.row.names = make.unique(rownames(reordered.matrix), sep = " ")
  plot.col.names = make.unique(colnames(reordered.matrix), sep = " ")

  # Create a data frame for ggplot
  heatmap.data = expand.grid(Row = factor(plot.row.names, levels = rev(plot.row.names)),
                             Col = factor(plot.col.names, levels = plot.col.names))
  heatmap.data$Value = as.vector(t(reordered.matrix))

  ## Create the Plot
  p = ggplot(heatmap.data, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = if (is.na(border.color)) NA else border.color, linewidth = 0.5) +
    labs(title = title, x = "", y = "", fill = "Value") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = fontsize),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = if (legend) "right" else "none",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
    ) +
    scale_y_discrete(position = "right")

  # Add theme elements for axis labels
  p = p + theme(
    axis.text.x = if (show.col.names) element_text(angle = col.name.angle, hjust = 1, size = col.name.size, color = label.color) else element_blank(),
    axis.text.y = if (show.row.names) element_text(size = row.name.size, color = label.color) else element_blank()
  )

  ## Apply Color Scale and Legend Options
  # Inbuilt Palettes
  inbuilt_palettes = list(
    # Sequential
    "YlOrRd" = c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026"),
    "YlOrBr" = c("#FFFFD4", "#FED98E", "#FE9929", "#EC7014", "#CC4C02", "#993404", "#662506"),
    "YlGnBu" = c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58"),
    "YlGn" = c("#FFFFE5", "#F7FCB9", "#D9F0A3", "#ADDD8E", "#78C679", "#41AB5D", "#238443", "#006837", "#004529"),
    "Reds" = c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FA6542", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"),
    "RdPu" = c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E", "#7A0177", "#49006A"),
    "Purples" = c("#FBF7FF", "#E9E3F8", "#D2C2EC", "#B8A2D3", "#9B74BE", "#7D4E9F", "#5F2783", "#431464", "#2B0B3C"),
    "PuBu" = c("#F1EEF6", "#D7B5D8", "#DF65B0", "#DD1C77", "#980043"),
    "PuBuGn" = c("#FFF7FB", "#ECE7F2", "#D0D1E6", "#A6BDDB", "#67A9CF", "#3690C0", "#02818A", "#016C59", "#014636"),
    "OrRd" = c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"),
    "Oranges" = c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704"),
    "Greys" = c("#FFFFFF", "#F0F0F0", "#D9D9D9", "#BDBDBD", "#969696", "#737373", "#525252", "#252525", "#000000"),
    "Greens" = c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"),
    "GnBu" = c("#F7FCF0", "#E0F3DB", "#CCEBC5", "#A8DDB5", "#7BCCC4", "#4EB3D3", "#2B8CBE", "#0868AC", "#084081"),
    "BuPu" = c("#F7F4F9", "#E7E1EF", "#D4B9DA", "#C994C7", "#DF65B0", "#E7298A", "#CE1256", "#990042"),
    "BuGn" = c("#F7FCF0", "#E0F5D5", "#CCEBC5", "#A8DDB5", "#7BCCC4", "#4EB3D3", "#2B8CBE", "#0868AC", "#084081"),
    "Blues" = c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")
  )

  # Determine the final color palette
  if (is.character(colors) && length(colors) == 1 && colors %in% names(inbuilt_palettes)) {
    final_colors_raw = inbuilt_palettes[[colors]]
  } else if (is.character(colors) && length(colors) > 1) {
    final_colors_raw = colors
  } else {
    stop("Invalid 'colors' input. Must be a character string for a predefined palette or a vector of colors.")
  }

  color_interpolator = colorRampPalette(final_colors_raw)
  if (is.null(breaks)) {
    # If breaks are not specified, use the number of colors in the palette
    final_breaks_count = length(final_colors_raw)
  } else if (is.numeric(breaks)) {
    # If a numeric break count is provided, use that
    final_breaks_count = breaks
  } else {
    stop("Invalid 'breaks' input. Must be a numeric value or NULL.")
  }

  color_palette = color_interpolator(final_breaks_count)

  # Determine the data range for the breaks
  value_range = range(heatmap.data$Value, na.rm = TRUE)

  # Create a sequence of breaks that are equally spaced
  breaks_seq = seq(from = value_range[1], to = value_range[2], length.out = length(color_palette) + 1)


  p = p + scale_fill_stepsn(colors = color_palette,
                            breaks = breaks_seq,
                            labels = function(x) round(x, 2))

  if (!show.legend.labels) {
    p = p + guides(fill = guide_legend(label = FALSE))
  } else if (!is.null(legend.labels)) {
    p = p + guides(fill = guide_legend(override.aes = list(fill = color_palette), labels = legend.labels))
  }

  ## Add Numbers to Cells
  if (display.numbers) {
    p = p + geom_text(
      mapping = aes(x = Col, y = Row, label = round(Value, round.digits)),
      color = numbers.color,
      size = final.numbers.fontsize * .pt,
      inherit.aes = FALSE,
      data = heatmap.data
    )
  }

  ## Set Cell Aspect Ratio
  if (!is.null(cell.ratio)) {
    p = p + coord_fixed(ratio = cell.ratio, expand = FALSE)
  } else {
    p = p + coord_fixed(expand = FALSE)
  }

  return(suppressWarnings(p))
}
