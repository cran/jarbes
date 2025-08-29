#' Plot Data Function
#'
#' This function plots the treatment effect (TE) against specified covariates or
#' creates a forest plot of treatment effects, arranging up to four plots in a grid.
#'
#' @param data A data frame with columns `TE`, `seTE`, and the covariates.
#' @param formula A one-sided formula (e.g., `~ x1 + x2`) for the x-axis variables.
#'                If not provided or `~1`, a forest plot is generated.
#' @param study.label A character string specifying the name of the column in `data`
#'                    that contains the study labels for the forest plot. If `NULL`,
#'                    default labels ("Study 1", "Study 2", etc.) are created.
#' @param x.lab A character vector for x-axis labels. Defaults to variable names.
#' @param y.lab Text for the y-axis label. Defaults to "TE" or "Study" for forest plots.
#' @param title.plot Text for the plot title. Defaults to `NULL`.
#' @param errbar.width The width of the error bars for regular plots. Defaults to `0.25`.
#' @param label.size A numeric value for the font size of the study labels in the
#'                   forest plot. Defaults to `6`.
#' @param forest.color The color of the points and error bars in the forest plot. Defaults to "blue".
#' @param forest.shape The shape of the points in the forest plot. Defaults to `19` (filled circle).
#' @param forest.size The size of the points in the forest plot. Defaults to `3`.
#' @param xlim A numeric vector of length 2 for the x-axis limits.
#' @param ylim A numeric vector of length 2 for the y-axis limits.
#' @param ... Additional arguments passed to `geom_point()`.
#'
#' @return A `gridExtra` plot object or a single `ggplot` object.
#'
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @export
plot_data = function(data,
                     formula = ~1,
                     study.label = NULL,
                     x.lab = NULL,
                     y.lab = "TE",
                     title.plot = NULL,
                     errbar.width = 0.25,
                     label.size = 6,
                     forest.color = "blue",
                     forest.shape = 19,
                     forest.size = 3,
                     xlim = NULL,
                     ylim = NULL,
                     ...) {

  seTE = TE = NULL

  # Check for forest plot condition
  is_forest_plot = identical(formula, ~1)

  if (is_forest_plot) {
    # If study.label is not provided, create default labels
    if (is.null(study.label)) {
      data$study_label_default = paste("Study", 1:nrow(data))
      study_label_col = "study_label_default"
    } else {
      if (!(study.label %in% names(data))) {
        stop(paste0("Column '", study.label, "' not found in the data frame."))
      }
      study_label_col = study.label
    }

    # Ensure y-axis labels are ordered correctly from top to bottom
    data[[study_label_col]] = factor(data[[study_label_col]], levels = unique(data[[study_label_col]]))

    p = ggplot(data = data, aes(x = TE, y = .data[[study_label_col]])) +
      geom_point(shape = forest.shape, size = forest.size, color = forest.color, ...) +
      geom_errorbarh(aes(xmin = TE - 2 * seTE, xmax = TE + 2 * seTE), height = 0.1, color = forest.color) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      theme_bw() +
      labs(
        x = if (is.null(x.lab)) "Estimate" else x.lab,
        y = if (is.null(y.lab) || y.lab == "TE") "Study" else y.lab,
        title = title.plot
      ) +
      theme(
        axis.text.y = element_text(size = label.size)
      )

    if (!is.null(xlim)) {
      p = p + coord_cartesian(xlim = xlim)
    }

    return(suppressWarnings(p))
  }

  # Plotting data for covariates
  x_vars = all.vars(formula)
  if (length(x_vars) == 0) {
    stop("Formula must specify at least one variable.")
  }

  if (length(x_vars) > 4) {
    warning("More than 4 variables specified. Only the first 4 will be plotted.")
    x_vars = x_vars[1:4]
  }

  plot_list = list()

  # Build plot for each variable
  for (i in 1:length(x_vars)) {
    x_var = x_vars[i]

    current_x_lab = if (!is.null(x.lab) && length(x.lab) >= i) {
      x.lab[i]
    } else {
      x_var
    }

    p = ggplot(data = data, aes(x = .data[[x_var]], y = TE)) +
      theme_bw() +
      labs(
        x = current_x_lab,
        y = y.lab,
        title = if (length(x_vars) > 1) {
          paste0(title.plot, " (", x_var, ")")
        } else {
          title.plot
        }
      ) +
      geom_errorbar(
        aes(ymin = TE - 2 * seTE, ymax = TE + 2 * seTE),
        width = errbar.width
      ) +
      geom_point(
        aes(size = 1/seTE),
        fill = "royalblue",
        shape = 21,
        ...
      ) +
      geom_hline(yintercept = 0, colour = "black", lty = 2)

    # Apply axis limits if provided
    if (!is.null(xlim)) p = p + xlim(xlim)
    if (!is.null(ylim)) p = p + ylim(ylim)

    plot_list[[x_var]] = p
  }

  # Return plot or grid
  if (length(x_vars) == 1 && length(plot_list) == 1) {
    return(suppressWarnings(plot_list[[1]]))
  }

  if (length(plot_list) > 0) {
    return(suppressWarnings(gridExtra::grid.arrange(grobs = plot_list)))
  } else {
    invisible(NULL)
  }
}
