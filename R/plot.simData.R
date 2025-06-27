#' Forest Plot for simData Objects
#'
#' Generates a forest plot of the simulated meta-analysis dataset,
#' showing the observed effect sizes (TE) with 95% confidence intervals.
#'
#' @name plot.simData
#' @method plot simData
#' @param x A `simData` object created by `simData()`.
#' @param x.lim Numeric vector of length 2 specifying the x-axis limits.
#' @param x.lab Text with the label of the x-axis.
#' @param y.lab Text with the label of the y-axis.
#' @param title.plot Text for the title of the plot.
#' @param bias_colors Named character vector specifying colors for bias categories.
#' @param bias_legend_title Text label for the bias legend.
#' @param ref_line_color Color for the vertical reference line corresponding to the mean.
#' @param ... Additional arguments (currently not used).
#'
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom grid grid.newpage grid.draw unit

# Declare global variables
utils::globalVariables(c("B.flag", "CI_lower", "CI_upper", "TE"))
#' @export

plot.simData <- function(x,
                         x.lim = NULL,
                         x.lab = "Observed Treatment Effect (TE)",
                         y.lab = "Study",
                         title.plot = "Forest Plot of Simulated Meta-Analysis",
                         bias_colors = c(
                           "No B" = "#1f77b4",       # Blue for unbiased
                           "Mild B" = "#fcae91",     # Light red
                           "Large B" = "#fb6a4a",    # Medium red
                           "Extreme B" = "#cb181d"   # Dark red
                         ),
                         bias_legend_title = "Bias",
                         ref_line_color = "darkgreen",
                         ...) {


  # Extract overall treatment effect
  mu <- attr(x, "mu")

  # Calculate 95% confidence intervals
  x$CI_lower <- x$TE - 1.96 * x$seTE
  x$CI_upper <- x$TE + 1.96 * x$seTE

  if (is.null(x.lim)) {
    x.lim <- range(c(x$CI_lower, x$CI_upper, mu), na.rm = TRUE)
  } else {
    x.lim <- range(c(x.lim, mu), na.rm = TRUE)
  }

  # Generate forest plot
  plot_forest <- ggplot(x, aes(
    x = TE,
    y = factor(rownames(x), levels = rev(rownames(x))),
    color = B.flag
  )) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.3, size = 1) +
    geom_vline(xintercept = mu, color = ref_line_color, size = 1) +
    scale_color_manual(values = bias_colors, name = bias_legend_title) +
    scale_x_continuous(
      limits = x.lim,
      breaks = sort(unique(c(pretty(x.lim), mu))),
      labels = function(b) {
        ifelse(abs(b - mu) < 1e-8, round(mu, 3), b)
      }
    ) +
    labs(
      title = title.plot,
      x = x.lab,
      y = y.lab
    ) +
    theme_bw(base_size = 8) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 6),
      plot.margin = unit(c(1, 1, 2, 1), "lines")
    )

  # Dynamically adjust height based on number of studies
  n_studies <- nrow(x)
  plot_height <- max(4, n_studies * 0.15)

  # Display the plot
  grid.newpage()
  grid.draw(grid.arrange(plot_forest, heights = unit(plot_height, "inches")))
}
