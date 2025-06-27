#' Compare Posterior Estimates from Bayesian Models
#'
#' Generates a caterpillar-style plot (forest plot) for visualizing and comparing
#' posterior parameter estimates from one or two Bayesian models.
#' This function is designed for use within the `jarbes` package and
#' supports models fitted via MCMC. It allows custom labels, credible intervals,
#' and styling for visual model comparison, particularly in meta-analytic and
#' hierarchical modeling contexts.
#'
#'
#' @param model1 An object containing MCMC draws. Various formats (e.g., arrays, matrices, data frames, `posterior::draws` objects) are accepted.
#' @param model2 Optional object containing MCMC draws. Accepted formats are the same as for `model1`.
#' @param pars Character vector of parameter names to include in the plot.
#' @param plotmath.labels Optional character vector for y-axis labels. If provided in R's plotmath syntax (e.g., for Greek letters or mathematical symbols), these labels will be displayed on the plot.
#' @param model1.name Text for the label of the first model.
#' @param model2.name Text for the label of the second model.
#' @param model.legend.title Text for the title of the model legend.
#' @param ref.lines Numeric value indicating vertical reference lines.
#' @param colors Character vector specifying the colors for models.
#' @param point.size Numeric value for the size of points in the plot.
#' @param point.shapes Numeric or character vector specifying the shapes for points, one for each model.
#' @param prob Numeric value for the probability mass to include in the inner interval.
#' @param prob.outer Numeric value for the probability mass to include in the outer interval.
#' @param point.est Text specifying the type of point estimate to show. Either `"median"` (the default), `"mean"`, or `"none"`.
#' @param x.lab Text with the label of the x-axis.
#' @param y.lab Text with the label of the y-axis.
#' @param inner.line.thickness Numeric value for the thickness of the inner interval line.
#' @param outer.line.thickness Numeric value for the thickness of the outer interval line.
#' @param ... \dots
#'
#' @import ggplot2
#' @importFrom bayesplot mcmc_intervals_data
#' @export
caterplot_compare = function(model1,
                             model2 = NULL,
                             pars,
                             plotmath.labels = NULL,
                             model1.name = "Model 1",
                             model2.name = "Model 2",
                             model.legend.title = "Model",
                             ref.lines = c(0),
                             colors = c("blue", "red"),
                             point.size = 3,
                             point.shapes = c(16, 17),
                             prob = 0.5,
                             prob.outer = 0.9,
                             point.est = "median",
                             x.lab = "Estimate",
                             y.lab= NULL,
                             inner.line.thickness = 2,
                             outer.line.thickness = 0.8, ...) {

  # Declare global variables
  m=y.pos=ll=hh=model=l=h=NULL

  # --- Data extraction and preparation ---
  df1 = mcmc_intervals_data(model1, pars = pars,
                            prob = prob, prob_outer = prob.outer, point_est = point.est)

  if (!is.null(model2)) {
    # --- Case: Comparing two models ---
    current.model.names = c(model1.name, model2.name)
    df1$model = current.model.names[1]
    df2 = mcmc_intervals_data(model2, pars = pars,
                              prob = prob, prob_outer = prob.outer, point_est = point.est)
    df2$model = current.model.names[2]
    df = rbind(df1, df2)

    # Assigns numeric y-axis positions to parameters, ensuring they are ordered according to the 'pars' argument.
    param.levels = intersect(pars, unique(df$parameter))
    df$y.pos = as.numeric(factor(df$parameter, levels = param.levels))
    # Apply vertical offset to separate estimates for each model on the y-axis, preventing overlap.
    df$y.pos = df$y.pos + ifelse(df$model == current.model.names[1], -0.15, 0.15)

    plot.colors = colors
    model.legend.labels = current.model.names
    plot.point.shapes = point.shapes # Use the new shapes argument

  } else {
    # --- Case: Plotting a single model ---
    current.model.names = model1.name
    df1$model = current.model.names
    df = df1

    # Assigns numeric y-axis positions to parameters, ensuring they are ordered according to the 'pars' argument.
    param.levels = intersect(pars, unique(df$parameter))
    df$y.pos = as.numeric(factor(df$parameter, levels = param.levels))

    plot.colors = colors[1]
    model.legend.labels = current.model.names
    plot.point.shapes = point.shapes[1] # Use the first shape for single model
  }

  # --- Prepare y-axis labels ---
  if (!is.null(plotmath.labels)) {
    y.labels = plotmath.labels[match(param.levels, pars)]
    # Convert character labels to R expressions for plotmath rendering.
    if(is.character(y.labels)) {
      y.labels = parse(text = y.labels)
    }

  } else {
    y.labels = param.levels
  }

  # --- Plot construction ---
  p = ggplot(df, aes(x = m, y = y.pos, color = model)) +
    geom_point(aes(shape = model), size = point.size)

  # Add outer credible interval bars.
  if (outer.line.thickness > 0) {
    p = p + geom_errorbarh(aes(xmin = ll, xmax = hh), height = 0, size = outer.line.thickness)
  }

  # Add inner credible interval bars.
  if (inner.line.thickness > 0) {
    p = p + geom_errorbarh(aes(xmin = l, xmax = h), height = 0, size = inner.line.thickness)
  }

  p = p +
    scale_y_continuous(
      breaks = 1:length(param.levels),
      labels = y.labels,
      trans = "reverse"
    ) +
    geom_vline(xintercept = ref.lines, linetype = "dashed", color = "grey40", linewidth = 0.8) +
    scale_color_manual(values = plot.colors, name = model.legend.title, labels = model.legend.labels) +
    # Use scale_shape_manual to control point shapes
    scale_shape_manual(values = plot.point.shapes, name = model.legend.title, labels = model.legend.labels) +
    labs(x = x.lab, y = y.lab) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8)) +
    geom_hline(yintercept = length(param.levels) + 0.5, color = "black", linetype = "solid", linewidth = 0.5)


  return(suppressWarnings(p))
}
