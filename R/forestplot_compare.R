#' @title Forest Plot Comparison for One or Two Models
#'
#' @description
#' Creates a forest plot comparing parameter estimates from either:
#' - a single model with one parameter set.
#' - a single model with two different parameter sets.
#' - two different models with matching parameters.
#'
#' @param model1 First model object. Can be an `mcmc.list` object, a model with BUGS output, a matrix, or a data frame of MCMC draws.
#' @param model2 Optional second model object. If omitted, the function can compare two parameter sets within model1.
#'   Accepted types are the same as `model1`.
#' @param pars Character vector of parameter names for the first model or first parameter set
#' @param pars2 Optional character vector of parameter names for the second parameter set in model1. Ignored if model2 is provided.
#' @param plotmath.labels Optional vector of plotmath expressions as strings to label parameters in the first model or set.
#' @param plotmath.labels2 Optional vector of plotmath expressions for parameters in the second parameter set in model1.
#' @param study.labels Optional character vector of labels to use for studies on the y-axis. If provided,
#'   these labels will replace parameter names, and each study label will correspond to a pair of estimates.
#' @param consolidate.param.labels Logical. If `TRUE` and `study.labels` is `NULL` for comparison plots,
#'   displays a single consolidated label for the parameter pair.
#'   If `FALSE` (default), displays a label for each individual parameter.
#'   Ignored for single model plots.
#' @param model1.name Label for the first model in the plot legend.
#' @param model2.name Label for the second model or parameter set in the legend.
#' @param model.legend.title Title for the model legend.
#' @param ref.lines Numeric vector of vertical reference lines to display.
#' @param colors Vector of colors for each model or parameter set.
#' @param point.size Size of point estimate markers.
#' @param point.shapes Shapes used for points corresponding to each model or set
#' @param prob Width of the inner credible interval (e.g. 0.5 for 50 \% interval).
#' @param prob.outer Width of the outer credible interval (e.g. 0.9 for 90\% interval).
#' @param point.est Which point estimate to use: "median" or "mean", default is "median".
#' @param x.lab Label for the x-axis.
#' @param y.lab Label for the y-axis.
#' @param inner.line.thickness Thickness of inner credible interval lines.
#' @param outer.line.thickness Thickness of outer credible interval lines.
#' @param x.lim Optional numeric vector of length two to fix x-axis limits.
#' @param param.distance Numeric, controls the vertical distance between points for different models/parameter sets
#'   within the same study. Default is 0.3.
#' @param group.spacing Numeric, controls the vertical spacing between different studies/group of parameters. Default is 1.
#' @param flip.coords Logical, if TRUE, flips the x and y axes, creating a horizontal plot. Default is FALSE.
#' @param x.label.angle Numeric. Angle for x-axis labels when `flip.coords = TRUE`. Default is 0.
#' @param study.label.size Numeric. Text size for study or parameter labels on the y-axis (or x-axis if flipped).
#' @param ... \dots
#'
#' @import ggplot2
#' @export
forestplot_compare = function(
    model1,
    model2 = NULL,
    pars,
    pars2 = NULL,
    plotmath.labels = NULL,
    plotmath.labels2 = NULL,
    study.labels = NULL,
    consolidate.param.labels = FALSE,
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
    y.lab = NULL,
    inner.line.thickness = 1.5,
    outer.line.thickness = .5,
    x.lim = NULL,
    param.distance = 0.3,
    group.spacing = 1,
    flip.coords = FALSE,
    x.label.angle = 0,
    study.label.size = 11,
    ...
) {

  # This function converts the MCMC output to a matrix ........................
  convert_model_to_matrix = function(model_obj)
  { # Handle mcmc lists
    if (inherits(model_obj, "mcmc.list")) {
      return(do.call(rbind, lapply(model_obj, as.matrix)))
    }
    # Handle Matrix or dataframe input
    else if (is.matrix(model_obj) || is.data.frame(model_obj)) {
      return(as.matrix(model_obj))
    }
    # Handle model inputs which have $BUGSoutput$sims.matrix
    else if (is.list(model_obj) && !is.null(model_obj$BUGSoutput) && !is.null(model_obj$BUGSoutput$sims.matrix)) {
      return(model_obj$BUGSoutput$sims.matrix)
    }
    else {
      stop("Unsupported model input type. Must be an 'mcmc.list' object, a matrix, a data.frame, or an object containing 'BUGSoutput$sims.matrix' (e.g., from rjags/jags).", call. = FALSE)
    }
  }
  # ............................................................................

  # This function calculates the posterior intervals ..........................
  calculate_intervals_df = function(draws_matrix, params_to_extract, prob, prob.outer, point.est)
  {

    # Compute quantile probabilities
    probs = c(
      0.5 - prob.outer / 2,
      0.5 - prob / 2,
      0.5 + prob / 2,
      0.5 + prob.outer / 2
    )

    available_params = intersect(params_to_extract, colnames(draws_matrix))

    results = lapply(available_params, function(param) {
      vals = draws_matrix[, param, drop=TRUE]

      p_est = if (point.est == "median") {
        median(vals, na.rm=TRUE)
      } else if (point.est == "mean") {
        mean(vals, na.rm=TRUE)
      } else {
        NA
      }

      qtiles = quantile(vals, probs, na.rm=TRUE)

      res = data.frame(
        parameter = param,
        m = p_est,
        l = qtiles[2],
        h = qtiles[3],
        ll = qtiles[1],
        hh = qtiles[4]
      )

    })

    do.call(rbind, results)
  }
  # ............................................................................

  # Model 1 and if Model 2 are converted to matrix ............................
  model1_converted = convert_model_to_matrix(model1)
  if (!is.null(model2)) model2_converted = convert_model_to_matrix(model2)

  # One model two parameters sets..............................................
  if (!is.null(pars2) && is.null(model2)) {
    if (length(pars) != length(pars2)) stop("'pars' and 'pars2' must have the same length.")
    df1 = calculate_intervals_df(model1_converted, pars, prob, prob.outer, point.est)
    df2 = calculate_intervals_df(model1_converted, pars2, prob, prob.outer, point.est)

    df_list = list()
    ordered_pars_with_data <- c()

    for (i in seq_along(pars)) {
      p1 = pars[i]
      p2 = pars2[i]
      row1 = df1[df1$parameter==p1, ]
      row2 = df2[df2$parameter==p2, ]

      if (nrow(row1)>0 && nrow(row2)>0) {
        row1$model = model1.name
        row2$model = model2.name

        base_y_pos = i * group.spacing
        row1$y.pos = base_y_pos - param.distance / 2
        row2$y.pos = base_y_pos + param.distance / 2

        row1$label_text = if (!is.null(plotmath.labels)) plotmath.labels[match(p1, pars)] else p1
        row2$label_text = if (!is.null(plotmath.labels2)) plotmath.labels2[match(p2, pars2)] else p2

        df_list = c(df_list, list(row1, row2))
        ordered_pars_with_data <- c(ordered_pars_with_data, pars[i])
      }
    }

    df = do.call(rbind, df_list)
    df$model = factor(df$model, levels = c(model1.name, model2.name))

    if (!is.null(study.labels)) {
      y.breaks_values = (seq_along(ordered_pars_with_data)) * group.spacing
      y.breaks = sort(y.breaks_values)
      original_indices_for_labels = match(y.breaks / group.spacing, seq_along(pars))
      y.labels = study.labels[original_indices_for_labels]
    } else if (!consolidate.param.labels) {
      df_labels_temp = df[, c("y.pos", "label_text")]
      unique_rows_idx = !duplicated(df_labels_temp)
      df_labels_unique = df_labels_temp[unique_rows_idx, ]
      df_labels = df_labels_unique[order(df_labels_unique$y.pos), ]
      y.breaks = df_labels$y.pos
      y.labels = parse(text = df_labels$label_text)
    } else {
      y.breaks = (seq_along(ordered_pars_with_data)) * group.spacing
      y.breaks = sort(y.breaks)
      original_indices_for_labels = match(ordered_pars_with_data, pars)
      y.labels = if (!is.null(plotmath.labels)) parse(text=plotmath.labels[original_indices_for_labels]) else parse(text=ordered_pars_with_data)
    }

    # Two models one set of parameters...........................................
  } else if (!is.null(model2)) {
    df1 = calculate_intervals_df(model1_converted, pars, prob, prob.outer, point.est)
    df2 = calculate_intervals_df(model2_converted, pars, prob, prob.outer, point.est)
    common = intersect(df1$parameter, df2$parameter)

    df1 = df1[df1$parameter %in% common, ]
    df2 = df2[df2$parameter %in% common, ]

    df_list = list()
    ordered_common_with_data <- c()

    for (i in seq_along(common)) {
      p = common[i]
      row1 = df1[df1$parameter == p, ]
      row2 = df2[df2$parameter == p, ]

      if (nrow(row1) > 0 && nrow(row2) > 0) {
        row1$model = model1.name
        row2$model = model2.name

        base_y_pos = i * group.spacing
        row1$y.pos = base_y_pos - param.distance / 2
        row2$y.pos = base_y_pos + param.distance / 2

        original_idx_p = match(p, pars)
        row1$label_text = if (!is.null(plotmath.labels)) plotmath.labels[original_idx_p] else p
        row2$label_text = if (!is.null(plotmath.labels)) plotmath.labels[original_idx_p] else p

        df_list = c(df_list, list(row1, row2))
        ordered_common_with_data <- c(ordered_common_with_data, p)
      }
    }

    df = do.call(rbind, df_list)
    df$model = factor(df$model, levels = c(model1.name, model2.name))

    if (!is.null(study.labels)) {
      y.breaks = sort((seq_along(ordered_common_with_data)) * group.spacing)
      original_indices_for_labels = match(y.breaks / group.spacing, seq_along(pars))
      y.labels = study.labels[original_indices_for_labels]
    } else if (!consolidate.param.labels) {
      df_labels_temp = df[, c("y.pos", "label_text")]
      df_labels_unique = df_labels_temp[!duplicated(df_labels_temp), ]
      df_labels = df_labels_unique[order(df_labels_unique$y.pos), ]
      y.breaks = df_labels$y.pos
      y.labels = parse(text = df_labels$label_text)
    } else {
      y.breaks = sort((seq_along(ordered_common_with_data)) * group.spacing)
      original_indices_for_labels = match(ordered_common_with_data, pars)
      y.labels = if (!is.null(plotmath.labels)) parse(text = plotmath.labels[original_indices_for_labels]) else parse(text = ordered_common_with_data)
    }

    # One model one set of parameters ............................................
  } else {
    df = calculate_intervals_df(model1_converted, pars, prob, prob.outer, point.est)
    df$model = model1.name
    df$model = factor(df$model, levels = c(model1.name))
    df$y.pos = seq_len(nrow(df)) * group.spacing

    y.labels = if (!is.null(study.labels)) study.labels else if (!is.null(plotmath.labels)) parse(text=plotmath.labels) else df$parameter
    y.breaks = seq_len(nrow(df)) * group.spacing
  }

  # Base plot ..................................................................
  p = ggplot(df, aes(x=m, y=y.pos, color=model)) +
    geom_point(aes(shape=model), size=point.size)+
    geom_errorbarh(aes(xmin=ll, xmax=hh), height=0, linewidth=outer.line.thickness)+
    geom_errorbarh(aes(xmin=l, xmax=h), height=0, linewidth=inner.line.thickness)

  if (!is.null(x.lim)) {
    p = p + coord_cartesian(xlim = x.lim)
  }

  p = p +
    geom_vline(xintercept=ref.lines, linetype="dashed", color="grey40") +
    scale_color_manual(values=colors, name=model.legend.title) +
    scale_shape_manual(values=point.shapes, name=model.legend.title) +
    labs(x=x.lab, y=y.lab) +
    theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.text.y = element_text(size = if (!flip.coords) study.label.size else NULL)
    )

  if (flip.coords) {
    p = p +
      scale_y_continuous(breaks=y.breaks, labels=y.labels) +
      coord_flip() +
      theme(axis.text.x = element_text(angle = x.label.angle, size = study.label.size))

    if (is.null(y.lab)) p = p + labs(y = x.lab)
    if (x.lab == "Estimate" && is.null(y.lab)) p = p + labs(x = NULL)
  } else {
    p = p +
      scale_y_continuous(breaks=y.breaks, labels=y.labels, trans="reverse")
  }

  return(suppressWarnings(p))
}

# Prevent R CMD check notes
m = y.pos = ll = hh = model = l = h = NULL
