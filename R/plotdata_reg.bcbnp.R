#' Plot Meta-Regression Results from a BCBNP Model
#'
#' This function visualizes the results of a bias-corrected Bayesian nonparametric
#' meta-regression. It plots the raw data with error bars, colors points based on
#' the posterior probability of bias (I > 0.5), and adds the posterior mean
#' regression line, the credible interval for the true effects, and the
#' prediction interval for new observations, including its median line.
#'
#' @param object The output object from the bcbnp function.
#' @param covariate A character string specifying the single covariate
#' from the data to be plotted (e.g., `"baseline_va"`).
#' @param title.plot The title for the plot.
#' @param ci.color The color of the credible interval lines. Defaults to "green".
#' @param pi.color The color of the prediction interval lines. Defaults to "purple".
#' @param mean.line.color The color of the mean regression line. Defaults to "magenta".
#' @param pi.mean.line.color The color of the prediction interval median line. Defaults to "gray".
#' @param point.color.unbiased The color of unbiased points. Defaults to "blue".
#' @param point.color.biased The color of biased points. Defaults to "red".
#' @param show_ci A logical value indicating whether to show the credible interval. Defaults to `TRUE`.
#' @param show_pi A logical value indicating whether to show the prediction interval. Defaults to `TRUE`.
#' @param errorbar_width The width of the error bars. Defaults to `0.25`.
#' @param x.lab A character string for the x-axis label. Defaults to the covariate name.
#' @param y.lab A character string for the y-axis label. Defaults to "Treatment Effect (TE)".
#' @param legend.name A character string for the color legend title. Defaults to "Posterior Bias".
#' @param size.legend.name A character string for the size legend title. Defaults to "1/seTE".
#' @param label.unbiased A character string for the label of the unbiased group. Defaults to "Unbiased".
#' @param label.biased A character string for the label of the biased group. Defaults to "Biased".
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
plotdata_reg.bcbnp = function(object, covariate, title.plot = NULL,
                              ci.color = "green",
                              pi.color = "purple",
                              mean.line.color = "magenta",
                              pi.mean.line.color = "gray",
                              point.color.unbiased = "blue",
                              point.color.biased = "red",
                              show_ci = TRUE,
                              show_pi = TRUE,
                              errorbar_width = 0.25,
                              x.lab = covariate,
                              y.lab = "Treatment Effect (TE)",
                              legend.name = "Posterior Bias",
                              size.legend.name = "1/seTE",
                              label.unbiased = "Unbiased",
                              label.biased = "Biased") {

  seTE = x = lower_pi = upper_pi = median_pi = lower_ci = upper_ci = mean_ci = NULL

  if (!covariate %in% names(object$data)) {
    stop(paste("Covariate '", covariate, "' not found in the original data.", sep = ""))
  }

  plot_data = object$data
  plot_data$I = object$BUGSoutput$mean$I[1:nrow(plot_data)]
  plot_data$col = ifelse(plot_data$I > 0.5, "biased", "unbiased")

  gamma_name = paste0("gamma_", covariate)
  mcmc_samples = object$BUGSoutput$sims.matrix
  mu_theta_samples = mcmc_samples[, "mu.theta"]
  gamma_samples = mcmc_samples[, gamma_name]
  tau_theta_samples = mcmc_samples[, "tau.theta"]

  # --- Build grid of x values for continuous lines ---
  xgrid = seq(min(plot_data[[covariate]]), max(plot_data[[covariate]]), length.out = 100)

  # Calculate fitted values for each posterior draw on the xgrid
  fit_matrix_ci = sapply(1:length(mu_theta_samples), function(i) {
    mu_theta_samples[i] + gamma_samples[i] * xgrid
  })

  # Take quantiles for the Credible Interval (CI)
  fit_quantiles_ci = apply(fit_matrix_ci, 1, quantile, probs = c(0.025, 0.5, 0.975))
  ci_df = data.frame(
    x = xgrid,
    lower_ci = fit_quantiles_ci[1,],
    mean_ci = fit_quantiles_ci[2,],
    upper_ci = fit_quantiles_ci[3,]
  )

  # Generate posterior predictive samples for each point on the grid
  y_new_matrix = sapply(xgrid, function(x) {
    # For each x-grid point, calculate the mean of the normal distribution
    mean_pred = mu_theta_samples + gamma_samples * x

    # Draw a new observation (y_new) for each MCMC sample with a zero-mean error term, then add it
    errors = rnorm(n = length(mean_pred), mean = 0, sd = tau_theta_samples)
    return(mean_pred + errors)
  })

  # Calculate quantiles for the Prediction Interval (PI), including the median
  fit_quantiles_pi = apply(y_new_matrix, 2, quantile, probs = c(0.025, 0.5, 0.975))
  pi_df = data.frame(
    x = xgrid,
    lower_pi = fit_quantiles_pi[1,],
    median_pi = fit_quantiles_pi[2,],
    upper_pi = fit_quantiles_pi[3,]
  )

  # Start plotting
  p = ggplot(plot_data, aes_string(x = covariate, y = "TE")) +
    theme_bw() +
    labs(x = x.lab, y = y.lab, title = title.plot) +
    geom_hline(yintercept = 0, colour = "black", size = 0.5)

  # Add Prediction Interval (PI) lines
  if (show_pi) {
    p = p +
      geom_line(data = pi_df, aes(x = x, y = lower_pi), color = pi.color, linetype = 1) +
      geom_line(data = pi_df, aes(x = x, y = upper_pi), color = pi.color, linetype = 1) +
      geom_line(data = pi_df, aes(x = x, y = median_pi), color = pi.mean.line.color, linetype = 2, lwd = 0.8) # New median line
  }

  # Add Credible Interval (CI) lines
  if (show_ci) {
    p = p +
      geom_line(data = ci_df, aes(x = x, y = lower_ci), color = ci.color, linetype = 1) +
      geom_line(data = ci_df, aes(x = x, y = upper_ci), color = ci.color, linetype = 1)
  }

  # Add posterior mean regression line (from CI data frame)
  p = p + geom_line(data = ci_df, aes(x = x, y = mean_ci), colour = mean.line.color, lwd = 1, linetype = 1)

  # Add observed data
  p = p +
    geom_errorbar(aes(ymin = TE - 2*seTE, ymax = TE + 2*seTE), width = errorbar_width) +
    geom_point(aes_string(x = covariate, y = "TE", size = "1/seTE", fill = "col"), shape = 21) +
    scale_fill_manual(values = c("biased" = point.color.biased, "unbiased" = point.color.unbiased),
                      breaks = c("unbiased", "biased"),
                      labels = c(label.unbiased, label.biased),
                      name = legend.name) +
    guides(size = guide_legend(title = size.legend.name))

  return(p)
}
