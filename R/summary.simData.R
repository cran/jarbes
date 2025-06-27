#' Summary method for simData objects
#'
#' Provides a summary of the simulated meta-analysis dataset using the known simulation parameters.
#'
#' @param object A simData object created by simData().
#' @param digits The number of significant digits printed. The default value is 3.
#' @param ... Additional arguments (currently not used).
#'
#' @return A formatted summary of the meta-analysis simulation.
#' @export
summary.simData <- function(object, digits = 3, ...) {
  
  # Extract key statistics
  N <- nrow(object)
  mu <- attr(object, "mu")
  sigma <- attr(object, "sigma")
  tau <- attr(object, "tau")
  mu.beta.1 <- attr(object, "mu.beta.1")
  mu.beta.2 <- attr(object, "mu.beta.2")
  mu.beta.3 <- attr(object, "mu.beta.3")
  n.B.1 <- attr(object, "n.B.1")
  n.B.2 <- attr(object, "n.B.2")
  n.B.3 <- attr(object, "n.B.3")
  
  # Sample size statistics
  n.total <- object$n.total
  min_n <- min(n.total)
  median_n <- median(n.total)
  max_n <- max(n.total)
  
  # Bias summary
  n_bias_low <- sum(object$B.flag == "Mild B")
  n_bias_mid <- sum(object$B.flag == "Large B")
  n_bias_high <- sum(object$B.flag == "Extreme B")
  
  # Summary
  cat("The data is simulated using the following parameters:\n")
  cat("====================================\n")
  cat("   Simulated Meta-Analysis Summary  \n")
  cat("====================================\n")
  cat(sprintf("%-45s %d\n", "Number of studies:", N))
  cat(sprintf("%-45s %.*f\n", "Mean observed effect (mu):", digits, mu))
  cat(sprintf("%-45s %.*f\n", "Within-study SD (sigma):", digits, sigma))
  cat(sprintf("%-45s %.*f\n", "Between-study SD (tau):", digits, tau))
  cat("====================================\n")
  cat(" Bias Parameters\n")
  cat("------------------------------------\n")
  cat(sprintf("Mild Bias:    %d studies\n", n.B.1))
  cat(sprintf("Large Bias:   %d studies\n", n.B.2))
  cat(sprintf("Extreme Bias: %d studies\n", n.B.3))
  cat("------------------------------------\n")
  cat(" Sample Size Summary \n")
  cat("------------------------------------\n")
  cat(sprintf("Min sample size:     %d\n", min_n))
  cat(sprintf("Median sample size:  %d\n", median_n))
  cat(sprintf("Max sample size:     %d\n", max_n))
  cat("====================================\n")
  
}