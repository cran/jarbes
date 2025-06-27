#' @title Simulate Data for Meta-Analysis
#'
#' @description This function simulates aggregated data for a meta-analysis, introducing biased studies at different levels.
#'
#' @param mu Scalar with the true pooled effect value.
#' @param sigma Scalar with the true intra-study standard deviation.
#' @param n.total A vector with the sample sizes of the studies.
#' @param tau Scalar with the between-studies standard deviation.
#' @param N Scalar with the total number of studies in the meta-analysis.
#' @param mu.beta.1 Scalar with the mean bias of studies in the mild bias class.
#' @param mu.beta.2 Scalar with the mean bias of studies in the large bias class.
#' @param mu.beta.3 Scalar with the mean bias of studies in the extreme bias class.
#' @param n.B.1 Scalar with the number of studies in the mild bias class.
#' @param n.B.2 Scalar with the number of studies in the large bias class.
#' @param n.B.3 Scalar with the number of studies in the extreme bias class.
#'
#' @return A dataframe with columns:
#' 	\item{TE}{Observed study's effect.}
#' 	\item{seTE}{Standard error of the study's effect.}
#' 	\item{theta}{True study's effect.}
#' 	\item{n.total}{Sample size of the study.}
#' 	\item{B.flag}{Bias category: "No B", "Mild B", "Large B", "Extreme B".}
#'
#' @examples
#'
#' set.seed(123)
#' simData(mu = 0, sigma = 1, n.total = rep(100, 10), tau = 0.5, N = 10,
#'         mu.beta.1 = 0.2, mu.beta.2 = 0.5, mu.beta.3 = 1,
#'         n.B.1 = 2, n.B.2 = 2, n.B.3 = 2)
#'
#' @export
simData <- function(mu, sigma, n.total, tau, N,
                    mu.beta.1, mu.beta.2, mu.beta.3,
                    n.B.1, n.B.2, n.B.3) {
  UseMethod("simData")
}


#' @export
# Define simData function
simData.default <- function(mu, sigma, n.total, tau, N, mu.beta.1, mu.beta.2, mu.beta.3, n.B.1, n.B.2, n.B.3) {
  theta <- rnorm(N, mean = mu, sd = tau)
  SE.y <- sigma * sqrt(1 / n.total)
  TE <- rnorm(N, mean = theta, sd = SE.y)
  Study.id <- paste("Study", 1:N, sep = "-")

  # Initialize bias flag
  B.flag <- rep("No B", N)

  # Apply biases if present
  if(any(c(n.B.1, n.B.2, n.B.3)>0)){
    case.B = which(c(n.B.1, n.B.2, n.B.3)>0)
    case.B = paste0(case.B, collapse = "")


    switch(case.B,
           "1" = {
             set.pick <- 1:N
             pick.outlier.1 <- sample(set.pick, n.B.1, replace = FALSE)
             TE[pick.outlier.1] <- TE[pick.outlier.1] + mu.beta.1
             B.flag[pick.outlier.1] <- "Mild B"
           },
           "2" = {
             set.pick <- 1:N
             pick.outlier.2 <- sample(set.pick, n.B.2, replace = FALSE)
             TE[pick.outlier.2] <- TE[pick.outlier.2] + mu.beta.2
             B.flag[pick.outlier.2] <- "Large B"
           },
           "3" = {
             set.pick <- 1:N
             pick.outlier.3 <- sample(set.pick, n.B.3, replace = FALSE)
             TE[pick.outlier.3] <- TE[pick.outlier.3] + mu.beta.3
             B.flag[pick.outlier.3] <- "Extreme B"
           },
           "12" = {
             set.pick <- 1:N
             pick.outlier.1 <- sample(set.pick, n.B.1, replace = FALSE)
             set.pick <- setdiff(set.pick, pick.outlier.1)
             pick.outlier.2 <- sample(set.pick, n.B.2, replace = FALSE)
             TE[pick.outlier.1] <- TE[pick.outlier.1] + mu.beta.1
             B.flag[pick.outlier.1] <- "Mild B"
             TE[pick.outlier.2] <- TE[pick.outlier.2] + mu.beta.2
             B.flag[pick.outlier.2] <- "Large B"
           },
           "23" = {
             set.pick <- 1:N
             pick.outlier.2 <- sample(set.pick, n.B.2, replace = FALSE)
             set.pick <- setdiff(set.pick, pick.outlier.2)
             pick.outlier.3 <- sample(set.pick, n.B.3, replace = FALSE)
             TE[pick.outlier.2] <- TE[pick.outlier.2] + mu.beta.2
             B.flag[pick.outlier.2] <- "Large B"
             TE[pick.outlier.3] <- TE[pick.outlier.3] + mu.beta.3
             B.flag[pick.outlier.3] <- "Extreme B"
           },
           "13" = {
             set.pick <- 1:N
             pick.outlier.1 <- sample(set.pick, n.B.1, replace = FALSE)
             set.pick <- setdiff(set.pick, pick.outlier.1)
             pick.outlier.3 <- sample(set.pick, n.B.3, replace = FALSE)
             pick.outlier <- c(pick.outlier.1, pick.outlier.3)

             TE[pick.outlier.1] <- TE[pick.outlier.1] + mu.beta.1
             B.flag[pick.outlier.1] <- "Mild B"

             TE[pick.outlier.3] <- TE[pick.outlier.3] + mu.beta.3
             B.flag[pick.outlier.3] <- "Extreme B"
           },
           "123" = {
             set.pick <- 1:N
             pick.outlier.1 <- sample(set.pick, n.B.1, replace = FALSE)
             pick.outlier.2 <- sample(set.pick[-pick.outlier.1], n.B.2, replace = FALSE)
             pick.outlier.3 <- sample(set.pick[-c(pick.outlier.1, pick.outlier.2)], n.B.3, replace = FALSE)

             TE[pick.outlier.1] <- TE[pick.outlier.1] + mu.beta.1
             B.flag[pick.outlier.1] <- "Mild B"

             TE[pick.outlier.2] <- TE[pick.outlier.2] + mu.beta.2
             B.flag[pick.outlier.2] <- "Large B"

             TE[pick.outlier.3] <- TE[pick.outlier.3] + mu.beta.3
             B.flag[pick.outlier.3] <- "Extreme B"
           }
    )
  } else {
    cat("No biased studies in this simulation \n")
  }

  # Create dataframe
  sims <- data.frame(TE, seTE = SE.y, theta, n.total, B.flag)
  rownames(sims) <- Study.id

  # Assign class "simData"
  class(sims) <- c("simData", class(sims))
  attr(sims, "mu") <- mu
  attr(sims, "sigma") <- sigma
  attr(sims, "tau") <- tau
  attr(sims, "N") <- N
  attr(sims, "mu.beta.1") <- mu.beta.1
  attr(sims, "mu.beta.2") <- mu.beta.2
  attr(sims, "mu.beta.3") <- mu.beta.3
  attr(sims, "n.B.1") <- n.B.1
  attr(sims, "n.B.2") <- n.B.2
  attr(sims, "n.B.3") <- n.B.3

  return(sims)
}
#' Generic print function for simData object.
#'
#'
#' @param x An object of class \code{simData} generated by the \code{simData} function.
#' @param ... Additional arguments passed to \code{print.data.frame}.
#'
#' @export
print.simData <- function(x, ...) {
  print.data.frame(x)
}
