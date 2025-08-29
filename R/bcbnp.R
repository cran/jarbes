#' @title Bias Corrected Meta-Analysis and Meta-Regression with Dirichlet Mixture Process of normal distributions
#'
#' @description This function performs a Bias Corrected Bayesian Nonparametric meta-analysis and meta-regression using a DPM of Normal distributions to the bias process
#'
#' @param data A data frame with at least two columns with the following names:
#' 1) TE = treatment effect,
#' 2) seTE = the standard error of the treatment effect.
#'
#' @param formula A formula for meta-regression (e.g., ~ x1 + x2). Use ~ 1 for no covariates.
#'
#' @param mean.mu.theta Prior mean of mu.theta, the default value is mean.mu.0 = 0.
#' @param sd.mu.theta Prior standard deviation of mu.theta, the default value is 10^-6.
#'
#' @param scale.sigma.between Prior for tau.theta: Scale parameter for scale.gamma distribution for the precision 1/tau.theta^2 between studies.
#' The default value is 0.5.
#' @param df.scale.between Prior for tau.theta: Degrees of freedom of the scale.gamma distribution for the precision 1/tau.theta^2 between studies.
#' The default value is 1, which results in a Half Cauchy distribution for the standard
#' deviation between studies tau.theta. Larger values e.g. 30 corresponds to a Half Normal distribution.
#'
#' @param mean.mu.k Prior mean for the mean mu.k bias cluster components. The default value is 0.
#' @param sd.mu.k Prior standard deviation for mu.k bias cluster component. The default is 10.
#'
#' @param g.0 Shape parameter the prior of the precision tau.k. The default value is 0.5 (tau.k~ dgamma(shape, rate))
#' @param g.1 Rate parameter the prior of the precision tau.k. The default value is 0.5 (tau.k~ dgamma(shape, rate))
#'
#' @param a.0 Parameter for the prior Beta distribution for the probability of bias. Default value is a0 = 0.5.
#' @param a.1 Parameter for the prior Beta distribution for the probability of bias. Default value is a1 = 1.
#'
#' @param alpha.0 Lower bound of the uniform prior for the concentration parameter for the DP,
#' the default value is 0.5.
#' @param alpha.1 Upper bound of the uniform prior for the concentration parameter for the DP,
#' the default value depends on the sample size, see the example below. We give as
#' working value alpha.1 = 2
#'
#' @param K Maximum number of clusters in the DPM, the default value depends on alpha.1, see the
#' example below. We give as working value K = 10.
#'
#' @param nr.chains Number of chains for the MCMC computations, default 2.
#' @param nr.iterations Number of iterations after adapting the MCMC, default is 10000. Some models may need more iterations.
#' @param nr.adapt Number of iterations in the adaptation process, default is 1000. Some models may need more iterations during adptation.
#' @param nr.burnin Number of iteration discard for burn-in period, default is 1000. Some models may need a longer burnin period.
#' @param nr.thin Thinning rate, it must be a positive integer, the default value 1.
#' @param parallel NULL -> jags, 'jags.parallel' -> jags.parallel execution
#'
#' @return This function returns an object of the class "bcbnp". This object contains the MCMC
#' output of each parameter and hyper-parameter in the model and
#' the data frame used for fitting the model.
#'
#'
#' @details The results of the object of the class bcbnp can be extracted with R2jags or with rjags. In addition a summary, a print and a plot functions are
#' implemented for this type of object.
#'
#'
#' @references Verde, P.E., and Rosner, G.L. (2025a), A Bias-Corrected Bayesian Nonparametric Model for Combining Studies With Varying Quality in Meta-Analysis. Biometrical Journal., 67: e70034. https://doi.org/10.1002/bimj.70034
#' @references Verde, P.E, and Rosner, G.L. (2025b), jarbes: an R package for bias corrected meta-analysis models
#'
#' @examples
#' \dontrun{
#' library(jarbes)
#'
#'
#' # Example: Stemcells
#'
#' data("stemcells")
#' stemcells$TE = stemcells$effect.size
#' stemcells$seTE = stemcells$se.effect
#'
#' # Beta(0.5, 1)
#' a.0 = 0.5
#' a.1 = 1
#'
#' # alpha.max
#' N = dim(stemcells)[1]
#' alpha.max = 1/5 *((N-1)*a.0 - a.1)/(a.0 + a.1)
#'
#' alpha.max
#'
#'# K.max
#' K.max = 1 + 5*alpha.max
#' K.max = round(K.max)
#'
#' K.max
#'
#' set.seed(20233)
#'
#' bcbnp.stemcell = bcbnp(stemcells
#' alpha.0 = 0.5,
#' alpha.1 = alpha.max,
#' a.0 = a.0,
#' a.1 = a.1,
#' K = K.max,
#' nr.chains = 4,
#' nr.iterations = 50000,
#' nr.adapt = 1000,
#' nr.burnin = 10000,
#' nr.thin = 4)
#'
#'
#' }
#'
#' @import R2jags
#' @import rjags
#'
#' @export
bcbnp <- function(
    data,
    formula = ~ 1,
    # Hyperpriors parameters............................................
    mean.mu.theta = 0,
    sd.mu.theta = 10,
    scale.sigma.between = 0.5,
    df.scale.between = 1,
    mean.mu.k = 0,
    sd.mu.k = 10,
    g.0 = 0.5,
    g.1 = 0.5,
    a.0 = 0.5,
    a.1 = 1,
    alpha.0 = 0.5,
    alpha.1 = 2,
    K = 10,
    # MCMC setup........................................................
    nr.chains = 2,
    nr.iterations = 10000,
    nr.adapt = 1000,
    nr.burnin = 1000,
    nr.thin = 1,
    parallel = NULL
) {
  # Coerce data to a data.frame to ensure compatibility
  if (inherits(data, "tbl_df")) {
    data <- as.data.frame(data)
  }

  # --- Checks ----------------------------------------------------------
  if (!is.null(parallel) && parallel != "jags.parallel")
    stop("The parallel option must be NULL or 'jags.parallel'")

  if (!all(c("TE", "seTE") %in% names(data)))
    stop("data must include 'TE' and 'seTE' columns.")

  # --- Outcome data with NA trick for prediction ------------------------
  y <- c(data$TE, NA)
  se.y <- c(data$seTE, 1)
  N <- length(y)

  # --- Design matrix and covariate parsing -----------------------------
  covars <- all.vars(formula)

  if (length(covars) == 0) {
    p.x <- 0
  } else {
    p.x <- length(covars)
  }

  has_covariates <- p.x > 0

  # Initialize JAGS model string components
  reg_term_list <- c("mu.theta")
  gamma_priors_block_list <- c()

  if (has_covariates) {
    for (v in covars) {
      if (is.numeric(data[[v]])) {
        reg_term_list <- c(reg_term_list, paste0("gamma_", v, " * ", v, "[i]"))
        gamma_priors_block_list <- c(gamma_priors_block_list, paste0("  # Prior for continuous covariate '", v, "'\n  gamma_", v, " ~ dnorm(0, 0.01)\n"))
      } else if (is.factor(data[[v]])) {
        reg_term_list <- c(reg_term_list, paste0("gamma_", v, "[", v, "[i]]"))
        gamma_priors_block_list <- c(gamma_priors_block_list,
                                     paste0(
                                       "  # Priors for factor '", v, "'\n",
                                       "  gamma_", v, "[1] <- 0\n",
                                       "  for (j in 2:p_", v, ") {\n",
                                       "    gamma_", v, "[j] ~ dnorm(0, 0.01)\n",
                                       "  }\n"
                                     ))
      } else {
        stop("Unsupported covariate type. Covariates must be numeric or factors.")
      }
    }
  }

  # Combine components into final JAGS model strings
  reg_term <- paste(reg_term_list, collapse = " + ")
  gamma_priors_block <- paste(gamma_priors_block_list, collapse = "\n")

  # --- Build model string dynamically ---------------------------------
  model.bugs <- paste0("
    model {
      for (i in 1:N) {
        # Likelihood with bias
        y[i] ~ dnorm(theta.bias[i], 1 / se.y[i]^2)
        theta.bias[i] <- theta[i] + (I[i] * beta[i])

        # Meta-(re)gression
        theta[i] ~ dnorm(", reg_term, ", inv.var.theta)

        # Bias model: spike & DP mixture
        I[i] ~ dbern(pi.B)
        beta[i]  ~ dnorm(mu.beta[group[i]], tau.beta[group[i]])
        group[i] ~ dcat(pi[1:K])

        # Study-to-cluster indicators
        for (j in 1:K) {
          gind[i, j] <- equals(j, group[i])
        }
      }

      # Stick-breaking weights
      q[1] ~ dbeta(1, alpha)
      p[1] <- q[1]
      for (k in 2:K) {
        q[k] ~ dbeta(1, alpha)
        p[k] <- q[k] * prod(1 - q[1:(k-1)])
      }

      # Normalized cluster probs
      for (k in 1:K) {
        pi[k] <- p[k] / sum(p[1:K])
      }

      # Priors for cluster-specific bias parameters
      for (k in 1:K) {
        mu.beta[k]  ~ dnorm(mean.mu.k, 1 / sd.mu.k^2)
        tau.beta[k] ~ dgamma(g.0, g.1)      # precision
        sigma.beta[k] <- pow(tau.beta[k], -0.5)
      }

    ", gamma_priors_block, "

      # Priors for overall mean and heterogeneity
      mu.theta ~ dnorm(mean.mu.theta, 1 / sd.mu.theta^2)
      inv.var.theta ~ dscaled.gamma(scale.sigma.between, df.scale.between)
      sigma.theta <- 1 / sqrt(inv.var.theta)
      tau.theta  <- 1 / sqrt(inv.var.theta)  # alias to match monitored name

      # Prior for bias probability and DP concentration
      pi.B  ~ dbeta(a.0, a.1)
      alpha ~ dunif(alpha.0, alpha.1)

      # Post-processing: co-clustering matrices
      for (i in 1:N) {
        new.group[i] <- (1 - I[i]) + I[i] * (group[i] + 1)
      }
      for (i in 1:N) {
        for (j in 1:N) {
          equalsmatrix.bias.1[i,j] <- equals(I[i], I[j])
          equalsmatrix.bias.2[i,j] <- equals(new.group[i], new.group[j])
        }
      }
      # Predictions for each x
      for (i in 1:N) {
        y.new.bias[i] ~ dnorm(theta.bias[i], 1/(se.y[i]^2))
        y.new[i] ~ dnorm(theta[i], 1/(se.y[i]^2))
      }
      # Count used clusters
      K.hat <- sum(cl[])
      for (j in 1:K) {
        sumind[j] <- sum(gind[, j])
        cl[j] <- step(sumind[j] - 1 + 0.01)
      }
    }
    ")

  # --- Data list for JAGS ---------------------------------------------
  data.bcbnp <- list(
    y = y,
    se.y = se.y,
    N = N,
    # Hyperparameters for the model of interest
    mean.mu.theta = mean.mu.theta,
    sd.mu.theta = sd.mu.theta,
    scale.sigma.between = scale.sigma.between,
    df.scale.between = df.scale.between,
    # Hyperparameters for the bias model
    mean.mu.k = mean.mu.k,
    sd.mu.k = sd.mu.k,
    g.0 = g.0,
    g.1 = g.1,
    a.0 = a.0,
    a.1 = a.1,
    alpha.0 = alpha.0,
    alpha.1 = alpha.1,
    K = K
  )

  # --- Parameters to monitor ------------------------------------------
  par.bcbnp <- c(
    "mu.theta","tau.theta","theta",
    "theta.bias","beta","mu.beta","sigma.beta","I","pi.B","pi","alpha",
    "K.hat","group","gind","equalsmatrix.bias.1","equalsmatrix.bias.2","new.group"
  )

  if (has_covariates) {
    for (v in covars) {
      if (is.numeric(data[[v]])) {
        data.bcbnp[[v]] <- c(data[[v]], mean(data[[v]], na.rm = TRUE))
        par.bcbnp <- c(par.bcbnp, paste0("gamma_", v))
      } else if (is.factor(data[[v]])) {
        mode_level <- names(which.max(table(data[[v]])))
        imputed_mode_value <- which(levels(data[[v]]) == mode_level)
        data.bcbnp[[v]] <- c(as.numeric(data[[v]]), imputed_mode_value)
        data.bcbnp[[paste0("p_", v)]] <- length(levels(data[[v]]))
        par.bcbnp <- c(par.bcbnp, paste0("gamma_", v))
      }
    }
    par.bcbnp <- c(par.bcbnp, "y.new", "y.new.bias")
  }

  # --- Fit with R2jags -------------------------------------------------
  if (is.null(parallel)) {
    model.bugs.connection <- textConnection(model.bugs)
    on.exit(close(model.bugs.connection), add = TRUE)
    results <- R2jags::jags(
      data = data.bcbnp,
      parameters.to.save = par.bcbnp,
      model.file = model.bugs.connection,
      n.chains = nr.chains,
      n.iter = nr.iterations,
      n.burnin = nr.burnin,
      n.thin = nr.thin,
      DIC = TRUE,
      pD = TRUE
    )
  } else if (parallel == "jags.parallel") {
    writeLines(model.bugs, "model.bugs")
    on.exit(unlink("model.bugs"), add = TRUE)
    results <- R2jags::jags.parallel(
      data = data.bcbnp,
      parameters.to.save = par.bcbnp,
      model.file = "model.bugs",
      n.chains = nr.chains,
      n.iter = nr.iterations,
      n.burnin = nr.burnin,
      n.thin = nr.thin,
      DIC = TRUE
    )
    results$BUGSoutput$pD <- results$BUGSoutput$DIC - results$BUGSoutput$mean$deviance
  }

  # --- Attach extras ---------------------------------------------------
  results$data <- data
  results$formula <- formula
  results$prior <- list(
    mean.mu = mean.mu.theta, sd.mu = sd.mu.theta,
    scale.sigma.between = scale.sigma.between, df.scale.between = df.scale.between,
    mean.mu.k = mean.mu.k, sd.mu.k = sd.mu.k, g.0 = g.0, g.1 = g.1, a.0 = a.0, a.1 = a.1,
    alpha.0 = alpha.0, alpha.1 = alpha.1, K = K
  )
  results$N <- N
  results$has_covariates <- has_covariates
  results$model_string <- model.bugs

  # --- Return results --------------------------------------------------
  class(results) <- "bcbnp"
  return(results)
}

#' @export
print.bcbnp <- function(x, digits, ...) {
  print(x$BUGSoutput, ...)
}
