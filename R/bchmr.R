#' Bias-Corrected Bayesian Nonparametric Model to combine aggregated and individual participant data for cross design synthesis.
#'
#'
#' This function performers a Bayesian cross design synthesis. The function fits a
#' hierarchical meta-regression model based on a BC-BNP model
#'
#'  The model is experimental and under construction for the version 2.2.5 (March 2025)
#'
#' @param data            Aggregated data results: a data frame where the first four columns containing the number of events in
#'                        the control group (yc), the number of patients in the control group (nc),
#'                        the number of events in the treatment group (yt) and the number of patients in the
#'                        treatment group (nt). If two.by.two = TRUE a data frame where each line
#'                        contains the trial results with column names: yc, nc, yt, nt.
#'
#' @param two.by.two      If TRUE indicates that the trial results are with names: yc, nc, yt, nt.
#'
#' @param dataIPD         Individual participant data: a data frame where
#'                        the first column is the outcome variable and the
#'                        other columns represent individual participant
#'                        charachteristics.
#'
#'
#' @param re              Random effects distribution for the resulting model. Possible
#'                        values are \emph{normal} for bivariate random effects and \emph{sm}
#'                        for scale mixtures.
#'
#'
#' @param mean.mu.1       Prior mean of baseline risk, default value is 0.
#'
#' @param mean.mu.2       Prior mean of treatment effect, default value is 0.
#'
#' @param mean.mu.phi     Prior mean of the bias parameter which measures the difference
#'                        between the baseline mean mu.1 and the intercept parameter of
#'                        the logistic regression of the individual participant data.
#'                        The defalut vaule is 0.
#'
#' @param sd.mu.1         Prior standard deviation of mu.1, default value is 1. The default prior of mu.1 is a
#'                        logistic distribution with mean 0 and dispersion 1. The implicit prior for mu.1 in
#'                        the probability scale is a uniform between 0 and 1.
#'
#' @param sd.mu.2         Prior standard deviation of mu.2, default value is 1. The default prior of mu.2 is a
#'                        logistic distribution with mean 0 and dispersion 1. The implicit prior for mu.2 in
#'                        the probability scale is a uniform between 0 and 1.
#'
#' @param sd.mu.phi       Prior standard deviation of mu.phi, default value is 1.
#'
#' @param sigma.1.upper   Upper bound of the uniform prior of sigma.1, default value is 5.
#'
#' @param sigma.2.upper   Upper bound of the uniform prior of sigma.2, default value is 5.
#'
#' @param sigma.beta.upper  Upper bound of the uniform prior of sigma.beta, default value is 5.
#'
#' @param mean.Fisher.rho Mean of rho in the Fisher scale, default value is 0.
#'
#' @param sd.Fisher.rho   Standard deviation of rho in the Fisher scale, default value is 1/sqrt(2).
#'
#'
#' @param nr.chains       Number of chains for the MCMC computations, default 5.
#'
#' @param nr.iterations   Number of iterations after adapting the MCMC, default is 10000. Some models may need more iterations.
#'
#' @param nr.adapt        Number of iterations in the adaptation process, default is 1000. Some models may need more iterations during adptation.
#'
#' @param nr.burnin       Number of iteration discarded for burnin period, default is 1000. Some models may need a longer burnin period.
#'
#' @param nr.thin         Thinning rate, it must be a positive integer, the default value 1.
#'
#' @return This function returns an object of the class "bchmr". This object contains the MCMC output of
#' each parameter and hyper-parameter in the model, the data frame used for fitting the model, and further model outputs.
#'
#' The results of the object of the class hmr can be extracted with R2jags. In addition
#' a summary, a print and a plot function are implemented for this type of object.
#'
#
#'
#' @references Verde, P. E. (2019) Learning from Clinical Evidence: The Hierarchical Meta-Regression Approach. Biometrical Journal. Biometrical Journal; 1-23.
#'
#' @references Verde, P.E., and Rosner, G.L. (2025), A Bias-Corrected Bayesian Nonparametric Model for Combining Studies With Varying Quality in Meta-Analysis. Biometrical Journal., 67: e70034. https://doi.org/10.1002/bimj.70034
#'
#'
#' @examples
#'
#'\dontrun{
#'library(jarbes)
#'
#' data("healing")
#' AD <- healing[, c("y_c", "n_c", "y_t", "n_t")]
#'
#' data("healingipd")
#'
#' IPD <- healingipd[, c("healing.without.amp", "PAD", "neuropathy",
#' "first.ever.lesion", "no.continuous.care", "male", "diab.typ2",
#' "insulin", "HOCHD", "HOS", "CRF", "dialysis", "DNOAP", "smoking.ever",
#' "diabdur", "wagner.class")]
#'
#' mx1 <- bchmr(AD, two.by.two = FALSE,
#'            dataIPD = IPD,
#'            re = "normal",
#'            sd.mu.1 = 2,
#'            sd.mu.2 = 2,
#'            sd.mu.phi = 2,
#'            sigma.1.upper = 5,
#'            sigma.2.upper = 5,
#'            sigma.beta.upper = 5,
#'            sd.Fisher.rho = 1.25,
#'            df.estimate = FALSE,
#'            df.lower = 3,
#'            df.upper = 10,
#'            nr.chains = 1,
#'            nr.iterations = 1500,
#'            nr.adapt = 100,
#'            nr.thin = 1)
#'
#' print(mx1)
#'
#'
#' # End of the examples.
#'
#' }
#'
#'
#' @import R2jags
#' @import rjags
#' @import graphics
#' @import stats
#' @import graphics
#'
#' @export
#'

bchmr <- function(data,
                two.by.two = TRUE,
                dataIPD,
                # Arguments for the model:
                re              = "normal",
                # Hyperpriors parameters............................................
                mean.mu.1       = 0,
                mean.mu.2       = 0,
                mean.mu.phi     = 0,
                sd.mu.1         = 1,
                sd.mu.2         = 1,
                sd.mu.phi       = 1,
                sigma.1.upper   = 5,
                sigma.2.upper   = 5,
                sigma.beta.upper = 5,
                mean.Fisher.rho = 0,
                sd.Fisher.rho   = 1/sqrt(2),
               # MCMC setup........................................................
                nr.chains       = 2,
                nr.iterations   = 10000,
                nr.adapt        = 1000,
                nr.burnin       = 1000,
                nr.thin         = 1)UseMethod("bchmr")

#' @export
#'
bchmr.default <- function(
    data,
    two.by.two = TRUE,
    dataIPD,
    # Arguments for the model:
    re              = "normal",
    # Hyperpriors parameters............................................
    mean.mu.1       = 0,
    mean.mu.2       = 0,
    mean.mu.phi     = 0,
    sd.mu.1         = 1,
    sd.mu.2         = 1,
    sd.mu.phi       = 1,
    sigma.1.upper   = 5,
    sigma.2.upper   = 5,
    sigma.beta.upper = 5,
    mean.Fisher.rho = 0,
    sd.Fisher.rho   = 1/sqrt(2),
    # MCMC setup........................................................
    nr.chains       = 2,
    nr.iterations   = 10000,
    nr.adapt        = 1000,
    nr.burnin       = 1000,
    nr.thin         = 1)
{



  # Setting up hyperparameters ...
  pre.mu.1 <- 1/(sd.mu.1*sd.mu.1)
  pre.mu.2 <- 1/(sd.mu.2*sd.mu.2)
  pre.mu.phi <- 1/ (sd.mu.phi * sd.mu.phi)
  pre.Fisher.rho <- 1/(sd.Fisher.rho * sd.Fisher.rho)

  # Setting up data nodes ...

  if(two.by.two == FALSE)
  {
    yc <- data[,1]
    nc <- data[,2]
    yt <- data[,3]
    nt <- data[,4]
  } else
  {
    yc <- data$yc
    nc <- data$nc
    yt <- data$yt
    nt <- data$nt
  }

  N <- length(nc)
  # Increase N for one observational study

  N <- N + 1


  # Setup IPD data for regression ................................................
  dataIPD$y <- dataIPD[,1]
  dataIPD <- dataIPD[,-1]       # <- corregido! Se generó un bug version > 1.7.0

  model.1 <- model.frame(y ~ ., data = dataIPD)
  X.IPD <- model.matrix(model.1, data = dataIPD)
  X.IPD <- X.IPD[ , dimnames(X.IPD)[[2]] != "(Intercept)", drop=FALSE]
  y.0.IPD <- model.response(model.1)

  if(is.factor(y.0.IPD)){y.0.IPD <- unclass(y.0.IPD) - 1}

  K <- dim(X.IPD)[2]
  M <- dim(X.IPD)[1]

  #.............................................................................
  # Data, initial values and parameters  .......................................
  data.model <-
    list(M = M,
         K = K,
         X.IPD = X.IPD,
         y.0.IPD = y.0.IPD,
         N = N,
         yc = yc,
         nc = nc,
         yt = yt,
         nt = nt,
         mean.mu.1 = mean.mu.1,
         mean.mu.2 = mean.mu.2,
         mean.mu.phi = mean.mu.phi,
         pre.mu.1 = pre.mu.1,
         pre.mu.2 = pre.mu.2,
         pre.mu.phi = pre.mu.phi,
         sigma.1.upper = sigma.1.upper,
         sigma.2.upper = sigma.2.upper,
         sigma.beta.upper = sigma.beta.upper,
         mean.Fisher.rho = mean.Fisher.rho,
         pre.Fisher.rho = pre.Fisher.rho
    )

   # Parameters to monitor ....................................................................
  parameters.model <-
    c("mu.1",
      "mu.2",
      "mu.phi",
      "sigma.1",
      "sigma.2",
      "rho",
      "Odds.pool",
      #    "Odds.new",
      "P_control.pool",
      #    "P_control.new",
      "alpha.0",
      "alpha.1",
      "beta.IPD",
      "sigma.beta",
      "x.subgroup",
      "eta.subgroup" # Verde's formula ....
    )

  # Model BUGS script

  # Block for data model ......................................................................
    dm <-
      "
    model
    {
    for(i in 1:(N-1))
    {
    yc[i] ~ dbin(pc[i], nc[i])
    yt[i] ~ dbin(pt[i], nt[i])
    "


   #----
    # Block for the link function................................................................
    link.logit <-
      "
    # Random effects model
    logit(pc[i]) <- theta.1[i]
    logit(pt[i]) <- theta.2[i] + logit(pc[i])
    "

    # Block for structural distribution .........................................................
    re.normal <-
      "
    theta.1[i] ~ dnorm(mu.1, pre.theta.1)
    theta.2[i] ~ dnorm(mu.2.1[i], pre.theta.2.1)

    #Conditional mean
    mu.2.1[i] <- mu.2 + rho * sigma.2 / sigma.1 * (theta.1[i] - mu.1)

    }

    # Hyper priors
    mu.1 ~  dlogis(mean.mu.1, pre.mu.1)
    mu.2 ~  dlogis(mean.mu.2, pre.mu.2)

    # Dispersion parameters
    sigma.1 ~ dunif(0, sigma.1.upper)
    sigma.2 ~ dunif(0, sigma.2.upper)
    pre.theta.1 <- 1/(sigma.1*sigma.1)
    pre.theta.2 <- 1/(sigma.2*sigma.2)

    # Conditional precision
    pre.theta.2.1 <- pre.theta.2 / (1 - rho*rho)

    # Correlation
    z ~ dnorm(mean.Fisher.rho, pre.Fisher.rho)
    rho <- 2*exp(z)/(1+exp(z)) - 1

    # Predictions ...
    mu.new[1] <- mu.1
    mu.new[2] <- mu.2

    Sigma.new[1, 1] <- pow(sigma.1, 2)
    Sigma.new[2, 2] <- pow(sigma.2, 2)
    Sigma.new[1, 2] <- rho * sigma.1 * sigma.2
    Sigma.new[2, 1] <- Sigma.new[1, 2]

    Omega.new[1:2, 1:2] <- inverse(Sigma.new[1:2, 1:2])
    theta.new[1:2] ~ dmnorm(mu.new[1:2], Omega.new[1:2 ,1:2])

    # Functional parameters
    alpha.0 <- (mu.2 - alpha.1 * mu.1)
    alpha.1 <- rho * sigma.2/sigma.1

    # Block for: link = logistic, re = normal and split.w = FALSE
    # Model for individual data ..................................................
    # Random effect for the observational study ..................................

    # Introduced mu.phi...........................................................
    mu.phi  ~ dlogis(mean.mu.phi, pre.mu.phi) # non-informative

    mu.obs <- mu.1 + mu.phi

    # Change mu.1 to mu.obs = mu.1 + mu.phi .......................................
    mu.2.1[N] <- mu.2 + rho * sigma.2 / sigma.1 * (theta.1[N] - mu.obs)
    theta.1[N] ~ dnorm(mu.obs, pre.theta.1)
    theta.2[N] ~ dnorm(mu.2.1[N], pre.theta.2.1)

    # Priors for beta.IPD .........................................................
    for( i in 1:K)
    {
    beta.IPD[i] ~ dnorm(0, pre.beta.IPD)
    }

    pre.beta.IPD <- 1/(sigma.beta*sigma.beta)
    sigma.beta ~ dunif(0, sigma.beta.upper)

    "


    #...........................................................................
    # Block of parameters of interest depending on the links ......................
    par.logit <- "
#..............................................................................
for(i in 1:M){
    y.0.IPD[i] ~ dbern(p0.IPD[i])
    logit(p0.IPD[i]) <- theta.1[N] + inprod(beta.IPD[1:K], X.IPD[i,1:K])
  }

# Verde's formula to estimate relative treatment effect in a subgroup .........
for(i in 1:K){
   x.subgroup[i] = mean(mu.phi + beta.IPD[i])
eta.subgroup[i]  = alpha.0 + alpha.1*(x.subgroup[i] - mu.1)
}

  # Parameters of interest
  # Pooled summaries ...
       Odds.pool <- exp(alpha.0)      # Here was a bug in version < 2.0.0
  P_control.pool <- ilogit(mu.1)

  }"

# List of possible models ...

# normal random effects
model.bugs <- paste(dm, link.logit, re.normal, par.logit)
#----

model.bugs.connection <- textConnection(model.bugs)

# Use R2jags as interface for JAGS ...
results <- jags(              data = data.model,
                              parameters.to.save = parameters.model,
                              model.file = model.bugs.connection,
                              n.chains = nr.chains,
                              n.iter = nr.iterations,
                              n.burnin = nr.burnin,
                              n.thin = nr.thin,
                              pD = TRUE)

# Close text conecction
close(model.bugs.connection)

# Extra outputs that are linked with other functions ...
IPD.names <- names(dataIPD)
IPD.names <- IPD.names[IPD.names!="y"]
results$beta.names <- paste("beta.", IPD.names, sep = "")

results$re <- re
results$data <- data
results$two.by.two <- two.by.two

# Priors ...

results$prior$mean.mu.1       = mean.mu.1
results$prior$mean.mu.1       = mean.mu.1
results$prior$mean.mu.phi     = mean.mu.phi
results$prior$sd.mu.1         = sd.mu.1
results$prior$sd.mu.2         = sd.mu.2
results$prior$sd.mu.phi       = sd.mu.phi
results$prior$sigma.1.upper   = sigma.1.upper
results$prior$sigma.2.upper   = sigma.2.upper
results$prior$mean.Fisher.rho = mean.Fisher.rho
results$prior$sd.Fisher.rho   = sd.Fisher.rho

class(results) <- c("bchmr")

return(results)
}




#' Generic print function for hmr object in jarbes.
#' @param x The object generated by the function bchmr.
#'
#' @param digits The number of significant digits printed. The default value is 3.
#'
#' @param intervals A numeric vector of probabilities with values in [0,1]. The default value is
#'                  intervals = c(0.025, 0.5, 0.975).
#' @param ... \dots
#'
#' @export


print.bchmr <- function(x, digits = 3,
                      intervals = c(0.025, 0.25, 0.5, 0.75, 0.975), ...)
{
  m <- x
  x <- x$BUGSoutput
  sims.matrix <- x$sims.matrix
  mu.vect <- apply(sims.matrix, 2, mean)
  sd.vect <- apply(sims.matrix, 2, sd)
  int.matrix <- apply(sims.matrix, 2, quantile, intervals)
  if (x$n.chains>1) {
    n.eff <- x$summary[,"n.eff"]
    Rhat <- x$summary[,"Rhat"]
  } else {
    n.eff <- Rhat <- NULL
  }
  summaryMatrix <- t(rbind(mu.vect, sd.vect, int.matrix, Rhat, n.eff))
  rownameMatrix <- rownames(summaryMatrix)

  # Names of the regression coefficients ...
  ind.names <- grep("beta.IPD", rownameMatrix)
  rownameMatrix[ind.names] <- m$beta.names

  rownames(summaryMatrix) <- rownameMatrix


  dev.idx <- match("deviance", rownameMatrix)
  if(any(!is.na(dev.idx))){
    summaryMatrix <- rbind(summaryMatrix[-dev.idx,], summaryMatrix[dev.idx,])
    rownames(summaryMatrix) <- c(rownameMatrix[-dev.idx], rownameMatrix[dev.idx])
  }




  if (!is.null(m$model.file))
    cat("Inference for Bugs model at \"", m$model.file, "\", ",
        sep = "")
  if (!is.null(m$program))
    cat("fit using ", m$program, ",", sep = "")
  cat("\n ", m$n.chains, " chains, each with ", m$n.iter, " iterations (first ",
      x$n.burnin, " discarded)", sep = "")
  if (x$n.thin > 1)
    cat(", n.thin =", x$n.thin)
  cat("\n n.sims =", x$n.sims, "iterations saved\n")
  print(round(summaryMatrix, digits), ...)
  if (x$n.chains > 1) {
    cat("\nFor each parameter, n.eff is a crude measure of effective sample size,")
    cat("\nand Rhat is the potential scale reduction factor (at convergence, Rhat=1).\n")
  }
  if (x$isDIC) {
    msgDICRule <- ifelse(x$DICbyR, "(using the rule, pD = var(deviance)/2)",
                         "(using the rule, pD = Dbar-Dhat)")
    cat(paste("\nDIC info ", msgDICRule, "\n", sep = ""))
    if (length(x$DIC) == 1) {
      cat("pD =", fround(x$pD, 1), "and DIC =", fround(x$DIC,
                                                       1))
    }
    else if (length(x$DIC) > 1) {
      print(round(x$DIC, 1))
    }
    cat("\nDIC is an estimate of expected predictive error (lower deviance is better).\n")
  }
  invisible(x)

}

fround <- function (x, digits) {
  format (round (x, digits), nsmall=digits)
}

pfround <- function (x, digits) {
  print (fround (x, digits), quote=FALSE)
}









