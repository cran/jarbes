#' @title Bayesian Meta-Analysis for Combining Studies
#'
#' @description This function performers a Bayesian meta-analysis
#'
#'
#'
#' @param data                A data frame with at least three columns with the following names:
#'                             1) TE = treatment effect,
#'                             2) seTE = the standard error of the treatment effect.
#'                             3) design = indicates study type or clustering subgroup.
#'
#' @param mean.mu.0             Prior mean of the overall mean parameter mu.0 (mean across designs), default value is 0.
#'
#' @param sd.mu.0               Prior standard deviation of mu.0 (mean across designs), the default value is 10.
#'
#' @param scale.sigma.between Prior scale parameter for scale gamma distribution for the
#'                            precision between study types. The default value is 0.5.
#'
#' @param df.scale.between    Degrees of freedom of the scale gamma distribution for the precision between study types.
#'                            The default value is 1, which results in a Half Cauchy distribution for the standard
#'                            deviation between studies. Larger values e.g. 30 corresponds to a Half Normal distribution.

#' @param scale.sigma.within Prior scale parameter for scale gamma distribution for the
#'                            precision within study types. The default value is 0.5.
#'
#' @param df.scale.within    Degrees of freedom of the scale gamma distribution for the precision within study types.
#'                            The default value is 1, which results in a Half Cauchy distribution for the standard
#'                            deviation between studies. Larger values e.g. 30 corresponds to a Half Normal distribution.
#'
#' @param nr.chains           Number of chains for the MCMC computations, default 2.
#' @param nr.iterations       Number of iterations after adapting the MCMC, default is 10000. Some models may need more iterations.
#' @param nr.adapt            Number of iterations in the adaptation process, default is 1000. Some models may need more iterations during adptation.
#' @param nr.burnin           Number of iteration discard for burn-in period, default is 1000. Some models may need a longer burnin period.
#' @param nr.thin             Thinning rate, it must be a positive integer, the default value 1.
#' @param be.quiet            Do not print warning message if the model does not adapt. The default value is FALSE. If you are not sure about the adaptation period choose be.quiet=TRUE.
#' @param r2jags              Which interface is used to link R to JAGS (rjags and R2jags), default value is R2Jags=TRUE.
#'
#' @return                    This function returns an object of the class "bmeta". This object contains the MCMC
#'                            output of each parameter and hyper-parameter in the model and
#'                            the data frame used for fitting the model.
#'
#'
#' @details The results of the object of the class bcmeta can be extracted with R2jags or with rjags. In addition a summary, a print and a plot functions are
#' implemented for this type of object.
#'
#'
#' @references Verde, P.E. (2021) A Bias-Corrected Meta-Analysis Model for Combining Studies of Different Types and Quality. Biometrical Journal; 1–17.
#'
#'
#' @examples
#' \dontrun{
#' library(jarbes)
#'
# data(ppvipd)
# b3lm1 = b3lmeta(ppvipd)
# summary(b3lm1, digits = 3, y.lim= c(0,2.2))
# plot(b3lm1, y.lim = c(0, 2.3))
#
# diagnostic(b3lm1, post.p.value.cut = 0.1, shape.forest = 4, lwd.forest = 0.75)
#
#'
#' }
#'
#' @import R2jags
#' @import rjags
#'
#'
#' @export
b3lmeta = function(
  data,
  # Hyperpriors parameters............................................
  mean.mu.0     = 0,
  sd.mu.0       = 10,
  scale.sigma.between = 0.5,
  df.scale.between    = 1,
  scale.sigma.within = 0.5,
  df.scale.within    = 1,

  # MCMC setup........................................................
  nr.chains       = 2,
  nr.iterations   = 10000,
  nr.adapt        = 1000,
  nr.burnin       = 1000,
  nr.thin         = 1,
  # Further options to link jags and R ...............................
  be.quiet        = FALSE,
  r2jags          = TRUE
)UseMethod("b3lmeta")



#' @export
b3lmeta.default = function(
  data,
  # Hyperpriors parameters............................................
  mean.mu.0     = 0,
  sd.mu.0       = 10,
  scale.sigma.between = 0.5,
  df.scale.between   = 1,
  scale.sigma.within = 0.5,
  df.scale.within    = 1,

  # MCMC setup........................................................
  nr.chains       = 2,
  nr.iterations   = 10000,
  nr.adapt        = 1000,
  nr.burnin       = 1000,
  nr.thin         = 1,

  # Further options to link jags and R ...............................
  be.quiet        = FALSE,
  r2jags          = TRUE
)
{

  # Data setup......................................................
  y = sort(data$TE)
  se.y = data$seTE[order(data$TE)]

  N = length(y)

 # if(N<3)stop("Low number of studies in the meta-analysis!")

  x = data$design[order(data$TE)]
  x = as.numeric(x)
  Ndesign <- length(table(x))

  # This list describes the data used by the BUGS script.
  data.b3lmeta <- list ("y", "se.y", "x", "N", "Ndesign",
                      "mean.mu.0", "sd.mu.0",
                      "scale.sigma.between",
                      "df.scale.between",
                      "scale.sigma.within",
                      "df.scale.within")

  # List of parameters
  par.b3lmeta  <- c("theta",
                    "mu.0",
                    "mu.0.new",
                    "mu",
                    "y.ghost",            # Approx Cross-validation...
                    "tau.theta.between",
                    "tau.theta.within",
                    "tau.theta.total")

  # Model in BUGS

  model.bugs <-
    "
  model
  {
   for( i in 1 : N) {
  # Likelihood of theta[i] ..........................................

            y[i] ~ dnorm(theta[i], pre.y[i])
        pre.y[i] <- pow(se.y[i], -2)

  # Grouped random effects ..........................................
    theta[i] ~ dnorm(mu[x[i]], pre.theta[x[i]])
      }

  # Between designs variability ....................................
    for(i in 1:Ndesign){
      mu[i] ~ dnorm(mu.0, pre.between)
    }


  # Priors for hyper-parameters ................................
    tau.theta.between <- 1/sqrt(pre.between)
    pre.between ~ dscaled.gamma(scale.sigma.between, df.scale.between)

    pre.mu.0 <- 1/sd.mu.0^2
    mu.0 ~ dnorm(mean.mu.0, pre.mu.0)

    for(i in 1:Ndesign){
      tau.theta.within[i] <- 1/sqrt(pre.theta[i])
      pre.theta[i] ~ dscaled.gamma(scale.sigma.within, df.scale.within)
      }

  # Predictive mu..............................................

  # Total variability

  tau.theta.total <- sqrt(sum(tau.theta.within[1:Ndesign]^2) + tau.theta.between^2)

  pre.total <- 1/tau.theta.total^2
  mu.0.new ~ dnorm(mu.0, pre.total)

# Approximate Bayesian Cross-Validation ......................
theta.ghost ~ dnorm(mu.0.new, pre.total)

for(i in 1:N)
{
    y.ghost[i] ~ dnorm(theta.ghost, pre.y[i])
}


  }
  "

model.bugs.connection <- textConnection(model.bugs)

if(r2jags == TRUE){
  # Use R2jags as interface for JAGS ...
  results <- jags(              data = data.b3lmeta,
                                parameters.to.save = par.b3lmeta,
                                model.file = model.bugs.connection,
                                n.chains = nr.chains,
                                n.iter = nr.iterations,
                                n.burnin = nr.burnin,
                                n.thin = nr.thin)
}
else {
  # Use rjags as interface for JAGS ...
  # Send the model to JAGS, check syntax, run ...
  jm <- jags.model(file     = model.bugs.connection,
                   data     = data.b3lmeta,
                   n.chains = nr.chains,
                   n.adapt  = nr.adapt,
                   quiet    = be.quiet)

  results <- coda.samples(jm,
                          variable.names = par.b3lmeta,
                          n.iter         = nr.iterations)
}

if(r2jags == FALSE)
{cat("You are using the package rjags as interface to JAGS.", "\n")
  cat("The plot functions for output analysis are not implemented in this jarbes version", "\n")
}

# Close text connection
close(model.bugs.connection)

# Extra outputs that are linked with other functions ...
results$data = data

# Hyperpriors parameters............................................
results$prior$mean.mu.0     = mean.mu.0
results$prior$sd.mu.0       = sd.mu.0
results$prior$scale.sigma.between = scale.sigma.between
results$prior$df.scale.between    = df.scale.between
results$prior$scale.sigma.within  = scale.sigma.within
results$prior$df.scale.within     = df.scale.within

results$N = N
results$x = x
results$Ndesign = Ndesign

class(results) = c("b3lmeta")

return(results)
}



#' Generic print function for bcmeta object in jarbes.
#'
#' @param x The object generated by the function bcmeta.
#'
#' @param digits The number of significant digits printed. The default value is 3.
#'
#' @param ... \dots
#'
#' @export
print.b3lmeta <- function(x, digits, ...)
{
  print(x$BUGSoutput,...)
}




#' Generic plot function for bcmeta object in jarbes.
#'
#' @param x The object generated by the bcmeta function.
#'
#' @param ... \dots
#'
#' @export
plot.b3lmeta <- function(x, ...)
{
  plot(x)
}


