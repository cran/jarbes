#' @title Bayesian Meta-Analysis for Combining Studies
#'
#' @description This function performers a Bayesian meta-analysis
#'
#'
#'
#' @param data                A data frame with at least two columns with the following names:
#'                             1) TE = treatment effect,
#'                             2) seTE = the standard error of the treatment effect.
#'
#' @param mean.mu             Prior mean of the overall mean parameter mu, default value is 0.
#'
#' @param sd.mu               Prior standard deviation of mu, the default value is 10.
#'
#' @param scale.sigma.between Prior scale parameter for scale gamma distribution for the
#'                            precision between studies. The default value is 0.5.
#'
#' @param df.scale.between    Degrees of freedom of the scale gamma distribution for the precision between studies.
#'                            The default value is 1, which results in a Half Cauchy distribution for the standard
#'                            deviation between studies. Larger values e.g. 30 corresponds to a Half Normal distribution.
#'
#'
#' @param nr.chains           Number of chains for the MCMC computations, default 2.
#' @param nr.iterations       Number of iterations after adapting the MCMC, default is 10000. Some models may need more iterations.
#' @param nr.adapt            Number of iterations in the adaptation process, default is 1000. Some models may need more iterations during adptation.
#' @param nr.burnin           Number of iteration discard for burn-in period, default is 1000. Some models may need a longer burnin period.
#' @param nr.thin             Thinning rate, it must be a positive integer, the default value 1.
#' @param be.quiet            Do not print warning message if the model does not adapt. The default value is FALSE. If you are not sure about the adaptation period choose be.quiet=TRUE.
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
#' @examples
#' \dontrun{
#' library(jarbes)
#'
#' #Example: ppvipd
#'
#' data(ppvipd)
#' bm1 = bmeta(ppvipd)
#'
#' summary(bm1)
#' plot(bm1, x.lim = c(-3, 1), y.lim = c(0, 3))
#'
#' diagnostic(bm1, study.names = ppvipd$name, post.p.value.cut = 0.1,
#'            lwd.forest = 1, shape.forest = 4)
#'
#' # Example: Stemcells
#'
#' data("stemcells")
#' stemcells$TE = stemcells$effect.size
#' stemcells$seTE = stemcells$se.effect
#'
#' bm2 = bmeta(stemcells)
#' summary(bm2)
#' plot(bm2, x.lim = c(-1, 7), y.lim = c(0, 1))
#'
#' diagnostic(bm2, study.names = stemcells$trial,
#'            post.p.value.cut = 0.05,
#'            lwd.forest = 0.5, shape.forest = 4)
#'
#' diagnostic(bm2, post.p.value.cut = 0.05,
#'            lwd.forest = 0.5, shape.forest = 4)
#' }
#'
#' @import R2jags
#' @import rjags
#'
#'
#' @export
bmeta = function(
  data,
  # Hyperpriors parameters............................................
  mean.mu     = 0,
  sd.mu       = 10,
  scale.sigma.between = 0.5,
  df.scale.between    = 1,

  # MCMC setup........................................................
  nr.chains       = 2,
  nr.iterations   = 10000,
  nr.adapt        = 1000,
  nr.burnin       = 1000,
  nr.thin         = 1,
  # Further options to link jags and R ...............................
  be.quiet        = FALSE
          )UseMethod("bmeta")



#' @export
bmeta.default = function(
  data,
  # Hyperpriors parameters............................................
  mean.mu     = 0,
  sd.mu       = 10,
  scale.sigma.between = 0.5,
  df.scale.between    = 1,

  # MCMC setup........................................................
  nr.chains       = 2,
  nr.iterations   = 10000,
  nr.adapt        = 1000,
  nr.burnin       = 1000,
  nr.thin         = 1,

  # Further options to link jags and R ...............................
  be.quiet        = FALSE
)
{

    # Data
     y = data$TE
  se.y = data$seTE
     N = length(y)

     if(N<3)warning("You have less than 3 studies. This is a low number
                    of studies for this meta-analysis!")
     # Approximate Bayesian Cross-Validation
     #y.ghost = rep(NA, N)

     # This list describes the data used by the BUGS script.
     data.bmeta <- list ("y", "se.y", "N",
                          "mean.mu",
                          "sd.mu",
                          "scale.sigma.between",
                          "df.scale.between")

     # List of parameters
     par.bmeta  <- c("theta",
                        "mu",
                     "mu.new",
                      "tau",
                     "y.ghost")

  # Model in BUGS

  model.bugs <-
  "
  model
  {
   for( i in 1 : N ) {
  # Likelihood of theta[i] ....................................

            y[i] ~ dnorm(theta[i], pre.y[i])
        pre.y[i] <- pow(se.y[i], -2)

# Random effects ..........................

          theta[i] ~ dnorm(mu, inv.var)
}


# Approximate Bayesian Cross-Validation ......................
    theta.ghost ~ dnorm(mu, inv.var)

for(i in 1:N)
{
    y.ghost[i] ~ dnorm(theta.ghost, pre.y[i])
}

# Priors for hyper-parameters ................................


# Variance components ...
   tau <- 1/sqrt(inv.var)
inv.var ~ dscaled.gamma(scale.sigma.between,
                          df.scale.between)

# Prior for mu
inv.var.mu <- pow(sd.mu, -2)
       mu   ~ dnorm(mean.mu, inv.var.mu)

# Predictive mu
mu.new ~ dnorm(mu, inv.var)

  }
  "

  model.bugs.connection <- textConnection(model.bugs)

  # Use R2jags as interface for JAGS ...
    results <- jags(              data = data.bmeta,
                                  parameters.to.save = par.bmeta,
                                  model.file = model.bugs.connection,
                                  n.chains = nr.chains,
                                  n.iter = nr.iterations,
                                  n.burnin = nr.burnin,
                                  n.thin = nr.thin)



  # Close text connection
  close(model.bugs.connection)

  # Extra outputs that are linked with other functions ...
  results$data = data

  # Hyperpriors parameters............................................
  results$prior$mean.mu     = mean.mu
  results$prior$sd.mu       = sd.mu
  results$prior$scale.sigma.between = scale.sigma.between
  results$prior$df.scale.between    = df.scale.between
  results$N

  class(results) = c("bmeta")

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
print.bmeta <- function(x, digits, ...)
{
  print(x$BUGSoutput,...)
}



