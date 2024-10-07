#' @title Bias Corrected Meta-Analysis with Dirichlet Process Priors
#'
#' @description This function performers a Bayesian meta-analysis with DP as random effects
#'
#'
#'
#' @param data                A data frame with at least two columns with the following names:
#'                             1) TE = treatment effect,
#'                             2) seTE = the standard error of the treatment effect.
#'
#' @param mean.mu.0            Prior mean of the mean of the base distribution default value is mean.mu.0 = 0.
#'
#' @param sd.mu.0              Prior standard deviation of the base distribution, the default value is 10.
#'
#' @param scale.sigma.between Prior scale parameter for scale gamma distribution for the
#'                            precision between studies. The default value is 0.5.
#'
#' @param df.scale.between    Degrees of freedom of the scale gamma distribution for the precision between studies.
#'                            The default value is 1, which results in a Half Cauchy distribution for the standard
#'                            deviation between studies. Larger values e.g. 30 corresponds to a Half Normal distribution.
#' @param B.lower             Lower bound of the bias parameter B, the default value is 0.
#' @param B.upper             Upper bound of the bias parameter B, the default value is 10.
#'
#' @param a.0                 Parameter for the prior Beta distribution for the probability of bias. Default value is a0 = 1.
#' @param a.1                 Parameter for the prior Beta distribution for the probability of bias. Default value is a1 = 1.
#'
#'
#' @param alpha.0             Lower bound of the uniform prior for the concentration parameter for the DPM,
#'                            default value is alpha.0 = 0.03.
#' @param alpha.1             Upper bound of the uniform prior for the concentration parameter for the DPM,
#'                            default value is alpha.1 = 10.
#'
#' @param K                   Maximum number of clusters in the DPM, default value is K = 30.
#'
#'
#' @param nr.chains           Number of chains for the MCMC computations, default 2.
#' @param nr.iterations       Number of iterations after adapting the MCMC, default is 10000. Some models may need more iterations.
#' @param nr.adapt            Number of iterations in the adaptation process, default is 1000. Some models may need more iterations during adptation.
#' @param nr.burnin           Number of iteration discard for burn-in period, default is 1000. Some models may need a longer burnin period.
#' @param nr.thin             Thinning rate, it must be a positive integer, the default value 1.
#' @param be.quiet            Do not print warning message if the model does not adapt. The default value is FALSE. If you are not sure about the adaptation period choose be.quiet=TRUE.
#'
#' @return                    This function returns an object of the class "bcdpmeta". This object contains the MCMC
#'                            output of each parameter and hyper-parameter in the model and
#'                            the data frame used for fitting the model.
#'
#'
#' @details The results of the object of the class bcdpmeta can be extracted with R2jags or with rjags. In addition a summary, a print and a plot functions are
#' implemented for this type of object.
#'
#'
#' @references Verde, P.E. (2021) A Bias-Corrected Meta-Analysis Model for Combining Studies of Different Types and Quality. Biometrical Journal; 1–17.
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
#' bm1 = bcdpmeta(stemcells)
#' summary(bm1)
#'
#' }
#'
#' @import R2jags
#' @import rjags
#'
#'
#' @export
bcdpmeta = function(
    data,
    # Hyperpriors parameters............................................
    mean.mu.0     = 0,
    sd.mu.0       = 10,
    scale.sigma.between = 0.5,
    df.scale.between    = 1,
    B.lower             = 0,
    B.upper             = 10,
    a.0                 = 1,
    a.1                 = 1,
    alpha.0             = 0.03,
    alpha.1             = 10,
    K                   = 30,
    # MCMC setup..........................................................
    nr.chains       = 2,
    nr.iterations   = 10000,
    nr.adapt        = 1000,
    nr.burnin       = 1000,
    nr.thin         = 1,
    # Further options to link jags and R ...............................
    be.quiet        = FALSE
)UseMethod("bcdpmeta")


#' @export
bcdpmeta.default = function(
    data,
    # Hyperpriors parameters............................................
    mean.mu.0           = 0,
    sd.mu.0             = 10,
    scale.sigma.between = 0.5,
    df.scale.between    = 1,
    B.lower             = 0,
    B.upper             = 10,
    a.0                 = 1,
    a.1                 = 1,
    alpha.0             = 0.03,
    alpha.1             = 10,
    K                   = 5,
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

  # Load module mix
  load.module("mix")

  # Data sort to facilitate the mixture setup
  y = sort(data$TE)
  se.y = data$seTE[order(data$TE)]
  N = length(y)

  if(N<5)stop("Low number of studies in the meta-analysis!")

  # Note: prepare a regression equation ...
  #x = data$design[order(data$TE)]
  #x = as.numeric(x)


  # T is the index that allocates each study to one of the two components
  T = rep(NA, N)
  T[1] = 1
  T[N] = 2

  # This list describes the data used by the BUGS script.
  data.bcdpmeta = list ("y", "se.y", "N", "T",
                        "mean.mu.0",
                        "sd.mu.0",
                        "scale.sigma.between",
                        "df.scale.between",
                        "B.lower",
                        "B.upper",
                        "a.0",
                        "a.1",
                        "alpha.0",
                        "alpha.1",
                        "K")

  # List of parameters
  par.bcdpmeta  <- c("mu.k",
                     #"tau.k",
                     "mu.0",
                     "sd.0",
                    # "p",
                     "alpha.I.0",
                     "alpha.I.1",
                     "p.bias",
                     "B")

  # Model in BUGS
  model.bugs.bcdp =
    "
  model
  {
   for( i in 1 : N ) {
  # Likelihood of theta.B[i] ............................
            y[i] ~ dnorm(theta.B[i], pre.y[i])
        pre.y[i] <- pow(se.y[i], -2)

  # Dirichlet Process Random effects ..................
          theta.B[i] <- mu.k[T[i], group[i]]     # Locations combine mixture and DP
                T[i] ~ dcat(p.bias[1:2])
            group[i] ~ dcat(pi[T[i], ])
           }

# Stick-Breaking process and sampling from G0 for I=0 i.e. T=1...............
             q.I.0[1] ~ dbeta(1, alpha.I.0)
             p.I.0[1] <- q.I.0[1]
              pi[1, 1] <- p.I.0[1]

    for (j in 2:K) {
             q.I.0[j] ~ dbeta(1, alpha.I.0)
             p.I.0[j] <- q.I.0[j]*(1 - q.I.0[j-1])*p.I.0[j-1]/q.I.0[j-1]
              pi[1, j] <- p.I.0[j]
           }

# Stick-Breaking process and sampling from G0 for I=1 i.e. T=2...............
             q.I.1[1] ~ dbeta(1, alpha.I.1)
             p.I.1[1] <- q.I.1[1]
              pi[2, 1] <- p.I.1[1]

    for (j in 2:K) {
             q.I.1[j] ~ dbeta(1, alpha.I.1)
             p.I.1[j] <- q.I.1[j]*(1 - q.I.1[j-1])*p.I.1[j-1]/q.I.1[j-1]
             pi[2, j] <- p.I.1[j]
           }


  # Priors for the clusters' parameters ......................

  for(k in 1:K){
       mu.k[1, k] ~ dnorm(mu.mix[1], inv.var.0)
       mu.k[2, k] ~ dnorm(mu.mix[2], inv.var.0)
     }

  # Parameters for dnormmix ..................................

  mu.mix[1] <- mu.0
  mu.mix[2] <- mu.mix[1] + B
           B ~ dunif(B.lower, B.upper) # Bias paramter

   p.bias[2] ~ dbeta(a.0, a.1)         # Probability of bias
   p.bias[1] <- 1 - p.bias[2]

  inv.var.mix[1] <- inv.var.0
  inv.var.mix[2] <- inv.var.0

  # Prior for mu.0 ..........................................
       mu.0 ~ dnorm(mean.mu.0, inv.var.mu.0)
       inv.var.mu.0 <- pow(sd.mu.0, -2)

  # Prior for inv.var.0 .....................................
  inv.var.0 ~ dscaled.gamma(scale.sigma.between,
                               df.scale.between)

  sd.0 <- pow(inv.var.0, -0.5)

  # Prior for alpha .........................................
  alpha.I.0 ~ dunif(alpha.0, alpha.1)
  alpha.I.1 ~ dunif(alpha.0, alpha.1)
  }
  "


  model.bugs.connection = textConnection(model.bugs.bcdp)

  # Use R2jags as interface for JAGS ...

  results <- jags( data = data.bcdpmeta,
                   parameters.to.save = par.bcdpmeta,
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
  results$prior$mean.mu             = mean.mu.0
  results$prior$sd.mu               = sd.mu.0
  results$prior$scale.sigma.between = scale.sigma.between
  results$prior$df.scale.between    = df.scale.between
  results$prior$B.lower             = B.lower
  results$prior$B.upper             = B.upper
  results$prior$a.0                 = a.0
  results$prior$a.1                 = a.1
  results$prior$alpha.0             = alpha.0
  results$prior$alpha.1             = alpha.1
  results$prior$K                   = K
  results$N

  class(results) = c("bcdpmeta")

  return(results)
}


#' Generic print function for bcdpmeta object in jarbes.
#'
#' @param x The object generated by the function dpmmeta.
#'
#' @param digits The number of significant digits printed. The default value is 3.
#'
#' @param ... \dots
#'
#' @export
print.bcdpmeta <- function(x, digits, ...)
{
  print(x$BUGSoutput,...)
}




