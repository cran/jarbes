#' @title Bayesian Meta-Analysis with Dirichlet Process Priors
#'
#' @description This function performers a Bayesian meta-analysis with DP as random effects
#'
#'
#'
#' @param data                A data frame with at least two columns with the following names:
#'                             1) TE = treatment effect,
#'                             2) seTE = the standard error of the treatment effect.
#'
#' @param x                   a covariate to perform meta-regression.
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
#'
#' @return                    This function returns an object of the class "dpmetareg". This object contains the MCMC
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
#'
#' # Example: Stemcells
#'
#' data("stemcells")
#' stemcells$TE = stemcells$effect.size
#' stemcells$seTE = stemcells$se.effect
#'
#' bm1 = dpmmeta(stemcells)
#' summary(bm1)
#' plot(bm1, x.lim = c(-1, 7), y.lim = c(0, 1))
#'
#' diagnostic(bm1, study.names = stemcells$trial,
#'            post.p.value.cut = 0.05,
#'            lwd.forest = 0.5, shape.forest = 4)
#'
#' diagnostic(bm1, post.p.value.cut = 0.05,
#'            lwd.forest = 0.5, shape.forest = 4)
#' }
#'
#' @import R2jags
#' @import rjags
#'
#'
#' @export
dpmetareg = function(
    data,
    x,
    # Hyperpriors parameters............................................
    # Hyperpriors parameters............................................
    mean.mu.0     = 0,
    sd.mu.0       = 10,
    scale.sigma.between = 0.5,
    df.scale.between    = 1,
    alpha.0             = 0.03,
    alpha.1             = 10,
    K                   = 30,
    # MCMC setup........................................................
    nr.chains       = 2,
    nr.iterations   = 10000,
    nr.adapt        = 1000,
    nr.burnin       = 1000,
    nr.thin         = 1)UseMethod("dpmetareg")


#' @export
dpmetareg.default = function(
    data,
    x,
    # Hyper-priors parameters............................................
    mean.mu.0     = 0,
    sd.mu.0       = 10,
    scale.sigma.between = 0.5,
    df.scale.between    = 1,
    alpha.0             = 0.03,
    alpha.1             = 10,
    K                   = 5,
    # MCMC setup........................................................
    nr.chains       = 2,
    nr.iterations   = 10000,
    nr.adapt        = 1000,
    nr.burnin       = 1000,
    nr.thin         = 1)
{

  # Data
  y = data$TE
  se.y = data$seTE
  N = length(y)
  x = x

  if(N<3)stop("Low number of studies in the meta-analysis!")
  # Approximate Bayesian Cross-Validation
  #y.ghost = rep(NA, N)


  # This list describes the data used by the BUGS script.
  data.dp = list ("y", "se.y", "x", "N",
                   "mean.mu.0",
                   "sd.mu.0",
                   "scale.sigma.between",
                   "df.scale.between",
                   "alpha.0",
                   "alpha.1",
                   "K")

  # List of parameters
  par.dp  <- c("mu.k",
               "beta.k",
                "mu.theta",
                "theta",
                "mu.0",
                "sd.0",
                "pi",
                "alpha",
                "gind")

  # Model in BUGS
  model.bugs <-
    "
  model
  {
   for( i in 1 : N )
   {
  # Likelihood of theta[i] ............................
            y[i] ~ dnorm(theta[i], pre.y[i])
        pre.y[i] <- pow(se.y[i], -2)

  # Dirichlet Process Random effects ..................
       #  theta[i] <- mu.k[group[i]] + beta.k[group[i]]*x[i]

          theta[i] ~ dnorm(mu.theta[i], inv.var.0)

       #   mu.theta[i] <- mu.k[group[i]] + beta.k*x[i] # It works well

           mu.theta[i] <- mu.k[group[i]] + beta.k[group[i]]*x[i]
           group[i] ~ dcat(pi[])



  #Frequency of study i belong to cluster j
  for(j in 1:K)
  {gind[i, j] <- equals(j, group[i])}
           }

 # Priors for the clusters' regression parameters ......................
  for(k in 1:K){
         mu.k[k] ~ dnorm(mu.0, inv.var.0)
       beta.k[k] ~ dnorm(mu.0, inv.var.0)
      # theta.mu[k] <- mu.k[k] + beta.k[k]*x[]
       }

  # Prior for mu.0
       mu.0 ~ dnorm(mean.mu.0, inv.var.mu.0)
       inv.var.mu.0 <- pow(sd.mu.0, -2)

     #  beta.k ~ dnorm(0, 0.001)

  # Prior for inv.var.0
  inv.var.0 ~ dscaled.gamma(scale.sigma.between,
                            df.scale.between)

  sd.0 <- pow(inv.var.0, -0.5)

  # Prior for alpha
  alpha ~ dunif(alpha.0, alpha.1)


  # Stick-Breaking process and sampling from G0...............
             q[1] ~ dbeta(1, alpha)
             p[1] <- q[1]
            pi[1] <- p[1]

    for (j in 2:K) {
             q[j] ~ dbeta(1, alpha)
             p[j] <- q[j]*(1 - q[j-1])*p[j-1]/q[j-1]
            pi[j] <- p[j]/sum(p[])                    # Make sure that pi[] adds to 1
    }
  }
  "

  model.bugs.connection <- textConnection(model.bugs)

  # Use R2jags as interface for JAGS ...

  results <- jags( data = data.dp,
                   parameters.to.save = par.dp,
                   model.file = model.bugs.connection,
                   n.chains = nr.chains,
                   n.iter = nr.iterations,
                   n.burnin = nr.burnin,
                   n.thin = nr.thin,
                   pD = TRUE)

  # Close text connection
  close(model.bugs.connection)

  # Extra outputs that are linked with other functions ...
  results$data = data

  # Hyperpriors parameters............................................
  results$prior$mean.mu     = mean.mu.0
  results$prior$sd.mu       = sd.mu.0
  results$prior$scale.sigma.between = scale.sigma.between
  results$prior$df.scale.between    = df.scale.between
  results$prior$alpha.0 = alpha.0
  results$prior$alpha.1 = alpha.1
  results$prior$K = K
  results$N

  class(results) = c("dpmetareg")

  return(results)
}



#' Generic print function for dpmeta object in jarbes.
#'
#' @param x The object generated by the function dpmmeta.
#'
#' @param digits The number of significant digits printed. The default value is 3.
#'
#' @param ... \dots
#'
#' @export
print.dpmetareg <- function(x, digits, ...)
{
  print(x$BUGSoutput,...)
}





