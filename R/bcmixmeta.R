#' @title Bias Corrected Meta-Analysis with Dirichlet Mixture Process Priors for the biased component
#'
#' @description This function performers a Bayesian meta-analysis with DPM as random effects
#'
#'
#'
#' @param data                A data frame with at least two columns with the following names:
#'                             1) TE = treatment effect,
#'                             2) seTE = the standard error of the treatment effect.
#'
#' @param mean.mu.0            Prior mean of the mean of the base distribution default value is mean.mu.0 = 0.
#'
#' @param sd.mu.0              Prior standard deviation of the base distribution, the default value is 10^-6.
#'
#' @param scale.sigma.between Prior scale parameter for scale gamma distribution for the
#'                            precision between studies. The default value is 0.5.
#'
#' @param df.scale.between    Degrees of freedom of the scale gamma distribution for the precision between studies.
#'                            The default value is 1, which results in a Half Cauchy distribution for the standard
#'                            deviation between studies. Larger values e.g. 30 corresponds to a Half Normal distribution.
#' @param scale.sigma.beta    Prior scale parameter for the scale.gamma distribution for the
#'                            precision between study biases.
#'
#' @param df.scale.beta      Degrees of freedom of the scale gamma distribution for the precision between
#'                            study biases. The default value is 1, which results in a Half Cauchy distribution
#'                            for the standard deviation between biases.
#'
#' @param B.lower             Lower bound of the bias parameter B, the default value is -15.
#' @param B.upper             Upper bound of the bias parameter B, the default value is 15.
#'
#' @param a.0                 Parameter for the prior Beta distribution for the probability of bias. Default value is a0 = 0.5.
#' @param a.1                 Parameter for the prior Beta distribution for the probability of bias. Default value is a1 = 1.
#'
#'
#' @param alpha.0             Lower bound of the uniform prior for the concentration parameter for the DP,
#'                            the default value is 0.5.
#' @param alpha.1             Upper bound of the uniform prior for the concentration parameter for the DP,
#'                            the default value depends on the sample size, see the example below. We give as
#'                            working value alpha.1 = 2
#'
#' @param K                   Maximum number of clusters in the DP, the default value depends on alpha.1, see the
#'                            example below. We give as working vaule K = 10.
#'
#'
#' @param sort.priors         Experimental option, indicates if a weak information regarding the means and the variances are used.
#'                            If sort.priors == TRUE then the Delta parameter is not used
#'                            and only the order of the means and variances are restricted.
#'
#'
#' @param nr.chains           Number of chains for the MCMC computations, default 2.
#' @param nr.iterations       Number of iterations after adapting the MCMC, default is 10000. Some models may need more iterations.
#' @param nr.adapt            Number of iterations in the adaptation process, default is 1000. Some models may need more iterations during adptation.
#' @param nr.burnin           Number of iteration discard for burn-in period, default is 1000. Some models may need a longer burnin period.
#' @param nr.thin             Thinning rate, it must be a positive integer, the default value 1.
#' @param be.quiet            Do not print warning message if the model does not adapt. The default value is FALSE. If you are not sure about the adaptation period choose be.quiet=TRUE.
#'
#' @return                    This function returns an object of the class "bcmixmeta". This object contains the MCMC
#'                            output of each parameter and hyper-parameter in the model and
#'                            the data frame used for fitting the model.
#'
#'
#' @details The results of the object of the class bcmixmeta can be extracted with R2jags or with rjags. In addition a summary, a print and a plot functions are
#' implemented for this type of object.
#'
#'
#' @references Verde, P.E. and Rosner, G. L. (2024) A Bias-Corrected Bayesian Nonparamteric Model for Combining
#'              Studies with Varying Quality in Meta-Analysis. Biometrical Journal; (under revision).
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
#'  N = dim(stemcells)[1]
#'  alpha.max = 1/5 *((N-1)*a.0 - a.1)/(a.0 + a.1)
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
#' bcmix.2.stemcell = bcmixmeta(stemcells,
#'                             mean.mu.0=0, sd.mu.0=100,
#'                             B.lower = -15,
#'                             B.upper = 15,
#'                             alpha.0 = 0.5,
#'                             alpha.1 = alpha.max,
#'                             a.0 = a.0,
#'                             a.1 = a.1,
#'                             K = K.max,
#'                             sort.priors = FALSE,
#'                             df.scale.between = 1,
#'                             scale.sigma.between = 0.5,
#'                             nr.chains = 4,
#'                             nr.iterations = 50000,
#'                             nr.adapt = 1000,
#'                             nr.burnin = 10000,
#'                             nr.thin = 4)
#'
#'
#'  diagnostic(bcmix.2.stemcell, y.lim = c(-1, 15), title.plot = "Default priors")
#'
#'
#'  bcmix.2.stemcell.mcmc <- as.mcmc(bcmix.1.stemcell$BUGSoutput$sims.matrix)
#'
#'
#'theta.names <- paste(paste("theta[",1:31, sep=""),"]", sep="")
#'theta.b.names <- paste(paste("theta.bias[",1:31, sep=""),"]", sep="")
#'
# Greeks ...
#'theta.b.greek.names <- paste(paste("theta[",1:31, sep=""),"]^B", sep="")
#'theta.greek.names <- paste(paste("theta[",1:31, sep=""),"]", sep="")
#'
#'
#'caterplot(bcmix.2.stemcell.mcmc,
#'          parms = theta.names,               # theta
#'          labels = theta.greek.names,
#'          greek = T,
#'          labels.loc="axis", cex =0.7,
#'          col = "black",
#'          style = "plain",
#'          reorder = F,
#'          val.lim =c(-6, 16),
#'          quantiles = list(outer=c(0.05,0.95),inner=c(0.16,0.84)),
#'          x.lab = "Effect: mean difference"
#')
#'title( "95% posterior intervals of studies' effects")
#'caterplot(bcmix.2.stemcell.mcmc,
#'          parms = theta.b.names,             # theta.bias
#'          labels = theta.greek.names,
#'          greek = T,
#'          labels.loc="no",
#'          cex = 0.7,
#'          col = "grey",
#'          style = "plain", reorder = F,
#'          val.lim =c(-6, 16),
#'          quantiles = list(outer=c(0.025,0.975),inner=c(0.16,0.84)),
#'          add = TRUE,
#'          collapse=TRUE, cat.shift= -0.5,
#')
#'
#'attach.jags(bcmix.2.stemcell, overwrite = TRUE)
#'abline(v=mean(mu.0), lwd =2, lty =2)
#'
#'legend(9, 20, legend = c("bias corrected", "biased"),
#'      lty = c(1,1), lwd = c(2,2), col = c("black", "grey"))
#'
#'
#' }
#'
#' @import R2jags
#' @import rjags
#'
#'
#' @export
bcmixmeta = function(
    data,
    # Hyperpriors parameters............................................
    # Hyperpriors parameters............................................
    mean.mu.0     = 0,
    sd.mu.0       = 10,
    scale.sigma.between = 0.5,
    df.scale.between    = 1,
    scale.sigma.beta    = 0.5,
    df.scale.beta       = 1,
    B.lower             = -15,
    B.upper             = 15,
    a.0                 = 0.5,
    a.1                 = 1,
    alpha.0             = 0.03,
    alpha.1             = 2,
    K                   = 10,
    sort.priors         = FALSE,
    # MCMC setup..........................................................
    nr.chains       = 2,
    nr.iterations   = 10000,
    nr.adapt        = 1000,
    nr.burnin       = 1000,
    nr.thin         = 1,
    # Further options to link jags and R ...............................
    be.quiet        = FALSE
)UseMethod("bcmixmeta")


#' @export
bcmixmeta.default = function(
    data,
    # Hyperpriors parameters............................................
    mean.mu.0           = 0,
    sd.mu.0             = 10,
    scale.sigma.between = 0.5,
    df.scale.between    = 1,
    scale.sigma.beta    = 0.5,
    df.scale.beta       = 1,
    B.lower             = -15,
    B.upper             = 15,
    a.0                 = 0.5,
    a.1                 = 1,
    alpha.0             = 0.03,
    alpha.1             = 2,
    K                   = 10,
    sort.priors         = FALSE,
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


  # Data sort to facilitate the mixture setup
  # y = sort(data$TE)
  # se.y = data$seTE[order(data$TE)]

  # Avoid order of the effects
  y = data$TE
  se.y = data$seTE

  N = length(y)

  if(N<5)stop("Low number of studies in the meta-analysis!")

  # T is the index that allocates each study to one of the two components
  T = rep(NA, N)

  #T[1] = 1
  #T[N] = 2

  # Index of the max and min to inform T[], but this does not work for bilateral bias
  min.index = which(y == min(y))
  max.index = which(y == max(y))
  T[min.index] = 1
  T[max.index] = 2

  # This list describes the data used by the BUGS script.
  data.bcmixmeta.delta = list ("y", "se.y", "N", "T",
                               "mean.mu.0",
                               "sd.mu.0",
                               "scale.sigma.between",
                               "df.scale.between",
                               "scale.sigma.beta",
                               "df.scale.beta",
                               "B.lower",
                               "B.upper",
                               "a.0",
                               "a.1",
                               "alpha.0",
                               "alpha.1",
                               "K")

  # This list describes the data used by the BUGS script.
  data.bcmixmeta.sort = list ("y", "se.y", "N", "T",
                              "mean.mu.0",
                              "sd.mu.0",
                              "scale.sigma.between",
                              "df.scale.between",
                              "a.0",
                              "a.1",
                              "alpha.0",
                              "alpha.1",
                              "K",
                              "B.upper")


  if(sort.priors == TRUE)
  {data.bcmixmeta = data.bcmixmeta.sort}
  else
  {data.bcmixmeta = data.bcmixmeta.delta}



  # List of parameters
  par.bcmixmeta  <- c("mu.k",
                      "mu.0",
                      "mu.new",
                      "sd.0",
                      "p",
                      "alpha",
                      "p.bias",
                      "B",
                      "sd.beta",
                      "K.hat",
                      "theta.bias",
                      "theta",
                      "beta",
                       "I")

  # Model in BUGS without Delta, but using sort() function in JAGS ...
  model.bugs.bcmixmeta.sort =
    "
  model
  {
   for( i in 1 : N ) {
  # Likelihood of theta.B[i] ...................................................
            y[i] ~ dnorm(theta.bias[i], pre.y[i])
        pre.y[i] <- pow(se.y[i], -2)

  # Mixture Process: Normal and Dirichlet process ..............................

# Aca tengo que probar..........................................................................
#theta.bc[i] ~ dnorm(mu[i], inv.var.0)
#mu[i] <- mu.bias[T[i], group[i]]
#T[i] ~  dcat(p.bias[1:2])
# ......

#mu.bias[1] <- mu.0
#mu.bias[2,1:K] <- DP ... algo asi...
#..............................................................................................

# theta.b[i] <- theta[i]*(1-I[i]) + beta[i]*I[i]
# theta.bc[i] <- theta[i]*(1-I[i]) + (beta[i]+theta[i])*I[i]

 theta.bias[i] <- theta[i]+ (beta[i]*I[i])

              I[i] <- T[i] - 1
              T[i] ~  dcat(p.bias[1:2])
          theta[i] ~  dnorm(mu.0, inv.var.0)

  # Dirichlet Process for biased component ............................
           beta[i] <- mu.k[group[i]]
           group[i] ~ dcat(pi[])

  for(j in 1:K)
  {
      gind[i, j] <- equals(j, group[i])
      }
   }


  # Counting the number of clusters ...
  K.hat <- sum(cl[])
  for(j in 1:K)
  {
  sumind[j] <- sum(gind[,j])
      cl[j] <- step(sumind[j] -1+0.01) # cluster j used in this iteration.
  }

  # Prior for the probability of bias
     p.bias[2] ~ dbeta(a.0, a.1)
     p.bias[1] <- 1 - p.bias[2]


  # Prior for the means and variances to be ordered .....
  # Important: JAGS starts the MCMC of mu.s for mu.s[1] = m.s[2], which makes
  # no sence, because mu.0 < mu.beta for the biased case.
  #
  # In addition, the MCMC mixing takes longer if we start at mu.s[1] = m.s[2]
  #
  # Therefore, mu.s[2] has a prior N(mean.mu.0 + B/2, ...)

  inv.var.mu.0 <- pow(sd.mu.0, -2)


        mu.s[1] ~ dnorm(mean.mu.0, inv.var.mu.0)
   inv.var.s[1] ~ dscaled.gamma(scale.sigma.between, df.scale.between)

#mean.mu.bias <- mean.mu.0 + B/2
        mu.s[2] ~ dnorm(mean.mu.0 + B.upper/2, inv.var.mu.0)
   inv.var.s[2] ~ dscaled.gamma(scale.sigma.between, df.scale.between)

  # Sort means ....
  mu.sorted <- sort(mu.s)
       mu.0 <- mu.sorted[1]
    mu.beta <- mu.sorted[2]

  # Mean bias
   B <- mu.beta - mu.0

  # Sort variances ...
  inv.var.sorted <- sort(inv.var.s)
       inv.var.0 <- inv.var.sorted[2]
    inv.var.beta <- inv.var.sorted[1]
            sd.0 <- pow(inv.var.0, -0.5)
         sd.beta <- pow(inv.var.beta, -0.5)

# Posterior predictive
mu.new ~ dnorm(mu.0,  inv.var.0)


  # Priors for the clusters' parameters biased component .......................

  for(k in 1:K){
       mu.k[k] ~ dnorm(mu.beta, inv.var.beta)
       }


 # Prior for alpha .............................................................
  alpha ~ dunif(alpha.0, alpha.1)

 # Stick-Breaking process and sampling from G0..................................

             q[1] ~ dbeta(1, alpha)
             p[1] <- q[1]
            pi[1] <- p[1]

    for (j in 2:K) {
             q[j] ~ dbeta(1, alpha)
             p[j] <- q[j]*(1 - q[j-1])*p[j-1]/q[j-1]
            pi[j] <- p[j]
           }

  }
  "




  # Model in BUGS with Delta ....
  model.bugs.bcmixmeta.delta =
    "
  model
  {
   for( i in 1 : N ) {
  # Likelihood of theta.B[i] ...................................................
            y[i] ~ dnorm(theta.bias[i], pre.y[i])
        pre.y[i] <- pow(se.y[i], -2)

  # Mixture Process: Normal and Dirichlet process ..............................

 #theta.bias[i] <- theta[i]*(1-I[i]) + beta[i]*I[i]

  theta.bias[i] <- theta[i] + (beta[i]*I[i])

              I[i] <- T[i] - 1
              T[i] ~  dcat(p.bias[1:2])
          theta[i] ~  dnorm(mu.0, inv.var.0)

  # Dirichlet Process for biased component ............................
           beta[i] <- mu.k[group[i]]
           group[i] ~ dcat(pi[])


  for(j in 1:K)
  {
      gind[i, j] <- equals(j, group[i])
      }
   }

  # Counting the number of clusters ...
  K.hat <- sum(cl[])
  for(j in 1:K)
  {
  sumind[j] <- sum(gind[,j])
      cl[j] <- step(sumind[j] -1+0.01) # cluster j used in this iteration.
  }


  # Prior for the probability of bias
     p.bias[2] ~ dbeta(a.0, a.1)
     p.bias[1] <- 1 - p.bias[2]

  # Priors for model of interest ..............................................

          mu.0 ~ dnorm(mean.mu.0, inv.var.mu.0)
          inv.var.mu.0 <- pow(sd.mu.0, -2)

     inv.var.0 ~ dscaled.gamma(scale.sigma.between, df.scale.between)

         sd.0 <- pow(inv.var.0, -0.5)

# Posterior predictive
mu.new ~ dnorm(mu.0,  inv.var.0)


  # Priors for the clusters' parameters biased component .......................

         B ~ dunif(B.lower, B.upper) # Bias parameter
  mu.bias <- mu.0 + B                # Here mu.bias = mu.0 + Bias
  for(k in 1:K){
      # mu.k[k] ~ dnorm(mu.0 + B, inv.var.0 + inv.var.beta) #
        mu.k[k] ~ dnorm(B, inv.var.beta) #                    DP on the betas_i
        }

  inv.var.beta ~ dscaled.gamma(scale.sigma.beta, df.scale.beta)
  sd.beta = pow(inv.var.beta, -0.5)


 # Prior for alpha .............................................................
  alpha ~ dunif(alpha.0, alpha.1)

 # Stick-Breaking process and sampling from G0..................................

             q[1] ~ dbeta(1, alpha)
             p[1] <- q[1]
            pi[1] <- p[1]

           for (j in 2:K) {
             q[j] ~ dbeta(1, alpha)
             p[j] <- q[j]*(1 - q[j-1])*p[j-1]/q[j-1]
            pi[j] <- p[j]
           }

  }
  "

  if(sort.priors == TRUE)
  {model.bugs.connection = textConnection(model.bugs.bcmixmeta.sort)}
  else
  {model.bugs.connection = textConnection(model.bugs.bcmixmeta.delta)}



  # Use R2jags as interface for JAGS ...

  results <- jags( data = data.bcmixmeta,
                   parameters.to.save = par.bcmixmeta,
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

  class(results) = c("bcmixmeta")

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
print.bcmixmeta <- function(x, digits, ...)
{
  print(x$BUGSoutput,...)
}



