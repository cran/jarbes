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
#' @param x                   a covariate to perform meta-regression.
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
#'                            example below. We give as working value K = 10.
#'
#'
#' @param bilateral.bias      Experimental option, which indicates if bias could be to the left and
#'                            to the right of the model of interest. If bilateral.bias==TRUE,
#'                            then the function generates three mean and sorts the means in two
#'                            groups: mean_bias_left, mean_theta, mean_bias_right.
#'
#'
#'
#' @param nr.chains           Number of chains for the MCMC computations, default 2.
#' @param nr.iterations       Number of iterations after adapting the MCMC, default is 10000. Some models may need more iterations.
#' @param nr.adapt            Number of iterations in the adaptation process, default is 1000. Some models may need more iterations during adptation.
#' @param nr.burnin           Number of iteration discard for burn-in period, default is 1000. Some models may need a longer burnin period.
#' @param nr.thin             Thinning rate, it must be a positive integer, the default value 1.
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
    x = NULL,
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
    bilateral.bias      = FALSE,
    # MCMC setup..........................................................
    nr.chains       = 2,
    nr.iterations   = 10000,
    nr.adapt        = 1000,
    nr.burnin       = 1000,
    nr.thin         = 1)UseMethod("bcmixmeta")


#' @export
bcmixmeta.default = function(
    data,
    x = NULL,
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
    bilateral.bias      = FALSE,
    # MCMC setup........................................................
    nr.chains       = 2,
    nr.iterations   = 10000,
    nr.adapt        = 1000,
    nr.burnin       = 1000,
    nr.thin         = 1
)
{


  # Data sort to facilitate the mixture setup
  # y = sort(data$TE)
  # se.y = data$seTE[order(data$TE)]

  # Avoid order of the effects
  y = data$TE
  se.y = data$seTE
  x = x            # a single covariate for meta-regression
  N = length(y)

  if(N<5)stop("Low number of studies in the meta-analysis!")

  # T is the index that allocates each study to one of the two components
  T = rep(NA, N)

  #T[1] = 1
  #T[N] = 2

  # Index of the max and min to inform T[], for bilateral bias it is adapted.

  if(bilateral.bias==TRUE){
    min.index = which(y == min(y, na.rm = TRUE))
    max.index = which(y == max(y, na.rm = TRUE))
    T[min.index] = 1
    T[max.index] = 3
  }
  else
  {
    min.index = which(y == min(y, na.rm = TRUE))
    max.index = which(y == max(y, na.rm = TRUE))
    T[min.index] = 1
    T[max.index] = 2
  }

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

  data.bcmixmeta.x = list ("y", "se.y", "x",
                           "N", "T",
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
  data.bcmixmeta.bilateral.bias = list ("y", "se.y", "N", "T",
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


  if(bilateral.bias == TRUE)
  {data.bcmixmeta = data.bcmixmeta.bilateral.bias}
  else
  {data.bcmixmeta = data.bcmixmeta.delta}

  if(!is.null(x)) data.bcmixmeta = data.bcmixmeta.x

  # List of parameters
  par.bcmixmeta  <- c("mu.k",
                      "mu.0",
                      "mu.new",
                      "sd.0",
                      "p.cluster",
                      "alpha",
                      "p.bias",
                      "B",
                      "sd.beta",
                      "K.hat",
                      "theta.bias",
                      "theta",
                      "beta",
                       "I",
                      "group",
                      "gind",
                      "equalsmatrix.bias.1",
                      "equalsmatrix.bias.2",
                      "new.group")
  # List of parameters
  par.bcmixmeta.x  <- c("mu.k",
                      "mu.0",
                      "mu.new",
                      "sd.0",
                      "p.cluster",
                      "alpha",
                      "p.bias",
                      "B",
                      "sd.beta",
                      "K.hat",
                      "theta.bias",
                      "theta",
                      "beta",
                      "I",
                      "group",
                      "gind",
                      "equalsmatrix",
                      "beta.x",
                      "mu.x",
                      "mu.x.pred",
                      "alpha.0.bias",
                      "alpha.1.bias")

  # List of parameters
  par.bcmixmeta.bilateral  <- c(
                      "beta.k",
                      "mu.0",
                      "mu.beta.left",
                      "mu.beta.right",
                      "mu.new",
                      "sd.0",
                      "sd.beta.left",
                      "sd.beta.right",
                      "p.cluster",
                      "alpha",
                      "p.bias",
                      "B.left",
                      "B.right",
                      "K.hat",
                      "theta.bias",
                      "theta",
                      "beta",
                      "I",
                      "gind")

  if(bilateral.bias == TRUE)
  {par.bcmixmeta = par.bcmixmeta.bilateral }
  else
  {par.bcmixmeta = par.bcmixmeta}

  if(!is.null(x)) par.bcmixmeta = par.bcmixmeta.x




  # Model in BUGS with mu_beta ....
  model.bugs.bcmixmeta.delta =
    "
  model
  {
   for( i in 1 : N ) {
  # Likelihood of theta.B[i] ...................................................
            y[i] ~ dnorm(theta.bias[i], pre.y[i])
        pre.y[i] <- pow(se.y[i], -2)

  # Mixture Process: Normal and Dirichlet process ..............................

  theta.bias[i] <- theta[i] + (beta[i]*I[i])

              I[i] <- T[i] - 1
              T[i] ~  dcat(p.bias[1:2])
          theta[i] ~  dnorm(mu.0, inv.var.0)

  # Dirichlet Process for biased component ............................
           beta[i] <- mu.k[group[i]]
           group[i] ~ dcat(p.cluster[])

  for(j in 1:K)
  {gind[i, j] <- equals(j, group[i])}   #Frequency of study i belong to cluster j

   }

  # Number of studies in the same bias cluster
  for(i in 1:N)
  {new.group[i] <- (1-I[i]) + I[i]*(group[i]+1)}

  for(i in 1:N){
    for(j in 1:N){
    equalsmatrix.bias.1[i,j] <- equals(I[i], I[j])
    equalsmatrix.bias.2[i,j] <- equals(new.group[i], new.group[j])
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

  #q.beta <- 0.01  # quality constant to increase the variability when I=1

  for(k in 1:K){
        mu.k[k] ~ dnorm(B, inv.var.beta) # DP on the betas_i
    #   mu.k[k] ~ dnorm(B, inv.var.0*q.beta)
        }

  inv.var.beta ~ dscaled.gamma(scale.sigma.beta, df.scale.beta)
  sd.beta = pow(inv.var.beta, -0.5)


 # Prior for alpha .............................................................
  alpha ~ dunif(alpha.0, alpha.1)

  #alpha <- alpha.1

 # Stick-Breaking process and sampling from G0..................................

             q[1] ~ dbeta(1, alpha)
             p[1] <- q[1]
     p.cluster[1] <- p[1]

           for (j in 2:(K-1)) {
             q[j] ~ dbeta(1, alpha)
             p[j] <- q[j]*(1 - q[j-1])*p[j-1]/q[j-1]
     p.cluster[j] <- p[j]/sum(p[])                    # Make sure that pi[] adds to 1
           }


  }
  "

# Model in BUGS with mu_beta ....
  model.bugs.bcmixmeta.x =
    "
  model
  {
   for( i in 1 : N ) {
  # Likelihood of theta.B[i] ...................................................
            y[i] ~ dnorm(theta.bias[i], pre.y[i])
        pre.y[i] <- pow(se.y[i], -2)

  # Mixture Process: Normal and Dirichlet process ..............................

  theta.bias[i] <- theta[i] + (beta[i]*I[i])

              I[i] <- T[i] - 1
              T[i] ~  dcat(p.bias[i,1:2])

              logit(p.bias[i,2]) = alpha.0.bias + alpha.1.bias*x[i]
              p.bias[i,1] = 1-p.bias[i,2]

              theta[i] ~  dnorm(mu.0, inv.var.0)

       #   theta[i] ~  dnorm(mu.x[i], inv.var.0)
       #   mu.x[i] <- mu.0 + beta.x*x[i]
       #  mu.x.pred[i] ~  dnorm(mu.x[i], inv.var.0)

  # Dirichlet Process for biased component ............................
           beta[i] <- mu.k[group[i]]
           group[i] ~ dcat(p.cluster[])

 for(j in 1:K) {gind[i, j] <- equals(j, group[i])}   #Frequency of study i belong to cluster j
   }

  # Number of studies in the same bias cluster
  for(i in 1:N){
    for(j in 1:N){
    equalsmatrix[i,j] <- equals(group[i], group[j])
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
   #  p.bias[2] ~ dbeta(a.0, a.1)
   #  p.bias[1] <- 1 - p.bias[2]
   alpha.0.bias ~ dnorm(0, 0.01)
   alpha.1.bias ~ dnorm(0, 0.01)

  # Priors for model of interest ..............................................

          mu.0 ~ dnorm(mean.mu.0, inv.var.mu.0)
          inv.var.mu.0 <- pow(sd.mu.0, -2)

          beta.x ~ dnorm(mean.mu.0, inv.var.mu.0)

     inv.var.0 ~ dscaled.gamma(scale.sigma.between, df.scale.between)

         sd.0 <- pow(inv.var.0, -0.5)

# Posterior predictive
mu.new ~ dnorm(mu.0,  inv.var.0)

# Priors for the clusters' parameters biased component .......................

         B ~ dunif(B.lower, B.upper) # Bias parameter
  mu.bias <- mu.0 + B                # Here mu.bias = mu.0 + Bias
  for(k in 1:K){
        mu.k[k] ~ dnorm(B, inv.var.beta) #                    DP on the betas_i
        }

  inv.var.beta ~ dscaled.gamma(scale.sigma.beta, df.scale.beta)
  sd.beta = pow(inv.var.beta, -0.5)


 # Prior for alpha .............................................................
  alpha ~ dunif(alpha.0, alpha.1)

  #alpha <- alpha.1

 # Stick-Breaking process and sampling from G0..................................

             q[1] ~ dbeta(1, alpha)
             p[1] <- q[1]
            p.cluster[1] <- p[1]

           for (j in 2:K) {
             q[j] ~ dbeta(1, alpha)
             p[j] <- q[j]*(1 - q[j-1])*p[j-1]/q[j-1]
           p.cluster[j] <- p[j]/sum(p[])              # Make sure that pi[] adds to 1
           }

  }
  "





  # Model in BUGS using bilateral bias ...
  model.bugs.bcmixmeta.bilateral.bias =
    "
  model
  {
   for( i in 1 : N ) {
  # Likelihood of theta.bias[i] ................................................
            y[i] ~ dnorm(theta.bias[i], pre.y[i])
        pre.y[i] <- pow(se.y[i], -2)

# Three groups mixtures: Normal and left/right Dirichlet process................

     theta.bias[i] <- theta[i]+ (beta[i]*I[i])
              I[i] <- T[i] - 2                 # I is -1, 0, 1
              T[i] ~  dcat(p.bias[1:3])
          theta[i] ~  dnorm(mu.0, inv.var.0)

  # Dirichlet Process for biased component .....................................
           beta[i] <- beta.k[T[i], group[i]] #
          group[i] ~  dcat(pi[T[i], ])       # group[i] is the index of the cluster within the process
                                             # pi[1, i] is the probability of the cluster i left
                                             # pi[2, i] = 0
                                             # pi[3, i] is the probability of the cluster i right

  for(j in 1:K)
  {
      gind[i, j] <- equals(j, group[i])
      }
   }


  # Stick-Breaking process on the left, i.e.  T= 1...............
              q.I.1[1] ~ dbeta(1, alpha)
             p.I.1[1] <- q.I.1[1]
             pi[1, 1] <- p.I.1[1]

    for (j in 2:K) {
             q.I.1[j] ~ dbeta(1, alpha)
             p.I.1[j] <- q.I.1[j]*(1 - q.I.1[j-1])*p.I.1[j-1]/q.I.1[j-1]
             pi[1, j] <- p.I.1[j]
    }

  # Stick-Breaking process on the right, i.e. T = 3...............
             q.I.3[1] ~ dbeta(1, alpha)
             p.I.3[1] <- q.I.3[1]
             pi[3, 1] <- p.I.3[1]

    for (j in 2:K) {
             q.I.3[j] ~ dbeta(1, alpha)
             p.I.3[j] <- q.I.3[j]*(1 - q.I.3[j-1])*p.I.3[j-1]/q.I.3[j-1]
             pi[3, j] <- p.I.3[j]
           }

  # For the model of interest, i.e. T = 2, pi[2, i] = 0
    for(i in 1:K)
      {
        pi[2, i] <- 0.01
      }

  # Counting the number of clusters ...
  K.hat <- sum(cl[])
  for(j in 1:K)
  {
  sumind[j] <- sum(gind[,j])
      cl[j] <- step(sumind[j]-1+0.01) # cluster j used in this iteration.
  }

  # Prior for the probability of bias
     p.bias[1] ~ dbeta(a.0, a.1)                  # Left bias
     p.bias[2] <- 1 - (p.bias[1]+p.bias[3])
     p.bias[3] ~ dbeta(a.0, a.1)                  # Right bias

  # Important: JAGS starts the MCMC of mu.s for mu.s[1] = m.s[2] = m.s[3], which makes
  # no sence for the biased case.
  #
  # In addition, the MCMC mixing takes longer if we start at mu.s[1] = m.s[2] = m.s[3]
  #
  # Therefore, mu.s[1] has a prior N(mean.mu.0 + B.lower, ...)
  #            mu.s[3] has a prior N(mean.mu.0 + B.upper, ...)

  inv.var.mu.0 <- pow(sd.mu.0, -2)

# mean.mu.bias <- mean.mu.0 + B.lower
        mu.s[1] ~ dnorm(mean.mu.0 + B.lower, inv.var.mu.0)
   inv.var.beta.left ~ dscaled.gamma(scale.sigma.between, df.scale.between)

        mu.s[2] ~ dnorm(mean.mu.0, inv.var.mu.0)
      inv.var.0 ~ dscaled.gamma(scale.sigma.between, df.scale.between)

# mean.mu.bias <- mean.mu.0 + B.upper
        mu.s[3] ~ dnorm(mean.mu.0 + B.upper, inv.var.mu.0)
   inv.var.beta.right ~ dscaled.gamma(scale.sigma.between, df.scale.between)

  # Sort means ....
  mu.sorted <- mu.s   #sort(mu.s)
     mu.beta.left <- mu.sorted[1]
             mu.0 <- mu.sorted[2]
    mu.beta.right <- mu.sorted[3]

  # sigmas ....
      sd.0 <- pow(inv.var.0, -0.5)
   sd.left <- pow(inv.var.beta.left, -0.5)
  sd.right <- pow(inv.var.beta.right, -0.5)

  # Mean bias ...
    B.left <- mu.beta.left - mu.0
   B.right <- mu.beta.right - mu.0

# Priors for the clusters' parameters biased component .......................
  for(k in 1:K){
      beta.k[1, k] ~ dnorm(mu.beta.left, inv.var.beta.left)
      beta.k[2, k] <- 0
      beta.k[3, k] ~ dnorm(mu.beta.right, inv.var.beta.right)
       }

# Prior for alpha .............................................................
  alpha ~ dunif(alpha.0, alpha.1)

# Posterior predictive
  mu.new ~ dnorm(mu.0,  inv.var.0)

  }
  "


  if(bilateral.bias == TRUE)
  {model.bugs.connection = textConnection(model.bugs.bcmixmeta.bilateral.bias)}
  else
  {model.bugs.connection = textConnection(model.bugs.bcmixmeta.delta)}


  if(!is.null(x)) {model.bugs.connection = textConnection(model.bugs.bcmixmeta.x)}


  # Use R2jags as interface for JAGS ...

  results <- jags( data = data.bcmixmeta,
                   parameters.to.save = par.bcmixmeta,
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


#' Generic print function for bcmixmeta object in jarbes.
#'
#' @param x The object generated by the function bcmixmeta.
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



