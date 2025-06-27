#' @title Bias-Corrected Meta-Analysis for Combining Studies of Different Types and Quality
#'
#' @description This function performers a Bayesian meta-analysis to jointly
#' combine different types of studies. The random-effects follows a finite
#' mixture of normal distributions.
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
#' @param B.lower             Lower bound of the bias parameter B, the default value is 0.
#' @param B.upper             Upper bound of the bias parameter B, the default value is 10.
#'
#' @param a.0                 Parameter for the prior Beta distribution for the probability of bias. Default value is a0 = 1.
#' @param a.1                 Parameter for the prior Beta distribution for the probability of bias. Default value is a1 = 1.
#'
#' @param nu                  Parameter for the Beta distribution for the quality weights. The default value is nu = 0.5.
#'
#' @param nu.estimate         If TRUE, then we estimate nu from the data.
#'
#' @param b.0                 If nu.estimate = TRUE, this parameter is the shape parameter of the prior Gamma distribution for nu.
#' @param b.1                 If nu.estimate = TRUE, this parameter is the rate parameter of the prior Gamma distribution for nu.
#'                            Note that E(nu) = b.0/b.1 and we need to choose b.0 << b.1.
#'
#' @param nr.chains           Number of chains for the MCMC computations, default 2.
#' @param nr.iterations       Number of iterations after adapting the MCMC, default is 10000. Some models may need more iterations.
#' @param nr.adapt            Number of iterations in the adaptation process, defualt is 1000. Some models may need more iterations during adptation.
#' @param nr.burnin           Number of iteration discared for burnin period, default is 1000. Some models may need a longer burnin period.
#' @param nr.thin             Thinning rate, it must be a positive integer, the default value 1.
#' @param parallel            NULL -> jags, 'jags.parallel' -> jags.parallel execution
#'
#' @return                    This function returns an object of the class "bcmeta". This object contains the MCMC
#'                            output of each parameter and hyper-parameter in the model and
#'                            the data frame used for fitting the model.
#'
#'
#' @details The results of the object of the class bcmeta can be extracted with R2jags or with rjags. In addition a summary, a print and a plot functions are
#' implemented for this type of object.
#'
#'
#' @references Verde, P. E. (2017) Two Examples of Bayesian Evidence Synthesis with the Hierarchical Meta-Regression Approach. Chap.9, pag 189-206. Bayesian Inference, ed. Tejedor, Javier Prieto. InTech.
#'
#' @references Verde, P.E. (2021) A Bias-Corrected Meta-Analysis Model for Combining Studies of Different Types and Quality. Biometrical Journal; 1–17.
#'
#'
#' @examples
#' \dontrun{
#' library(jarbes)
#'
#' # Example ppvipd data
#'
#' data(ppvipd)
#'
# bcm1 = bcmeta(ppvipd, a.0 = 8, a.1 = 3) # 8 OS and 3 RCT...
# summary(bcm1)
# plot(bcm1)
#
# diagnostic(bcm1, study.names = ppvipd$name,
#            post.p.value.cut = 0.1, shape.forest = 4, lwd.forest = 0.75)
#
#
# diagnostic(bcm1, post.p.value.cut = 0.1,
#            shape.forest = 4, lwd.forest = 0.75,
#            y.lim = c(0, 3),
#            color.data.points = "blue",
#            bias.plot = TRUE,
#            S = 10000)
#
# # Example with stem cells
#
# data("stemcells")
# stemcells$TE = stemcells$effect.size
# stemcells$seTE = stemcells$se.effect
#
# bcm2 = bcmeta(stemcells, mean.mu = 0, sd.mu = 1000,
#               nr.iterations = 50000,
#               nr.adapt = 20000,
#               nr.thin = 2)
#
# summary(bcm2)
#
# plot(bcm2, x.lim = c(-5, 15), y.lim = c(0, .8))
#
# diagnostic(bcm2, study.names = stemcells$trial,
#            post.p.value.cut = 0.05,
#            shape.forest = 4, lwd.forest = 0.75)
#
# # only bias plot
# diagnostic(bcm2, study.names = stemcells$trial,
#            post.p.value.cut = 0.05,
#            cross.val.plot = FALSE,
#            shape.forest = 4, lwd.forest = 0.75)
#
# # only CV plot
# diagnostic(bcm2, study.names = stemcells$trial,
#            post.p.value.cut = 0.05,
#            cross.val.plot = TRUE,
#            bias.plot = FALSE,
#            shape.forest = 4, lwd.forest = 0.75)
#
#'
#' }
#'

#' @import R2jags
#' @import rjags
#'
#'
#' @export
bcmeta = function(
  data,
  # Hyperpriors parameters............................................
  mean.mu     = 0,
  sd.mu       = 10,
  scale.sigma.between = 0.5,
  df.scale.between    = 1,
  B.lower     = 0,
  B.upper     = 10,

  a.0         = 1,
  a.1         = 1,
  nu          = 0.5,
  nu.estimate = FALSE,
  b.0 = 1,
  b.1 = 2,

  # MCMC setup........................................................
  nr.chains       = 2,
  nr.iterations   = 10000,
  nr.adapt        = 1000,
  nr.burnin       = 1000,
  nr.thin         = 1,
  parallel        = NULL
           )UseMethod("bcmeta")



#' @export
bcmeta.default = function(
  data,
  # Hyperpriors parameters............................................
  mean.mu     = 0,
  sd.mu       = 10,
  scale.sigma.between = 0.5,
  df.scale.between    = 1,
  B.lower     = 0,
  B.upper     = 10,
  a.0         = 1,
  a.1         = 1,
  nu          = 0.5,
  nu.estimate = FALSE,
  b.0         = 1,
  b.1         = 2,


  # MCMC setup........................................................
  nr.chains       = 2,
  nr.iterations   = 10000,
  nr.adapt        = 1000,
  nr.burnin       = 1000,
  nr.thin         = 1,
  parallel        = NULL)
{
  if(!is.null(parallel) && parallel != "jags.parallel") stop("The parallel option must be NULL or 'jags.parallel'")
  # Mixture of Normals Random effects meta-analysis
  #    y = sort(data$TE, na.last = TRUE)
  # se.y = data$seTE[order(data$TE, na.last = TRUE)]
  #    N = length(y)

     # Avoid order of the effects
     y = data$TE
     se.y = data$seTE

     N = length(y)

     if(N<5)stop("Low number of studies in the meta-analysis!")

     # Note: prepare a regression equation ...
     #x = data$design[order(data$TE)]
     #x = as.numeric(x)


     #T is the index that allocates each study to one of the two components
     T = rep(NA, N)
     min.index = which(y == min(y, na.rm = TRUE))
     max.index = which(y == max(y, na.rm = TRUE))
     T[min.index] = 1
     T[max.index] = 2

     # This list describes the data used by the BUGS script.
     data.bcmeta <-
       list(y = y,
            se.y = se.y,
            N = N,
            T = T,
            mean.mu = mean.mu,
            sd.mu = sd.mu,
            scale.sigma.between = scale.sigma.between,
            df.scale.between = df.scale.between,
            B.lower = B.lower,
            B.upper = B.upper,
            a.0 = a.0,
            a.1 = a.1,
            nu = nu
       )


     if(nu.estimate==TRUE) {
       data.bcmeta = data.bcmeta[-13]              # Remove "nu"
       data.bcmeta = c(data.bcmeta, "b.0", "b.1")} # Add hyper constants for "nu

     # List of parameters
     par.bcmeta  <- c("theta.bc",
                      "theta",
                      "theta.bias",
                      "mu",
                      "mu.new",  # BC component
                      "y.ghost", # Approx Cross-Validation
                      "tau",
                      "p.bias",
                      "I",
                      "w",
                      "w.diag", # Diagnostic weights
                      "B")

     if(nu.estimate==TRUE) {par.bcmeta = c(par.bcmeta, "nu" )}


  # Model in BUGS

  model.bugs <-
  "
  model
  {
   for( i in 1 : N ) {
  # Likelihood of theta[i] ..........................................

            y[i] ~ dnorm(theta.bc[i], pre.y[i])
       pre.y[i] <- pow(se.y[i], -2)

  # Mixture heterocedastic random effects ..........................

      #theta.bc[i] <- theta[i]*(1-I[i]) + theta.bias[i]*I[i]

      theta.bc[i] <- theta[i]+ (theta.bias[i]*I[i])


             I[i] <- T[i] - 1
          theta[i] ~ dnorm(mu[1], prec.tau[i])
     theta.bias[i] ~ dnorm(mu[2], prec.tau[i])

            T[i] ~ dcat(p.bias[1:2])
     prec.tau[i] <- inv.var[T[i]] * w[T[i],i]  #Slash parametrization
          w[1,i] <- 1

   # Slash distribution ...
      w[2,i] ~ dbeta(nu, 1)

   # Diagnostic weights
   w.diag[i] <- 1/w[T[i],i]
}

# Posterior predictive based on the Bias Corrected component...

mu.new.1 ~ dnorm(mu[1], inv.var[1])
#mu.new.2 ~ dnorm(mu[2], inv.var[2])
#y.bias.new ~ dbern(p.bias[1])

#mu.new <- (1-y.bias.new)*mu.new.1 + y.bias.new*mu.new.2
#inv.var.new <- (1-y.bias.new)*inv.var[1] + y.bias.new*inv.var[2]

mu.new <- mu.new.1
inv.var.new <- inv.var[1]

# Approximate Bayesian Cross-Validation ......................
theta.ghost ~ dnorm(mu.new, inv.var.new)

for(i in 1:N)
{
    y.ghost[i] ~ dnorm(theta.ghost, pre.y[i])
}



# Priors for hyper-parameters .................................

# Prior of probability classes

p.bias[2] ~ dbeta(a.0, a.1)
p.bias[1] <- 1 - p.bias[2]

# Variance components ...

tau <- 1/sqrt(inv.var[1])

inv.var[1] ~ dscaled.gamma(scale.sigma.between,
                          df.scale.between)

inv.var[2] <- inv.var[1]

# Prior for mu
inv.var.mu <- pow(sd.mu, -2)
      mu[1] ~ dnorm(mean.mu, inv.var.mu)
          B ~ dunif(B.lower, B.upper)
    # mu[2] <- mu[1] + B
    mu[2] <- B
  }
  "

if (is.null(parallel)) { #execute R2jags
  model.bugs.connection <- textConnection(model.bugs)

  # Use R2jags as interface for JAGS ...

  results <- jags( data = data.bcmeta,
                   parameters.to.save = par.bcmeta,
                   model.file = model.bugs.connection,
                   n.chains = nr.chains,
                   n.iter = nr.iterations,
                   n.burnin = nr.burnin,
                   n.thin = nr.thin,
                   DIC = TRUE,
                   pD=TRUE)

  # Close text connection
  close(model.bugs.connection)
}else if(parallel == "jags.parallel"){
  writeLines(model.bugs, "model.bugs")
  results <- jags.parallel(     data = data.bcmeta,
                                parameters.to.save = par.bcmeta,
                                model.file = "model.bugs",
                                n.chains = nr.chains,
                                n.iter = nr.iterations,
                                n.burnin = nr.burnin,
                                n.thin = nr.thin,
                                DIC=TRUE)

  #Compute pD from result
  results$BUGSoutput$pD = results$BUGSoutput$DIC - results$BUGSoutput$mean$deviance

  # Delete model.bugs on exit ...
  unlink("model.bugs")
}


  # Extra outputs that are linked with other functions ...
  results$data = data

  # Hyperpriors parameters............................................
  results$prior$mean.mu     = mean.mu
  results$prior$sd.mu       = sd.mu
  results$prior$scale.sigma.between = scale.sigma.between
  results$prior$df.scale.between    = df.scale.between
  results$prior$B.lower     = B.lower
  results$prior$B.upper     = B.upper

  results$prior$a.0         = a.0
  results$prior$a.1         = a.1
  results$prior$nu          = nu
  results$prior$nu.estimate = nu.estimate
  results$prior$b.0 = b.0
  results$prior$b.1 = b.1

  class(results) = c("bcmeta")

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
print.bcmeta <- function(x, digits, ...)
{
  print(x$BUGSoutput,...)
}





