#' @title Generalized Evidence Synthesis with Direct Penalization of Observational Studies
#'
#'
#' @description This function fits a hierarchical meta-regression model based on
#' a three levels random-effects model with direct penalization of observational studies.
#' 
#'
##' @param data            A data frame where the first colum is the treatment effect TE, the second column
##'                        is the standard error of the treatment effect seTE and the third column is a factor
##'                        indicating the statistical design of the study.
##' 
##' @param EmBi           Indicates if Empirical Bias Penalization is performed. The default is FALSE.
##' @param mu.bias        Empirical mean bias. The default value is -1.25 (BRANDO study, Savovic et al. 2012).
##' @param sigma.bias     Empirical bias sd. The default value is 0.1 (BRANDO study, Savovic et al. 2012).
##' 
##' @param ExPe           Indicates if Explicit Bias Penalization is performed.  The default is FALSE. 
##' @param K              Prior expected value of the weights used for observational studies, E(w_OS)=K. 
##'                       The default value is 0.5.
##' @param df             The degrees of freedom for the scale mixture distribution. The default value is 4.                      
##' 
##' @param mean.mu.0       Prior mean of the overall mean parameter mu.0. The default value is 0.
##' @param sd.mu.0         Prior standard deviation of mu.0
##' @param scale.sigma.within Prior scale parameter for the standard deviation within study type. The default value is 0.5.
##' @param scale.sigma.between Prior scale parameter for the standard deviation between study type. The default value is 0.5.
##' @param df.scale.within Degrees of freedom of the standard deviation within study type. The default value is 1.
##' @param df.scale.between Degrees of freedom of the standard deviation between study type. The default value is 1.
##' 
##' @param nr.chains       Number of chains for the MCMC computations, default 2.
##' @param nr.iterations   Number of iterations after adapting the MCMC, default is 10000. Some models may need more iterations.
##' @param nr.adapt        Number of iterations in the adaptation process, defualt is 1000. Some models may need more iterations during adptation.
##' @param nr.burnin       Number of iteration discared for burnin period, default is 1000. Some models may need a longer burnin period.
##' @param nr.thin         Thinning rate, it must be a positive integer, the default value 1.
##' @param be.quiet        Do not print warning message if the model does not adapt. The default value is FALSE. If you are not sure about the adaptation period choose be.quiet=TRUE.
##' @param r2jags          Which interface is used to link R to JAGS (rjags and R2jags), default value is R2Jags=TRUE.
##'
##' @return This function returns an object of the class "ges". This object contains the MCMC output of
##' each parameter and hyper-parameter in the model, the data frame used for fitting the model and 
##' type of random effects distribution
##'
##' The results of the object of the class ges can be extracted with R2jags or with rjags. In addition
##' a summary, a print and a plot functions are implemented for this type of object.
##' 
##' @references Verde, P.E. and Curcio, D. (2017) Hierarchical Meta-Regression Modelling: The Case of The Pneumococcal Polysaccharide Vaccine. Technical Report. 
##'
##' @examples
##' 
##' \dontrun{
##' 
##'  library(jarbes)
##'  data(ppvipv)
##'  
##'  # Three levels hierarchical model ....
##'  m1.ges = ges(ppvipv)
##'  summary(m1.ges)
##'  
##'  # Three levels hierarchical mode with empirical bias penalization...
##'  m2.ges = ges(ppvipv, EmBi = TRUE)
##'  summary(m2.ges)
##'  
##'  # Three levels hierarchical model with explicit bias penalization ...
##'  m3.ges = ges(ppvipv, ExPe = TRUE, K = 0.5, df = 10)
##'  summary(m3.ges)
##' }
#'
#' 
#' @import R2jags
#' @import rjags
#'
#' @export
#'

ges <- function(
  data,
  # Empirical Bias information .......................................
  EmBi = FALSE,
  mu.bias = -1.25,
  sigma.bias = 0.1,
  # Explicit Bias Penalization .......................................
  ExPe            = FALSE,
  # E(w_OS) = K 
  K               = 0.5,
  df              = 4,
  # Hyperpriors parameters............................................
  mean.mu.0       = 0,
  sd.mu.0         = 10,
  scale.sigma.within = 0.5,
  scale.sigma.between = 0.5,
  df.scale.within = 1,
  df.scale.between = 1,
  # MCMC setup........................................................
  nr.chains       = 2,
  nr.iterations   = 20000,
  nr.adapt        = 1000,
  nr.burnin       = 1000,
  nr.thin         = 1,
  # Further options to link jags and R ...............................
  be.quiet        = FALSE,
  r2jags          = TRUE)UseMethod("ges")

#' @export
#' 
ges.default <- function(
  data,
  # Empirical Bias information .......................................
  EmBi = FALSE,
  mu.bias = -1.25,
  sigma.bias = 0.1,
  # Explicit Bias Penalization .......................................
  ExPe            = FALSE,
  # E(w_OS) = K 
  K               = 0.5,
  df              = 4,
  # Hyperpriors parameters............................................
  mean.mu.0       = 0,
  sd.mu.0         = 10,
  scale.sigma.within = 0.5,
  scale.sigma.between = 0.5,
  df.scale.within = 1,
  df.scale.between = 1,
  # MCMC setup........................................................
  nr.chains       = 2,
  nr.iterations   = 20000,
  nr.adapt        = 1000,
  nr.burnin       = 1000,
  nr.thin         = 1,
  # Further options to link jags and R ...............................
  be.quiet        = FALSE,
  r2jags          = TRUE)
{
  # Model errors checking....
  if(EmBi == TRUE & ExPe == TRUE)stop("This random effects distribution is not implemented")
  
  # Load module "glm" in jags....
  load.module("glm")
  
  # Data setup......................................................
  y <- data$TE
  se.y <- data$seTE
  N <- length(y)
  x <- data$design
  x <- as.numeric(x)
  Ndesign <- length(table(x))
  NRCTs <- table(x)[1]
  
  # This list describes the data used by the BUGS script.
  if(EmBi == TRUE & ExPe == FALSE)
    data.model <- list( y = y,
                        se.y = se.y,
                        x = x,
                        N = N,
                        Ndesign = Ndesign,
                   #     NRCTs = NRCTs,
                        mean.mu.0 = mean.mu.0,
                        sd.mu.0 = sd.mu.0,
                        scale.sigma.within = scale.sigma.within,
                        scale.sigma.between = scale.sigma.between,
                        df.scale.within = df.scale.within,
                        df.scale.between = df.scale.between,
                        mu.bias = mu.bias,
                        sigma.bias = sigma.bias)
  else
    if(ExPe == TRUE & EmBi == FALSE)
      data.model <- list( y = y,
                          se.y = se.y,
                          x = x,
                          N = N,
                          Ndesign = Ndesign,
                          NRCTs = NRCTs,
                          mean.mu.0 = mean.mu.0,
                          sd.mu.0 = sd.mu.0,
                          scale.sigma.within = scale.sigma.within,
                          scale.sigma.between = scale.sigma.between,
                          df.scale.within = df.scale.within,
                          df.scale.between = df.scale.between, 
                          df  = df,
                          K  = K)
  else
      data.model <- list( y = y,
                          se.y = se.y,
                          x = x,
                          N = N,
                          Ndesign = Ndesign,
                          mean.mu.0 = mean.mu.0,
                          sd.mu.0 = sd.mu.0,
                          scale.sigma.within = scale.sigma.within,
                          scale.sigma.between = scale.sigma.between,
                          df.scale.within = df.scale.within,
                          df.scale.between = df.scale.between) 
  
  # Parameters to monitor ............................................
  parameters.model <-
    c("theta",
      "sigma.theta",
      "sigma.theta.between",
      "mu.0",  
      "mu",
      "OR",
      "OR.pool",
      "mu.new",
      "OR.new",
      "OR.pool.new")
  
  # Parameters for different models ...
  if(EmBi == TRUE)parameters.model <- c(parameters.model, "mu.design")
  
  if(ExPe == TRUE)parameters.model <- c(parameters.model, "w")
  
  # Model construction
  blueprint <- function(EmBi = FALSE, ExPe = FALSE)
  {
    
    # Model 1:
    # Bayesian three levels hierarchical model (3LHM)
    model.1 <- 
    "model
    {
    for( i in 1 : N) {
    # Likelihood of theta[i] ...........................................
    y[i] ~ dnorm(theta[i], pre.y[i])
    pre.y[i] <- pow(se.y[i], -2)
    
    # Group random effects ...........................................
    theta[i] ~ dnorm(mu[x[i]], pre.theta[x[i]])
    }
    
    # Between designs variability ....................................
    for(i in 1:Ndesign){
    mu[i] ~ dnorm(mu.0, pre.between)
    }
    
    # Priors.........................................................
    sigma.theta.between <- 1/sqrt(pre.between)
    pre.between ~ dscaled.gamma(scale.sigma.between, df.scale.between)
    pre.mu.0 <- 1/sd.mu.0^2
    mu.0 ~ dnorm(mean.mu.0, pre.mu.0)
    
    for(i in 1:Ndesign){
    sigma.theta[i] <- 1/sqrt(pre.theta[i])
    pre.theta[i] ~ dscaled.gamma(scale.sigma.within, df.scale.within) 
    }
    
    # Functional parameters of interest..............................  
    
    OR.pool <- exp(mu.0)
    mu.pool.new ~ dnorm(mu.0, pre.between)T(-4, 4)
    OR.pool.new <- exp(mu.pool.new) 
    
    for(i in 1:Ndesign){
    mu.new[i] ~ dnorm(mu[i], pre.theta[i])T(-4, 4)
    OR[i] <- exp(mu[i])
    OR.new[i] <- exp(mu.new[i])
    }
    }"

    # Model 2: 
    # Bayesian 3LHM with Empirical Bias...
    # The empirical bias model is based on the mixture model derived from a 
    # meta-analysis of empirical bias evidence.
    # I have to implement the mixture stuff here!
    model.2 <-
    "model
    {
    for( i in 1 : N) {
    # Likelihood of theta[i] ...........................................
    y[i] ~ dnorm(theta[i], pre.y[i])
    pre.y[i] <- pow(se.y[i], -2)
    
    # Group random effects ...........................................
    theta[i] ~ dnorm(mu.design[i], pre.theta[x[i]])
    mu.design[i] <- mu[x[i]] + bias[x[i]]
    }
    
    # Between designs variability ....................................
    for(i in 1:Ndesign){
    mu[i] ~ dnorm(mu.0, pre.between)
    }
    
    # Bias correction for observational studies ......................
    bias[1] <- 0 # RCTs
    
    # Observational studies ...
    # Default values:
    #  mu.bias = -0.2
    #  sigma.bias = 0.1

    pre.bias <- pow(sigma.bias, -2)
    
    for(k in 2:Ndesign)
    {
    bias[k] ~ dnorm(mu.bias, pre.bias)
    }
    
    # Priors.........................................................
    
    sigma.theta.between <- 1/sqrt(pre.between)
    pre.between ~ dscaled.gamma(scale.sigma.between, df.scale.between)
    
    pre.mu.0 <- 1/sd.mu.0^2
    mu.0 ~ dnorm(mean.mu.0, pre.mu.0)
    
    for(i in 1:Ndesign){
    sigma.theta[i] <- 1/sqrt(pre.theta[i])
    pre.theta[i] ~ dscaled.gamma(scale.sigma.within, df.scale.within) 
    }
    
    # Functional parameters of interest..............................  
    
    OR.pool <- exp(mu.0)
    mu.pool.new ~ dnorm(mu.0, pre.between)T(-4, 4)
    OR.pool.new <- exp(mu.pool.new) 
    
    for(i in 1:Ndesign){
    mu.new[i] ~ dnorm(mu.design[i], pre.theta[i])T(-4, 4)
    OR[i] <- exp(mu.design[i])
    OR.new[i] <- exp(mu.new[i])
    bias.mu[i] <- mu.0 - mu[i]
    }
    }"

    # Model 3: This is the Verde (2016, 2017) explicit bias penalization ....
    # Bayesian 3LHM with scale mixture random effects 
    # Mean OS fixed at K with a dispersion of degrees of freedom fixed at df.
    model.3 <- 
    "model
    {
    for( i in 1 : N) {
    # Likelihood of theta[i] ...........................................
    y[i] ~ dnorm(theta[i], pre.y[i])
    pre.y[i] <- pow(se.y[i], -2)
    
    # Group random effects ...........................................
    theta[i] ~ dnorm(mu[x[i]], pre.theta[i])
    pre.theta[i] <- w[i]*inv.sigma.2[x[i]]
    }
    
    for(i in 1:Ndesign)
    {
    inv.sigma.2[i] ~ dscaled.gamma(scale.sigma.within, df.scale.within)
    sigma.theta[i] <- 1/sqrt(inv.sigma.2[i])
    }
    
    # Within variability for observational studies and penalization ..........
    
    for(i in 1:NRCTs)
    {
    w[i] <- 1
    }
    
    # E(w_OS) = K
    alpha <- df / 2
    delta <- df*(1-K)/K
    beta <- (df + delta)/2
    
    #startOS <- NRCTs+1
    for(i in (NRCTs+1):N)  
    {
    w[i] ~ dgamma(alpha, beta)
    }
    
    # Between designs variability ....................................
    for(i in 1:Ndesign)
    {
    mu[i] ~ dnorm(mu.0, pre.between)
    }
    
    # Priors.........................................................
    sigma.theta.between <- 1/sqrt(pre.between)
    pre.between ~ dscaled.gamma(scale.sigma.between, df.scale.between)
    
    pre.mu.0 <- 1/sd.mu.0^2
    mu.0 ~ dnorm(mean.mu.0, pre.mu.0)
    
    # Functional parameters of interest..............................  
    OR.pool <- exp(mu.0)
    mu.pool.new ~ dnorm(mu.0, pre.between)T(-4, 4)
    OR.pool.new <- exp(mu.pool.new) 
    
    for(i in 1:Ndesign)
    {
    w.new[i] ~ dgamma(alpha, beta)
    pre.theta.new[i] <- w.new[i]*inv.sigma.2[x[i]]
    mu.new[i] ~ dnorm(mu[i], pre.theta.new[i])T(-4, 4)
    OR[i] <- exp(mu[i])
    OR.new[i] <- exp(mu.new[i])
    bias.mu[i] <- mu.0 - mu[i]
    }
    }"

    if(EmBi == FALSE & ExPe == FALSE) return(model.1)
    else if(EmBi == TRUE) return(model.2)
    else return(model.3)
  }
  
  model.bugs <- blueprint(EmBi, # Empirical Bias modeling ...
                          ExPe) # Explicit Penalization ...

  if(r2jags == TRUE){
    # Use R2jags as interface for JAGS ...
    results <- jags(              data = data.model,
                                  parameters.to.save = parameters.model,
                                  #   inits = inits.model,
                                  model.file = textConnection(model.bugs),
                                  n.chains = nr.chains,
                                  n.iter = nr.iterations,
                                  n.burnin = nr.burnin,
                                  n.thin = nr.thin
    )
  }
  else {
    # Use rjags as interface for JAGS ...
    # Send the model to JAGS, check syntax, run ...
    jm <- jags.model(file     = textConnection(model.bugs),
                     data     = data.model,
                     #  inits    = inits.model,
                     n.chains = nr.chains,
                     n.adapt  = nr.adapt,
                     quiet    = be.quiet)
    
    results <- coda.samples(jm,
                            variable.names = parameters.model,
                            n.iter         = nr.iterations)
  }
  
  if(r2jags == FALSE)
  {cat("You are using the package rjags as interface to JAGS.", "\n")
    cat("The plot functions for output analysis are not implemented in this version", "\n")
  }
  
  # Extra outputs that are linked with other functions
  
  results$data <- data
  results$EmBi <- EmBi
  results$ExPe <- ExPe
  class(results) <- c("ges")
  return(results)
  }




#' Generic print function for ges object in jarbes.
#' 
#' @param x The object generated by the function ges.
#'
#' @param digits The number of significant digits printed. The default value is 3.
#'
#' @param ... \dots
#'
#' @export

print.ges <- function(x, digits, ...)
{
  print(x$BUGSoutput,...)
  
}





#' Generic summary function for ges object in jarbes
#' @param object The object generated by the ges function.
#'
#' @param digits The number of significant digits printed. The default value is 3.
#'
#' @param intervals A numeric vector of probabilities with values in [0,1]. The default value is
#'                  intervals = c(0.025, 0.5, 0.975).
#'
#' @param ... \dots
#'
#' @export
#'
summary.ges <- function(object, digits = 3,  intervals = c(0.025, 0.5, 0.975), ...)
{
  
  print(object$BUGSoutput, digits= digits, intervals = intervals, ...)
  
}



#' Generic plot function for ges object in jarbes.
#' 
#' @param x The object generated by the ges function.
#' 
#' @param ... \dots
#' 
#' @export
#' 
plot.ges <- function(x, ...)
{
  plot(x)
}
















