#' @title Generalized Evidence Synthesis with Indirect Penalization of Observational Studies
#'
#' @description This function performers a Bayesian meta-analysis to jointly
#' combine different types of studies. The random-effects follows a finite
#' mixture of normals.
#' 
#'
#'
#' @param data            A data frame with at least three variables with the following names:
#'                        1) TE = treatment effect, 
#'                        2) seTE = the standard error of the treatment effect
#'                        3) design = the study type. The other columns of the data 
#'                        frame may correspond to the study name and other covariates that
#'                        characterize the study and could influence internal validity bias
#'                        (e.g. evaluation of risk of bias).
#' 
#' @param mean.mu        Prior mean of the overall mean parameter mu, default value is 0.
#' @param sd.mu          Prior standard deviation of mu, the default value is 10^-6.
#' @param scale.sigma    Prior scale parameter for the standard deviation between studies, 
#'                       the default value is 0.5.
#'                       
#' @param K.lower        Lower bound of the bias parameter K, the default value is -10.
#' @param K.upper        Upper bound of the bias parameter K, the default value is 0.                       
#' 
#' 
#' @param nr.chains       Number of chains for the MCMC computations, default 2.
#' @param nr.iterations   Number of iterations after adapting the MCMC, default is 10000. Some models may need more iterations.
#' @param nr.adapt        Number of iterations in the adaptation process, defualt is 1000. Some models may need more iterations during adptation.
#' @param nr.burnin       Number of iteration discared for burnin period, default is 1000. Some models may need a longer burnin period.
#' @param nr.thin         Thinning rate, it must be a positive integer, the default value 1.
#' @param be.quiet        Do not print warning message if the model does not adapt. The default value is FALSE. If you are not sure about the adaptation period choose be.quiet=TRUE.
#' @param r2jags          Which interface is used to link R to JAGS (rjags and R2jags), default value is R2Jags=TRUE.
#'
#' @return This function returns an object of the class "gesmix". 
#' This object contains the MCMC output of
#' each parameter and hyper-parameter in the model and 
#' the data frame used for fitting the model.
#' 
#'  
#' @details The results of the object of the class gesmix can be extracted with R2jags
#' or with rjags. In addition a summary, a print and a plot functions are 
#' implemented for this type of object.
#'
#'
#' @references Verde, P. E. (2017) Two Examples of Bayesian Evidence Synthesis with the Hierarchical Meta-Regression Approach. Chap.9, pag 189-206. Bayesian Inference, ed. Tejedor, Javier Prieto. InTech.
#' 
#' @references Verde, P.E. and Curcio, D. (2019) Hierarchical Meta-Regression Modelling: The Case of The Pneumococcal Polysaccharide Vaccine. Technical Report. 
#
#' @examples
#' \dontrun{
#' 
#'  library(jarbes)
#'  
#' out = gesmix(data = ppvipd)
#' attach.jags(out, overwrite = TRUE)
#'
#' par(mfrow = c(1,2))
#' hist(lambda[,1],  breaks = 100, probability = TRUE,
#'     ylim = c(0, 3),
#'     xlim = c(-2, 0.85), main = "",
#'     xlab = "Treatment Effect: log(Odds Ratio)")
#'
#' lines(density(lambda[,2]), lwd = 3, col = "red")
#' lines(density(lambda[,1]), lwd = 3, col = "blue")
#' abline(v = 0, lty = 2)
#' legend(-2, 3, legend = c("RCTs' Component",
#'                         "Obsevationals' Component"),
#'       col = c("blue", "red"),
#'       lty = rep(1, 2),
#'       lwd = rep(2, 2))
#' hist(c(lambda[,1],lambda[,2]+0.337), 
#'     breaks = 150,
#'     probability = TRUE,
#'     ylim = c(0, 3.5),
#'     xlim = c(-2.1, 0.2), main = "",
#'     xlab = "Treatment Effect: log(Odds Ratio)")
#' lines(density(c(lambda[,1],lambda[,2]+0.337)), col = "blue", lwd =3)
#' abline(v = quantile(c(lambda[,1],lambda[,2]+0.337), 
#'                    probs = c(0.025,0.5,0.975)), col ="blue",
#'       lty = 2)
#' m.pool <- log(0.43)
#' sd.pool <- (log(0.54) - log(0.34))/(2*1.96)
#' curve(dnorm(x, m = m.pool, s = sd.pool), 
#'      from = -3, to = 0.5, add = TRUE, col = "black", lwd =3)
#' legend(-2.1, 3, 
#'       legend = c("Bias-adjusted",
#'                  "Biased result"),
#'       col = c("blue", "black"),
#'       lty = rep(1, 2),
#'       lwd = rep(2, 2))
#' abline(v = 0, lty = 2)
#' arrows(-1.5, 1.1, -1, 1.55, lwd =2, length=0.1, angle=20)
#' text(-1.8, 1, "Ignoring study types")
#' par(mfrow = c(1,1))
#' }
#' @import R2jags
#' @import rjags
#' @export
#'
#'
gesmix <- function(
  data,
  # Hyperpriors parameters............................................
  mean.mu     = 0,
  sd.mu       = 10^6,
  scale.sigma = 0.5,
  K.lower     = -10,     
  K.upper     = 0,
  
  # MCMC setup........................................................
  nr.chains       = 2,
  nr.iterations   = 10000,
  nr.adapt        = 1000,
  nr.burnin       = 1000,
  nr.thin         = 1,
  # Further options to link jags and R ...............................
  be.quiet        = FALSE,
  r2jags          = TRUE
)UseMethod("gesmix")


#'
#' @export
#' 

gesmix.default <- function(
  data,
  # Hyperpriors parameters............................................
  mean.mu     = 0,
  sd.mu       = 10^6,
  scale.sigma = 0.5,
  K.lower     = -10,     
  K.upper     = 0,
  
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
  
  # Mixture of Normals Random effects meta-analysis
  y = sort(data$TE)
  se.y = data$seTE[order(data$TE)]
  N = length(y)
  x = data$design[order(data$TE)]
  x = as.numeric(x)
  
  # T is the index that allocates each observation to one of the
  # two normals.
  T = rep(NA, N)
  T[1] = 1
  T[N] = 2
  
  # This list describes the data used by the BUGS script.
  data.gesmix = list(y = y, se.y = se.y, N = N, T = T, x = x)
  
  # List of parameters
  par.gesmix <- c("theta",
                    "theta.corrected",
                    "tau",
                    "P",
                    "K",
                    "lambda",
                    "alpha")
  
  # Model in BUGS
  
  model.bugs <- "model
  {
  for( i in 1 : N ) {
  # Likelihood of theta[i] ...........................................
  y[i] ~ dnorm(theta[i], pre.y[i])
  pre.y[i] <- pow(se.y[i], -2)
  
  # Mixture random effects ...........................................
  theta[i] ~ dnorm(mu[i], prec.tau)
  mu[i] <- lambda[T[i]]    # mean of component i=1,2
  T[i] ~ dcat(P[i, 1:2])        # random index
  
  logit(P[i, 2]) <- alpha.0 + alpha[x[i]]
  P[i, 1] <- 1 - P[i, 2]
  
  # Bias correction for theta .......................................
  
  theta.corrected[i] <- theta[i]*(1-P[i,2]) + (theta[i]-K)*P[i,2]
  
  # theta corrected is theta if Pr bias is 0
  # theta is theta - K if Pr bias is 1
  }
  
  K ~ dunif(-10, 0.0)            # prior for K
  lambda[2] <- lambda[1] + K     # calculation of mean lambda[1]
  lambda[1] ~ dnorm(0.0, 0.001)  # prior for lambda[1]
  
  prec.tau <- pow(tau, -2)
  tau ~ dunif(0, 100)     # prior for sigma
  
  alpha.0 ~ dnorm(0.0, 0.001)
  alpha[1] <- 0
  alpha[2] ~ dnorm(0.0, 0.001)
  alpha[3] ~ dnorm(0.0, 0.001)
  }"

  model.bugs.connection <- textConnection(model.bugs)
  
  if(r2jags == TRUE){
    # Use R2jags as interface for JAGS ...
    results <- jags(              data = data.gesmix,
                                  parameters.to.save = par.gesmix,
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
                     data     = data.gesmix,
                     #  inits    = inits.model,
                     n.chains = nr.chains,
                     n.adapt  = nr.adapt,
                     quiet    = be.quiet)
    
    results <- coda.samples(jm,
                            variable.names = par.gesmix,
                            n.iter         = nr.iterations)
  }
  
  if(r2jags == FALSE)
  {cat("You are using the package rjags as interface to JAGS.", "\n")
    cat("The plot functions for output analysis are not implemented in this jarbes version", "\n")
  }
  
  # Close text conection
  close(model.bugs.connection)
  
  # Extra outputs that are linked with other functions ...
  
  results$data = data
  class(results) = c("gesmix")
  return(results)
}



#' Generic print function for metamixt object in jarbes.
#' 
#' @param x The object generated by the function metamixt.
#'
#' @param digits The number of significant digits printed. The default value is 3.
#'
#' @param ... \dots
#'
#' @export

print.gesmix <- function(x, digits, ...)
{
  print(x$BUGSoutput,...)
  
}





#' Generic summary function for metamixt object in jarbes
#' @param object The object generated by the metamixt function.
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
summary.gesmix <- function(object, digits = 3,  intervals = c(0.025, 0.5, 0.975), ...)
{
  
  print(object$BUGSoutput, digits= digits, intervals = intervals, ...)
  
}



#' Generic plot function for metamixt object in jarbes.
#' 
#' @param x The object generated by the metamixt function.
#' 
#' @param ... \dots
#' 
#' @export
#' 
plot.gesmix <- function(x, ...)
{
  plot(x)
}






