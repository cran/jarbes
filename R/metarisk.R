#' Bayesian meta-analysis to investigate the relationship between treatment 
#' effect and underlying risk.
#'
#' This function performers a Bayesian meta-analysis to analyse heterogeneity of
#' the treatment effect as a function of the baseline risk. The function fits a 
#' hierarchical meta-regression model based on a bivariate random effects model. 
#' 
#' The number of events in the control and treated group are modeled with two
#' conditional Binomial distributions and the random-effects are based on a 
#' bivariate scale mixture of Normals.
#' 
#' The function calculates the implicit hierarchical meta-regression, where the 
#' treatment effect is regressed to the baseline risk (rate of events in the control
#' group). The scale mixture weights are used to adjust for internal validity and 
#' structural outliers identification.
#' 
#' Computations are done by calling JAGS (Just Another Gibbs Sampler) to perform
#' MCMC (Markov Chain Monte Carlo) sampling and returning an object of the
#' class \emph{mcmc.list}.
#'
#' Installation of JAGS: It is important to note that R 3.3.0 introduced a major change in the
#' use of toolchain for Windows. This new toolchain is incompatible with older packages written in C++.
#' As a consequence, if the installed version of JAGS does not match the R installation, then the rjags
#' package will spontaneously crash. Therefore, if a user works with R version >= 3.3.0, then JAGS must
#' be installed with the installation program JAGS-4.2.0-Rtools33.exe. For users who continue using R 3.2.4 or
#' an earlier version, the installation program for JAGS is the default installer JAGS-4.2.0.exe.
#'
#'
#' @param data            A data frame where the first four columns containing the number of events in
#'                        the control group (yc), the number of patients in the control group (nc), 
#'                        the number of events in the treatment group (yt) and the number of patients in the
#'                        treatment group (nt). If two.by.two = TRUE a data frame where each line 
#'                        contains the trial results with column names: yc, nc, yt, nt.
#'
#' @param two.by.two      If TRUE indicates that the trial results are with names: yc, nc, yt, nt.
#'
#' @param re              Random effects distribution for the resulting model. Possible
#'                        values are \emph{normal} for bivariate random effects and \emph{sm} 
#'                        for scale mixtures.
#'
#' @param link            The link function used in the model. Possible values are
#'                        \emph{logit}, \emph{cloglog} \emph{probit}.
#'
#' @param mean.mu.1       Prior mean of baseline risk, default value is 0.
#'
#' @param mean.mu.2       Prior mean of treatment effect, default value is 0.
#'
#' @param sd.mu.1         Prior standard deviation of mu.1, default value is 1. The default prior of mu.1 is a 
#'                        logistic distribution with mean 0 and dispersion 1. The implicit prior for mu.1 in 
#'                        the probability scale is a uniform between 0 and 1.
#'
#' @param sd.mu.2         Prior standard deviation of mu.2, default value is 1. The default prior of mu.2 is a 
#'                        logistic distribution with mean 0 and dispersion 1. The implicit prior for mu.2 in 
#'                        the probability scale is a uniform between 0 and 1.
#'
#' @param sigma.1.upper   Upper bound of the uniform prior of sigma.1, default value is 5.
#'
#' @param sigma.2.upper   Upper bound of the uniform prior of sigma.2, default value is 5.
#'
#' @param mean.Fisher.rho Mean of rho in the Fisher scale default value is 0.
#'
#' @param sd.Fisher.rho   Standard deviation of rho in the Fisher scale, default value is 1/sqrt(2).
#'
#' @param df              If de.estimate = FALSE, then df is the degrees of freedom for the scale mixture distribution, default value is 4.
#'
#' @param df.estimate     Estimate the posterior of df. The default value is FALSE.
#'
#' @param df.lower        Lower bound of the prior of df. The default value is 3.
#'
#' @param df.upper        Upper bound of the prior of df. The default value is 30.
#'
#' @param split.w         Split the w parameter in two independent weights one for each random effect. 
#'                        The default value is FALSE.
#'
#' @param nr.chains       Number of chains for the MCMC computations, default 5.
#'
#' @param nr.iterations   Number of iterations after adapting the MCMC, default is 10000. Some models may need more iterations.
#'
#' @param nr.adapt        Number of iterations in the adaptation process, defualt is 1000. Some models may need more iterations during adptation.
#'
#' @param nr.burnin       Number of iteration discared for burnin period, default is 1000. Some models may need a longer burnin period.
#'
#' @param nr.thin         Thinning rate, it must be a positive integer, the default value is 1.
#'
#' @param be.quiet        Do not print warning message if the model does not adapt default value is FALSE. If you are not sure about the adaptation period choose be.quiet=TRUE.
#'
#' @param r2jags          Which interface is used to link R to JAGS (rjags and R2jags) default value is R2Jags=TRUE.
#'
#' @return This function returns an object of the class "metarisk". This object contains the MCMC output of
#' each parameter and hyper-parameter in the model, the data frame used for fitting the model, the link function,
#' type of random effects distribution and the splitting information for conflict of evidence analysis.
#'
#' The results of the object of the class metadiag can be extracted with R2jags or with rjags. In addition
#' a summary, a print and a plot functions are implemented for this type of object.
#'
#' @references Verde, P.E, Ohmann, C., Icks, A. and Morbach, S. (2016) Bayesian evidence synthesis and combining randomized and nonrandomized results: a case study in diabetes. Statistics in Medicine. Volume 35, Issue 10, 10 May 2016, Pages: 1654 to 1675.
#'
#' @references Verde, P.E. and Curcio, D. (2017) Hierarchical Meta-Regression Modelling: The Case of The Pneumococcal Polysaccharide Vaccine. Technical Report. 
#' 
#' @references Verde, P.E. (2018) The Hierarchical Meta-Regression Approach and Learning from Clinical Evidence. Technical report.
#'
#' @examples
#' 
#' 
#' library(jarbes)
#' 
#' # This example is used to test the function and it runs in about 5 seconds.
#' # In a real application you must increase the number of MCMC interations. 
#' # For example use: nr.burnin = 5000 and nr.iterations = 20000
#' 
#' data(ppvcap)
#' dat <- ppvcap[, c("yt","nt","yc","nc")]
#' res.hmr.1 <- metarisk(dat,             # Data frame
#'                       two.by.two = FALSE,    # Data is given as: (tp, n1, fp, n2)
#'                       re = "sm",             # Random effects distribution
#'                       link = "logit",        # Link function
#'                       sd.Fisher.rho   = 1.7, # Prior standard deviation of correlation
#'                       split.w = TRUE,        # Split the studies' weights
#'                       df.estimate = TRUE,    # Estimate the posterior of df parameter.
#'                       nr.burnin = 100,       # Iterations for burnin
#'                       nr.iterations = 1000,  # Total iterations
#'                       nr.chains = 2,         # Number of chains
#'                       r2jags = TRUE)         # Use r2jags as interface to jags
#'                       
#' print(res.hmr.1, digits = 3)
#' 
#'
#' 
#' \donttest{
#' 
#' # The following examples corresponds to Section 4 in Verde (2018).
#' # These are simulated examples to investigate internal and external validity
#' # bias in meta-analysi.
#' 
#'
#'library(MASS)
#'library(ggplot2)
#'library(jarbes)
#'library(gridExtra)
#'library(mcmcplots)
#'
#' #Experiment 1: External validity bias
#'
#'set.seed(2018)
#' # mean control 
#' pc <- 0.7
#' # mean treatment 
#'pt <- 0.4
#' # reduction of "amputations" odds ratio
#'OR <- (pt/(1-pt)) /(pc/(1-pc))
#'OR
#'
#'# mu_2
#'log(OR)
#'mu.2.true <- log(OR)
#' #sigma_2
#' sigma.2.true <- 0.5
#' # mu_1
#'mu.1.true <- log(pc/(1-pc))
#'mu.1.true
#'#sigma_1
#'sigma.1.true <- 1
#'# rho
#'rho.true <- -0.5
#'Sigma <- matrix(c(sigma.1.true^2, sigma.1.true*sigma.2.true*rho.true,
#'                  sigma.1.true*sigma.2.true*rho.true, sigma.2.true^2), 
#'                  byrow = TRUE, ncol = 2)
#'Sigma
#'
#'theta <- mvrnorm(35, mu = c(mu.1.true, mu.2.true),
#'                 Sigma = Sigma )
#'                 
#'                 
#'plot(theta, xlim = c(-2,3))
#'abline(v=mu.1.true, lty = 2)
#'abline(h=mu.2.true, lty = 2)
#'abline(a = mu.2.true, b=sigma.2.true/sigma.1.true * rho.true, col = "red")
#'abline(lm(theta[,2]~theta[,1]), col = "blue")
#'
#'# Target group
#'mu.T <- mu.1.true + 2 * sigma.1.true
#'abline(v=mu.T, lwd = 3, col = "blue")
#'eta.true <- mu.2.true + sigma.2.true/sigma.1.true*rho.true* mu.T
#'eta.true
#'exp(eta.true)
#'abline(h = eta.true, lty =3, col = "blue")
#'# Simulation of each primary study:
#'n.c <- round(runif(35, min = 30, max = 60),0)
#'n.t <- round(runif(35, min = 30, max = 60),0)
#'y.c <- y.t <- rep(0, 35)
#'p.c <- exp(theta[,1])/(1+exp(theta[,1]))
#'p.t <- exp(theta[,2]+theta[,1])/(1+exp(theta[,2]+theta[,1]))
#'for(i in 1:35)
#'{
#'  y.c[i] <- rbinom(1, n.c[i], prob = p.c[i])
#'  y.t[i] <- rbinom(1, n.t[i], prob = p.t[i])
#'}
#'
#'AD.s1 <- data.frame(yc=y.c, nc=n.c, yt=y.t, nt=n.t)
#'
#'##########################################################
#'incr.e <- 0.05
#'incr.c <- 0.05
#'n11 <- AD.s1$yt 
#'n12 <- AD.s1$yc 
#'n21 <- AD.s1$nt - AD.s1$yt 
#'n22 <- AD.s1$nc - AD.s1$yc
#'AD.s1$TE <- log(((n11 + incr.e) * (n22 + incr.c))/((n12 + incr.e) * (n21 + incr.c)))
#'AD.s1$seTE <- sqrt((1/(n11 + incr.e) + 1/(n12 + incr.e) + 
#'                      1/(n21 + incr.c) + 1/(n22 + incr.c)))
#'
#'Pc <- ((n12 + incr.c)/(AD.s1$nc + 2*incr.c))
#'
#'AD.s1$logitPc <- log(Pc/(1-Pc))
#'
#'AD.s1$Ntotal <- AD.s1$nc + AD.s1$nt
#' rm(list=c("Pc", "n11","n12","n21","n22","incr.c", "incr.e"))
#'
#'
#'dat.points <- data.frame(TE = AD.s1$TE, logitPc = AD.s1$logitPc, N.total = AD.s1$Ntotal)
#'###################################################################
#'
#'res.s1 <- metarisk(AD.s1, two.by.two = FALSE, sigma.1.upper = 5,
#'                   sigma.2.upper = 5, 
#'                   sd.Fisher.rho = 1.5)
#'
#'print(res.s1, digits = 4)
#'
#'
#'library(R2jags)
#'attach.jags(res.s1)
#'eta.hat <- mu.2 + rho*sigma.2/sigma.1*(mu.T - mu.1)
#'mean(eta.hat)
#'sd(eta.hat)
#'
#'OR.eta.hat <- exp(eta.hat)
#'mean(OR.eta.hat)
#'sd(OR.eta.hat)
#'quantile(OR.eta.hat, probs = c(0.025, 0.5, 0.975))
#'
#'ind.random <- sample(1:18000, size = 80, replace = F)
#'
#'##############################################################               
#' p1 <- ggplot(dat.points, aes(x = logitPc, y = TE, size = N.total)) +
#'       xlab("logit Baseline Risk")+
#'       ylab("log(Odds Ratio)")+
#'       geom_point(shape = 21, colour = "blue") + scale_size_area(max_size=10)+
#'       scale_x_continuous(name= "Rate of The Control Group (logit scale)", 
#'                        limits=c(-2, 5)) +
#'      geom_vline(xintercept = mu.T, colour = "blue", size = 1, lty = 1) +
#'        geom_hline(yintercept = eta.true, colour = "blue", size = 1, lty = 1)+
#'          geom_abline(intercept=beta.0[ind.random], 
#'                    slope=beta.1[ind.random],alpha=0.3,
#'                    colour = "grey", size = 1.3, lty = 2)+
#'          geom_abline(intercept = mean(beta.0[ind.random]), 
#'          slope = mean(beta.1[ind.random]),
#'          colour = "black", size = 1.3, lty = 2)+
#'      geom_abline(intercept = mu.2.true, slope = sigma.2.true/sigma.1.true * rho.true,
#'      colour = "blue", size = 1.2)+ theme_bw() 
#'      
#'      
#'# Posterior of eta.hat
#'
#'eta.df <- data.frame(x = OR.eta.hat)
#'
#'
#'p2 <- ggplot(eta.df, aes(x = x)) + 
#'  xlab("Odds Ratio") +
#'  ylab("Posterior distribution")+ 
#'  geom_histogram(fill = "royalblue", colour = "black", alpha= 0.4, bins=80) +
#'  geom_vline(xintercept = exp(eta.true), colour = "black", size = 1.7, lty = 1)+
#'  geom_vline(xintercept = c(0.28, 0.22, 0.35), colour = "black", size = 1, lty = 2)+
#'  theme_bw()
#'
#'
#' grid.arrange(p1, p2, ncol = 2, nrow = 1)



#'
#' #Experiment 2: Internal validity bias and assesing conflict of evidence between the RCTs.
#' 
#' 
#' set.seed(2018)
#' ind <- sample(1:35, size=5, replace = FALSE)
#' ind
#'AD.s4.contaminated <- AD.s1[ind,1:4]
#'AD.s4.contaminated$yc <- AD.s1$yt[ind]
#'AD.s4.contaminated$yt <- AD.s1$yc[ind]
#'AD.s4.contaminated$nc <- AD.s1$nt[ind]
#'AD.s4.contaminated$nt <- AD.s1$nc[ind]
#'AD.s4.contaminated <- rbind(AD.s4.contaminated,
#'                          AD.s1[-ind,1:4])
#'                          
#' res.s4 <- metarisk(AD.s4.contaminated, 
#'                    two.by.two = FALSE,
#'                    re = "sm",
#'                    sigma.1.upper = 3,
#'                    sigma.2.upper = 3, 
#'                    sd.Fisher.rho = 1.5,
#'                    df.estimate = TRUE)
#'                    
#'print(res.s4, digits = 4)
#'
#'attach.jags(res.s4)
#'
#'w.s <- apply(w, 2, median)
#'w.u <- apply(w, 2, quantile, prob = 0.75)
#'w.l <- apply(w, 2, quantile, prob = 0.25)
#'
#'studies <- c(ind,c(1,3,4,5,6,8,9,10,11,13,14,17,18,19,20:35))
#'
#'
#'dat.weights <- data.frame(x = studies,
#'                          y = w.s,
#'                         ylo  = w.l,
#'                         yhi  = w.u)
#'                         
#'# Outliers:
#'w.col <- studies %in% ind
#'w.col.plot <- ifelse(w.col, "black", "grey")
#'w.col.plot[c(9,17, 19,27,34,35)] <- "black"
#'
#' w.plot <- function(d){
#'   # d is a data frame with 4 columns
#'   # d$x gives variable names
#'   # d$y gives center point
#'   # d$ylo gives lower limits
#'   # d$yhi gives upper limits
#'   
#'   p <- ggplot(d, aes(x=x, y=y, ymin=ylo, ymax=yhi) )+ 
#'        geom_pointrange( colour=w.col.plot, lwd =0.8)+
#'        coord_flip() + geom_hline(yintercept = 1, lty=2)+ 
#'        xlab("Study ID") +
#'        ylab("Scale mixture weights") + theme_bw()
#'        return(p)}
#'        
#'w.plot(dat.weights)
#'
#' #List of other possible statistical models:
#' #    1) Different link functions: logit, cloglog and probit
#'
#' #    2) Different random effects distributions:
#' #       "normal" or "sm = scale mixtures".
#'
#' #    3) For the scale mixture random effects:
#' #       split.w = TRUE => "split the weights".
#'
#' #    4) For the scale mixture random effects:
#' #       df.estimate = TRUE => "estimate the degrees of freedom".
#'
#' #    5) For the scale mixture random effects:
#' #       df.estimate = TRUE => "estimate the degrees of freedom".
#'
#' #    6) For the scale mixture random effects:
#' #       df = 4 => "fix the degrees of freedom to a particual value".
#' #       Note that df = 1 fits a Cauchy bivariate distribution to 
#' #       the random effects.
#' #End of the examples
#' 
#' }
#' 
#' 
#' @import R2jags
#' @import rjags
#'
#' @export
#'

metarisk <- function(data,
                     two.by.two = TRUE,
                     # Arguments for the model:
                     re              = "normal",
                     link            = "logit",
                     # Hyperpriors parameters............................................
                     mean.mu.1       = 0,
                     mean.mu.2       = 0,
                     sd.mu.1         = 1,
                     sd.mu.2         = 1,
                     sigma.1.upper   = 5,
                     sigma.2.upper   = 5,
                     mean.Fisher.rho = 0,
                     sd.Fisher.rho   = 1/sqrt(2),
                     df              = 4,
                     df.estimate     = FALSE,
                     df.lower        = 3,
                     df.upper        = 20,
                     # Split weights
                     split.w         = FALSE,
                     # MCMC setup........................................................
                     nr.chains       = 2,
                     nr.iterations   = 10000,
                     nr.adapt        = 1000,
                     nr.burnin       = 1000,
                     nr.thin         = 1,
                     # Further options to link jags and R ...............................
                     be.quiet        = FALSE,
                     r2jags          = TRUE
                     )UseMethod("metarisk")

#' @export
#' 
metarisk.default <- function(
          # Data
           data,
           two.by.two = TRUE,
          # Arguments for the model:
          re              = "normal",
          # Parametrization...................................................
          link            = "logit",
					# Hyperpriors parameters............................................
					mean.mu.1       = 0,
					mean.mu.2       = 0,
					sd.mu.1         = 1,
					sd.mu.2         = 1,
					sigma.1.upper   = 3,
					sigma.2.upper   = 3,
					mean.Fisher.rho = 0,
					sd.Fisher.rho   = 1/sqrt(2),
					df              = 4,
					df.estimate     = FALSE,
					df.lower        = 3,
					df.upper        = 20,
          # Split weights
          split.w         = FALSE,
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

# Model errors checking-----

if(re=="normal" & split.w==TRUE)stop("Normal random effects and splitting weights are not compatible options")

re.test <- re %in% c("normal", "sm")

if(!re.test)stop("This random effects distribution is not implemented")

link.test <- link %in% c("logit", "cloglog", "probit")

if(!link.test)stop("This link function is not implemented")

	# Setting up hyperparameters ...
	      pre.mu.1 <- 1/(sd.mu.1*sd.mu.1)
	      pre.mu.2 <- 1/(sd.mu.2*sd.mu.2)
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

  # Data errors
  if(yc>nc || yt>nt)stop("the data is inconsistent")
  #if(missing(data))stop("NAs are not alow in this function")

#-----
# Data, initial values and parameters .....................................................
	data.model <-
            list(N = N,
	              yc = yc,
	              nc = nc,
	              yt = yt,
	              nt = nt,
	       mean.mu.1 = mean.mu.1,
	       mean.mu.2 = mean.mu.2,
	        pre.mu.1 = pre.mu.1,
	        pre.mu.2 = pre.mu.2,
	   sigma.1.upper = sigma.1.upper,
	   sigma.2.upper = sigma.2.upper,
	 mean.Fisher.rho = mean.Fisher.rho,
	  pre.Fisher.rho = pre.Fisher.rho
	 )

	if(re == "sm" & df.estimate == TRUE){data.model$df.lower <- df.lower
	                                     data.model$df.upper <- df.upper}
	if(re == "sm" & df.estimate == FALSE){data.model$df <- df}

  # Parameters to monitor ....................................................................
	parameters.model <-
    c("mu.1",
      "mu.2",
      "sigma.1",
      "sigma.2",
      "rho",
      "Odds.pool",
      "Odds.new",
      "P_control.pool",
      "P_control.new",
      "beta.0",
      "beta.1")

# This take the weights for the scale mixture random effects model, etc...

        if(re == "sm" & split.w == TRUE & df.estimate == TRUE)
          parameters.model <- c(parameters.model, "w1", "w2", "p.w1", "p.w2", "df")
        else
          if(re == "sm" & split.w == TRUE & df.estimate == FALSE)
            parameters.model <- c(parameters.model, "w1", "w2", "p.w1", "p.w2")
	      else
          if(re=="sm" & split.w == FALSE & df.estimate == TRUE)
            parameters.model <- c(parameters.model, "w", "p.w", "df")
	     else
	       if(re=="sm" & split.w == FALSE & df.estimate == FALSE)
	         parameters.model <- c(parameters.model, "w", "p.w")


# Model BUGS script

#----
blueprint <- function(link = "logit", re = "normal", split.w = FALSE, df.estimate = FALSE)
{

   if(split.w == FALSE & df.estimate == TRUE ) re <- "sm.df"
   if(split.w == TRUE  & df.estimate == FALSE) re <- "sm.split"
   if(split.w == TRUE  & df.estimate == TRUE ) re <- "sm.split.df"

#----
  # Block for data model ......................................................................
  dm <-
    "
  model
  {
  for(i in 1:N)
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

  link.cloglog <-
    "
  # Random effects model 
  cloglog(pc[i]) <- theta.1[i]
  cloglog(pt[i]) <- theta.2[i] + cloglog(pc[i])
  "

  link.probit <-
    "
  # Random effects model 
  probit(pc[i]) <- theta.1[i]
  probit(pt[i]) <- theta.2[i] + probit(pc[i])
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
  beta.0 <- (mu.2 - beta.1 * mu.1)
	beta.1 <- rho * sigma.2/sigma.1

  "

#----
  re.sm <-
    "
  theta.1[i] ~ dnorm(mu.1, pre.theta.1[i])
  theta.2[i] ~ dnorm(mu.2.1[i], pre.theta.2.1[i])
	
	#Conditional mean
	mu.2.1[i] <- mu.2 + rho * sigma.2 / sigma.1 * (theta.1[i] - mu.1)
	
	# Conditional precision
	pre.theta.2.1[i] <- pre.theta.2[i] / (1 - rho * rho)
	
	lambda[i] ~ dchisqr(df)
    	w[i] <- df/lambda[i]
 
	pre.theta.1[i] <- pre.1 / w[i]
	pre.theta.2[i] <- pre.2 / w[i]

  # Probability of outlier	
	p.w[i] <- step(w[i]-0.99)
  }

  # Hyper priors
  mu.1 ~  dlogis(mean.mu.1, pre.mu.1)
	mu.2 ~  dlogis(mean.mu.2, pre.mu.2)
	
	# Dispersion parameters
	sigma.1 ~ dunif(0, sigma.1.upper)
	sigma.2 ~ dunif(0, sigma.2.upper)
	pre.1 <- 1/(sigma.1*sigma.1)
	pre.2 <- 1/(sigma.2*sigma.2)

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

  lambda.new ~ dchisqr(df)
	w.new <- df/lambda.new
	Omega.new[1:2,1:2] <- inverse(Sigma.new[1:2, 1:2]) / w.new
	
	theta.new[1:2] ~ dmnorm(mu.new[1:2], Omega.new[1:2 ,1:2])

  # Functional parameters
  beta.0 <- (mu.2 - beta.1 * mu.1)
  beta.1 <- rho * sigma.2/sigma.1
  "

#--------------------------------------------------

  re.sm.split <-
    "
  theta.1[i] ~ dnorm(mu.1, pre.theta.1[i])
  theta.2[i] ~ dnorm(mu.2.1[i], pre.theta.2.1[i])
	
	#Conditional mean
	mu.2.1[i] <- mu.2 + rho * sigma.2 / sigma.1 * (theta.1[i] - mu.1)
	
	# Conditional precision
	pre.theta.2.1[i] <- pre.theta.2[i] / (1 - rho * rho)

	pre.theta.1[i] <- pre.1 / w1[i]
	pre.theta.2[i] <- pre.2 / w2[i]

  lambda1[i] ~ dchisqr(df)
  w1[i] <- df/lambda1[i]
	
	lambda2[i] ~ dchisqr(df)
	w2[i] <- df/lambda2[i]

  # Probability of outlier
  p.w1[i] <- step(w1[i]-0.99)
	p.w2[i] <- step(w2[i]-0.99)
}

 # Hyper priors
 mu.1 ~  dlogis(mean.mu.1, pre.mu.1)
 mu.2 ~  dlogis(mean.mu.2, pre.mu.2)

 # Dispersion parameters
 sigma.1 ~ dunif(0, sigma.1.upper)
 sigma.2 ~ dunif(0, sigma.2.upper)
 pre.1 <- 1/(sigma.1*sigma.1)
 pre.2 <- 1/(sigma.2*sigma.2)

 # Correlation
 z ~ dnorm(mean.Fisher.rho, pre.Fisher.rho)
 rho <- 2*exp(z)/(1+exp(z)) - 1

 # Predictions ...
 mu.new[1] <- mu.1
 mu.new[2] <- mu.2

 lambda1.new ~ dchisqr(df)
 lambda2.new ~ dchisqr(df)

 w1.new <- df/lambda1.new
 w2.new <- df/lambda2.new

 Sigma.new[1, 1] <- pow(sigma.1, 2)/w1.new
 Sigma.new[2, 2] <- pow(sigma.2, 2)/w2.new
 Sigma.new[1, 2] <- rho * sigma.1 * sigma.2 / sqrt(w1.new*w2.new)
 Sigma.new[2, 1] <- Sigma.new[1, 2]
 Omega.new[1:2,1:2] <- inverse(Sigma.new[1:2, 1:2])

 theta.new[1:2] ~ dmnorm(mu.new[1:2], Omega.new[1:2 ,1:2])

  # Functional parameters
  beta.0 <- (mu.2 - beta.1 * mu.1)
  beta.1 <- rho * sigma.2/sigma.1
"

#------------------------------------------------------------------------

prior.of.df <- "
# Degrees of freedom
  a <- 1/df.upper
  b <- 1/df.lower
 df <- 1/U
  U ~ dunif(a, b)
"

re.sm.df <- paste(re.sm, prior.of.df)

#------------------------------------------------------------------------

re.sm.split.df <- paste(re.sm.split, prior.of.df)

#------------------------------------------------------------------------

#----

# Block of parameters of interest depending on the links ......................
  par.logit <- "
  # Parameters of interest
  # Pooled summaries ...
       Odds.pool <- ilogit(beta.0)
  P_control.pool <- ilogit(mu.1)

  # Predictive summaries ...
       Odds.new <- ilogit(theta.new[2])
  P_control.new <- ilogit(theta.new[1])
}"

  par.cloglog <- "
  # Parameters of interest
  # Pooled summaries ...
       Odds.pool <- icloglog(beta.0)
  P_control.pool <- icloglog(mu.1)

  # Predictive summaries ...
       Odds.new <- icloglog(theta.new[2])
  P_control.new <- icloglog(theta.new[1])
}"

  par.probit <-
  "
  # Parameters of interest
  # Pooled summaries ...
       Odds.pool <- phi(beta.0)
  P_control.pool <- phi(mu.1)

  # Predictive summaries ...
        Odds.new <- phi(theta.new[2])
    P_control.new <- phi(theta.new[1])
  }
  "

# List of possible models for hte ...

# normal random effects 
m1 <- paste(dm, link.logit,   re.normal, par.logit)
m2 <- paste(dm, link.cloglog, re.normal, par.cloglog)
m3 <- paste(dm, link.probit,  re.normal, par.probit)

# sm random effects
m4 <- paste(dm, link.logit,   re.sm, par.logit)
m5 <- paste(dm, link.cloglog, re.sm, par.cloglog)
m6 <- paste(dm, link.probit,  re.sm, par.probit)

# sm random effects with two t-distributions one for each random effect
m7 <- paste(dm, link.logit,   re.sm.split,  par.logit)
m8 <- paste(dm, link.cloglog, re.sm.split,  par.cloglog)
m9 <- paste(dm, link.probit,  re.sm.split,  par.probit)

# sm random effects
m10 <- paste(dm, link.logit,   re.sm.df, par.logit)
m11 <- paste(dm, link.cloglog, re.sm.df, par.cloglog)
m12 <- paste(dm, link.probit,  re.sm.df, par.probit)

# sm random effects with two t-distributions one for each random effect
m13 <- paste(dm, link.logit,   re.sm.split.df,  par.logit)
m14 <- paste(dm, link.cloglog, re.sm.split.df,  par.cloglog)
m15 <- paste(dm, link.probit,  re.sm.split.df,  par.probit)


#----
  switch(re,
         normal = switch(link,
                           logit = return(m1),
                         cloglog = return(m2),
                          probit = return(m3),
                         stop("The model you requested is not implemented.")
                         ),
             sm = switch(link,
                           logit = return(m4),
                         cloglog = return(m5),
                          probit = return(m6),
                         stop("The model you requested is not implemented.")
                         ),
         sm.split = switch(link,
                           logit = return(m7),
                         cloglog = return(m8),
                          probit = return(m9),
                         stop("The model you requested is not implemented.")
                         ),
           sm.df = switch(link,
                           logit = return(m10),
                         cloglog = return(m11),
                          probit = return(m12),
                         stop("The model you requested is not implemented.")
         ),
         sm.split.df = switch(link,
                           logit = return(m13),
                         cloglog = return(m14),
                          probit = return(m15),
                         stop("The model you requested is not implemented.")
         ),
        stop("The model you requested is not implemented.")
        )
}

#----
model.bugs <- blueprint(link, re, split.w, df.estimate)
model.bugs.connection <- textConnection(model.bugs)

if(r2jags == TRUE){
  # Use R2jags as interface for JAGS ...
  results <- jags(              data = data.model,
                  parameters.to.save = parameters.model,
                            #   inits = inits.model,
                          model.file = model.bugs.connection,
                            n.chains = nr.chains,
                              n.iter = nr.iterations,
                            n.burnin = nr.burnin,
                              n.thin = nr.thin
                       )
  }
  else {
  # Use rjags as interface for JAGS ...
  # Send the model to JAGS, check syntax, run ...
	jm <- jags.model(file     = model.bugs.connection,
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
   cat("The plot functions for output analysis are not implemented in this jarbes version", "\n")
}

# Close text conection
close(model.bugs.connection)

# Extra outputs that are linked with other functions
results$link <- link
results$re <- re
results$data <- data
results$two.by.two <- two.by.two
#results$r2jags <- r2jags
results$split.w <- split.w
#results$df.estimate <- df.estimate

class(results) <- c("metarisk")

return(results)
}



#' Generic print function for metarisk object in jarbes.
#' 
#' @param x The object generated by the function metarisk.
#'
#' @param digits The number of significant digits printed. The default value is 3.
#'
#' @param ... \dots
#'
#' @export

print.metarisk <- function(x, digits, ...)
{
  print(x$BUGSoutput,...)
  
}





#' Generic summary function for metarisk object in jarbes
#' @param object The object generated by the metarisk function.
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
summary.metarisk <- function(object, digits = 3,  intervals = c(0.025, 0.5, 0.975), ...)
{
  
  print(object$BUGSoutput, digits= digits, intervals = intervals, ...)
  
}



#' Generic plot function for metarisk object in jarbes.
#' 
#' @param x The object generated by the metarisk function.
#' 
#' @param ... \dots
#' 
#' @export
#' 
plot.metarisk <- function(x, ...)
{
  plot(x)
}
















