#' @title Bayesian meta-analysis for cross-design synthesis.
#'
#'
#' @description This function performers a Bayesian cross-design synthesis. The function fits a
#' hierarchical meta-regression model based on a bivariate random effects model.
#' The number of events in the control and treated group are modeled with two
#' conditional Binomial distributions and the random-effects are based on a
#' bivariate scale mixture of Normals. The individual participant data is modeled
#' as a Bayesian logistic regression for participants in the control group.
#' Coefficients in the regression are modeled as exchangeables.
#'
#'
#'
#' @param data            Aggregated data results: a data frame where the first four columns containing the number of events in
#'                        the control group (yc), the number of patients in the control group (nc),
#'                        the number of events in the treatment group (yt) and the number of patients in the
#'                        treatment group (nt). If two.by.two = TRUE a data frame where each line
#'                        contains the trial results with column names: yc, nc, yt, nt.
#'
#' @param two.by.two      If TRUE indicates that the trial results are with names: yc, nc, yt, nt.
#'
#' @param dataIPD         Individual participant data: a data frame where
#'                        the first column is the outcome variable and the
#'                        other columns represent individual participant
#'                        charachteristics.
#'
#'
#' @param re              Random effects distribution for the resulting model. Possible
#'                        values are \emph{normal} for bivariate random effects and \emph{sm}
#'                        for scale mixtures.
#'
#' @param link            The link function used in the model. Possible values are
#'                        \emph{logit}, \emph{cloglog} \emph{probit}.
#'
#'
#' @param mean.mu.1       Prior mean of baseline risk, default value is 0.
#'
#' @param mean.mu.2       Prior mean of treatment effect, default value is 0.
#'
#' @param mean.mu.phi     Prior mean of the bias parameter which measures the difference
#'                        between the baseline mean mu.1 and the intercept parameter of
#'                        the logistic regression of the individual participant data.
#'                        The defalut vaule is 0.
#'
#' @param sd.mu.1         Prior standard deviation of mu.1, default value is 1. The default prior of mu.1 is a
#'                        logistic distribution with mean 0 and dispersion 1. The implicit prior for mu.1 in
#'                        the probability scale is a uniform between 0 and 1.
#'
#' @param sd.mu.2         Prior standard deviation of mu.2, default value is 1. The default prior of mu.2 is a
#'                        logistic distribution with mean 0 and dispersion 1. The implicit prior for mu.2 in
#'                        the probability scale is a uniform between 0 and 1.
#'
#' @param sd.mu.phi       Prior standard deviation of mu.phi, default value is 1.
#'
#' @param sigma.1.upper   Upper bound of the uniform prior of sigma.1, default value is 5.
#'
#' @param sigma.2.upper   Upper bound of the uniform prior of sigma.2, default value is 5.
#'
#' @param sigma.beta.upper  Upper bound of the uniform prior of sigma.beta, default value is 5.
#'
#' @param mean.Fisher.rho Mean of rho in the Fisher scale, default value is 0.
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
#'
#'
#'
#' @param nr.chains       Number of chains for the MCMC computations, default 5.
#'
#' @param nr.iterations   Number of iterations after adapting the MCMC, default is 10000. Some models may need more iterations.
#'
#' @param nr.adapt        Number of iterations in the adaptation process, default is 1000. Some models may need more iterations during adptation.
#'
#' @param nr.burnin       Number of iteration discarded for burnin period, default is 1000. Some models may need a longer burnin period.
#'
#' @param nr.thin         Thinning rate, it must be a positive integer, the default value 1.
#'
#' @param be.quiet        Do not print warning message if the model does not adapt default value is FALSE. If you are not sure about the adaptation period choose be.quiet=TRUE.
#'
#' @param r2jags          Which interface is used to link R to JAGS (rjags and R2jags) default value is R2Jags TRUE.
#'
#' @return This function returns an object of the class "hmr". This object contains the MCMC output of
#' each parameter and hyper-parameter in the model, the data frame used for fitting the model, the link function,
#' type of random effects distribution and the splitting information for conflict of evidence analysis.
#' The results of the object of the class metadiag can be extracted with R2jags or with rjags. In addition
#' a summary, a print and a plot function are implemented for this type of object.
#'
#' @details The function calculates the implicit hierarchical meta-regression, where the
#' treatment effect is regressed to the baseline risk (rate of events in the control
#' group). The scale mixture weights are used to adjust for internal validity and
#' structural outliers identification. This is used to predict the treatment effect for
#' subgroups of individual participant data.
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
#' @references Verde, P.E, Ohmann, C., Icks, A. and Morbach, S. (2016) Bayesian evidence synthesis and combining randomized and nonrandomized results: a case study in diabetes. Statistics in Medicine. Volume 35, Issue 10, 10 May 2016, Pages: 1654 to 1675.
#'
#' @references Verde, P.E. (2019) The hierarchical meta-regression approach and learning from clinical evidence.
#' Biometrical Journal. 1 - 23.
#'
#'
#' @examples
#' \dontrun{
#'
#' # Example: from Verde 2019, Section 5
#'
#' library(jarbes)
#'
#' data("healing")
#' AD <- healing[, c("y_c", "n_c", "y_t", "n_t")]
#'
#' data("healingipd")
#'
#' IPD <- healingipd[, c("healing.without.amp", "PAD", "neuropathy",
#' "first.ever.lesion", "no.continuous.care", "male", "diab.typ2",
#' "insulin", "HOCHD", "HOS", "CRF", "dialysis", "DNOAP", "smoking.ever",
#' "diabdur", "wagner.class")]
#'
#' mx2 <- hmr(AD, two.by.two = FALSE,
#'            dataIPD = IPD,
#'            re = "sm",
#'            link = "logit",
#'            sd.mu.1 = 2,
#'            sd.mu.2 = 2,
#'            sd.mu.phi = 2,
#'            sigma.1.upper = 5,
#'            sigma.2.upper = 5,
#'            sigma.beta.upper = 5,
#'            sd.Fisher.rho = 1.25,
#'            df.estimate = FALSE,
#'            df.lower = 3,
#'            df.upper = 10,
#'            nr.chains = 1,
#'            nr.iterations = 1500,
#'            nr.adapt = 100,
#'            nr.thin = 1)
#'
#' print(mx2)
#'
#'# This experiment corresponds to Section 4 in Verde (2018).
#'#
#'# Experiment: Combining aggretated data from RCTs and a single
#'# observational study with individual participant data.
#'#
#'# In this experiment we assess conflict of evidence between the RCTs
#'# and the observational study with a partially identified parameter
#'# mu.phi.
#'#
#'# We run two simulated data: 1) mu.phi = 0.5 which is diffucult to
#'# identify. 2) mu.phi = 2 which can be identify. The simulations are
#'# used to see if the hmr() function can recover mu.phi.
#'
#'
#'
#'library(MASS)
#'library(ggplot2)
#'library(jarbes)
#'library(gridExtra)
#'library(mcmcplots)
#'library(R2jags)
#'
#'
#'# Simulation of the IPD data
#'
#'invlogit <- function (x)
#'{
#'  1/(1 + exp(-x))
#'}
#'
#'# Data set for mu.phi = 0.5 .........................................
#'
#' # Parameters Aggregated Data:
#'
#' # mean control
#' pc <- 0.7
#' # mean treatment
#' pt <- 0.4
#' # reduction of "amputations" odds ratio
#' OR <- (pt/(1-pt)) /(pc/(1-pc))
#' # mu_2: treatment effect ...
#' log(OR)
#' mu.2.true <- log(OR)
#' # mu_1
#' mu.1.true <- log(pc/(1-pc)) # Baseline risk
#' mu.1.true
#' #sigma_1 # Between studies variability
#' sigma.1.true <- 1
#' #sigma_2
#' sigma.2.true <- 0.5
#' # rho: correlation between treatment effect and baseline risk
#' rho.true <- -0.5
#'
#' Sigma <- matrix(c(sigma.1.true^2, sigma.1.true*sigma.2.true*rho.true,
#'                   sigma.1.true*sigma.2.true*rho.true, sigma.2.true^2),
#'                   byrow = TRUE, ncol = 2)
#'
#' # Parameters values IPD
#'
#'# Parameters values
#'mu.phi.true <- 0.5
#'beta0 <- mu.1.true + mu.phi.true
#'beta1 <- 2.5
#'beta2 <- 2
#'
#'# Regression variables
#'
#'x1 <- rnorm(200)
#'x2 <- rbinom(200, 1, 0.5)
#'
#'# Binary outcome as a function of "b0 + b1 * x1 + b2 * x2"
#'
#'y <- rbinom(200, 1,
#'          invlogit(beta0 + beta1 * x1 + beta2 * x2))
#'
#'
#'# Preparing the plot to visualize the data
#'jitter.binary <- function(a, jitt = 0.05)
#'
#'  ifelse(a==0, runif(length(a), 0, jitt),
#'         runif(length(a), 1-jitt, 1))
#'
#'
#'
#'plot(x1, jitter.binary(y), xlab = "x1",
#'     ylab = "Success probability")
#'
#'curve(invlogit(beta0 + beta1*x),
#'      from = -2.5, to = 2.5, add = TRUE, col ="blue", lwd = 2)
#'curve(invlogit(beta0 + beta1*x + beta2),
#'      from = -2.5, to = 2.5, add = TRUE, col ="red", lwd =2)
#'legend("bottomright", c("b2 = 0", "b2 = 2"),
#'       col = c("blue", "red"), lwd = 2, lty = 1)
#'
#'noise <- rnorm(100*20)
#'dim(noise) <- c(100, 20)
#'n.names <- paste(rep("x", 20), seq(3, 22), sep="")
#'colnames(noise) <- n.names
#'
#'data.IPD <- data.frame(y, x1, x2, noise)
#'
#' # Aggregated Data
#'
#' # Experiment 1: External validity bias
#
#' theta <- mvrnorm(35, mu = c(mu.1.true, mu.2.true),
#'                   Sigma = Sigma )
#'
#'  plot(theta, xlim = c(-2,3))
#'  abline(v=mu.1.true, lty = 2)
#'  abline(h=mu.2.true, lty = 2)
#'  abline(a = mu.2.true, b=sigma.2.true/sigma.1.true * rho.true, col = "red")
#'  abline(lm(theta[,2]~theta[,1]), col = "blue")
#
#'  # Target group
#'  mu.T <- mu.1.true + 2 * sigma.1.true
#'  abline(v=mu.T, lwd = 3, col = "blue")
#'  eta.true <- mu.2.true + sigma.2.true/sigma.1.true*rho.true* mu.T
#'  eta.true
#'  exp(eta.true)
#'  abline(h = eta.true, lty =3, col = "blue")
#'  # Simulation of each primary study:
#'  n.c <- round(runif(35, min = 30, max = 60),0)
#'  n.t <- round(runif(35, min = 30, max = 60),0)
#'  y.c <- y.t <- rep(0, 35)
#'  p.c <- exp(theta[,1])/(1+exp(theta[,1]))
#'  p.t <- exp(theta[,2]+theta[,1])/(1+exp(theta[,2]+theta[,1]))
#'  for(i in 1:35)
#'  {
#'    y.c[i] <- rbinom(1, n.c[i], prob = p.c[i])
#'    y.t[i] <- rbinom(1, n.t[i], prob = p.t[i])
#'  }
#'
#'  AD.s1 <- data.frame(yc=y.c, nc=n.c, yt=y.t, nt=n.t)
#'
#'
#'  incr.e <- 0.05
#'  incr.c <- 0.05
#'  n11 <- AD.s1$yt
#'  n12 <- AD.s1$yc
#'  n21 <- AD.s1$nt - AD.s1$yt
#'  n22 <- AD.s1$nc - AD.s1$yc
#'  AD.s1$TE <- log(((n11 + incr.e) * (n22 + incr.c))/((n12 + incr.e) * (n21 + incr.c)))
#'  AD.s1$seTE <- sqrt((1/(n11 + incr.e) + 1/(n12 + incr.e) +
#'                        1/(n21 + incr.c) + 1/(n22 + incr.c)))
#'
#'  Pc <- ((n12 + incr.c)/(AD.s1$nc + 2*incr.c))
#'
#'  AD.s1$logitPc <- log(Pc/(1-Pc))
#'
#'  AD.s1$Ntotal <- AD.s1$nc + AD.s1$nt
#'   rm(list=c("Pc", "n11","n12","n21","n22","incr.c", "incr.e"))
#'
#'# Application of HMR .......................................
#'
#'res.s2 <- hmr(AD.s1, two.by.two = FALSE,
#'              dataIPD = data.IPD,
#'              sd.mu.1 = 2,
#'              sd.mu.2 = 2,
#'              sd.mu.phi = 2,
#'              sigma.1.upper = 5,
#'              sigma.2.upper = 5,
#'              sd.Fisher.rho = 1.5)
#'
#'
#'print(res.s2)
#'
#'# Data set for mu.phi = 2 ..................................
#'# Parameters values
#'
#'mu.phi.true <- 2
#'beta0 <- mu.1.true + mu.phi.true
#'beta1 <- 2.5
#'beta2 <- 2
#'
#'# Regression variables
#'x1 <- rnorm(200)
#'x2 <- rbinom(200, 1, 0.5)
#'# Binary outcome as a function of "b0 + b1 * x1 + b2 * x2"
#'y <- rbinom(200, 1,
#'            invlogit(beta0 + beta1 * x1 + beta2 * x2))
#'
#'# Preparing the plot to visualize the data
#'jitter.binary <- function(a, jitt = 0.05)
#'
#'  ifelse(a==0, runif(length(a), 0, jitt),
#'         runif(length(a), 1-jitt, 1))
#'
#'plot(x1, jitter.binary(y), xlab = "x1",
#'     ylab = "Success probability")
#'
#'curve(invlogit(beta0 + beta1*x),
#'      from = -2.5, to = 2.5, add = TRUE, col ="blue", lwd = 2)
#'curve(invlogit(beta0 + beta1*x + beta2),
#'      from = -2.5, to = 2.5, add = TRUE, col ="red", lwd =2)
#'legend("bottomright", c("b2 = 0", "b2 = 2"),
#'       col = c("blue", "red"), lwd = 2, lty = 1)
#'
#'noise <- rnorm(100*20)
#'dim(noise) <- c(100, 20)
#'n.names <- paste(rep("x", 20), seq(3, 22), sep="")
#'colnames(noise) <- n.names
#'
#'data.IPD <- data.frame(y, x1, x2, noise)
#'
#'# Application of HMR ................................................
#'
#'res.s3 <- hmr(AD.s1, two.by.two = FALSE,
#'              dataIPD = data.IPD,
#'              sd.mu.1 = 2,
#'              sd.mu.2 = 2,
#'              sd.mu.phi = 2,
#'              sigma.1.upper = 5,
#'              sigma.2.upper = 5,
#'              sd.Fisher.rho = 1.5
#')
#'
#'print(res.s3)
#'
#'# Posteriors for mu.phi ............................
#' attach.jags(res.s2)
#' mu.phi.0.5 <- mu.phi
#' df.phi.05 <- data.frame(x = mu.phi.0.5)
#'
#' attach.jags(res.s3)
#' mu.phi.1 <- mu.phi
#' df.phi.1 <- data.frame(x = mu.phi.1)
#'
#'
#' p1 <- ggplot(df.phi.05, aes(x=x))+
#'  xlab(expression(mu[phi])) +
#'  ylab("Posterior distribution")+
#'  xlim(c(-7,7))+
#'  geom_histogram(aes(y=..density..),fill = "royalblue",
#'               colour = "black", alpha= 0.4, bins=60) +
#'  geom_vline(xintercept = 0.64, colour = "black", size = 1.7, lty = 2)+
#'  geom_vline(xintercept = 0.5, colour = "black", size = 1.7, lty = 1)+
#'  stat_function(fun = dlogis,
#'                n = 101,
#'                args = list(location = 0, scale = 1), size = 1.5) + theme_bw()
#'
#' p2 <- ggplot(df.phi.1, aes(x=x))+
#'  xlab(expression(mu[phi])) +
#'  ylab("Posterior distribution")+
#'  xlim(c(-7,7))+
#'  geom_histogram(aes(y=..density..),fill = "royalblue",
#'                 colour = "black", alpha= 0.4, bins=60) +
#'  geom_vline(xintercept = 2.2, colour = "black", size = 1.7, lty = 2)+
#'  geom_vline(xintercept = 2, colour = "black", size = 1.7, lty = 1)+
#'  stat_function(fun = dlogis,
#'                n = 101,
#'                args = list(location = 0, scale = 1), size = 1.5) + theme_bw()
#'
#' grid.arrange(p1, p2, ncol = 2, nrow = 1)
#'
#'
#' # Catter plots for regression coefficients ...........................
#' library(rjags)
#' var.names <- names(data.IPD[-1])
#' v <- paste("beta", names(data.IPD[-1]), sep = ".")
#' mcmc.x.2 <- as.mcmc.rjags(res.s2)
#' mcmc.x.3 <- as.mcmc.rjags(res.s3)
#'
#' greek.names <- paste(paste("beta[",1:22, sep=""),"]", sep="")
#' par.names <- paste(paste("beta.IPD[",1:22, sep=""),"]", sep="")
#'
#' caterplot(mcmc.x.2,
#'          parms = par.names,
#'          col = "black", lty = 1,
#'          labels = greek.names,
#'          greek = TRUE,
#'          labels.loc="axis", cex =0.7,
#'          style = "plain",reorder = FALSE,
#'          denstrip = FALSE)
#'
#' caterplot(mcmc.x.3,
#'          parms = par.names,
#'          col = "grey", lty = 2,
#'          labels = greek.names,
#'          greek = TRUE,
#'          labels.loc="axis", cex =0.7,
#'          style = "plain", reorder = FALSE,
#'          denstrip = FALSE,
#'          add = TRUE,
#'          collapse=TRUE, cat.shift=-0.5)
#'
#' abline(v=0, lty = 2, lwd = 2)
#' abline(v =2, lty = 2, lwd = 2)
#' abline(v =2.5, lty = 2, lwd = 2)
#'
#' # End of the examples.
#'
#' }
#'
#' @import R2jags
#' @import rjags
#' @import graphics
#' @import stats
#' @import graphics
#'
#' @export
#'

hmr <- function(data,
                     two.by.two = TRUE,
                     dataIPD,
                     # Arguments for the model:
                     re              = "normal",
                     link            = "logit",
                     # Hyperpriors parameters............................................
                     mean.mu.1       = 0,
                     mean.mu.2       = 0,
                     mean.mu.phi     = 0,
                     sd.mu.1         = 1,
                     sd.mu.2         = 1,
                     sd.mu.phi       = 1,
                     sigma.1.upper   = 5,
                     sigma.2.upper   = 5,
                    sigma.beta.upper = 5,
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
)UseMethod("hmr")

#' @export
#'
hmr.default <- function(
  data,
  two.by.two = TRUE,
  dataIPD,
  # Arguments for the model:
  re              = "normal",
  link            = "logit",
  # Hyperpriors parameters............................................
  mean.mu.1       = 0,
  mean.mu.2       = 0,
  mean.mu.phi     = 0,
  sd.mu.1         = 1,
  sd.mu.2         = 1,
  sd.mu.phi       = 1,
  sigma.1.upper   = 5,
  sigma.2.upper   = 5,
  sigma.beta.upper = 5,
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
  check.model.1 = re=="normal" & split.w==TRUE
  if(check.model.1)stop("Normal random effects and splitting weights are not compatible options")

  re.test <- re %in% c("normal", "sm")
  if(!re.test)stop("This random effects distribution is not implemented")

  link.test <- link %in% c("logit", "cloglog", "probit")
  if(!link.test)stop("This link function is not implemented")

  # Setting up hyperparameters ...
  pre.mu.1 <- 1/(sd.mu.1*sd.mu.1)
  pre.mu.2 <- 1/(sd.mu.2*sd.mu.2)
  pre.mu.phi <- 1/ (sd.mu.phi * sd.mu.phi)
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
  # Increase N for one observacional study

  N <- N + 1

  # Data errors
  # fix  --- failure: length > 1 in coercion to logical ---

  check.dat <- any(yc>nc || yt>nt)
  if(check.dat)stop("the data is inconsistent")
  #if(any(missing(data)))stop("NAs are not alow in this function")

  # Setup IPD data for regression ................................................
  dataIPD$y <- dataIPD[,1]
  #dataIPD <- dataIPD[,-1]

  model.1 <- model.frame(y ~ ., data = dataIPD)
  X.IPD <- model.matrix(model.1, data = dataIPD)
  X.IPD <- X.IPD[ , dimnames(X.IPD)[[2]] != "(Intercept)", drop=FALSE]
  y.0.IPD <- model.response(model.1)

  check.y.IPD <- is.factor(y.0.IPD)
  if(check.y.IPD){y.0.IPD <- unclass(y.0.IPD) - 1}

  K <- dim(X.IPD)[2]
  M <- dim(X.IPD)[1]

  #.............................................................................
  # Data, initial values and parameters  .......................................
  data.model <-
    list(M = M,
         K = K,
         X.IPD = X.IPD,
         y.0.IPD = y.0.IPD,
         N = N,
         yc = yc,
         nc = nc,
         yt = yt,
         nt = nt,
         mean.mu.1 = mean.mu.1,
         mean.mu.2 = mean.mu.2,
         mean.mu.phi = mean.mu.phi,
         pre.mu.1 = pre.mu.1,
         pre.mu.2 = pre.mu.2,
         pre.mu.phi = pre.mu.phi,
         sigma.1.upper = sigma.1.upper,
         sigma.2.upper = sigma.2.upper,
         sigma.beta.upper = sigma.beta.upper,
         mean.Fisher.rho = mean.Fisher.rho,
         pre.Fisher.rho = pre.Fisher.rho
    )

  check.model.2 = re == "sm" & df.estimate == TRUE
  if(check.model.2){
  data.model$df.lower <- df.lower
  data.model$df.upper <- df.upper}

  check.model.3 = re == "sm" & df.estimate == FALSE
  if(check.model.3){data.model$df <- df}

  # Parameters to monitor ....................................................................
  parameters.model <-
    c("mu.1",
      "mu.2",
      "mu.phi",
      "sigma.1",
      "sigma.2",
      "rho",
      "Odds.pool",
      "Odds.new",
      "P_control.pool",
      "P_control.new",
      "beta.0",
      "beta.1",
      "beta.IPD",
      "sigma.beta")

  # This take the weights for the scale mixture random effects model, etc...

  # Model parameters

  check.q1 = re == "sm" & split.w == TRUE & df.estimate == TRUE
  check.q2 = re == "sm" & split.w == TRUE & df.estimate == FALSE
  check.q3 = re=="sm" & split.w == FALSE & df.estimate == TRUE
  check.q4 = re=="sm" & split.w == FALSE & df.estimate == FALSE

  if(check.q1)
    parameters.model <- c(parameters.model, "w1", "w2", "p.w1", "p.w2", "df")
  else
    if(check.q2)
      parameters.model <- c(parameters.model, "w1", "w2", "p.w1", "p.w2")
  else
    if(check.q3)
      parameters.model <- c(parameters.model, "w", "p.w", "df")
  else
    if(check.q4)
      parameters.model <- c(parameters.model, "w", "p.w")


  # Model BUGS script

  #----
  blueprint <- function(link = "logit", re = "normal", split.w = FALSE, df.estimate = FALSE)
  {

    model.q1 = split.w == FALSE & df.estimate == TRUE
    model.q2 = split.w == TRUE  & df.estimate == FALSE
    model.q3 = split.w == TRUE  & df.estimate == TRUE

    if(model.q1) re <- "sm.df"
    if(model.q2) re <- "sm.split"
    if(model.q3) re <- "sm.split.df"

    #----
    # Block for data model ......................................................................
    dm <-
      "
    model
    {
    for(i in 1:(N-1))
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


    # Block for: link = logistic, re = normal and split.w = FALSE
    # Model for individual data ..................................................
    # Random effect for the observational study ..................................


    # Introduced mu.phi...........................................................
    mu.phi  ~ dlogis(mean.mu.phi, pre.mu.phi) # non-informative

    mu.obs <- mu.1 + mu.phi

    # Change mu.1 to mu.obs = mu.1 + mu.phi .......................................
    mu.2.1[N] <- mu.2 + rho * sigma.2 / sigma.1 * (theta.1[N] - mu.obs)
    theta.1[N] ~ dnorm(mu.obs, pre.theta.1)
    theta.2[N] ~ dnorm(mu.2.1[N], pre.theta.2.1)

    # Priors for beta.IPD .........................................................
    for( i in 1:K)
    {
    beta.IPD[i] ~ dnorm(0, pre.beta.IPD)
    }

    pre.beta.IPD <- 1/(sigma.beta*sigma.beta)
    sigma.beta ~ dunif(0, sigma.beta.upper)

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

# Block for: link = logistic, re = sm and split.w = FALSE
# Model for individual data ..................................................
    # Random effect for the observational study ..................................

    # Introduced mu.phi...........................................................
    mu.phi  ~ dlogis(mean.mu.phi, pre.mu.phi) # non-informative


    mu.obs <- mu.1 + mu.phi

    # Structure of w[N]...........................................................
    # Case: random weight same as RCTs
    lambda[N] ~ dchisqr(df)
    w[N] <- df/lambda[N]

    pre.theta.1[N] <- pre.1 / w[N]
    pre.theta.2[N] <- pre.2 / w[N]
    pre.theta.2.1[N] <- pre.theta.2[N] / (1-rho*rho)

    # Probability of outlier.......................................................
    p.w[N] <- step(w[N]-0.99)

    # Change mu.1 to mu.obs = mu.1 + mu.phi .......................................
    mu.2.1[N] <- mu.2 + rho * sigma.2 / sigma.1 * (theta.1[N] - mu.obs)
    theta.1[N] ~ dnorm(mu.obs, pre.theta.1[N])
    theta.2[N] ~ dnorm(mu.2.1[N], pre.theta.2.1[N])

    # Priors for beta.IPD .........................................................
    for( i in 1:K)
    {
    beta.IPD[i] ~ dnorm(0, pre.beta.IPD)
    }

    pre.beta.IPD <- 1/(sigma.beta*sigma.beta)
    sigma.beta ~ dunif(0, sigma.beta.upper)
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

    # Introduced mu.phi...........................................................
    mu.phi  ~ dlogis(mean.mu.phi, pre.mu.phi) # non-informative

    mu.obs <- mu.1 + mu.phi

    # Structure of w[N]...........................................................
    # Case: random weight same as RCTs
    lambda1[N] ~ dchisqr(df)
    w1[N] <- df/lambda1[N]

    lambda2[N] ~ dchisqr(df)
    w2[N] <- df/lambda2[N]

    pre.theta.1[N] <- pre.1 / w1[N]
    pre.theta.2[N] <- pre.2 / w2[N]
    pre.theta.2.1[N] <- pre.theta.2[N] / (1-rho*rho)

    # Probability of outlier.......................................................
    p.w1[N] <- step(w1[N]-0.99)

    # Change mu.1 to mu.obs = mu.1 + mu.phi .......................................
    mu.2.1[N] <- mu.2 + rho * sigma.2 / sigma.1 * (theta.1[N] - mu.obs)
    theta.1[N] ~ dnorm(mu.obs, pre.theta.1[N])
    theta.2[N] ~ dnorm(mu.2.1[N], pre.theta.2.1[N])

    # Priors for beta.IPD .........................................................
    for( i in 1:K)
    {
    beta.IPD[i] ~ dnorm(0, pre.beta.IPD)
    }

    pre.beta.IPD <- 1/(sigma.beta*sigma.beta)
    sigma.beta ~ dunif(0, sigma.beta.upper)
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


    # Block of parameters of interest depending on the links ......................
    par.logit <- "
#..............................................................................
for(i in 1:M){
    y.0.IPD[i] ~ dbern(p0.IPD[i])
    logit(p0.IPD[i]) <- theta.1[N] + inprod(beta.IPD[1:K], X.IPD[i,1:K])
  }

    # Parameters of interest
    # Pooled summaries ...
    Odds.pool <- ilogit(beta.0)
    P_control.pool <- ilogit(mu.1)

    # Predictive summaries ...
    Odds.new <- ilogit(theta.new[2])
    P_control.new <- ilogit(theta.new[1])
    }"

    par.cloglog <- "

#..............................................................................
for( i in 1:M){
    y.0.IPD[i] ~ dbern(p0.IPD[i])
    cloglog(p0.IPD[i]) <- theta.1[N] + inprod(beta.IPD[1:K], X.IPD[i,1:K])
  }

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
#..............................................................................
for( i in 1:M){
    y.0.IPD[i] ~ dbern(p0.IPD[i])
    probit(p0.IPD[i]) <- theta.1[N] + inprod(beta.IPD[1:K], X.IPD[i,1:K])
}

    # Parameters of interest
    # Pooled summaries ...
    Odds.pool <- phi(beta.0)
    P_control.pool <- phi(mu.1)

    # Predictive summaries ...
    Odds.new <- phi(theta.new[2])
    P_control.new <- phi(theta.new[1])
    }
    "

    # List of possible models ...

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

    #...
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

  # Extra outputs that are linked with other functions ...
  IPD.names <- names(dataIPD)
  IPD.names <- IPD.names[IPD.names!="y"]
  results$beta.names <- paste("beta.", IPD.names, sep = "")

  results$link <- link
  results$re <- re
  results$data <- data
  results$two.by.two <- two.by.two
  #results$r2jags <- r2jags
  results$split.w <- split.w
  #results$df.estimate <- df.estimate

  class(results) <- c("hmr")

  return(results)
  }




#' Generic print function for hmr object in jarbes.
#' @param x The object generated by the function hmr.
#'
#' @param digits The number of significant digits printed. The default value is 3.
#'
#' @param intervals A numeric vector of probabilities with values in [0,1]. The default value is
#'                  intervals = c(0.025, 0.5, 0.975).
#' @param ... \dots
#'
#' @export


print.hmr <- function(x, digits = 3,
                        intervals = c(0.025, 0.25, 0.5, 0.75, 0.975), ...)
{
  m <- x
  x <- x$BUGSoutput
  sims.matrix <- x$sims.matrix
  mu.vect <- apply(sims.matrix, 2, mean)
  sd.vect <- apply(sims.matrix, 2, sd)
  int.matrix <- apply(sims.matrix, 2, quantile, intervals)
  if (x$n.chains>1) {
    n.eff <- x$summary[,"n.eff"]
    Rhat <- x$summary[,"Rhat"]
  } else {
    n.eff <- Rhat <- NULL
  }
  summaryMatrix <- t(rbind(mu.vect, sd.vect, int.matrix, Rhat, n.eff))
  rownameMatrix <- rownames(summaryMatrix)

  # Names of the regression coefficients ...
  ind.names <- grep("beta.IPD", rownameMatrix)
  rownameMatrix[ind.names] <- m$beta.names

  rownames(summaryMatrix) <- rownameMatrix


  dev.idx <- match("deviance", rownameMatrix)
  if(any(!is.na(dev.idx))){
    summaryMatrix <- rbind(summaryMatrix[-dev.idx,], summaryMatrix[dev.idx,])
    rownames(summaryMatrix) <- c(rownameMatrix[-dev.idx], rownameMatrix[dev.idx])
  }




  if (!is.null(m$model.file))
    cat("Inference for Bugs model at \"", m$model.file, "\", ",
        sep = "")
  if (!is.null(m$program))
    cat("fit using ", m$program, ",", sep = "")
  cat("\n ", m$n.chains, " chains, each with ", m$n.iter, " iterations (first ",
      x$n.burnin, " discarded)", sep = "")
  if (x$n.thin > 1)
    cat(", n.thin =", x$n.thin)
  cat("\n n.sims =", x$n.sims, "iterations saved\n")
  print(round(summaryMatrix, digits), ...)
  if (x$n.chains > 1) {
    cat("\nFor each parameter, n.eff is a crude measure of effective sample size,")
    cat("\nand Rhat is the potential scale reduction factor (at convergence, Rhat=1).\n")
  }
  if (x$isDIC) {
    msgDICRule <- ifelse(x$DICbyR, "(using the rule, pD = var(deviance)/2)",
                         "(using the rule, pD = Dbar-Dhat)")
    cat(paste("\nDIC info ", msgDICRule, "\n", sep = ""))
    if (length(x$DIC) == 1) {
      cat("pD =", fround(x$pD, 1), "and DIC =", fround(x$DIC,
                                                       1))
    }
    else if (length(x$DIC) > 1) {
      print(round(x$DIC, 1))
    }
    cat("\nDIC is an estimate of expected predictive error (lower deviance is better).\n")
  }
  invisible(x)

}

fround <- function (x, digits) {
  format (round (x, digits), nsmall=digits)
}

pfround <- function (x, digits) {
  print (fround (x, digits), quote=FALSE)
}







#' Generic summary function for hmr object in jarbes
#' @param object The object generated by the hmr function.
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
summary.hmr <- function(object, digits = 3,  intervals = c(0.025, 0.5, 0.975), ...)
{

  print(object$BUGSoutput, digits= digits, intervals = intervals, ...)

}




#' Generic plot function for hmr object in jarbes.
#'
#' @param x The object generated by the hmr function.
#'
#' @param ... \dots
#'
#' @export
#'
plot.hmr <- function(x, ...)
{
  plot(x)
}






