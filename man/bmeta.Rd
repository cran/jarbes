% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmeta.R
\name{bmeta}
\alias{bmeta}
\title{Bayesian Meta-Analysis for Combining Studies}
\usage{
bmeta(
  data,
  mean.mu = 0,
  sd.mu = 10,
  scale.sigma.between = 0.5,
  df.scale.between = 1,
  nr.chains = 2,
  nr.iterations = 10000,
  nr.adapt = 1000,
  nr.burnin = 1000,
  nr.thin = 1,
  be.quiet = FALSE,
  parallel = NULL
)
}
\arguments{
\item{data}{A data frame with at least two columns with the following names:
1) TE = treatment effect,
2) seTE = the standard error of the treatment effect.}

\item{mean.mu}{Prior mean of the overall mean parameter mu, default value is 0.}

\item{sd.mu}{Prior standard deviation of mu, the default value is 10.}

\item{scale.sigma.between}{Prior scale parameter for scale gamma distribution for the
precision between studies. The default value is 0.5.}

\item{df.scale.between}{Degrees of freedom of the scale gamma distribution for the precision between studies.
The default value is 1, which results in a Half Cauchy distribution for the standard
deviation between studies. Larger values e.g. 30 corresponds to a Half Normal distribution.}

\item{nr.chains}{Number of chains for the MCMC computations, default 2.}

\item{nr.iterations}{Number of iterations after adapting the MCMC, default is 10000. Some models may need more iterations.}

\item{nr.adapt}{Number of iterations in the adaptation process, default is 1000. Some models may need more iterations during adptation.}

\item{nr.burnin}{Number of iteration discard for burn-in period, default is 1000. Some models may need a longer burnin period.}

\item{nr.thin}{Thinning rate, it must be a positive integer, the default value 1.}

\item{be.quiet}{Do not print warning message if the model does not adapt. The default value is FALSE. If you are not sure about the adaptation period choose be.quiet=TRUE.}

\item{parallel}{NULL -> jags, 'jags.parallel' -> jags.parallel execution}
}
\value{
This function returns an object of the class "bmeta". This object contains the MCMC
                           output of each parameter and hyper-parameter in the model and
                           the data frame used for fitting the model.
}
\description{
This function performers a Bayesian meta-analysis
}
\details{
The results of the object of the class bcmeta can be extracted with R2jags or with rjags. In addition a summary, a print and a plot functions are
implemented for this type of object.
}
\examples{
\dontrun{
library(jarbes)

#Example: ppvipd

data(ppvipd)
bm1 = bmeta(ppvipd)

summary(bm1)
plot(bm1, x.lim = c(-3, 1), y.lim = c(0, 3))

diagnostic(bm1, study.names = ppvipd$name, post.p.value.cut = 0.1,
           lwd.forest = 1, shape.forest = 4)

# Example: Stemcells

data("stemcells")
stemcells$TE = stemcells$effect.size
stemcells$seTE = stemcells$se.effect

bm2 = bmeta(stemcells)
summary(bm2)
plot(bm2, x.lim = c(-1, 7), y.lim = c(0, 1))

diagnostic(bm2, study.names = stemcells$trial,
           post.p.value.cut = 0.05,
           lwd.forest = 0.5, shape.forest = 4)

diagnostic(bm2, post.p.value.cut = 0.05,
           lwd.forest = 0.5, shape.forest = 4)
}

}
\references{
Verde, P.E. (2021) A Bias-Corrected Meta-Analysis Model for Combining Studies of Different Types and Quality. Biometrical Journal; 1–17.
}
