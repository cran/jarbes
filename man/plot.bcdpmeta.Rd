% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.bcdpmeta.R
\name{plot.bcdpmeta}
\alias{plot.bcdpmeta}
\title{Generic plot function for bcdpmeta object in jarbes.}
\usage{
\method{plot}{bcdpmeta}(
  x,
  x.lim = c(-3, 3),
  y.lim = c(0, 2),
  x.lab = "Treatment Effect: log(OR)",
  y.lab = "Posterior",
  title.plot.1 = "Model Components",
  title.plot.2 = "Bias Corrected Meta-Analysis",
  ...
)
}
\arguments{
\item{x}{The object generated by the bcmeta function.}

\item{x.lim}{Numeric vector of length 2 specifying the x-axis limits.}

\item{y.lim}{Numeric vector of length 2 specifying the y-axis limits.}

\item{x.lab}{Text with the label of the x-axis.}

\item{y.lab}{Text with the label of the y-axis.}

\item{title.plot.1}{Text for the posterior means by component (biased and bias corrected).}

\item{title.plot.2}{Text for the posterior mean (pooled and predictive).}

\item{...}{\dots}
}
\description{
Generic plot function for bcdpmeta object in jarbes.
}
