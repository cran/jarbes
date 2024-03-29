% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.hmr.R
\name{plot.hmr}
\alias{plot.hmr}
\title{Generic plot function for hmr object in jarbes.}
\usage{
\method{plot}{hmr}(
  x,
  x.lim = c(-5, 2.8),
  y.lim = c(-2, 1),
  x.lab = "Rate of The Control Group (logit scale)",
  y.lab = "No improvement <- Treatment effect -> Improvement",
  title.plot = "Treatment Effect Against Baseline Risk",
  names = NULL,
  name.side = NULL,
  ...
)
}
\arguments{
\item{x}{The object generated by the hmr function.}

\item{x.lim}{Numeric vector of length 2 specifying the x-axis limits.}

\item{y.lim}{Numeric vector of length 2 specifying the y-axis limits.}

\item{x.lab}{Text with the label of the x-axis.}

\item{y.lab}{Text with the label of the y-axis.}

\item{title.plot}{Text for setting a title in the plot.}

\item{names}{Add IPD names to the plot.}

\item{name.side}{Set the side of each name in the plot relative to the vertical line.}

\item{...}{\dots}
}
\description{
Generic plot function for hmr object in jarbes.
}
