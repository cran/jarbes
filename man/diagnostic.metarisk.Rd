% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostic.metarisk.R
\name{diagnostic.metarisk}
\alias{diagnostic.metarisk}
\title{Diagnostic function for metarisk object in jarbes}
\usage{
\method{diagnostic}{metarisk}(
  object,
  median.w = 1.5,
  study.names,
  size.forest = 0.4,
  lwd.forest = 0.2,
  shape.forest = 23,
  ...
)
}
\arguments{
\item{object}{The object generated by the function hmr.}

\item{median.w}{Change color if median of a weight > median.w. The default value is 1.5.}

\item{study.names}{Character vector containing names of the studies used.}

\item{size.forest}{Size of the center symbol mark in the forest-plot lines}

\item{lwd.forest}{Thickness of the lines in the forest-plot}

\item{shape.forest}{Type of symbol for the center mark in the forest-plot lines}

\item{...}{\dots}
}
\description{
This function performers a specially designed diagnostic for a metarisk object
}
